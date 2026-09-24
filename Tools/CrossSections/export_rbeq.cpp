/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCRBEQ.H"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

/** Export a linear-interpolation table using the sampler's host source
 * definition. */
int
main (int argc, char** argv)
{
    using namespace BackgroundMCCRBEQ;
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: export_rbeq N2|O2 "
                     "iaa_thesis_2023|elmolcs_b8643810 FILE [raw]\n";
        return 1;
    }
    try {
        std::string const target = argv[1];
        if (target != "N2" && target != "O2") {
            throw std::invalid_argument("Unknown target");
        }
        auto const model = parse(argv[2]);
        bool const positive = argc == 4;
        if (!positive && std::string(argv[4]) != "raw") {
            throw std::invalid_argument("The optional argument must be raw");
        }
        auto const parameters = shells(target == "N2", model);
        auto const minimum = parameters.back().binding_energy;
        constexpr double maximum = 1.0e9;
        constexpr double tolerance = 2.0e-4;
        auto const evaluate = [&] (double energy) {
            return crossSection(energy, parameters, positive);
        };
        std::map<double, double> knots;
        for (int i = 0; i <= 128; ++i) {
            auto const energy = minimum * std::pow(maximum / minimum, i / 128.0);
            knots[energy] = evaluate(energy);
        }
        knots[minimum] = 0;
        knots[maximum] = evaluate(maximum);
        for (auto const& shell : parameters) {
            knots[shell.binding_energy] = evaluate(shell.binding_energy);
            auto const threshold = findPositiveThreshold(
                shell, [&] (double e) { return rbeqTerms(e, shell).total > 0; });
            if (threshold > shell.binding_energy * (1 + 1.0e-6)) {
                knots[threshold] = evaluate(threshold);
            }
        }
        double peak = 0;
        for (auto const& entry : knots) {
            peak = std::max(peak, std::abs(entry.second));
        }
        for (auto left = knots.begin(); std::next(left) != knots.end();) {
            auto const right = std::next(left);
            bool split = false;
            for (auto const fraction : {.25, .5, .75}) {
                auto const energy = left->first + fraction * (right->first - left->first);
                auto const actual = evaluate(energy);
                auto const linear = (1 - fraction) * left->second + fraction * right->second;
                split = split || std::abs(actual - linear) >
                                     tolerance * std::max(std::abs(actual), peak * 1.0e-8);
            }
            if (split && right->first - left->first > 2.0e-6 * left->first) {
                auto const midpoint = .5 * (left->first + right->first);
                knots[midpoint] = evaluate(midpoint);
            } else {
                ++left;
            }
        }
        std::ofstream output(argv[3]);
        if (!output) {
            throw std::runtime_error("Cannot open output");
        }
        output << "# rbeq_model = " << name(model) << '\n'
               << "# rbeq_normalization = " << (positive ? "positive_part" : "raw_signed") << '\n'
               << "# rbeq_target = " << target << '\n'
               << "# E [eV], sigma [m^2]; linear interpolation; maximum source "
                  "energy 1 GeV\n"
               << "# Relative interpolation target 0.0002; floor 1e-8 of peak\n"
               << std::setprecision(17);
        for (auto const& entry : knots) {
            output << entry.first << ' ' << entry.second << '\n';
        }
        if (!output) {
            throw std::runtime_error("Cannot write output");
        }
        std::cout << target << ' ' << name(model) << ": " << knots.size() << " knots\n";
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
