/* Copyright 2021-2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Modern Electron, Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_SCATTERING_PROCESS_H_
#define WARPX_PARTICLES_COLLISION_SCATTERING_PROCESS_H_

#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_RandomEngine.H>
#include <AMReX_Vector.H>

#include <string>
#include <vector>

enum class ScatteringProcessType {
    INVALID,
    ELASTIC,
    TWOPRODUCT_REACTION,
    CHARGE_EXCHANGE,
    EXCITATION,
    IONIZATION,
    ATTACHMENT,
};

class ScatteringProcess
{
public:
    ScatteringProcess (
                const std::string& scattering_process,
                const std::string& cross_section_file,
                amrex::ParticleReal energy,
                ScatteringAngleModel scattering_angle_model
                );

    template <typename InputVector>
    ScatteringProcess (
                const std::string& scattering_process,
                const InputVector&& energies,
                const InputVector&& sigmas,
                amrex::ParticleReal energy,
                ScatteringAngleModel scattering_angle_model
                );

    ~ScatteringProcess() = default;

    ScatteringProcess (ScatteringProcess const&)            = delete;
    ScatteringProcess& operator= (ScatteringProcess const&) = delete;
    ScatteringProcess (ScatteringProcess &&)                = default;
    ScatteringProcess& operator= (ScatteringProcess &&)     = default;

    /** Read the given cross-section data file to memory.
     *
     * @param cross_section_file the path to the file containing the cross-
     *        section data
     * @param energies vector storing energy values in eV
     * @param sigmas vector storing cross-section values
     *
     */
    static
    void readCrossSectionFile (
                               const std::string& cross_section_file,
                               amrex::Vector<amrex::ParticleReal>& energies,
                               amrex::Gpu::HostVector<amrex::ParticleReal>& sigmas
                               );

    static
    void sanityCheckEnergyGrid (
                                const amrex::Vector<amrex::ParticleReal>& energies
                                );

    struct Executor {
        /** Get the collision cross-section using a simple linear interpolator. If
         * the energy value is lower (higher) than the given energy range the
         * first (last) cross-section value is used.
         *
         * @param E_coll collision energy in eV
         *
         */
        [[nodiscard]]
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal getCrossSection (amrex::ParticleReal E_coll) const
        {
            if (E_coll < m_energy_lo) {
                return m_sigma_lo;
            } else if (E_coll > m_energy_hi) {
                return m_sigma_hi;
            } else {
                // Locate the pair of energy grid points bracketing E_coll using a
                // bisection search. This does not assume an evenly spaced energy grid.
                // `amrex::bisect` returns idx_1 such that
                // m_energies_data[idx_1] <= E_coll < m_energies_data[idx_1 + 1].
                const int idx_1 = amrex::bisect(m_energies_data, 0, m_grid_size - 1, E_coll);
                const int idx_2 = idx_1 + 1;

                // linearly interpolate to the given energy value
                const amrex::ParticleReal temp =
                    (E_coll - m_energies_data[idx_1])
                    / (m_energies_data[idx_2] - m_energies_data[idx_1]);
                return m_sigmas_data[idx_1] + (m_sigmas_data[idx_2] - m_sigmas_data[idx_1]) * temp;
            }
        }

        amrex::ParticleReal* m_energies_data = nullptr;
        amrex::ParticleReal* m_sigmas_data = nullptr;
        int m_grid_size = 0;
        amrex::ParticleReal m_energy_lo, m_energy_hi, m_sigma_lo, m_sigma_hi, m_dE;
        amrex::ParticleReal m_energy_penalty;
        ScatteringProcessType m_type;
        ScatteringAngleModel m_scattering_angle_model;
        /** True if this process creates or removes particles. */
        bool m_produces_products = false;
    };

    [[nodiscard]]
    Executor const& executor () const {
#ifdef AMREX_USE_GPU
        return m_exe_d;
#else
        return m_exe_h;
#endif
    }

    [[nodiscard]] amrex::ParticleReal getCrossSection (amrex::ParticleReal E_coll) const
    {
        return m_exe_h.getCrossSection(E_coll);
    }

    /** Set a one-time multiplicative conversion from the tabulated values to m^2.
     *
     * The unscaled values are retained in double precision. This is important for
     * three-body attachment tables in m^5, whose entries can underflow in a
     * single-precision build before multiplication by the third-body density.
     */
    void setCrossSectionMultiplier (double multiplier);

    [[nodiscard]] amrex::ParticleReal getEnergyPenalty () const { return m_exe_h.m_energy_penalty; }
    [[nodiscard]] amrex::ParticleReal getMinEnergyInput () const { return m_exe_h.m_energy_lo; }
    [[nodiscard]] amrex::ParticleReal getMaxEnergyInput () const { return m_exe_h.m_energy_hi; }
    [[nodiscard]] amrex::ParticleReal getEnergyInputStep () const { return m_exe_h.m_dE; }

    [[nodiscard]] amrex::Vector<amrex::ParticleReal> const& getEnergyGrid () const
    {
        return m_energies;
    }

    [[nodiscard]] amrex::Gpu::HostVector<amrex::ParticleReal> const& getCrossSectionGrid () const
    {
        return m_sigmas_h;
    }

    [[nodiscard]] std::string const& name () const { return m_name; }
    [[nodiscard]] ScatteringProcessType type () const { return m_exe_h.m_type; }
    [[nodiscard]] ScatteringAngleModel scatteringAngleModel () const
    {
        return m_exe_h.m_scattering_angle_model;
    }

    static
    ScatteringProcessType parseProcessType(const std::string& process);

private:
    static void readCrossSectionFileRaw (
        const std::string& cross_section_file,
        amrex::Vector<amrex::ParticleReal>& energies,
        std::vector<double>& sigmas);

    void init (const std::string& scattering_process, amrex::ParticleReal energy,
               ScatteringAngleModel scattering_angle_model);

    std::string m_name;
    amrex::Vector<amrex::ParticleReal> m_energies;
    std::vector<double> m_sigmas_unscaled;

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_energies_d;
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sigmas_d;
    Executor m_exe_d;
#endif
    amrex::Gpu::HostVector<amrex::ParticleReal> m_sigmas_h;
    Executor m_exe_h;
};

#endif // WARPX_PARTICLES_COLLISION_SCATTERING_PROCESS_H_
