/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLENUMBER_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLENUMBER_H_

#include "ReducedDiags.H"

#include <string>

/**
 * Compute particle counts or physical charge for kinetic and fluid species.
 * Fluid physical populations are integrated from their density fields; their
 * macroparticle counts are zero. Charge mode includes kinetic ionization levels.
 */
class ParticleNumber : public ReducedDiags
{
public:

    /**
     * constructor
     * @param[in] rd_name reduced diags names
     * @param[in] charge output physical charge instead of macro- and physical-particle counts
     */
    explicit ParticleNumber(const std::string& rd_name, bool charge = false);

    /**
     * Compute total and per-species counts, or charge when charge mode is selected.
     *
     * @param[in] step current time step
     */
    void ComputeDiags(int step) final;

private:
    bool m_charge = false;
};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLENUMBER_H_
