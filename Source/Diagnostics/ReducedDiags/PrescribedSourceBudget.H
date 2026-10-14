/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PRESCRIBED_SOURCE_BUDGET_H_
#define WARPX_PRESCRIBED_SOURCE_BUDGET_H_

#include "ReducedDiags.H"

#include <string>
#include <vector>

/** Physical production and energy budgets for each prescribed beam collision source. */
class PrescribedSourceBudget final : public ReducedDiags
{
public:
    explicit PrescribedSourceBudget (std::string const& name);
    void ComputeDiags (int step) override;

private:
    std::vector<std::string> m_sources;
};

#endif
