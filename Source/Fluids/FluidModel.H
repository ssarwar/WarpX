/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FLUID_MODEL_H_
#define WARPX_FLUID_MODEL_H_

/** Evolution law for a fluid species. Prescribed species have no momentum mesh. */
enum class FluidModel { ColdRelativistic, RigidBeam, Immobile };

#endif
