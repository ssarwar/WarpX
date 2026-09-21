# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

fluids = Bucket("fluids", species_names=[])
fluids_list = []


def new_fluid_species(name, **parameters):
    """Register a fluid input bucket; parameters use the text-input names."""
    if name in fluids.species_names:
        raise ValueError(f"Fluid species {name!r} is already registered")
    species = Bucket(name, **parameters)
    fluids.species_names.append(name)
    fluids_list.append(species)
    return species
