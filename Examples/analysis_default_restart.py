#!/usr/bin/env python3

import argparse
import os

import numpy as np
import yt


def check_restart(filename, tolerance=1e-12):
    """
    Compare output data generated from initial run with output data generated after restart.

    Parameters
    ----------
    filename : str
        Name of the plotfile containing the output data generated after restart.
    tolerance : float, optional (default = 1e-12)
        Relative error between restart and original data must be smaller than tolerance.
    """
    # Load output data generated after restart
    ds_restart = yt.load(filename)

    # yt 4.0+ has rounding issues with our domain data:
    # RuntimeError: yt attempted to read outside the boundaries
    # of a non-periodic domain along dimension 0.
    if "force_periodicity" in dir(ds_restart):
        ds_restart.force_periodicity()

    ad_restart = ds_restart.covering_grid(
        level=0,
        left_edge=ds_restart.domain_left_edge,
        dims=ds_restart.domain_dimensions,
    )

    # Load output data generated from initial run
    benchmark = os.path.join(os.getcwd().replace("_restart", ""), filename)
    ds_benchmark = yt.load(benchmark)

    # yt 4.0+ has rounding issues with our domain data:
    # RuntimeError: yt attempted to read outside the boundaries
    # of a non-periodic domain along dimension 0.
    if "force_periodicity" in dir(ds_benchmark):
        ds_benchmark.force_periodicity()

    ad_benchmark = ds_benchmark.covering_grid(
        level=0,
        left_edge=ds_benchmark.domain_left_edge,
        dims=ds_benchmark.domain_dimensions,
    )

    # GPU redistribution and checkpoint reading need not preserve storage order.
    # Match every attribute by the persistent (creation CPU, particle ID) pair.
    # Check identities exactly before comparing physical values at the existing
    # tolerance, so missing or duplicated particles cannot be hidden by sorting.
    particle_orders = {}
    for species, attribute in ds_benchmark.field_list:
        if attribute != "particle_id":
            continue
        orders, identities = [], []
        for data in (ad_benchmark, ad_restart):
            ids = data[species, "particle_id"].v.reshape(-1)
            cpus = data[species, "particle_cpu"].v.reshape(-1)
            order = np.lexsort((ids, cpus))
            pairs = np.column_stack((cpus[order], ids[order]))
            assert not np.any(np.all(pairs[1:] == pairs[:-1], axis=1)), species
            orders.append(order)
            identities.append(pairs)
        np.testing.assert_array_equal(*identities, err_msg=species)
        particle_orders[species] = orders

    # Loop over all fields (all particle species, all particle attributes, all grid fields)
    # and compare output data generated from initial run with output data generated after restart
    print(f"\ntolerance = {tolerance}")
    print()
    for field in ds_benchmark.field_list:
        dr = ad_restart[field].squeeze().v
        db = ad_benchmark[field].squeeze().v
        if field[0] in particle_orders:
            benchmark_order, restart_order = particle_orders[field[0]]
            db = db.reshape(-1)[benchmark_order]
            dr = dr.reshape(-1)[restart_order]
        assert dr.shape == db.shape, field
        if db.size == 0:
            continue
        error = np.amax(np.abs(dr - db))
        if np.amax(np.abs(db)) != 0.0:
            error /= np.amax(np.abs(db))
        print(f"field: {field}; error = {error}")
        assert error < tolerance
    print()


if __name__ == "__main__":
    # define parser
    parser = argparse.ArgumentParser()

    # add arguments: output file path
    parser.add_argument(
        "--path",
        help="path to output file",
        type=str,
        required=True,
    )

    # add arguments: relative tolerance
    default_tolerance = 1e-12
    parser.add_argument(
        "--rtol",
        help="relative tolerance between restart and original",
        type=float,
        required=False,
        default=default_tolerance,
    )

    # parse arguments
    args = parser.parse_args()

    # compare restart results against original results
    check_restart(filename=args.path, tolerance=args.rtol)
