# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Exercise actual lookup paths and bounded majorants with real O2 source data."""

import argparse
import re
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--program", type=Path, required=True)
p.add_argument("--sampler", type=Path, required=True)
p.add_argument("--data-dir", type=Path, required=True)
p.add_argument("--source-dir", type=Path)
p.add_argument("--output", type=Path, required=True)
args = p.parse_args()
program, sampler = args.program.resolve(), args.sampler.resolve()
root = args.output.resolve()
root.mkdir(parents=True, exist_ok=True)
source = (args.source_dir or args.data_dir).resolve() / "O2/IAA/rotation_1_3.txt"
raw = np.loadtxt(source)
legacy = root / "legacy-source.txt"
legacy.write_text(
    "\n".join(
        line
        for line in source.read_text().splitlines()
        if not any(
            key in line
            for key in ["energy_min_eV", "energy_max_eV", "outside_energy_range"]
        )
    )
    + "\n"
)


def inputs(name, energy=1, table=source, reciprocal=False):
    directory = root / name
    directory.mkdir(exist_ok=True)
    proper = np.sqrt(energy * (energy + 2 * 510998.95069)) / 510998.95069
    text = f"""algo.maxwell_solver = Yee
amrex.the_arena_init_size = 8388608
algo.particle_shape = 1
amr.blocking_factor = 1
amr.max_grid_size = 8
amr.max_level = 0
amr.n_cell = 8
boundary.field_hi = periodic
boundary.field_lo = periodic
boundary.particle_hi = periodic
boundary.particle_lo = periodic
geometry.dims = 1
geometry.prob_hi = 10000
geometry.prob_lo = 0
particles.species_names = electrons
electrons.species_type = electron
electrons.injection_style = nuniformpercell
electrons.num_particles_per_cell_each_dim = 4
electrons.profile = constant
electrons.density = 1
electrons.momentum_distribution_type = gaussian
electrons.ux_m = 0
electrons.uy_m = 0
electrons.uz_m = {proper:.17g}
electrons.ux_th = 0
electrons.uy_th = 0
electrons.uz_th = 0
electrons.do_not_deposit = 1
electrons.do_not_gather = 1
electrons.initialize_self_fields = 0
collisions.collision_names = gas
gas.type = background_mcc
gas.species = electrons
gas.background_mass = {31.9988 * 1.66053906660e-27:.17g}
gas.background_density = 1
gas.background_temperature = 0
warpx.const_dt = 1e-14
warpx.verbose = 0
warpx.random_seed = 2026
max_step = 1
benchmark.steps = 1
benchmark.warmup = 1
"""
    if reciprocal:
        prefix = args.data_dir.resolve() / "O2/IAA"
        text += f'''gas.scattering_processes = elastic
gas.elastic_cross_section = "{source.parent / "elastic.txt"}"
gas.elastic_scattering_angle_model = IAA
gas.elastic_rotation_model = reciprocal_hybrid
gas.elastic_rotation_file = "{prefix / "reciprocal_hybrid_300K/thermal_rotation.rot"}"
gas.elastic_rotational_temperature = 300
'''
    else:
        text += f'''gas.scattering_processes = excitation_rotation
gas.excitation_rotation_cross_section = "{table}"
gas.excitation_rotation_energy = {raw[0, 0]:.17g}
gas.excitation_rotation_scattering_angle_model = isotropic
'''
    path = directory / "inputs"
    path.write_text(text)
    return path


def run(path, overrides=(), failure=False):
    result = subprocess.run(
        [str(program), str(path), *overrides],
        cwd=path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (path.parent / "run.log").write_text(result.stdout)
    if failure:
        assert result.returncode != 0 and re.search(
            "range|metadata|interval|limit", result.stdout, re.I
        ), result.stdout
    else:
        assert result.returncode == 0, result.stdout
    return result.stdout


bounds = []
for linear in [0, 1]:
    for is_legacy in [False, True]:
        path = inputs(
            f"bound-{linear}-{is_legacy}", table=legacy if is_legacy else source
        )
        output = run(
            path, [f"benchmark.linear_selector={linear}", "benchmark.print_majorant=1"]
        )
        value = float(re.search(r"MCC_MAJORANT gas ([0-9.eE+-]+)", output)[1])
        bounds.append(value)
    path = inputs(f"outside-{linear}", energy=40)
    run(
        path,
        [f"benchmark.linear_selector={linear}", "gas.background_density=1e30"],
        failure=True,
    )
    path = inputs(f"reciprocal-outside-{linear}", energy=1.01e9, reciprocal=True)
    run(
        path,
        [f"benchmark.linear_selector={linear}", "gas.background_density=1e30"],
        failure=True,
    )
direct = inputs("direct-outside")
run(direct, ["benchmark.direct_query=40"], failure=True)
for energy in [0, float(raw[0, 0]), 20]:
    path = inputs(f"direct-{energy:g}")
    run(path, [f"benchmark.direct_query={energy:.17g}"])
grid = np.unique(np.r_[raw[:, 0], np.geomspace(raw[0, 0], 20, 20000)])
sigma = np.interp(grid, raw[:, 0], raw[:, 1])
expected = np.max(
    sigma
    * 299792458
    * np.sqrt(grid * (grid + 2 * 510998.95069))
    / (grid + 510998.95069)
)
for bounded, unbounded in [(bounds[0], bounds[1]), (bounds[2], bounds[3])]:
    assert bounded >= expected * (1 - 2e-6) and bounded <= expected * (1 + 2e-5), (
        bounded,
        expected,
    )
    assert unbounded > 100 * bounded, (bounded, unbounded)
assert abs(bounds[0] / bounds[2] - 1) < 2e-6

for name, replacement in [
    ("inverted", "# energy_min_eV = 21"),
    ("unsupported", "# outside_energy_range = zero"),
    ("uncovered", "# energy_max_eV = 40"),
    ("missing", "# energy_max_eV ="),
    ("trailing", "# energy_max_eV = 20 garbage"),
    ("syntax", "# outside_energy_range error"),
]:
    table = root / (name + ".txt")
    key = replacement.split()[1]
    lines = [
        replacement if line.startswith("# " + key + " ") else line
        for line in source.read_text().splitlines()
    ]
    table.write_text("\n".join(lines) + "\n")
    run(inputs(name, table=table), failure=True)

bundle = args.data_dir.resolve() / "O2/IAA/reciprocal_hybrid_300K/thermal_rotation.rot"
base = bundle.read_text().splitlines()
count, parts, _ = map(int, base[3].split())
for name, mutate in [
    (
        "temperature",
        lambda rows: rows.__setitem__(2, "250 " + rows[2].split(" ", 1)[1]),
    ),
    (
        "energy",
        lambda rows: rows.__setitem__(2, rows[2].replace("1000000000", "2000000000")),
    ),
    ("duplicate", lambda rows: rows.__setitem__(5, rows[4])),
]:
    directory = root / ("bundle-" + name)
    directory.mkdir(exist_ok=True)
    lines = list(base)
    mutate(lines)
    index = directory / "thermal_rotation.rot"
    index.write_text("\n".join(lines) + "\n")
    for line in base[4 + count : 4 + count + parts]:
        filename = line.split()[0]
        destination = directory / filename
        if not destination.exists():
            destination.symlink_to(bundle.parent / filename)
    result = subprocess.run(
        [
            str(sampler),
            f"file={index}",
            f"cross_section={source.parent / 'elastic.txt'}",
            f"output={directory / 'samples.txt'}",
            "samples=16",
            "amrex.the_arena_init_size=8388608",
        ],
        cwd=directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (directory / "run.log").write_text(result.stdout)
    assert result.returncode != 0, name

# Only the metadata is changed in this reader fixture. It is not a physical
# 300.15 K dataset and must never be used for the thermal reference checks.
directory = root / "decimal-temperature-metadata"
directory.mkdir(exist_ok=True)
lines = list(base)
lines[2] = "300.15 " + lines[2].split(" ", 1)[1]
index = directory / "thermal_rotation.rot"
index.write_text("\n".join(lines) + "\n")
for line in base[4 + count : 4 + count + parts]:
    filename = line.split()[0]
    destination = directory / filename
    if not destination.exists():
        destination.symlink_to(bundle.parent / filename)
path = inputs("decimal-default", reciprocal=True)
text = path.read_text().replace("gas.elastic_rotational_temperature = 300\n", "")
text = text.replace(
    "gas.background_temperature = 0", "gas.background_temperature = 300.15"
)
text = text.replace(str(bundle), str(index))
path.write_text(text)
run(path)

lookup = next(
    (line for line in base[4 : 4 + count] if line.startswith("angular_lookup ")), None
)
if lookup is not None:
    directory = root / "invalid-lookup"
    directory.mkdir(exist_ok=True)
    (directory / "thermal_rotation.rot").write_text(bundle.read_text())
    offset = int(lookup.split()[3])
    part, position = divmod(offset, 32 * 1024 * 1024)
    filenames = [line.split()[0] for line in base[4 + count : 4 + count + parts]]
    for number, filename in enumerate(filenames):
        destination = directory / filename
        if destination.exists() or destination.is_symlink():
            destination.unlink()
        if number == part:
            shutil.copyfile(bundle.parent / filename, destination)
        else:
            destination.symlink_to(bundle.parent / filename)
    with (directory / filenames[part]).open("r+b") as stream:
        stream.seek(position)
        stream.write(struct.pack("<I", 0xFFFFFFFF))
    result = subprocess.run(
        [
            str(sampler),
            f"file={directory / 'thermal_rotation.rot'}",
            f"cross_section={source.parent / 'elastic.txt'}",
            f"output={directory / 'samples.txt'}",
            "samples=16",
            "amrex.the_arena_init_size=8388608",
        ],
        cwd=directory,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (directory / "run.log").write_text(result.stdout)
    assert result.returncode != 0 and "lookup" in result.stdout.lower()
print(
    "PASS: direct, cached and fallback bounds; finite majorants; malformed inputs",
    bounds,
    flush=True,
)
