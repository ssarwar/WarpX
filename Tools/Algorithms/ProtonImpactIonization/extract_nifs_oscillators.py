# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Extract the N2/O2 optical tables from an external NIFS-DATA-109 PDF.

This research utility needs pypdf; the production model does not. The source
contains total photoabsorption, not channel-resolved photoelectron spectra.
Keep the original cross section, wavelength, and mass-attenuation columns as
independent transcription checks. PDF page numbers below are one-based.
"""

import argparse
import hashlib
import json
import re
from pathlib import Path

SOURCE_URL = (
    "https://nifs-repository.repo.nii.ac.jp/record/11706/files/NIFS-DATA-109.pdf"
)
NUMBER = re.compile(r"[+-]?\d+\.\d+E[+-]\d{2}")


def rows(text, columns):
    """Read scientific-notation columns, including adjacent PDF text tokens."""
    values = [float(value) for value in NUMBER.findall(text)]
    if len(values) % columns:
        raise ValueError(f"Expected {columns} columns, found {len(values)} values")
    return [values[i : i + columns] for i in range(0, len(values), columns)]


def extract(path):
    from pypdf import PdfReader

    reader = PdfReader(path)
    if len(reader.pages) != 265:
        raise ValueError("Expected the 265-page NIFS-DATA-109 report")
    pages = {i: reader.pages[i - 1].extract_text() for i in range(47, 65)}
    n2_discrete = rows(pages[47].split("Table I.")[1].split("Table II.")[0], 3)
    n2_core = rows(pages[47].split("Table II.")[1].split("Table III.")[0], 3)
    n2_continuum = rows(
        pages[47].split("Table III.")[1] + "".join(pages[i] for i in range(48, 51)),
        5,
    )
    o2_line_numbers = NUMBER.findall(
        pages[53].split("Table I.")[1].split("Table II.")[0]
    )
    if len(o2_line_numbers) != 25:
        raise ValueError("Expected 22 band-wavelength endpoints and one O2 core line")
    o2_core = [[float(value) for value in o2_line_numbers[-3:]]]
    o2_continuum = rows(
        pages[53].split("Table II.", maxsplit=1)[1]
        + "".join(pages[i] for i in range(54, 64)),
        5,
    )
    # These are integrated band strengths, not pointwise oscillator densities.
    # Transcribed from PDF page 53, Table I; retain the interval endpoints.
    o2_bands = [
        [9.75, 10.17, 0.00833],
        [10.17, 10.44, 0.00707],
        [10.44, 10.62, 0.00077],
        [10.62, 10.71, 0.00066],
        [10.71, 10.84, 0.00140],
        [10.84, 10.98, 0.00081],
        [10.98, 11.17, 0.00076],
        [11.17, 11.33, 0.00050],
        [11.33, 11.52, 0.00147],
        [11.52, 11.59, 0.000419],
        [11.59, 12.07, 0.005496],
    ]
    result = {
        "source": {
            "title": "Oscillator strength spectra and related quantities of 9 atoms "
            "and 23 molecules over the entire energy region",
            "authors": "Sakamoto et al.",
            "report": "NIFS-DATA-109",
            "published": "2010-04-08",
            "retrieved": "2026-09-09",
            "url": SOURCE_URL,
            "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            "kind": "evaluated total photoabsorption, based on Berkowitz (2002)",
            "normalization": "selected using TRK and polarizability sum rules; "
            "these checks are not independent validation",
            "continuum_columns": [
                "loss_eV",
                "df_dloss_per_eV",
                "absorption_Mb",
                "mass_attenuation_cm2_per_g",
                "wavelength_angstrom",
            ],
            "line_columns": ["loss_eV", "oscillator_strength", "wavelength_angstrom"],
            "band_columns": ["lower_eV", "upper_eV", "integrated_strength"],
        },
        "N2": {
            "pages": [47, 48, 49, 50, 51, 52],
            "electrons": 14,
            "molar_mass": 28.0134,
            "first_threshold_eV": 15.58,
            "lines": n2_discrete + n2_core,
            "bands": [],
            "continuum": n2_continuum,
            "continuum_gaps_eV": [],
        },
        "O2": {
            "pages": list(range(53, 65)),
            "electrons": 16,
            "molar_mass": 31.9988,
            "first_threshold_eV": 12.07,
            "lines": o2_core,
            "bands": o2_bands,
            "continuum": o2_continuum,
            # The intervening strength is tabulated separately in Table I.
            "continuum_gaps_eV": [[9.75, 12.07]],
        },
    }
    for target in ("N2", "O2"):
        table = result[target]["continuum"]
        if any(b[0] <= a[0] for a, b in zip(table, table[1:])):
            raise ValueError(f"{target}: non-increasing continuum energies")
        if any(value < 0 for row in table for value in row):
            raise ValueError(f"{target}: negative tabulated quantity")
        for energy, density, cross_section, attenuation, wavelength in table:
            if abs(energy * wavelength / 12398.4 - 1) > 1e-4:
                raise ValueError(f"{target}: inconsistent wavelength at {energy} eV")
            if density and abs(cross_section / density / 109.76097 - 1) > 1e-4:
                raise ValueError(f"{target}: inconsistent optical units at {energy} eV")
            expected = (
                cross_section * 1e-18 * 6.02214179e23 / result[target]["molar_mass"]
            )
            if expected and abs(attenuation / expected - 1) > 1e-4:
                raise ValueError(f"{target}: inconsistent mass units at {energy} eV")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = extract(args.pdf)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    for target in ("N2", "O2"):
        data = result[target]
        print(
            target,
            len(data["continuum"]),
            "continuum rows,",
            len(data["lines"]),
            "lines",
        )
    print("Source SHA-256:", result["source"]["sha256"])


if __name__ == "__main__":
    main()
