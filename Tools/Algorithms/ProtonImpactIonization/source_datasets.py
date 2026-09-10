# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Source-aware research data readers; no production dependencies or downloads.

Keep measurement-derived fits, evaluated spectra, and ab initio calculations
distinct. In particular, Rudd (1983), Table V is calculated from Eq. (17),
and Cheng (1989), Fig. 1 has been adjusted to recommended total cross sections.
Neither supplies independent pointwise absolute normalization measurements.
See the proton_impact_ionization.rst theory documentation for source restrictions.
"""

import hashlib
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from reference import BETHE_CONSTANT, ELECTRON_REST_ENERGY, PROTON_REST_ENERGY

# sigma_photo [Mb] = PHOTO_FACTOR * df/dW [1/eV]. This is a conversion
# constant, not a multiplicity or normalization parameter.
PHOTO_FACTOR = 109.76097
HC_EV_NM = 1239.841984332

RUDD_1983_ENERGIES = 1e3 * np.array(
    [5, 7, 10, 15, 20, 30, 50, 70, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000]
)
# Table V, p. 3253, 10^-20 m^2 = 10^-16 cm^2. Exclude the printed
# 5000-keV row: it extrapolates beyond the experiment's 4000-keV upper limit.
RUDD_1983_FITTED_TOTALS = {
    "N2": np.array(
        [
            2.11e-16,
            2.72e-16,
            3.45e-16,
            4.32e-16,
            4.90e-16,
            5.57e-16,
            5.96e-16,
            5.86e-16,
            5.47e-16,
            4.78e-16,
            4.19e-16,
            3.36e-16,
            2.42e-16,
            1.91e-16,
            1.46e-16,
            1.06e-16,
            0.842e-16,
            0.603e-16,
        ]
    ),
    "O2": np.array(
        [
            1.67e-16,
            2.26e-16,
            3.01e-16,
            3.96e-16,
            4.63e-16,
            5.44e-16,
            5.98e-16,
            5.95e-16,
            5.60e-16,
            4.93e-16,
            4.36e-16,
            3.53e-16,
            2.58e-16,
            2.05e-16,
            1.58e-16,
            1.17e-16,
            0.933e-16,
            0.675e-16,
        ]
    ),
}


def rudd_1983_uncertainty(energy):
    """Published uncertainty of the fitted totals; interpolate in log energy.

    Quoted on p. 3254: 25/20/15/10/8 percent at 5/10/25/100/500 keV,
    and 8 percent above 500 keV. These are not independent standard errors.
    """
    energy = np.asarray(energy, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any((energy < 5e3) | (energy > 4e6)):
        raise ValueError("The Rudd (1983) experiment covers 5--4000 keV")
    return np.interp(
        np.log(energy),
        np.log([5e3, 10e3, 25e3, 100e3, 500e3]),
        [0.25, 0.20, 0.15, 0.10, 0.08],
    )


def source_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass
class OpticalTable:
    """Cross section as a function of photon/energy-loss energy, never T."""

    energy: np.ndarray
    cross_section_Mb: np.ndarray

    def __post_init__(self):
        self.energy = np.asarray(self.energy, dtype=float)
        self.cross_section_Mb = np.asarray(self.cross_section_Mb, dtype=float)
        if self.energy.ndim != 1 or self.cross_section_Mb.shape != self.energy.shape:
            raise ValueError("Expected matching one-dimensional optical columns")
        if len(self.energy) < 2 or np.any(np.diff(self.energy) <= 0):
            raise ValueError("Optical energies must increase strictly")
        if np.any(~np.isfinite(self.energy)) or np.any(self.energy <= 0):
            raise ValueError("Optical energies must be finite and positive")
        if np.any(~np.isfinite(self.cross_section_Mb)) or np.any(
            self.cross_section_Mb < 0
        ):
            raise ValueError("Optical cross sections must be finite and nonnegative")

    def __call__(self, energy, outside="raise"):
        energy = np.asarray(energy, dtype=float)
        if np.any(~np.isfinite(energy)):
            raise ValueError("Requested energies must be finite")
        if outside not in ("raise", "zero"):
            raise ValueError("Expected raise or zero for outside support")
        if outside == "raise" and np.any(
            (energy < self.energy[0]) | (energy > self.energy[-1])
        ):
            raise ValueError("Optical data do not cover the requested energy")
        return np.interp(energy, self.energy, self.cross_section_Mb, left=0, right=0)

    def integrate(self, lower, upper, moment=0):
        """Integrate sigma(W)*W**moment without extrapolating.

        Sixteen-point Gauss integration treats each piecewise-linear interval
        separately, preserving narrow resonances on the supplied grid.
        """
        if not self.energy[0] <= lower < upper <= self.energy[-1]:
            raise ValueError("Integration interval is outside tabulated support")
        x = np.r_[
            lower, self.energy[(self.energy > lower) & (self.energy < upper)], upper
        ]
        nodes, weights = np.polynomial.legendre.leggauss(16)
        energy = x[:-1, None] + np.diff(x)[:, None] * (nodes + 1) / 2
        return float(
            np.sum(self(energy) * energy**moment * weights * np.diff(x)[:, None] / 2)
        )


def load_leiden(path, check_closure=True):
    """Read the four-column wavelength table and check its additive channels.

    Wavelength is in nm; cross sections are in cm^2. Converting the abscissa
    to photon energy does not multiply the cross section by a Jacobian.
    """
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 4 or np.any(~np.isfinite(data)):
        raise ValueError("Expected four finite Leiden columns")
    if np.any(np.diff(data[:, 0]) <= 0) or np.any(data[:, :2] < 0):
        raise ValueError("Invalid wavelengths or cross sections")
    # The O2 dissociation column contains subtraction roundoff as small as
    # -1.23e-32 cm^2 where absorption equals ionization (~1e-17 cm^2).
    # Admit only roundoff relative to the local absorption, not a physical
    # negative cross section. Preserve the source file unchanged.
    roundoff = 100 * np.finfo(float).eps * data[:, 1, None]
    if np.any(data[:, 2:] < -roundoff):
        raise ValueError(
            "Negative channel cross section exceeds floating-point roundoff"
        )
    if np.any(data[:, 3] > data[:, 1] * (1 + 3e-5)):
        raise ValueError("Ionization exceeds absorption")
    if check_closure and not np.allclose(
        data[:, 1], data[:, 2] + data[:, 3], rtol=3e-5, atol=1e-35
    ):
        raise ValueError(
            "Leiden absorption does not equal dissociation plus ionization"
        )
    energy = HC_EV_NM / data[::-1, 0]
    return {
        "absorption": OpticalTable(energy, data[::-1, 1] * 1e18),
        "ionization": OpticalTable(energy, np.maximum(data[::-1, 3], 0) * 1e18),
    }


def load_mahla_o2(path):
    """Read the publisher's O2 tables without extracting arbitrary ZIP paths.

    DOI 10.60893/figshare.jcp.30757583.v1 supplies a total and *three grouped*
    partial cross sections, not seven separately tabulated ionic channels.
    Raw theoretical thresholds are retained; no energy shifts are applied.
    """
    tables = {}
    with ZipFile(path) as archive:
        for key, figure in (
            ("total", "1"),
            ("pi_g", "2a"),
            ("pi_u", "2b"),
            ("sigma_g", "2c"),
        ):
            name = f"Supplemental_Material_O2_O3/O2/Figure_{figure}.txt"
            text = archive.read(name).decode("utf-8-sig")
            data = np.loadtxt(StringIO(text), skiprows=4)
            if data.shape[1] != 2:
                raise ValueError(f"Expected two columns in {name}")
            tables[key] = OpticalTable(*data.T)
    return tables


def cheng_1989_figure1(pixel_shift=(0, 0)):
    """Approximate manual digitization of open-circle O2 SDCS markers only.

    Provenance: PDF page 2, Fig. 1, source SHA-256 and rendering command in
    the theory documentation. Coordinates refer to the 510x900 source-page crop.
    Retain the calibration and raw marker positions for independent review.
    Two-pixel reading sensitivity is separate from experimental uncertainty.
    Overlapping triangles/crosses and theoretical lines are not selected.
    A uniform (dx, dy) pixel shift permits a calibration sensitivity check;
    it is not an additional physical parameter or a statistical error model.
    """
    pixels = {
        7.5e3: [
            [97, 65],
            [140, 64],
            [173, 66],
            [199, 70],
            [219, 77],
            [244, 89],
            [262, 100],
        ],
        50e3: [
            [95, 344],
            [140, 343],
            [173, 340],
            [197, 337],
            [217, 331],
            [239, 329],
            [258, 334],
            [286, 336],
            [321, 355],
            [347, 385],
            [365, 416],
            [379, 439],
            [389, 458],
        ],
        150e3: [
            [94, 629],
            [137, 626],
            [171, 622],
            [198, 619],
            [240, 614],
            [285, 606],
            [345, 611],
            [362, 617],
            [377, 627],
            [407, 666],
        ],
    }
    calibrations = {
        7.5e3: (26, 275, 1.0),
        50e3: (310, 561, 10.0),
        150e3: (595, 846, 20.0),
    }
    rows = []
    for energy, points in pixels.items():
        xy = np.asarray(points, dtype=float) + np.asarray(pixel_shift)
        t = 10 ** ((xy[:, 0] - 26) / 147)
        top, bottom, top_value = calibrations[energy]
        y = top_value * 10 ** (-4 * (xy[:, 1] - top) / (bottom - top))
        equivalent = energy * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY
        # The paper's E=W+I1 uses I1=13.1 eV (Table I), not the 12.07-eV
        # adiabatic threshold. Y=T_equiv*(W+I1)^2*SDCS/(4*pi*a0^2*R^2).
        sdcs = y * BETHE_CONSTANT / equivalent / (t + 13.1) ** 2
        rows.extend(zip(np.full_like(t, energy), t, sdcs, strict=True))
    return np.asarray(rows)


def pjg_1976_figure5(pixel_shift=(0, 0)):
    """Measured N2 markers reproduced in PJG Fig. 5, not PJG's fitted curves.

    Coordinates refer to PDF page 8 rendered at 480 dpi with Poppler crop
    x=2050, y=3300, W=950, H=1750. Source SHA-256:
    79cf4e557b78c019a1fbdcf7d7d3a04dfd338ed2679a10f3ca2a383c1aa7c448.
    The spectra have separate vertically displaced log axes. Each has its
    own 10^-17 cm^2/eV anchor; applying the bottom axis to all is incorrect.
    Return E, T, SDCS, source index (0: Crooks/Rudd, 1: Toburen).
    A coherent pixel shift is a reading sensitivity, not a confidence interval.
    """
    markers = {
        (50e3, 0): [
            (242, 114),
            (350, 110),
            (408, 137),
            (458, 184),
            (488, 220),
            (530, 287),
            (565, 386),
            (596, 492),
        ],
        (100e3, 0): [
            (242, 318),
            (349, 308),
            (405, 335),
            (454, 389),
            (488, 416),
            (522, 456),
            (565, 527),
            (594, 592),
            (635, 728),
        ],
        (300e3, 0): [
            (242, 547),
            (352, 537),
            (409, 558),
            (454, 625),
            (488, 650),
            (529, 701),
            (565, 749),
            (594, 810),
            (634, 888),
            (653, 937),
            (670, 978),
        ],
        (300e3, 1): [
            (267, 533),
            (397, 561),
            (454, 610),
            (515, 681),
            (557, 744),
            (588, 799),
            (644, 902),
            (677, 970),
            (718, 1076),
        ],
        (1e6, 1): [
            (267, 785),
            (355, 781),
            (397, 803),
            (454, 875),
            (511, 948),
            (557, 1022),
            (583, 1070),
            (645, 1170),
            (676, 1229),
            (718, 1322),
            (744, 1378),
            (767, 1415),
        ],
    }
    anchors = {50e3: 171, 100e3: 357, 300e3: 542, 1e6: 729}
    rows = []
    for (energy, source), positions in markers.items():
        xy = np.asarray(positions, dtype=float) + np.asarray(pixel_shift)
        t = 10 ** ((xy[:, 0] - 210) / (565 / 3))
        s = 1e-17 * 10 ** ((anchors[energy] - xy[:, 1]) / 184)
        rows.extend(
            zip(np.full_like(t, energy), t, s, np.full_like(t, source), strict=True)
        )
    return np.asarray(rows)


def rudd_1979_figure7(pixel_shift=(0, 0)):
    """Sparse, unambiguous N2 measured markers, used only for validation.

    PDF page 6 at 300 dpi, Poppler crop x=350,y=405,W=860,H=1080.
    Source SHA-256:
    e67147178ae3ce85d6d35caa1e4b4113503ed872b4a2f0e8b655421911ebe275.
    Select the separated 5/20/70-keV curves; do not assign crowded near-zero
    markers to a guessed incident energy. The linear T axis and logarithmic
    SDCS axis have different transformations. The printed unit is m^2/eV.
    """
    markers = {
        5e3: [(119, 212), (141, 288), (163, 357), (205, 486), (289, 690), (395, 948)],
        20e3: [(289, 423), (395, 579), (502, 729), (629, 863), (758, 993)],
        70e3: [
            (164, 174),
            (207, 219),
            (291, 289),
            (397, 373),
            (503, 442),
            (630, 526),
            (758, 612),
        ],
    }
    rows = []
    for energy, positions in markers.items():
        xy = np.asarray(positions, dtype=float) + np.asarray(pixel_shift)
        t = (xy[:, 0] - 78) * 175 / (822 - 78)
        s = 1e4 * 10 ** (-21 - (xy[:, 1] - 167) * 4 / (1005 - 167))
        rows.extend(zip(np.full_like(t, energy), t, s, strict=True))
    return np.asarray(rows)
