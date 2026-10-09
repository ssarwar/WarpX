# Numerical inputs and source provenance

The September 2026 source search and assessments are recorded in the
[dataset audit](history/reports/DATASET_AUDIT.md), [final report](history/reports/FINAL_PJG.md)
and [theory chapter](../../../../Docs/source/theory/multiphysics/proton_impact_ionization.rst).
This archive freezes the files actually used, rather than replacing them
with whatever a mutable download URL supplies later. Retrieval dates are
not measurement dates. No new literature search or refit is claimed by this
archival update.

## Proton-impact data

The current [experimental.py](../experimental.py) and
[source_datasets.py](../source_datasets.py) contain the numerical proton
constraints, figure calibration, digitized coordinates, units and source
PDF hashes. The original research versions are also preserved in
`history/preproduction-source.tar.gz`.

| Source | Use and qualifications |
| --- | --- |
| [Porter, Jackman and Green (1976)](https://doi.org/10.1063/1.432812), Table III, Figs. 5--6 | Printed model parameters; Fig. 5 N2 measured markers; Fig. 6 inconsistency investigation |
| [Garvey, Porter and Green (1977), Ref. 26](https://doi.org/10.1063/1.323427) | Erratum defining the original-model baseline |
| [Crooks and Rudd (1971), Table I](https://doi.org/10.1103/PhysRevA.3.1628) | N2/O2 measured electron-production totals at 50--300 keV; shared normalization uncertainty |
| [Toburen (1971)](https://doi.org/10.1103/PhysRevA.3.216) | N2 measured spectral markers reproduced in PJG Fig. 5 |
| [Rudd (1979), Table I and Fig. 7](https://doi.org/10.1103/PhysRevA.20.787) | N2 totals, mean secondary energies and measured spectral comparisons; scale tied to the 1971 measurements |
| [Rudd et al. (1983), Table V](https://doi.org/10.1103/PhysRevA.28.3244) | Authors' fitted negative-charge yield, 18 rows through 3 MeV; experiment covers 5--4000 keV; extrapolated 5 MeV row excluded |
| [Rudd et al. (1985)](https://doi.org/10.1103/RevModPhys.57.965) | Recommended total curve for comparison and earlier fits, not independent raw measurements |
| [Rudd et al. (1992)](https://doi.org/10.1103/RevModPhys.64.441) | Recommended SDCS formula as a weak interpolation guide; orbital/optical comparisons and Bethe-limit discussion |
| [Cheng, Rudd and Hsu (1989), Fig. 1](https://doi.org/10.1103/PhysRevA.40.3599) | O2 digitized measured shapes at 7.5, 50 and 150 keV, with the authors' recommended-total normalization retained |

Tabulated numerical facts and our digitizations are included; full articles
and their page images are not. Electron yield is not positive-ion yield:
capture, fragmentation and multiplicity prevent an automatic conversion.
The final objective profiles correlated spectral normalizations and downweights
low-secondary-energy markers; it is not a chi-squared fit with invented
independent error bars. The earlier "2--100 eV" Rudd SDCS comparison was
against the analytic recommendation, not a newly obtained experimental table.
For optical comparisons the abscissa is energy loss W, not secondary energy T.

## Archived optical and stopping inputs

| File | Source, role and units |
| --- | --- |
| [nifs-optical.json](inputs/nifs-optical.json) | N2/O2 numerical extraction from [Sakamoto et al., NIFS-DATA-109 (2010)](https://nifs-repository.repo.nii.ac.jp/records/11706); 239/729 continuum rows, lines and bands; per-column units and PDF hash in `source` |
| [leiden-o2.txt](inputs/leiden-o2.txt) | [Leiden O2 file](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/O2/O2.txt), compiled 2020-03-10; wavelength nm, absorption/dissociation/ionization cm2; ionization column supplies the final O2 optical constraint |
| [leiden-n2-0.1nm.txt](inputs/leiden-n2-0.1nm.txt) | [Leiden coarse N2 file](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/N2/N2_0.1nm.txt); alternative optical comparison, not the final N2 anchor |
| [pstar_reference.json](inputs/pstar_reference.json) | Frozen copy of the source-rounded [NIST PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) nitrogen/oxygen gas tables retrieved 2026-09-07; seven named columns, 132 rows per gas; electronic stopping in MeV cm2/g |
| [mahla-o2-o3-supplement.zip](inputs/mahla-o2-o3-supplement.zip) | Mahla and Mehnen (2025), [numerical supplement v1](https://doi.org/10.60893/figshare.jcp.30757583.v1); O2/O3 theoretical photoionization tables, retained unchanged for the source audit |
| [mahla-supplement-metadata.json](inputs/mahla-supplement-metadata.json) and [listing](inputs/mahla-supplement-list.json) | Publisher/figshare metadata identifying version, authors, file digest and license |

The extracted NIFS numerical data retain their source attribution; this is
not a redistribution of the report's prose or page images. NIFS is evaluated
total absorption based on Berkowitz's compilation, not channel-resolved
inclusive electron production. Sum rules used to select the evaluation are
not independent confirmation of the fitted model. Only N2's 25--100 eV
absorption proxy is used in the final optical residual.

For Leiden, cite [Heays, Bosman and van Dishoeck (2017)](https://doi.org/10.1051/0004-6361/201628742)
and [Hrodmarsson and van Dishoeck (2023)](https://doi.org/10.1051/0004-6361/202346645),
as well as the underlying measurements named in each unchanged file header.
The fitted O2 range principally uses [Brion et al. (1979)](https://doi.org/10.1016/0368-2048(79)85032-X)
and [Holland et al. (1993)](https://doi.org/10.1016/0301-0104(93)80148-3),
including Holland's ionization efficiency. These third-party numerical files
are source snapshots with attribution, not data authored by WarpX; no new
license grant for them is asserted here.

The Mahla and Mehnen supplement, titled *cross sections*, is distributed
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as recorded
in its archived metadata. Authors: Sapna Mahla and Bilel Mehnen. Associated
paper: [Photoionization dynamics of O2 and O3 for atmospheric and astrophysical
applications](https://doi.org/10.1063/5.0307721). The ZIP is unmodified. Our
derived comparisons are not endorsed by the authors. Its raw O2 threshold
offset and grouped-partial nonclosure are reported, not silently repaired;
it is a comparison dataset, not the final absolute optical anchor.

The PSTAR comparison uses **electronic** stopping, not electronic plus nuclear
stopping. The source model's effective-pair loss integral contains emitted
electron kinetic energy plus its effective binding contribution, not all
nonionizing excitation loss. Equality is neither fitted nor required. The
source is ballistic, so this diagnostic is not energy removed from the beam.

## Assessed but not final numerical constraints

- [Gallagher et al. (1988)](https://doi.org/10.1063/1.555821) informed channel
  semantics and optical reliability. JILA Data Center Report No. 32 numerical
  partial-channel tables were unavailable; none were fabricated.
- [Sun et al. (2005)](https://doi.org/10.1088/1009-1963/14/7/019) and
  [Lin et al. (2013)](https://doi.org/10.1088/1674-1056/22/2/023404) supplied
  finite-q interpretation and normalization checks. There is no new numerical
  full-GOS fit in this archive. Lin's 0.23/0.91 values are q2, not q.
- [López-Patiño et al. (2016)](https://doi.org/10.1016/j.ijms.2016.05.014)
  was not fitted: the available abstract did not establish conversion of the
  measured ion/fragment channels to inclusive electron production.
- The provisional [NCAR GLOW](https://github.com/NCAR/GLOW/tree/6cb880d0e112810dbe9a79203f0183f9b07e3473)
  optical tables and UBC/LXCat downloads have redistribution restrictions and
  are not included. Their identifying hashes are in the external inventory.
  Derived exploratory results and our scripts are preserved. These inputs
  are unnecessary for reproducing the final model.
- The unused 57-MiB high-resolution Leiden N2 download is identified but not
  duplicated; the actual coarse comparison input is included.

The [external-source inventory](inputs/external-source-inventory.json) records
available local paper/download filenames, sizes and hashes without their
content. A hash identifies the consulted copy; it is not a permission to
redistribute it. Rendered paper pages, OCR, HTML boilerplate, browser cookies
and unrelated resources are excluded. The theory derivation also cites
[Bhabha (1938)](https://doi.org/10.1098/rspa.1938.0017) and
[Salvat and Heredia (2024)](https://doi.org/10.1016/j.nimb.2023.165157);
the independent free-collision checks are retained in the current tests.
