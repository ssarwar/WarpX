Sources and preprocessing
=========================

* Schmalzried, IAA thesis (2023): Eqs. 11.21b, 11.22, 11.24, 11.35,
  11.40; database construction and comparison chapters. The supplied elmolcs
  snapshot supplies the inclusive elastic rate and normalized elastic DCS.
* Kutz and Meyer, Physical Review A 51, 3819 (1995): the supplied elementary
  N2 integral tables. The generated numerical ``iaa`` N2 rotational tables
  with the incorrect resonance location are not used.
* Read and Andrick, Journal of Physics B 4, 911 (1971): short-lived, nearly
  pure resonance angular tensors. Their use for all vibrationally elastic
  scattering without a direct background is not justified.
* Jung et al., Journal of Physics B 15, 3535 (1982): Fig. 5 and Table 1,
  vibrationally elastic N2 at 2.22/2.47 eV and 500 +/- 30 K. The seven-angle
  digitization has finite reading uncertainty. Reverse branches in the
  experimental spectral fit were constrained; they are not independent
  tests of detailed balance.
* Morrison et al., Physical Review A 55, 2786 (1997): low-energy MERT and
  threshold-corrected calculations. The printed Eqs. 24--25 and Table III,
  evaluated literally in their stated eV units, do not reproduce Figs. 5--6.
  No coefficient or unit repair is guessed. The real elmolcs/Itikawa rates
  are retained.
* Gote and Ehrhardt, Journal of Physics B 28, 3957 (1995): Table 1,
  elementary transfer-rank contributions at 10--200 eV and 10--160 degrees.
  In the source JSON, -1 denotes a published <1% limit. The nominal fit uses
  the midpoint of that interval. Rows that cannot be imposed as normalized
  probabilities over the reported channels are replaced by interpolation
  between consistent neighboring rows. This is a stated data-completion
  choice, not a claim that the paper supplies exact normalized probabilities.
  For example, the printed 200 eV, 50 degree row already exceeds 100% in
  its first three entries.
* Bhattacharyya and Goswami, Physical Review A 28, 713 (1983): Tables II--III,
  potential B with cutoff radius 2 a0. Columns in the JSON are N=1->1,
  N=1->3, and the rotationally inclusive total, in a0^2. The momentum-transfer
  columns have the same ordering. This is a theoretical model without
  exchange and with a truncated potential expansion; it is not a complete
  experimental rotational DCS set.

The N2 Jung JSON contains nonnegative elementary tensor fractions inferred
from the digitized thermal branches by Clebsch--Gordan recoupling. It is used
with the explicit missing-angle correction and the Kutz--Meyer energy prior
in the reference code. These inferred fractions are not raw measurements.

The O2 completion preserves the published 1->3 integral and first angular
moment. The unresolved higher-rank distribution uses the spectator prior and
an entropy projection under the IAA inclusive angular marginal. This extra
closure is stated explicitly in the model; it is not attributed to the paper.

The canonical IAA rotational constants are retained: 1.998 cm^-1 for N2 and
1.438 cm^-1 for O2. The scattering distances are 2.068 a0 for the N2
Kutz--Meyer/spectator construction and 2.281 a0 for O2. They are source
scattering parameters; the rotational constant is not silently recomputed
from the scattering distance. Ground-relative energies are used for recoil
metadata, including the N=1 ground level of O2.
