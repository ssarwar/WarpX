# Overview

This explains how to generate and contribute to the WarpX documentation.
See [How to write documentation](source/developers/how_to_write_the_docs.rst)
for reStructuredText, Sphinx and Doxygen conventions.

## Beam and collision documentation

The added beam and neutral-gas features follow the existing manual's division
between usage, theory and development:

- [Beam-in-air setup](source/usage/beam_air.rst): species choices, PICMI,
  examples, diagnostics, restart and convergence.
- [Collision data](source/usage/collision_data.rst): production directories,
  file formats, units, model choices and manifests.
- [Collision theory](source/theory/multiphysics/collisions.rst): Background MCC,
  electron-air excitation/attachment, electron ionization, elastic/rotation,
  and proton ionization. Related derivations and figures live with their model
  under `source/theory/multiphysics/collisions/`.
- [Prescribed-fluid theory](source/theory/prescribed_fluids.rst): Gaussian
  projection, continuity, quiet sources and immobile charge.
- [Development and verification](source/developers/beam_air/index.rst): runtime
  implementation, data preparation, analytical tests and CPU/GPU evidence.

[The parameter list](source/usage/parameters.rst) is the single reference for
input names, defaults and restrictions. [The Python reference](source/usage/python.rst)
documents the corresponding API. Tool READMEs retain reproduction commands and
dated evidence; historical audit reports do not override the model specification.
Physical cross sections and their provenance belong in warpx-data, not Docs.

## Generating the documentation

### Installing the requirements

Install the Python requirements for compiling the documentation:
```
cd Docs/
python3 -m pip install -r requirements.txt
```

### Compiling the documentation

Still in the `Docs/` directory, type
```
make html
```
You can then open the file `build/html/index.html` with a standard web browser (e.g. Firefox), in order to visualize the results on your local computer.

### Cleaning the documentation

In order to remove all of the generated files, use:
```
make clean
```
