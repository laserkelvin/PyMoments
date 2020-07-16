# PyMoments

![Poetry CI build](https://github.com/laserkelvin/PyMoments/workflows/Poetry%20CI%20build/badge.svg)
[![DOI](https://zenodo.org/badge/280154706.svg)](https://zenodo.org/badge/latestdoi/280154706)

This is a Python library for calculating rotational constants from some standard forms of
input, for example XYZ and Z-matrix representations. Additionally, there is some high-level
functionality that will help streamline certain workflows involving rotational constants;
the most useful probably is the ability to automatically compute isotopic shifts for common
isotopologues.

## Instructions

This package is managed using `poetry`, and so the first and foremost requirement is to install it
using `pip install poetry`. Once you have `poetry` installed, clone this repository, and run `poetry install`.

After installing, you should have access to `pymoments` as both a library and a command-line call
accessible simply by running `pymoments`:

```
➜ pymoments --help
Usage: pymoments [OPTIONS] INPUT_FILE

  Read in internal coordinates from an input file, and print out the
  rotational constants, etc. formatted for reading.

  Takes the path to the file as input, as well as `format` argument that
  specifies which parser to use. Currently, three formats are supported:

  1. Legacy internal coordinate format

  2. Z-matrix/internal coordinate format

  3. XYZ format

  For the latter two formats, there is additional support for automatic
  generation of isotopologues for a molecule. The options `iso` and
  `abundance` control whether these routines are used, and in the latter
  case, what the minimum fractional abundance to consider. For example, the
  default value is 1e-4, which is sufficiently low to account for deuterium.

  Parameters ---------- filepath : str     Path to the legacy ZMAT file

  format : str     Format of the input file; this determines the parser
  used to create the `Molecule` object.

Options:
  --format [legacy|zmat|xyz]  Specifies the file format to parse.
  --iso                       Controls whether isotopologues are generated.
  --abundance FLOAT           Decimal natural abundance threshold for rare
                              isotopes.

  --help                      Show this message and exit.
```

The minimum input for this command is the filepath to an input file. Additional
options are available, which support different file formats. The preferred way
for `pymoments` is to use the `--format=zmat`, which uses the internal coordinate
representation that looks something like this:

```
O
H 1 0.79
H 1 0.79 2 108.0
```

Cartesian coordinates are also supported:

```
O 0.030541 0.042037 -0.000000
H -0.759459 0.042037 -0.000000
H 0.274665 -0.709298 0.000000
```

You can find these in the `tests` folder. As an example of how to run the command line
too, you can run `pymoments tests/h2o.xyz --format=xyz --iso` in the root directory.