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

If you use PyMoments for your work, consider citing the Zenodo DOI (provided in the badge above).

## Usage

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
too, you can run `pymoments tests/h2o.xyz --format=xyz --iso` in the root directory. This
is the kind of output you can expect:

```
===================== Primary input
Formula: OHH
Masses (AMU): [15.99491462  1.00782503  1.00782503]
Short masses (AMU): [16.  1.  1.]
===================== Parameters
Rotational constants (MHz): [1309350.48351127  613808.06588457  417901.00358681]
Inertial axis vectors:
[[ 8.09016994e-01 -5.87785252e-01  0.00000000e+00]
 [-5.87785252e-01 -8.09016994e-01 -6.12323400e-17]
 [-3.59914664e-17 -4.95380036e-17  1.00000000e+00]]
===================== Derived values
Asymmetry parameter: -0.5605
Inertial defect (amu A**2): 0.0000
Scaling factor: 1.0
```

### Rare isotopologues

When using the command line tool, you can specify the flag `--iso`, which will automatically
generate all isotopologues with natural abundance above a certain value specified with `--abundance=1e-4`.
This value is set by default, and includes deuterium isotopologues. Symmetry is taken into account
in a dumb way: before handing you the isotopologue data, we check that the rotational constants are
unique. In an ideal world, you can probably do a better job with internal coordinate or cartesian
coordinate transforms, but that is left to the reader.

### `pymoments` as a library

The command line tool provides a very quick way to calculate isotopologues, or just to do sanity checks
with the rotational constants and structure. There is additional functionality provided, however,
if you use `pymoments` as a library instead. For example, you can manipulate `Molecule` objects
with math operations and comparisons:

```python
water = Molecule.from_zmat_str("""O
    H 1 0.79
    H 1 0.79 2 108.0"""
    )
com, rotational_constants, pmm = water.orient()

isotopologues = water.generate_isotopologues()

# scalar scaling factor; can also be a NumPy array
scaling = 0.9729
# returns a new list of scaled values
scaled_isotopologues = [isotopologue * scaling for isotopologue in isotopologues]

# another way to obtain scaling factors
d2o = isotopologues[-1]
# this sets the scaling attribute of water; can access through water.scaling
water / d2o
# returns a new list of scaled values
scaled_isotopologues = [isotopologue * water for isotopologue in isotopologues]
```

### Contributing

Open to pull requests and contributions.
