# PyMoments

This is a Python library for calculating rotational constants from some standard forms of
input, for example XYZ and Z-matrix representations. Additionally, there is some high-level
functionality that will help streamline certain workflows involving rotational constants;
the most useful probably is the ability to automatically compute isotopic shifts for common
isotopologues.

## Instructions

This package is managed using `poetry`, and so the first and foremost requirement is to install it
using `pip install poetry`. Once you have `poetry` installed, clone this repository, and run `poetry install`.

After installing, you should have access to `pymoments` as both a library, as well as a command-line call
accessible simply by running `pymoments`:

```bash
➜ pymoments --help
Usage: pymoments [OPTIONS] INPUT_FILE

  Read in internal coordinates from an input file, and print out the
  rotational constants, etc. formatted for reading.

  Takes the path to the file as input, as well as `format` argument that
  specifies which parser to use. Currently, three formats are supported:

  1. Legacy internal coordinate format 2. Z-matrix/internal coordinate
  format 3. XYZ format

  Parameters ---------- filepath : str     Path to the legacy ZMAT file

  format : str     Format of the input file; this determines the parser
  used to create the `Molecule` object.

  Raises ------ FileNotFoundError     If the filepath specified does not
  exist

Options:
  --format [legacy|zmat|xyz]
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

