
import click

from pymoments.classes import Molecule


@click.command()
@click.argument(
    "input_file", type=click.File("r")
)
@click.option(
    "--format",
    default="zmat",
    type=click.Choice(["legacy", "zmat", "xyz"], case_sensitive=False),
)
def run_pymoments(input_file, format):
    """
    Read in internal coordinates from an input file, and print
    out the rotational constants, etc. formatted for reading.
    
    Takes the path to the file as input, as well as `format`
    argument that specifies which parser to use. Currently,
    three formats are supported:
    
    1. Legacy internal coordinate format
    
    2. Z-matrix/internal coordinate format
    
    3. XYZ format

    Parameters
    ----------
    filepath : str
        Path to the legacy ZMAT file
    
    format : str
        Format of the input file; this determines the parser
        used to create the `Molecule` object.
    """
    zmat_str = input_file.read()
    if format == "legacy":
        parser = Molecule.from_legacy_zmat
    elif format == "zmat":
        parser = Molecule.from_zmat_string
    elif format == "xyz":
        parser = Molecule.from_xyz
    else:
        raise NotImplementedError("Invalid option for the file format.")
    molecule = parser(zmat_str)
    _ = molecule.orient()
    click.echo(molecule.dump())
