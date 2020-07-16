
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
    help="Specifies the file format to parse."
)
@click.option(
    "--iso",
    default=False,
    is_flag=True,
    help="Controls whether isotopologues are generated."
)
@click.option(
    "--abundance",
    default=1e-4,
    type=float,
    help="Decimal natural abundance threshold for rare isotopes."
)
def run_pymoments(input_file, format, iso, abundance=1e-4):
    """
    Read in internal coordinates from an input file, and print
    out the rotational constants, etc. formatted for reading.
    
    Takes the path to the file as input, as well as `format`
    argument that specifies which parser to use. Currently,
    three formats are supported:
    
    1. Legacy internal coordinate format
    
    2. Z-matrix/internal coordinate format
    
    3. XYZ format
    
    For the latter two formats, there is additional support for automatic
    generation of isotopologues for a molecule. The options `iso` and `abundance`
    control whether these routines are used, and in the latter case, what
    the minimum fractional abundance to consider. For example, the default
    value is 1e-4, which is sufficiently low to account for deuterium.

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
    if iso and format != "legacy":
        isotopologues = molecule.generate_isotopologues(abundance)
        click.echo("\n".join(iso.dump() for iso in isotopologues))
    elif iso and format == "legacy":
        raise NotImplementedError(
            "Isotopologue generation is not supported for legacy ZMAT."
        )
