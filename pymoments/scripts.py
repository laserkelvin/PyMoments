
from pathlib import Path

import click

from pymoments.classes import Molecule


@click.command()
def moments_from_legacy_zmat(filepath: str):
    """
    Read in internal coordinates from a legacy ZMAT file, and print
    out the rotational constants, etc. in a formatted way.
    
    Takes the path to the file as input. An example of the expected
    structure for the legacy format looks like this:
    
    h2c3s calculations; CPL 326, 530 (2000)    
    H2CCCS 
    6	 1    1
    1	 0    0    0	0.0	        0.0        0.0    31.972070
    2	 1    0    0	1.594000    0.0        0.0    12.000000
    3	 2    1    0	1.275000  180.0        0.0    12.000000
    4    3    2    1    1.328000  180.0        0.0    12.000000
    5    4    3    2    1.096000  122.0        0.0     1.007825
    6    4    3    2    1.096000  122.0      180.0     1.007825
    0

    Parameters
    ----------
    filepath : str
        Path to the legacy ZMAT file

    Raises
    ------
    FileNotFoundError
        If the filepath specified does not exist
    """
    filepath = Path(filepath)
    if filepath.exists():
        zmat_str = filepath.read_text()
    else:
        raise FileNotFoundError("Filepath specified is bad!")
    molecule = Molecule.from_legacy_zmat(zmat_str)
    _ = molecule.orient()
    click.echo(molecule.dump())


def moments_from_zmat(filepath: str):
    filepath = Path(filepath)
    if filepath.exists():
        zmat_str = filepath.read_text()
    else:
        raise FileNotFoundError("Filepath specified is bad!")
    molecule = Molecule.from_legacy_zmat(zmat_str)
    _ = molecule.orient()
    print(molecule.dump())