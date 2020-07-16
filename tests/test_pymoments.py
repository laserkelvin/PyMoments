from pymoments import __version__
from pymoments import Molecule, Atom

import numpy as np


def test_version():
    assert __version__ == '0.1.0'


def test_com():
    """
    Test to make sure the center of mass calculation is working
    correctly; this also checks that the Z-matrix conversion is
    working as well.
    """
    zmat_str = """
    O
    H 1 0.79
    H 1 0.79 2 108.0
    """
    Molecule.from_zmat_string(zmat_str)
    test_com = np.array([2.53356346e-17, 1.58347716e-18, 0.00000000e+00])
    calc_com = Molecule.compute_com(shift=False)
    assert np.allclose(test_com, calc_com)
    