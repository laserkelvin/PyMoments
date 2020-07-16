from warnings import warn

from pymoments import __version__
from pymoments import Molecule, Atom

import numpy as np


def test_version():
    assert __version__ == "0.1.0"


def test_isotopologues():
    zmat_str = """
    S
    C 1 1.594
    C 2 1.275 1 180.0
    C 3 1.328 2 180.0 1 0.0
    H 4 1.096 3 122.0 2 0.0
    H 4 1.096 3 122.0 2 180.0
    """
    molecule = Molecule.from_zmat_string(zmat_str)
    assert len(molecule) == 6
    _ = molecule.orient()
    isotopologues = molecule.generate_isotopologues(min_abundance=1e-4)
    # we expect at the default level oxygen 18, so two isotopologues
    assert len(isotopologues) == 96
    # raise Exception([iso.rot_con for iso in isotopologues])
    raise Exception("\n".join([iso.dump() for iso in isotopologues]))


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
    molecule = Molecule.from_zmat_string(zmat_str)
    test_com = np.array([3.05458538e-02, 4.20427610e-02, 2.57437663e-18])
    calc_com = molecule.compute_com(False)
    assert np.allclose(test_com, calc_com)


def test_legacy_zmat():
    zmat_str = """
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
    """
    molecule = Molecule.from_legacy_zmat(zmat_str)
    assert len(molecule) == 6
    # make sure the masses are being read correctly
    ref_mass = [31.972070, 12.000000, 12.000000, 12.000000, 1.007825, 1.007825]
    assert np.allclose(molecule.get_masses(), ref_mass)
    # make sure the cartesian coordinates are the same
    ref_coords = np.array(
        [
            [0.00000000e00, 0.00000000e00, 0.00000000e00],
            [1.59400000e00, 0.00000000e00, 0.00000000e00],
            [2.86900000e00, -1.56142467e-16, 9.56096861e-33],
            [4.19700000e00, -4.81408657e-16, 2.94777785e-32],
            [4.77779151e00, -9.29460713e-01, 5.69130544e-17],
            [4.77779151e00, 9.29460713e-01, -1.70739163e-16],
        ]
    )
    assert np.allclose(molecule.get_coords(), ref_coords)
    # make sure the COM is correct
    com = molecule.compute_com()
    ref_com = np.array([1.62087129e00, -1.25094696e-16, -1.63635732e-18])
    assert np.abs(com - ref_com).sum() <= 1e-2
    # test the whole shebang, and compare rotational constants
    com, rotcon, pmm = molecule.orient()
    reference = np.array([290228.46266522,   2496.61473982,   2475.32143207])
    assert np.allclose(reference, rotcon)
