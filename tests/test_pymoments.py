
from warnings import warn

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
    ref_coords = np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
       [ 1.59400000e+00,  0.00000000e+00,  0.00000000e+00],
       [ 2.86900000e+00, -1.56142467e-16,  9.56096861e-33],
       [ 4.19700000e+00, -4.81408657e-16,  2.94777785e-32],
       [ 4.77779151e+00, -9.29460713e-01,  5.69130544e-17],
       [ 4.77779151e+00,  9.29460713e-01, -1.70739163e-16]])
    assert np.allclose(molecule.get_coords(), ref_coords)
    # make sure the COM is correct
    com = molecule.compute_com()
    ref_com = np.array([1.62087129e+00, -1.25094696e-16, -1.63635732e-18])
    assert np.abs(com - ref_com).sum() <= 1e-2
    com, rotcon, pmm = molecule.orient()
    reference = np.array([41792.97459405,  2349.96288366,  2224.86188517])
    assert np.allclose(reference, rotcon)
    # test the dump
    print(molecule.dump())