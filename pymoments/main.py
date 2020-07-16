
from dataclasses import dataclass
from collections import namedtuple
from typing import List

import numpy as np
from scipy import constants
from scipy.spatial.transform import Rotation as R


def compute_xyz(atom: "Atom", atom_list: List[str]):
    """
    This is an impure function that modifies the xyz coordinates of
    a given atom in place. The flow of this routine is determined by
    the `index` attribute of an `Atom`, which in turn dictates what 
    kind of connectivity is possible (i.e. `Bond`, `Angle`, `Dihedral`).

    Parameters
    ----------
    atom : Atom object
        Reference to an `Atom` object.
    atom_list : List[Atom]
        List holding all the `Atom` objects in the Z-matrix.
        
    Returns
    -------
    coords
        NumPy 1D array containing coordinates of the atom
    """
    index = atom.index
    # if it's the first atom, set as origin
    if index == 0:
        coords = np.zeros(3)
    # if it's the second atom, add as distance
    elif index == 1:
        coords = np.zeros(3)
        coords[0] += atom.bond.value
    # start dealing in angles
    elif index > 1:
        # get indices of the atoms
        i, j = atom.bond.i, atom.angle.j
        # retrieve coordinates for the two atoms defined in the
        # connectivity
        a_coords, b_coords = atom_list[i].xyz, atom_list[j].xyz
        dist, ang = atom.bond.value, atom.angle.value
        # for the third atom, constrain in-plane with the first two
        if index == 2:
            tor_angle = np.deg2rad(90.)
            c_coords = np.array([0., 0., 1.])
        # for all other atoms, go nuts
        else:
            tor_angle, tor_index = atom.dihedral.value, atom.dihedral.k
            c_coords = atom_list[tor_index].xyz
        
        # do a crap load of trigonometry
        v1 = a_coords - b_coords
        v2 = a_coords - c_coords

        n = np.cross(v1, v2)
        nn = np.cross(v1, n)

        n /= np.linalg.norm(n)
        nn /= np.linalg.norm(nn)

        n *= -np.sin(tor_angle)
        nn *= np.cos(tor_angle)

        v3 = n + nn
        v3 /= np.linalg.norm(v3)
        v3 *= dist * np.sin(ang)

        v1 /= np.linalg.norm(v1)
        v1 *= dist * np.cos(ang)
        coords = a_coords + v3 - v1
    return coords


def compute_angle(u: np.ndarray, v: np.ndarray) -> float:
    """
    Compute the angle between two vectors. The vectors do not 
    have to be normalized.

    Parameters
    ----------
    u, v : np.ndarray
        NumPy 1D arrays of the same dimensionality

    Returns
    -------
    float
        Angle between vectors `u` and `v`.
    """
    angle = np.arccos(np.dot(u, v) / np.dot(np.linalg.norm(u), np.linalg.norm(v)))
    return angle


def rotate_coordinates(coords: np.ndarray, axis_coords: np.ndarray) -> np.ndarray:
    """
    Given a set of coordinates, `coords`, and the eigenvectors of the principal
    moments of inertia tensor, use the scipy `Rotation` class to rotate the
    coordinates into the principal axis frame.

    Parameters
    ----------
    coords : np.ndarray
        NumPy 1D array containing xyz coordinates
    axis_coords : np.ndarray
        NumPy 2D array (shape 3x3) containing the principal axis
        vectors

    Returns
    -------
    np.ndarray
        NumPy 1D array containing the rotated coordinates.
    """
    # Create a Rotation object from the eigenvectors
    r_mat = R.from_matrix(axis_coords)
    # transform the coordinates into the principal axis
    return r_mat.apply(coords)

