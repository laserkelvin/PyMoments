
from typing import List
from collections import namedtuple
from copy import deepcopy

import numpy as np
from mendeleev import element
from scipy import constants

from pymoments.main import compute_xyz, compute_angle, rotate_coordinates

# These namedtuples are convenient for having some structure in
# how the different types of connectivity can be tracked
Bond = namedtuple("Bond", "i value")
Angle = namedtuple("Angle", "j value")
Dihedral = namedtuple("Dihedral", "k value")


class Atom:    
    def __init__(self, index: int, symbol: str, mass: float, **kwargs):
        self.index = index
        self.symbol = symbol
        self.mass = mass
        self.bond = None
        self.angle = None
        self.dihedral = None
        self.xyz = np.zeros(3)
        self.__dict__.update(kwargs)
    
    def __repr__(self):
        coords = [f"{value:.6f}" for value in self.xyz.tolist()]
        coords = " ".join(coords)
        return f"{self.symbol} {coords}"
    

class Molecule:
    """
    A high-level implementation of a `Molecule`, where the primary
    function is to store a collection of `Atom` objects as a list
    attribute.
    
    Additional functionality include convenience functions for computing
    rotational constants for a set of isotopologues, as well as determining
    scaling factors.
    """
    def __init__(self, atoms=None):
        self.atoms = atoms
        self.com = False
        self.inertial = False
        self.scaling = 1.
        self.rot_con = np.zeros(3)
        self.pmm = np.zeros((3, 3))
    
    def __repr__(self):
        return "\n".join([str(atom) for atom in self.atoms])
    
    def __truediv__(self, other):
        if type(other) == type(self):
            self.scaling = self.rot_con / other.rot_con
            return None
        elif type(other) == float or type(other) == np.ndarray:
            new_molecule = deepcopy(self)
            new_molecule.scaling = other
            new_molecule.rot_con / other
            return new_molecule
        else:
            raise NotImplementedError("Illegal division; can only divide with `Molecule`, float, or NumPy array objects!")
    
    def __mul__(self, other):
        if type(other) == type(self):
            new_molecule = deepcopy(self)
            new_molecule.rot_con *= other.scaling
            return new_molecule
        elif type(other) == float or type(other) == np.ndarray:
            new_molecule = deepcopy(self)
            new_molecule.rot_con *= other
            new_molecule.scaling = other
            return new_molecule
        else:
            raise NotImplementedError("Illegal mutiplication; can only multiply with `Molecule`, float, or NumPy array objects!")
    
    @classmethod
    def from_zmat_string(cls, zmat_str: str):
        """
        Parse a Gaussian/CFOUR format of Z-matrix, creating `Atom` objects
        for each line entry which provides a nice transparent, abstract method
        for tracking parameters of each row. The `Atom` objects are collected
        up at the end, and stored in an instance of the `Molecule` class.
        
        An example of the input string looks like this:
        
        zmat_str = "
        O
        H 1 0.79
        H 1 0.79 2 108.00
        "
        
        Each atom symbol is followed by the atom index (connectivity)
        and the corresponding parameter (bond, angle, dihedral).
        
        Rare isotopes can also be specified; for deuterium you can replace
        H with D, and for others you can use square brackets:
        zmat_str = "
        O[18]
        H 1 0.79
        D 1 0.79 2 108.0
        "

        Parameters
        ----------
        zmat : str
            String containing 

        Returns
        -------
        molecule
            Instance of a `Molecule` object
        """
        zmat = zmat.strip().split("\n")
        natoms = len(zmat)
        xyz = np.zeros((natoms, 3), dtype=float)
        atoms = list()
        for index, line in enumerate(zmat):
            split_line = line.split()
            symbol = split_line[0]
            mass = element(symbol).mass
            parameters = {key: None for key in ["bond", "angle", "dihedral"]}
            # first atom has nothing
            if index != 0:
                # some type conversions
                connectivity = [int(value) for value in split_line[1::2]]
                params = [float(value) for value in split_line[2::2]]
                # loop over the three types of connectivity for each atom. This loop
                # will terminate early if it runs out of connectivity/params
                for c, p, con_type in zip(connectivity, params, [Bond, Angle, Dihedral]):
                    name = con_type.__name__.lower()
                    # convert angles to radians
                    if name != "bond":
                        p = np.deg2rad(p)
                    parameters[name] = con_type(c - 1, p)
            # generate an Atom object with all the data
            atom = Atom(index, symbol, mass, **parameters)
            atoms.append(atom)
        
        # loop again, this time calculating the coordinates
        for atom in atoms:
            atom.xyz = compute_xyz(atom, atoms)
        # create a Molecule object
        molecule = cls(atoms)
        return molecule
    
    @classmethod
    def from_xyz(cls, xyz_str: str):
        xyz_str = xyz_str.strip().split("\n")
        natoms = len(xyz_str)
        atoms = list()
        # loop over each line of the XYZ file and parse out
        # the atomic symbol and coordinates
        for index, line in enumerate(xyz_str):
            split_line = line.split()
            symbol = split_line[0]
            coords = np.array([float(value) for value in line[1:]])
            mass = element(symbol).mass
            atoms.append(Atom(index, symbol, mass, **{"xyz": coords}))
        molecule = cls(atoms)
        return molecule
    
    @classmethod
    def from_legacy(cls, zmat_str: str):
        """
        The legacy input looks like this:
        "
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
        "

        Parameters
        ----------
        zmat_str : str
            [description]
        """
        zmat_str = zmat_str.strip().split("\n")
        # skip two comment lines
        zmat_str = zmat_str[2:]
        desc = zmat_str.pop(0).split()
        natoms = int(desc[0])
    
    def get_coords(self):
        return np.vstack([atom.xyz for atom in self.atoms])
    
    def get_masses(self):
        return np.array([atom.mass for atom in self.atoms])
    
    def get_symbols(self):
        return "".join([atom.symbol for atom in self.atoms])
    
    def compute_com(self, shift=False) -> np.ndarray:
        """
        Compute the center of mass coordinates for the `Molecule`.
        This function will more or less return zeros if the molecule
        has already been shifted to the center of mass.

        Parameters
        ----------
        shift : bool, optional
            Toggle whether to automatically move to a center of
            mass representation, by default False

        Returns
        -------
        np.ndarray
            NumPy 1D array containing center of mass xyz
        """
        coords = self.get_coords()
        masses = self.get_masses()
        # vectorized computation of the center of mass
        com = np.sum(masses[:,None] * coords, axis=0) / masses.sum()
        if shift and not self.com:
            self.com = True
            for atom in self.atoms:
                atom.xyz -= com
        return com
    
    def compute_inertia_tensor(self, shift=False):
        """
        Calculate the moments of inertia tensor, and diagonalize it to
        obtain the principal moments of inertia and axes. For transparency
        sake, we convert the mass and coordinates into SI units of meters
        and kg, and perform the final conversion of rotational constants
        into MHz using `scipy.constants` so that people can track the unit
        conversions appropriately.
        
        Keep in mind that the rotational constants returned are actually
        sorted: these may not correspond to the actual orientation of the
        principal axes, and you will have to make that judgement call yourself.

        Parameters
        ----------
        shift : bool, optional
            Toggles whether the `Molecule` is rotated into the principal
            axis orientation, by default False

        Returns
        -------
        np.ndarray, np.ndarray
            Rotational constants in MHz, and eigenvectors of the
            principal axis system.
        """
        coords = self.get_coords()
        masses = self.get_masses()
        # unit conversions; everything is better in SI
        coords *= 1e-9   # to meters
        masses *= constants.atomic_mass
        inertia_tensor = np.zeros((3, 3))
        # hard coded inertia matrix elements
        inertia_tensor[0,0] = np.sum((coords[:,1]**2. + coords[:,2]**2.) * masses[:,None])
        inertia_tensor[1,1] = np.sum((coords[:,0]**2. + coords[:,2]**2.) * masses[:,None])
        inertia_tensor[2,2] = np.sum((coords[:,0]**2. + coords[:,1]**2.) * masses[:,None])
        # off-diagonal elements
        inertia_tensor[0,1] = -np.sum(coords[:,0] * coords[:,1] * masses[:,None])
        inertia_tensor[0,2] = -np.sum(coords[:,0] * coords[:,2] * masses[:,None])
        inertia_tensor[1,2] = -np.sum(coords[:,1] * coords[:,2] * masses[:,None])
        inertia_tensor = np.maximum(inertia_tensor, inertia_tensor.T)
        # compute principal moments of inertia
        pmi, pmm = np.linalg.eig(inertia_tensor)
        # convert PMI to MHz
        rot_con = constants.h / (8 * (np.pi)**2 * (constants.c * 100.) * pmi) * constants.c / 10.
        # if we request for a shift, and we haven't already done so
        # we can rotate the atomic coordinates to the principal axis orientation
        if shift and not self.inertial:
            self.inertial = True
            for atom in self.atoms:
                atom.xyz = rotate_coordinates(atom.xyz, pmm)
        ordering = np.argsort(rot_con)[::-1]
        return rot_con[ordering], pmm[ordering]
    
    def final_frame(self):
        com = self.compute_com(True)
        (rotational_constants, inertial_vector) = self.compute_inertia_tensor(True)
        return com, rotational_constants, inertial_vector