from typing import List
from collections import namedtuple
from copy import deepcopy
from warnings import warn
from itertools import product

import numpy as np
from mendeleev import element
from scipy import constants

from pymoments.main import (
    compute_xyz,
    compute_angle,
    rotate_coordinates,
    kappa,
    inertial_defect,
)

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
        self.scaling = 1.0
        self.rot_con = np.zeros(3)
        self.pmm = np.zeros((3, 3))

    def __repr__(self):
        return "\n".join([str(atom) for atom in self.atoms])

    def __len__(self):
        return len(self.atoms)

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
            raise NotImplementedError(
                "Illegal division; can only divide with `Molecule`, float, or NumPy array objects!"
            )

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
            raise NotImplementedError(
                "Illegal mutiplication; can only multiply with `Molecule`, float, or NumPy array objects!"
            )

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
        zmat_str = zmat_str.strip().split("\n")
        natoms = len(zmat_str)
        xyz = np.zeros((natoms, 3), dtype=float)
        atoms = list()
        for index, line in enumerate(zmat_str):
            split_line = line.split()
            symbol = split_line[0]
            # this is a quick one-liner using list-comprehensions to get the most abundant mass
            isotopes = [
                isotope for isotope in element(symbol).isotopes if isotope.abundance
            ]
            mass = max(isotopes, key=lambda x: x.abundance).mass
            parameters = {key: None for key in ["bond", "angle", "dihedral"]}
            # first atom has nothing
            if index != 0:
                # some type conversions
                connectivity = [int(value) for value in split_line[1::2]]
                params = [float(value) for value in split_line[2::2]]
                # loop over the three types of connectivity for each atom. This loop
                # will terminate early if it runs out of connectivity/params
                for c, p, con_type in zip(
                    connectivity, params, [Bond, Angle, Dihedral]
                ):
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
        """
        Create a `Molecule` object from an XYZ string. This does not follow the
        standard .xyz file format, where the first two lines are the number of
        atoms and a comment line respectively; rather, this is the format where
        only the Atom X Y Z per line is required.
        
        For example:
        xyz_str = "
        O 0.030541 0.042037 -0.000000
        H -0.759459 0.042037 -0.000000
        H 0.274665 -0.709298 0.000000
        "

        Parameters
        ----------
        xyz_str : str
            String containing atom and XYZ specification, with each line
            corresponding to one atom.

        Returns
        -------
        `molecule`
            Instance of a `Molecule` object
        """
        xyz_str = xyz_str.strip().split("\n")
        natoms = len(xyz_str)
        atoms = list()
        # loop over each line of the XYZ file and parse out
        # the atomic symbol and coordinates
        for index, line in enumerate(xyz_str):
            split_line = line.split()
            symbol = split_line[0]
            coords = np.array([float(value) for value in split_line[1:]])
            # get the most abundant isotope
            isotopes = [
                isotope for isotope in element(symbol).isotopes if isotope.abundance
            ]
            mass = max(isotopes, key=lambda x: x.abundance).mass
            atoms.append(Atom(index, symbol, mass, **{"xyz": coords}))
        molecule = cls(atoms)
        return molecule

    @classmethod
    def from_legacy_zmat(cls, zmat_str: str):
        """
        Create a `Molecule` object using "legacy" input. This file format is
        not recommended as it is unable to take advantage of some of the newer
        functionality, but is supported just for backwards compatibility.
        
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
        zmat : str
            String containing 

        Returns
        -------
        molecule
            Instance of a `Molecule` object
        """
        zmat_str = zmat_str.strip().split("\n")
        # skip two comment lines
        zmat_str = zmat_str[2:]
        desc = zmat_str.pop(0).split()
        natoms = int(desc[0])
        atoms = list()
        for index, line in enumerate(zmat_str):
            if index == natoms:
                break
            else:
                parameters = {key: None for key in ["bond", "angle", "dihedral"]}
                split_line = line.split()
                mass = float(split_line[-1])
                # No symbols are defined for legacy ZMAT specification, and we
                # can't really just infer from mass
                symbol = "X"
                if index == 0:
                    pass
                else:
                    # Read in the connectivity and parameters
                    for offset, con_type in zip(range(1, 4), [Bond, Angle, Dihedral]):
                        name = con_type.__name__.lower()
                        # get the connectivity
                        connection = int(split_line[offset])
                        if connection == 0:
                            pass
                        else:
                            value = float(split_line[offset + 3])
                            if offset != 1:
                                # convert angles to radians
                                value = np.deg2rad(value)
                        parameters[name] = con_type(connection - 1, value)
                atom = Atom(index, symbol, mass, **parameters)
                atoms.append(atom)
        for atom in atoms:
            atom.xyz = compute_xyz(atom, atoms)
        molecule = cls(atoms)
        return molecule

    def get_coords(self):
        return np.vstack([atom.xyz for atom in self.atoms])

    def get_masses(self):
        return np.array([atom.mass for atom in self.atoms])

    def get_symbols(self):
        return "".join([atom.symbol for atom in self.atoms])

    def modify_atom_masses(self, masses: np.ndarray, copy=False):
        """
        Modify the atom masses of this molecule. This function can
        operate in two ways: if `copy=True`, then a new `Molecule`
        object is returned with the new masses. Otherwise, the masses
        are modified in-place.

        Parameters
        ----------
        masses : np.ndarray
            [description]
        """
        assert len(masses) == len(self)
        print(f"Old masses: {self.get_masses()}")
        if copy:
            new_molecule = deepcopy(self)
            new_molecule.modify_atom_masses(masses, copy=False)
            return new_molecule
        else:
            for atom, new_mass in zip(self.atoms, masses):
                atom.mass = new_mass
            return None

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
        com = np.sum(masses[:, None] * coords, axis=0) / masses.sum()
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
        coords *= 1e-9  # to meters
        masses *= constants.atomic_mass
        inertia_tensor = np.zeros((3, 3))
        # hard coded inertia matrix elements
        inertia_tensor[0, 0] = np.sum(
            (coords[:, 1] ** 2.0 + coords[:, 2] ** 2.0) * masses[:, None]
        )
        inertia_tensor[1, 1] = np.sum(
            (coords[:, 0] ** 2.0 + coords[:, 2] ** 2.0) * masses[:, None]
        )
        inertia_tensor[2, 2] = np.sum(
            (coords[:, 0] ** 2.0 + coords[:, 1] ** 2.0) * masses[:, None]
        )
        # off-diagonal elements
        inertia_tensor[0, 1] = -np.sum(coords[:, 0] * coords[:, 1] * masses[:, None])
        inertia_tensor[0, 2] = -np.sum(coords[:, 0] * coords[:, 2] * masses[:, None])
        inertia_tensor[1, 2] = -np.sum(coords[:, 1] * coords[:, 2] * masses[:, None])
        inertia_tensor = np.maximum(inertia_tensor, inertia_tensor.T)
        # compute principal moments of inertia
        pmi, pmm = np.linalg.eig(inertia_tensor)
        # convert PMI to MHz
        rot_con = (
            constants.h
            / (8 * (np.pi) ** 2 * (constants.c * 100.0) * pmi)
            * constants.c
            / 10.0
        )
        # if we request for a shift, and we haven't already done so
        # we can rotate the atomic coordinates to the principal axis orientation
        if shift and not self.inertial:
            self.inertial = True
            for atom in self.atoms:
                atom.xyz = rotate_coordinates(atom.xyz, pmm)
        # This sorts the rotational constants in order of A > B > C, and similarly
        # the principal axes vectors too (row order)
        ordering = np.argsort(rot_con)[::-1]
        return rot_con[ordering], pmm[ordering]

    def orient(self):
        """
        Shifts the molecular cartesian coordinates into the center of mass and
        principal axis frame sequentially. We compute the COM corrections first,
        followed by the inertial corrections.

        Returns
        -------
        np.ndarray, np.ndarray, np.ndarray
            Returns the COM, Rotational constants in MHz, and inertial
            axes vectors.
        """
        com = self.compute_com(True)
        (rotational_constants, inertial_vector) = self.compute_inertia_tensor(True)
        self.rot_con, self.pmm = rotational_constants, inertial_vector
        return com, rotational_constants, inertial_vector

    def compute_kappa(self):
        if not self.com or not self.inertial:
            warn("Not in center of mass or principal axis frame; not meaningful!")
        return kappa(*self.rot_con)

    def compute_inertial_defect(self):
        if not self.com or not self.inertial:
            warn("Not in center of mass or principal axis frame; not meaningful!")
        return inertial_defect(*self.rot_con)

    def dump(self):
        template = """===================== Primary input
Formula: {symbols}
Masses (AMU): {mass}
===================== Parameters
Rotational constants (MHz): {rot_con}
Inertial axis vectors:
{inertial_vector}
===================== Derived values
Asymmetry parameter: {kappa:.4f}
Inertial defect (amu A**2): {defect:.4f}
        """
        parameter_dict = {
            "symbols": self.get_symbols(),
            "mass": self.get_masses(),
            "rot_con": self.rot_con,
            "inertial_vector": self.pmm,
            "kappa": self.compute_kappa(),
            "defect": self.compute_inertial_defect(),
        }
        return template.format_map(parameter_dict)

    def generate_isotopologues(self, min_abundance=0.001):
        masses = list()
        for atom in self.atoms:
            isotopes = [isotope.mass for isotope in element(atom).isotopes if isotope.abundance]
            isotopes = filter(isotopes, lambda x: x.abundance >= min_abundance)
            masses.append([isotope.mass for isotope in isotopes])
        isotopologues = list()
        # iterate through every combination
        for iso_masses in product(*masses):
            iso =self.modify_atom_masses(iso_masses, copy=True)
            _ = iso.orient()
            isotopologues.append(iso)
        return isotopologues