"""
script to make supercells in CONQUEST coordinates format
"""

from src.conquest import *
from dataclasses import dataclass
import copy


class SUPERCELL:
    def __init__(self, Nx: int, Ny: int, Nz: int, coords: CONQUEST_COORDINATES_PROCESSOR) -> None:
        """
        Args:
            Nx (int): Number of repeats along "a" lattice vector
            Ny (int): Number of repeats along "b" lattice vector
            Nz (int): Number of repeats along "c" lattice vector
            coords (CONQUEST_COORDINATES): CONQUEST_COORDINATES instance for the cell to be supercell'd
        """
        if Nx < 0 or Ny < 0 or Nz < 0:
            raise ValueError("One of, or multiple of, Nx Ny, Nz was not at least 0.")
        if (not isinstance(Nx, int)) or (not isinstance(Ny, int)) or (not isinstance(Nz, int)):
            raise TypeError("One of, or multiple of, Nx Ny, Nz was not an integer.")
        
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.coords = coords
        # Create new CONQUEST_COORDINATES
        self.supercell_coords: CONQUEST_COORDINATES = CONQUEST_COORDINATES(
            self.coords.CONQUEST_input
        )
        self.scale_lattice_vectors()
        self.new_num_atoms()
        self.create_supercell()
        self.supercell_coords.assign_atom_labels()
        self.supercell_coords.index_to_atom_map()

    def scale_lattice_vectors(self) -> None:
        """
        CONQUEST deals with orthorhombic cells only
        """
        self.supercell_coords.lattice_vectors = [
            [self.coords.lattice_vectors[0][0] * (self.Nx + 1), 0, 0],
            [0, self.coords.lattice_vectors[1][1] * (self.Ny + 1), 0],
            [0, 0, self.coords.lattice_vectors[2][2] * (self.Nz + 1)],
        ]

    def new_num_atoms(self) -> int:
        natoms = int(self.coords.natoms) * (self.Nx + self.Ny + self.Nz)
        self.supercell_coords.natoms = str(natoms)
        return natoms

    def create_supercell(self) -> None:
        """
        In terms of fractional coordinates, we set new coords of the original atoms to be
        x' = x/Nx, y' = y/Ny, z' = z/Nz

        We then take these new coordinates, and create new atoms of the same elements with coords (x', y', z') for every Nx, Ny, Nz
        i.e. we add 1/Nx, 1/Ny, 1/Nz in the appropriate direction.

        Example:
        Consider a bcc crystal defined by A: (0,0,0) and B: (1/2, 1/2, 1/2). Suppose we want to make a 3x2x2 supercell, We expect 24 atoms in the final cell.
        i.e. duplicate 3 in the "a" direction, and twice in the other 2 directions. We first rescale the lattice parameters, then we rescale the fractional coordinates

        (0,0,0) is trivial, but B gets rescaled to (1/2 * 1/3, 1/2 * 1/2, 1/2 * 1/2) = (1/6, 1/4, 1/4)

        However, the original (0,0,0) now has duplicates in x,y,z, namely the new atoms at (0,0,0) + {(1/3, 0, 0), (0, 1/2, 0), (0,0,1/2), ...}
        """
        # Rescale original atoms
        if (self.Nx == 0 and self.Ny == 0 and self.Nz == 0): return
        for atom in self.coords.Atoms:
            self.supercell_coords.Atoms.append(
                Atom(
                    atom.species,
                    [atom.coords[0] / (self.Nx + 1), atom.coords[1] / (self.Ny + 1), atom.coords[2] / (self.Nz + 1)],
                    can_move=atom.can_move,
                    label=atom.label,
                )
            )
        # We will be cloning and displacing the scaled original atoms to place the new atoms
        # so we keep a deep copy to loop through, and avoid infinitely duplicating atoms
        initial_scaled_atoms = copy.deepcopy(self.supercell_coords.Atoms)
        for atom in initial_scaled_atoms:
            for l in range(1,self.Nx,1):
                for m in range(1,self.Ny,1):
                    for n in range(1,self.Nz,1):
                        new_atom = Atom(
                                atom.species,
                                [
                                    atom.coords[0] + (l / self.Nx),
                                    atom.coords[1] + (m / self.Ny),
                                    atom.coords[2] + (n / self.Nz),
                                ],
                                can_move=atom.can_move,
                                label=atom.label,
                            )
                        self.supercell_coords.Atoms.append(
                            new_atom
                        )