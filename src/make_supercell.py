"""
script to make supercells in CONQUEST coordinates format
"""

from conquest import *
from dataclasses import dataclass
import copy


class supercell:
    def __init__(self, Nx: int, Ny: int, Nz: int, coords: CONQUEST_COORDINATES_PROCESSOR) -> None:
        """_summary_

        Args:
            Nx (int): Number of repeats along "a" lattice vector
            Ny (int): Number of repeats along "b" lattice vector
            Nz (int): Number of repeats along "c" lattice vector
            coords (CONQUEST_COORDINATES): CONQUEST_COORDINATES instance for the cell to be super'd
        """
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

    def scale_lattice_vectors(self) -> None:
        self.supercell_coords.lattice_vectors = [
            [self.coords.lattice_vectors[0][0] * self.Nx, 0, 0],
            [0, self.coords.lattice_vectors[1][1] * self.Ny, 0],
            [0, 0, self.coords.lattice_vectors[2][2] * self.Nz],
        ]

    def new_num_atoms(self) -> None:
        self.supercell_coords.natoms = str(int(self.coords.natoms) * self.Nx * self.Ny * self.Nz)

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
        for atom in self.coords.Atoms:
            self.supercell_coords.Atoms.append(
                Atom(
                    atom.species,
                    [atom.coords[0] / self.Nx, atom.coords[1] / self.Ny, atom.coords[2] / self.Nz],
                    can_move=atom.can_move,
                    label=atom.label,
                )
            )
        for atom in self.supercell_coords.Atoms:
            for l in range(self.Nx):
                for m in range(self.Ny):
                    for n in range(self.Nz):
                        self.supercell_coords.Atoms.append(
                            Atom(
                                atom.species,
                                [
                                    atom.coords[0] + (l / self.Nx),
                                    atom.coords[1] + (m / self.Ny),
                                    atom.coords[2] + (n / self.Nz),
                                ],
                                can_move=atom.can_move,
                                label=atom.label,
                            )
                        )
