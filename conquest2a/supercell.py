import copy
from collections.abc import Sequence
from conquest2a.conquest import conquest_coordinates, conquest_coordinates_processor, Atom


class supercell:
    def __init__(
        self, repeats_x: int, repeats_y: int, repeats_z: int, coords: conquest_coordinates_processor
    ) -> None:
        """
        Args:
            Nx (int): Number of repeats along "a" lattice vector
            Ny (int): Number of repeats along "b" lattice vector
            Nz (int): Number of repeats along "c" lattice vector
            coords (CONQUEST_COORDINATES): CONQUEST_COORDINATES instance for the cell to be supercell'd
        """
        if repeats_x < 0 or repeats_y < 0 or repeats_z < 0:
            raise ValueError("One of, or multiple of, Nx Ny, Nz was not at least 0.")
        if (
            (not isinstance(repeats_x, int))
            or (not isinstance(repeats_y, int))
            or (not isinstance(repeats_z, int))
        ):
            raise TypeError("One of, or multiple of, Nx Ny, Nz was not an integer.")

        self.repeats_x = repeats_x
        self.repeats_y = repeats_y
        self.repeats_z = repeats_z
        self.coords = coords
        # Create new CONQUEST_COORDINATES
        self.supercell_coords: conquest_coordinates = conquest_coordinates(
            self.coords.conquest_input
        )
        self.scale_lattice_vectors()
        self.create_supercell()
        self.supercell_coords.natoms = str(self.new_num_atoms())
        self.supercell_coords.assign_atom_labels()
        self.supercell_coords.index_to_atom_map()

    def scale_lattice_vectors(self) -> None:
        """
        CONQUEST deals with orthorhombic cells only
        """
        self.supercell_coords.lattice_vectors = [
            [self.coords.lattice_vectors[0][0] * (self.repeats_x + 1), 0, 0],
            [0, self.coords.lattice_vectors[1][1] * (self.repeats_y + 1), 0],
            [0, 0, self.coords.lattice_vectors[2][2] * (self.repeats_z + 1)],
        ]

    def new_num_atoms(self) -> int:
        natoms = len(self.supercell_coords.Atoms)
        self.supercell_coords.natoms = str(natoms)
        return natoms

    def create_atom(
        self,
        species: str,
        can_move: Sequence[str],
        atom_number: int,
        label: str,
        coord_0: int | float,
        coord_1: int | float,
        coord_2: int | float,
        disp_0: int | float = 0.0,
        disp_1: int | float = 0.0,
        disp_2: int | float = 0.0,
    ) -> Atom:
        return Atom(
            species,
            [
                coord_0 + disp_0,
                coord_1 + disp_1,
                coord_2 + disp_2,
            ],
            can_move=can_move,
            label=label,
            number=atom_number
        )

    def range(self, upper_bound: int) -> list[int] | range:
        if upper_bound == 0:
            return [0]
        return range(0, upper_bound + 1, 1)

    def create_supercell(self) -> None:
        """
        In terms of fractional coordinates, we set new coords of the original atoms to be
        x' = x/Nx, y' = y/Ny, z' = z/Nz

        We then take these new coordinates, and create new atoms of the same elements,
        with coords (x', y', z') for every Nx, Ny, Nz
        i.e. we add 1/Nx, 1/Ny, 1/Nz in the appropriate direction.

        Example:
        Consider a bcc crystal defined by A: (0,0,0) and B: (1/2, 1/2, 1/2). 
        Suppose we want to make a 3x2x2 supercell: expect 24 atoms in the final cell.
        i.e. duplicate 3 in the "a" direction, twice in the other 2 directions.
        We first rescale the lattice parameters, then we rescale the fractional coordinates

        (0,0,0) is trivial, but B gets rescaled to (1/2 * 1/3, 1/2 * 1/2, 1/2 * 1/2) = (1/6, 1/4, 1/4)

        However, the original (0,0,0) now has duplicates in x,y,z, 
        namely the new atoms at (0,0,0) + {(1/3, 0, 0), (0, 1/2, 0), (0,0,1/2), ...}
        """
        # No repeats at all -> just return original crystal
        if self.repeats_x == 0 and self.repeats_y == 0 and self.repeats_z == 0:
            self.supercell_coords.Atoms = copy.deepcopy(self.coords.Atoms)
            self.supercell_coords.natoms = self.coords.natoms
            return
        for atom in self.coords.Atoms:
            for l in self.range(self.repeats_x):
                for m in self.range(self.repeats_y):
                    for n in self.range(self.repeats_z):
                        new_atom = Atom(
                            atom.species,
                            [
                                (atom.coords[0] + l) / (self.repeats_x + 1),
                                (atom.coords[1] + m) / (self.repeats_y + 1),
                                (atom.coords[2] + n) / (self.repeats_z + 1),
                            ],
                            can_move=atom.can_move,
                            label=atom.label,
                            number=len(self.supercell_coords.Atoms) + 1
                        )
                        self.supercell_coords.Atoms.append(new_atom)
