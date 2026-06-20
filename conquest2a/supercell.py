import copy
import numpy as np
from conquest2a.conquest import conquest_coordinates, conquest_coordinates_processor, Atom


class supercell:
    """Class to produce supercells. CONQUEST deals with orthorhombic cells only. The resulting supercell will therefore be orthorhombic.

    :param repeats_x: Number of repeats along the :math:`a` lattice vector
    :type repeats_x: ``int``
    :param repeats_y: Number of repeats along the :math:`b` lattice vector
    :type repeats_y: ``int``
    :param repeats_z: Number of repeats along the :math:`c` lattice vector
    :type repeats_z: ``int``
    :param coords_proc: The :class:`conquest_coordinates_processor` to use.
    :type coords_proc: conquest_coordinates_processor
    :raises ValueError: If any of the ``repeats_*`` is not a positive integer
    """

    def __init__(
        self,
        repeats_x: int,
        repeats_y: int,
        repeats_z: int,
        coords_proc: conquest_coordinates_processor,
    ) -> None:

        if repeats_x < 0 or repeats_y < 0 or repeats_z < 0:
            raise ValueError(
                "One of, or multiple of, repeats_x repeats_y, repeats_z was not at least 0."
            )

        self.repeats_x = repeats_x
        self.repeats_y = repeats_y
        self.repeats_z = repeats_z
        self.coords_proc = coords_proc
        # Create new CONQUEST_COORDINATES
        supercell_coords_instance = conquest_coordinates(self.coords_proc.coords.conquest_input)
        self.supercell_coords: conquest_coordinates = supercell_coords_instance
        self.scale_lattice_vectors()
        self.create_supercell()
        self.supercell_coords.natoms = str(self.new_num_atoms())
        self.supercell_coords.assign_atom_labels()
        self.supercell_coords.index_to_atom_map()

    def scale_lattice_vectors(self) -> None:
        """
        Produce the supercell's lattice vectors based on the repeats supplied.
        """
        repeat_matrix = np.array(
            [
                [(self.repeats_x + 1), 0, 0],
                [0, (self.repeats_y + 1), 0],
                [0, 0, (self.repeats_z + 1)],
            ]
        )
        self.supercell_coords.lattice_vectors = np.matmul(
            self.coords_proc.coords.lattice_vectors, repeat_matrix
        )

    def new_num_atoms(self) -> int:
        """Gets the number of atoms in the new supercell. Needed for the coordinates file.

        :return: The number of atoms in the new supercell
        :rtype: ``int``
        """
        natoms = len(self.supercell_coords.atoms)
        self.supercell_coords.natoms = str(natoms)
        return natoms

    def range(self, upper_bound: int) -> list[int] | range:
        if upper_bound == 0:
            return [0]
        return range(0, upper_bound + 1, 1)

    def create_supercell(self) -> None:
        """This method creates the :class:`Atom` s of the new supercell, looping through the original atoms, creating new ones corresponding to the repeats in each direction.

        In terms of fractional coordinates, we set new coords of the original atoms to be
        :math:`x' = x/N_x, y' = y/N_y, z' = z/N_z`

        We then take these new coordinates, and create new atoms of the same elements,
        with coords :math:`(x', y', z')` for every ``repeats_x,y,z`` respectively.
        i.e. we add :math:`1/N_x, 1/N_y, 1/N_z` in the appropriate direction.

        Example:
        Consider a bcc crystal defined by A: (0,0,0) and B: (1/2, 1/2, 1/2).
        Suppose we want to make a 3x2x2 supercell: expect 24 atoms in the final cell.
        i.e. duplicate 3 in the "a" direction, twice in the other 2 directions.
        We first rescale the lattice parameters, then we rescale the fractional coordinates

        (0,0,0) is trivial, but B gets rescaled to (1/2 * 1/3, 1/2 * 1/2, 1/2 * 1/2)
        = (1/6, 1/4, 1/4)

        However, the original (0,0,0) now has duplicates in x,y,z,
        namely the new atoms at (0,0,0) + {(1/3, 0, 0), (0, 1/2, 0), (0,0,1/2), ...}
        """
        # No repeats at all -> just return original crystal
        if self.repeats_x == 0 and self.repeats_y == 0 and self.repeats_z == 0:
            self.supercell_coords.atoms = copy.deepcopy(self.coords_proc.coords.atoms)
            self.supercell_coords.natoms = self.coords_proc.coords.natoms
            return
        for atom in self.coords_proc.coords.atoms:
            for l in self.range(self.repeats_x):
                for m in self.range(self.repeats_y):
                    for n in self.range(self.repeats_z):
                        new_atom = Atom(
                            atom.species,
                            np.array(
                                [
                                    (atom.coords[0] + l) / (self.repeats_x + 1),
                                    (atom.coords[1] + m) / (self.repeats_y + 1),
                                    (atom.coords[2] + n) / (self.repeats_z + 1),
                                ]
                            ),
                            can_move=atom.can_move,
                            label=atom.label,
                            number=len(self.supercell_coords.atoms) + 1,
                            spins=atom.spins,
                            forces=atom.forces,
                        )
                        self.supercell_coords.atoms.append(new_atom)
        self.supercell_coords.get_cartesian_positions()
