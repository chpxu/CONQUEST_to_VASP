from dataclasses import dataclass, field
from typing import Callable
import os
import re
import importlib.resources
from pathlib import Path
from collections.abc import Sequence
import numpy as np
import ase
from conquest2a.constants import BOHR_TO_ANGSTROM_VOLUME
import conquest2a._types as c2at


@dataclass
class Atom:
    """Class which holds Atom data.

    :param species: The integer referring to the species as defined in the ``Conquest_input`` file.
    :type species: ``int``
    :param coords: Fractional coordinates of the atom in the interval :math:`[0,1)`. Numbers outside this range are wrapped back into the range.
    :type coords: :ref:`REAL ARRAY <types>`
    :param can_move: Whether the atom is allowed to move in each Cartesian axis.
    :type can_move: ``Sequence[str]`` - list of 'T' or 'F'.
    :param number: The atom number, i.e. its position in the coordinates file.
    :type number: ``int``
    :param label: The atom element
    :type label: ``str``
    :param cart_coords: Cartesian coordinates of the atom.
    :type cart_coords: :ref:`REAL ARRAY <types>`
    :param forces: Force vector of the atom, defaults to ``np.array([0.0, 0.0, 0.0])``.
    :type forces: :ref:`REAL ARRAY <types>`
    :param spins: Spin moment on the atom,  defaults to ``np.array([0.0, 0.0, 0.0])``.
    :type spins: :ref:`REAL ARRAY <types>`
    """
    species: int
    coords: c2at.REAL_ARRAY
    can_move: Sequence[str]
    number: int
    label: str = ""
    cart_coords: c2at.REAL_ARRAY = field(init=False)
    forces: c2at.REAL_ARRAY = field(default_factory=lambda: np.array([0.0, 0.0, 0.0]))
    spins: c2at.REAL_ARRAY = field(default_factory=lambda: np.array([0.0, 0.0, 0.0]))

    def __str__(self) -> str:
        def fmt_array(arr: c2at.REAL_ARRAY) -> str:
            out_str: Callable[[c2at.REAL_NUMBER], str] = lambda x: f"{x:.5f}"
            return f"({', '.join(out_str(x) for x in arr)})"

        return (
            f"Atom {self.number} ({self.label})\n"
            f"  Species Index : {self.species}\n"
            f"  Frac coords   : {fmt_array(self.coords)}\n"
            f"  Cart coords   : {fmt_array(self.cart_coords)}\n"
            f"  Can move      : {', '.join(self.can_move)}\n"
            f"  Force         : {fmt_array(self.forces)}\n"
            f"  Spin          : {fmt_array(self.spins)}\n"
        )
    def to_ase(self) -> ase.Atom:
        """Method to return an equivalent ASE ```Atom`` <https://ase-lib.org/ase/atom.html>`__. 
        The ``spins`` field is set to the ``Atom.magmom`` whilst the ``forces`` field is set to ``Atom.momentum``.

        
        :returns: ase.Atom: Equivalent ASE ``Atom`` object.
        :rtype: ``ase.Atom``
        """
        return ase.Atom(symbol=self.label, position=self.coords, magmom=self.spins, momentum=self.forces)


class conquest_input:
    def __init__(self, species_dict: dict[int, str]) -> None:
        """`conquest_input` serves as the main entrypoint describing the species involved in the simulation.

        :param species_dict: hashmap species_index <-> element_label
                It is not completely reliable to directly read conquest_input
                E.g., for multiple spins, must duplicate an element and call it a
                new Conquest species, however labels can be any alphanumeric string
                Therefore we expect a dictionary to be passed in independently.
        :type species_dict: dict[int, str]
        :raises ValueError: If incorrect elements are passed in, the code will abort.
        """
        self.element_file = "elements.txt"
        self.species_dict = species_dict
        # elements_from_file = LIBRARY.joinpath(self.element_file)
        elements_from_file = importlib.resources.read_text("conquest2a", self.element_file)
        self.elements = [e.strip() for e in elements_from_file.split(",")]
        self.allowed_element_labels: list[str] = self.elements
        # elements_from_file.close()
        if not self.dict_contains_only_real_elements():
            raise ValueError("Provided species map contains fake chemical elements.")
        self.unique_elements = list(set(self.species_dict.values()))

    def dict_contains_only_real_elements(self) -> bool:
        species_dict_values = list(self.species_dict.values())
        return set(species_dict_values).issubset(self.allowed_element_labels)


class processor_base:
    def __init__(self, path: str, err_str: str | None = None) -> None:
        """Generic parent class for any processor classes

        :param path: Path to a file
        :type path: ``str``
        :param err_str: An error message to print if something goes wrong, defaults to None
        :type err_str: ``str | None``, optional
        """
        self.input_path = path.strip()
        self.abs_input_path: Path
        self.err_str = err_str
        self.re_float = re.compile(r"[-+]?\d*\.\d+")
        self.re_index = re.compile(r"\d+")

    def resolve_path(self, filename: str | None = None) -> None:
        """Checks for a file's existence and validity.

        :param filename: Path to a file, defaults to None, which will use the ``path`` passed into the class instance.
        :type filename: ``str | None``, optional
        :raises FileNotFoundError: Exits if the file is not found.
        :raises RuntimeError: Will exit if the file is corrupted
        """
        abs_coord_path: Path = Path("")
        if filename is None:
            abs_coord_path = Path(os.path.abspath(self.input_path))
        else:
            abs_coord_path = Path(os.path.abspath(filename))
        if not abs_coord_path.exists():
            raise FileNotFoundError(f"{abs_coord_path} not found.")
        if abs_coord_path.is_file() and os.stat(abs_coord_path).st_size <= 0:
            raise RuntimeError(f"{abs_coord_path} was an existing file, but has no file contents.")
        self.abs_input_path = abs_coord_path

    def open_file(self) -> None:
        pass

    def remove_blank_lines(self, lines: list[str]) -> list[str]:
        return [line for line in lines if line.strip()]


class conquest_coordinates:
    """Class which holds data about the system, including the :class:`Atom` s in the system.
    
    :param conquest_input: :class:`conquest_input` instance.
    :type conquest_input: ``conquest_input``
    """
    def __init__(
        self,
        conquest_input: conquest_input,
    ) -> None:
        
        self.atoms: list[Atom] = []
        self.conquest_input = conquest_input
        self.natoms: str
        self.element_map: dict[str, list[Atom]]
        self.lattice_vectors: c2at.REAL_ARRAY = np.array([])
        self.cart_position_vectors: c2at.REAL_ARRAY = np.array([])

    def get_cartesian_positions(self) -> c2at.REAL_ARRAY:
        """Returns the Cartesian position of all the atoms in the system and attaches to its :class:`Atom` field.

        :return: The 3D vector of the Cartesian position.
        :rtype: :ref:`REAL ARRAY <types>`
        """
        atom_frac_pos = np.vstack([atom.coords for atom in self.atoms])
        cart_coords = atom_frac_pos @ self.lattice_vectors.T
        for atom, cart_coord in zip(self.atoms, cart_coords):
            atom.cart_coords = cart_coord
        self.cart_position_vectors = cart_coords
        return cart_coords

    def assign_atom_labels(self) -> None:
        """Assign each Atom its label.
        If the ``conquest_input`` does not define labels for all species, these species labels will be silently skipped.
        """
        for atom in self.atoms:
            atom.label = self.conquest_input.species_dict[atom.species]

    def index_to_atom_map(self) -> None:
        """Form a dictionary with keys an element label, and values a list of all the :class:`Atom` s with that label. External file formats, such as ``.vasp``, require a count of the number of atoms per element."""
        ele_to_atom: dict[str, list[Atom]] = {}
        for element in self.conquest_input.unique_elements:
            # print(list(a for a in self.Atoms if a.label == element))
            ele_to_atom[element] = list(a for a in self.atoms if a.label == element)
        self.element_map = ele_to_atom

    def number_of_elements(self) -> dict[str, int]:
        """Function to get the number of atoms of each element. 

        :return: Returns a dictionary of the number of atoms per element
        :rtype: ``dict[str, int]``
        """
        num_eles: dict[str, int] = {}
        for element in list(self.element_map.keys()):
            num_eles[element] = len(self.element_map[element])
        return num_eles


class conquest_coordinates_processor(processor_base):
    """Class which extracts data from a CONQUEST coordinates file and populates a :class:`conquest_coordinates` instance.

        :param path: Path of the CONQUEST coordinates file to read.
        :type path: ``str``
        :param conquest_input: :class:`conquest_input` instance.
        :type conquest_input: conquest_input
    """
    def __init__(self, path: str, conquest_input: conquest_input) -> None:
        processor_base.__init__(
            self, path=path, err_str="Error opening specified CONQUEST coordinates file."
        )
        self.coords = conquest_coordinates(conquest_input=conquest_input)
        self.resolve_path()
        self.open_file()
        self.coords.get_cartesian_positions()
        self.coords.assign_atom_labels()
        self.coords.index_to_atom_map()
        self.volume_bohr: float = (
            self.coords.lattice_vectors[0][0] * self.coords.lattice_vectors[1][1]
            + self.coords.lattice_vectors[2][2]
        )
        self.volume_ang = self.volume_bohr * BOHR_TO_ANGSTROM_VOLUME

    def open_file(self) -> None:
        """This method reads a CONQUEST coordinate file.

        CONQUEST coords file split into 3 main chunks:
            * first 3 lines are lattice vectors
            * fourth line is the total number of atoms in the unit cell
            * the following lines describe each atom and look like
                <double> <double> <double> <int> <char> <char> <char>
        """
        with open(self.abs_input_path, "r", encoding="utf-8") as conquest_coord_file:
            conquest_lattice_data_str = [next(conquest_coord_file).strip() for _ in range(3)]
            cell_lattice_vectors: list[c2at.REAL_ARRAY] = []
            for lattice_vect in conquest_lattice_data_str:
                coords: c2at.REAL_ARRAY = np.fromstring(lattice_vect, sep=" ")
                cell_lattice_vectors.append(coords)
            self.coords.lattice_vectors = np.vstack(cell_lattice_vectors)
            self.coords.natoms = next(conquest_coord_file)
            atom_data = conquest_coord_file.readlines()
            atom_data_stripped = [atom for atom in atom_data if atom.strip()]
            atom_number = 1
            for atom in atom_data_stripped:
                split_atom_data = atom.strip().split()
                self.coords.atoms.append(
                    Atom(
                        species=int(split_atom_data[3]),
                        can_move=split_atom_data[4:],
                        coords=np.array(split_atom_data[:3]).astype(float),
                        number=atom_number,
                    )
                )
                atom_number += 1
        conquest_coord_file.close()
        self.coords.get_cartesian_positions()


class atom_charge(processor_base):
    """
    Class to process AtomCharge.dat from CONQUEST output files.

    It is assumed each row of AtomCharge.dat is arranged such that it is equivalent to the
    same CONQUEST input coordinates file.

    In particular, make use of the conquest_coordinates class to contain the list of Atoms
    
    :param coordinates: The :class:`conquest_coordinates` instance to use
    :type coordinates: :class:`conquest_coordinates`
    :param atom_charge_path: Path to the ``AtomCharge.dat`` file.
    :type atom_charge_path: ``str``
    """

    def __init__(self, coordinates: conquest_coordinates, atom_charge_path: str) -> None:
        processor_base.__init__(
            self,
            path=atom_charge_path,
            err_str="Error opening specified CONQUEST AtomCharge.dat file.",
        )
        self.coordinates = coordinates
        self.atom_charge_path = atom_charge_path
        self.abs_atom_charge_path: Path = Path(os.path.abspath(self.atom_charge_path))
        self.conquest_charge_data: list[c2at.REAL_ARRAY] = []
        self.resolve_path()
        self.open_file()
        self.assign_atom_charge()

    def open_file(self) -> None:
        """Reads in the AtomCharge.dat
        Format of the file is 3 columns: total, spin up spin down
        CONQUEST only deals with collinear spins
        """
        with open(self.abs_atom_charge_path, "r", encoding="utf-8") as conquest_charge_file:
            for line in conquest_charge_file:
                total_up_down = line.split()
                if total_up_down:
                    self.conquest_charge_data.append(np.array(total_up_down).astype(float))
                else:
                    continue

        conquest_charge_file.close()

    def assign_atom_charge(self) -> None:
        """
        Assign each Atom its spin values from the AtomCharge.dat file: up - down
        """
        for i, atom in enumerate(self.coordinates.atoms):
            split_charge_data = self.conquest_charge_data[i]
            atom.spins = np.array([0.0, 0.0, split_charge_data[1] - split_charge_data[2]])


class block_processor:
    """Generic class to process CONQUEST output files which are split into blocks via `&`

    Sometimes these blocks are categorised by spins and have comments
    # at the start of the section
    """

    def __init__(self) -> None:
        self.blocks: list[c2at.REAL_ARRAY] = []
        self.current_block: list[c2at.REAL_ARRAY] = []  # temp storage
        self.re_float = re.compile(r"[-+]?\d*\.\d+")
        self.re_index = re.compile(r"\d+")
        self.num_spins: int = 0

    def process_headers(self, line: str) -> None:
        """This function can be overridden by subclasses to process header lines starting with #
        Bandstructure and DOS/pDOS files have slightly different headers.
        """
        pass

    def process_block(self, line: str) -> None:
        """Extracts blocks, which are separated by ``&`` on its own newline, for data processing."""
        pass

    def read_file(self, filename: str) -> None:
        """Generic method which calls :func:`~process_headers` and :func:`~process_block`.

        :param filename: Path to file
        :type filename: ``str``
        """
        self.blocks = []
        with open(filename, "r", encoding="utf-8") as f:
            for line in f:
                stripped_line = line.strip()
                if not stripped_line:
                    continue
                if line.startswith("#"):
                    self.process_headers(line=stripped_line)
                else:
                    self.process_block(stripped_line)
