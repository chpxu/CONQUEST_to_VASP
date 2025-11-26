from dataclasses import dataclass, field
import os
from os.path import abspath
from pathlib import Path
from collections.abc import Sequence
@dataclass
class Atom:
    species: str
    coords: list[float | int]
    can_move: Sequence[str]
    spins: list[float | int ]  = field(default_factory=lambda: [0.0, 0.0, 0.0])
    label: str = ""

        
class conquest_input:
    def __init__(self, species_dict: dict[str, str]) -> None:
        """Constructor for conquest_input

        Args:
            species_dict (dict[int, str]): dict mapping the species index to an element label
                It is not completely reliable to directly read conquest_input
                E.g., for multiple spins, must duplicate an element and call it a new Conquest species,
                however labels can be any alphanumeric string
                Therefore we expect a dictionary to be passed in independently.
        """
        self.species_dict = species_dict
        self.allowed_element_labels: list[str] = [
            "H",
            "He",
            "Ne",
            "Ar",
            "Kr",
            "Xe",
            "Rn",
            "Og",
            "Li",
            "Na",
            "K",
            "Rb",
            "Cs",
            "Fr",
            "Be",
            "Mg",
            "Ca",
            "Sr",
            "Ba",
            "Ra",
            "B",
            "C",
            "N",
            "O",
            "F",
            "Al",
            "Si",
            "P",
            "S",
            "Cl",
            "Ga",
            "Ge",
            "As",
            "Se",
            "Br",
            "In",
            "Sn",
            "Sb",
            "Te",
            "I",
            "Tl",
            "Pb",
            "Bi",
            "Po",
            "At",
            "Nh",
            "Fl",
            "Mc",
            "Lv",
            "Ts",
            "Sc",
            "Ti",
            "V",
            "Cr",
            "Mn",
            "Fe",
            "Co",
            "Ni",
            "Cu",
            "Zn",
            "Y",
            "Zr",
            "Nb",
            "Mo",
            "Tc",
            "Ru",
            "Rh",
            "Pd",
            "Ag",
            "Cd",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Pm",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
            "La",
            "Th",
            "Pa",
            "U",
            "Np",
            "Pu",
            "Am",
            "Cm",
            "Bk",
            "Cf",
            "Es",
            "Fm",
            "Md",
            "No",
            "Lr",
            "Rf",
            "Db",
            "Sg",
            "Bh",
            "Hs",
            "Mt",
            "Ds",
            "Rg",
            "Cn",
        ]
        if not self.dict_contains_only_real_elements():
            raise Exception("Provided species map contains fake chemical elements.")
        self.unique_elements = list(set(self.species_dict.values()))

    def dict_contains_only_real_elements(self) -> bool:
        species_dict_values = list(self.species_dict.values())
        return set(species_dict_values).issubset(self.allowed_element_labels)

class processor_base:
    def __init__(
        self, path: Path | str, err_str: str | None = None
    ) -> None:
        self.input_path = path
        self.abs_input_path: str | Path
        self.err_str = err_str
    def resolve_path(self) -> None:
        try:
            abs_coord_path = Path(abspath(self.input_path))
            assert abs_coord_path.exists() is True
            assert abs_coord_path.is_file() is True
            assert os.stat(abs_coord_path).st_size > 0
            self.abs_input_path = abs_coord_path
        except FileNotFoundError as e:
            print(e)
            if self.err_str is not None:
                print(self.err_str)
    def open_file(self) -> None:
        pass
    def remove_blank_lines(self, lines: list[str]) -> list[str]:
        return [line for line in lines if line.strip()]
class conquest_coordinates:
    def __init__(
        self,
        conquest_input: conquest_input,
    ) -> None:
        self.Atoms: list[Atom] = []
        self.conquest_input = conquest_input
        self.lattice_vectors: list[list[float]] = []
        self.natoms: str
        self.element_map: dict[str, list[Atom]]
    def assign_atom_labels(self) -> None:
        """Assign each Atom its label"""
        for atom in self.Atoms:
            atom.label = self.conquest_input.species_dict[atom.species]

    def index_to_atom_map(self) -> None:
        """Every Atom now has its element label. External file formats require a count of the number of Atoms per element, so we now form a dict of elements to Atoms in preparation for writing"""
        ele_to_atom: dict[str, list[Atom]] = {}
        for element in self.conquest_input.unique_elements:
            print(list(a for a in self.Atoms if a.label == element))
            ele_to_atom[element] = list(a for a in self.Atoms if a.label == element)
        self.element_map = ele_to_atom

    def number_of_elements(self) -> dict[str, int]:
        num_eles: dict[str, int] = {}
        for element in list(self.element_map.keys()):
            num_eles[element] = len(self.element_map[element])
        return num_eles


class conquest_coordinates_processor(conquest_coordinates, processor_base):
    def __init__(self, path: Path | str, conquest_input: conquest_input) -> None:
        conquest_coordinates.__init__(self,conquest_input=conquest_input)
        processor_base.__init__(self, path=path, err_str="Error opening specified CONQUEST coordinates file.")
        # super().__init__(conquest_input=conquest_input, path=path)
        # self.input_coord_path = path
        # self.abs_coord_path: str | Path
        try:
            self.resolve_path()
            self.open_file()
            self.assign_atom_labels()
            self.index_to_atom_map()
        except FileNotFoundError as e:
            print(e)

    def open_file(self) -> None:
        """
        CONQUEST coords file split into 3 main chunks:
        first 3 lines are lattice vectors
        fourth line is the total number of atoms in the unit cell
        the following lines are of the for_summary_m:
          <double> <double> <double> <int> <char> <char> <char>
        """
        with open(self.abs_input_path, "r", encoding="utf-8") as conquest_coord_file:
            conquest_lattice_data_str = [next(conquest_coord_file).strip() for _ in range(3)]
            for lattice_vect in conquest_lattice_data_str:
                coords = lattice_vect.split()
                self.lattice_vectors.append([float(x) for x in coords])
            self.natoms = next(conquest_coord_file)
            atom_data = conquest_coord_file.readlines()
            atom_data_stripped = [atom for atom in atom_data if atom.strip()]
            for atom in atom_data_stripped:
                split_atom_data = atom.strip().split()
                self.Atoms.append(
                    Atom(
                        species=split_atom_data[3],
                        can_move=split_atom_data[4:],
                        coords=list(map(float, split_atom_data[:3])),
                    )
                )
        conquest_coord_file.close()
           
class atom_charge(processor_base):
    """
    Class to process AtomCharge.dat from CONQUEST output files.

    It is assumed each row of AtomCharge.dat is arranged such that it is equivalent to the CONQUEST input coordinates file.

    In particular, make use of the conquest_cordinates class to contain the list of Atoms 
    """
    def __init__(self, coordinates: conquest_coordinates, atom_charge_path: Path | str) -> None:
        processor_base.__init__(self, path=atom_charge_path, err_str="Error opening specified CONQUEST AtomCharge.dat file.")
        self.coordinates = coordinates
        self.atom_charge_path = atom_charge_path
        self.abs_atom_charge_path: str | Path
        self.conquest_charge_data: list[list[float | int]] = []

        try:
            self.resolve_path()
            self.open_file()
        except FileNotFoundError as e:
            print(e)
        self.assign_atom_charge()
    
    def open_file(self) -> None:
        """
        CONQUEST AtomCharge.dat
          <double> <double> <double>
          total, "up", "down" - CONQUEST only deals with collinear spins
        """
        with open(self.abs_atom_charge_path, "r", encoding="utf-8") as conquest_charge_file:
            for line in conquest_charge_file:
                total_up_down = line.split()
                if total_up_down:
                    total_up_down = [float(x) for x in total_up_down]
                    self.conquest_charge_data.append(total_up_down)
                else:
                    continue
            
        conquest_charge_file.close()
    def assign_atom_charge(self) -> None:
        """
        Assign each Atom its spin values from the AtomCharge.dat file: up - down
        """
        for i, atom in enumerate(self.coordinates.Atoms):
            split_charge_data = self.conquest_charge_data[i]
            atom.spins = [0.0,0.0, split_charge_data[1]-split_charge_data[2]]