from ctypes import Union
from dataclasses import dataclass
import os
from os.path import abspath
from pathlib import Path
from typing import Literal
 
@dataclass
class Atom:
    species: str
    coords: list[float]
    can_move: list[Literal['T', 'F']]
    label: str = ""


class CONQUEST_INPUT:
    def __init__(self, species_dict: dict[str, str]) -> None:
        """Constructor for CONQUEST_INPUT

        Args:
            species_dict (dict[int, str]): dict mapping the species index to an element label
            It is not completely reliable to directly read Conquest_input because in order to have multiple spins, you essentially must duplicate an element (and call it a new Conquest species) - however the labels are free to be named. Moreover, in the case isotopes are used, differentiating by elemental mass is not sufficient.

            Therefore we expect a dictionary to be passed in independently. This can be useful for large-scale post-processing too.
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
        assert self.dict_contains_only_real_elements()
        self.unique_elements = list(set(self.species_dict.values()))

    def dict_contains_only_real_elements(self) -> bool:
        species_dict_values = list(self.species_dict.values())
        return set(species_dict_values).issubset(self.allowed_element_labels)


class CONQUEST_COORDINATES:
    def __init__(
        self,
        CONQUEST_input: CONQUEST_INPUT,
    ) -> None:
        self.Atoms: list[Atom] = []
        self.CONQUEST_input = CONQUEST_input
        self.lattice_vectors: list[list[float]] = []
        self.natoms: str
        self.element_map: dict[str, list[Atom]]
    def assign_atom_labels(self) -> None:
        """Assign each Atom its label"""
        for atom in self.Atoms:
            atom.label = self.CONQUEST_input.species_dict[atom.species]

    def index_to_atom_map(self) -> None:
        """Every Atom now has its element label. External file formats require a count of the number of Atoms per element, so we now form a dict of elements to Atoms in preparation for writing"""
        ele_to_atom: dict[str, list[Atom]] = dict()
        for element in self.CONQUEST_input.unique_elements:
            print(list(a for a in self.Atoms if a.label == element))
            ele_to_atom[element] = list(a for a in self.Atoms if a.label == element)
        self.element_map = ele_to_atom

    def number_of_elements(self) -> dict[str, int]:
        num_eles: dict[str, int] = dict()
        for element in list(self.element_map.keys()):
            num_eles[element] = len(self.element_map[element])
        return num_eles


class CONQUEST_COORDINATES_PROCESSOR(CONQUEST_COORDINATES):
    def __init__(self, path: Path | str, CONQUEST_input: CONQUEST_INPUT) -> None:
        super().__init__(CONQUEST_input)
        self.input_coord_path = path
        self.open_file()
        self.assign_atom_labels()
        self.index_to_atom_map()

    def resolve_path(self) -> bool:
        try:
            abs_coord_path = Path(abspath(self.input_coord_path))
            assert abs_coord_path.exists() is True
            assert abs_coord_path.is_file() is True
            assert os.stat(abs_coord_path).st_size > 0
            self.abs_coord_path = abs_coord_path
            return True
        except Exception as e:
            print(e)
            return False

    def open_file(self) -> None:
        """
        CONQUEST coords file split into 3 main chunks:
        first 3 lines are lattice vectors
        fourth line is the total number of atoms in the unit cell
        the following lines are of the for_summary_m:
          <double> <double> <double> <int> <char> <char> <char>
        """
        if self.resolve_path():
            with open(self.abs_coord_path, "r", encoding="utf-8") as CONQUEST_coord_file:
                conquest_lattice_data_str = [next(CONQUEST_coord_file).strip() for _ in range(3)]
                print(conquest_lattice_data_str)
                for lattice_vect in conquest_lattice_data_str:
                    coords = lattice_vect.split()
                    self.lattice_vectors.append([float(x) for x in coords])
                self.natoms = next(CONQUEST_coord_file)
                atom_data = CONQUEST_coord_file.readlines()
                for atom in atom_data:
                    split_atom_data = atom.strip().split()
                    self.Atoms.append(
                        Atom(
                            species=split_atom_data[3],
                            can_move=split_atom_data[4:],
                            coords=list(map(float, split_atom_data[:3])),
                        )
                    )
            CONQUEST_coord_file.close()
            return
        print("Error opening specified CONQUEST coordinates file")

   