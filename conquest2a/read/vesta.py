"""
This module takes a VESTA file and converts it to CONQUEST coordinates.
Only reads in CELLP (cell params), STRUC (atom positions) and VECTR (arrows, treating as spin)
"""

import re
import numpy as np
import conquest2a._types as c2at
from conquest2a.constants import ANGSTROM_TO_BOHR
from conquest2a.conquest import processor_base, conquest_input, conquest_coordinates, Atom
from conquest2a.writers import conquest_writer


class vesta_to_conquest(processor_base):
    def __init__(self, vesta_path: str, output_path: str, conq_in: conquest_input) -> None:
        super().__init__(vesta_path)
        self.resolve_path()
        self.vesta_path = self.abs_input_path
        self.output_path = output_path
        self.content: str = ""  # holds block as giant string
        self.atom_to_vector: dict[int, int] = {}
        self.atom_to_spin: dict[int, int] = {}
        self.conquest_input = conq_in
        self.conq_coords = conquest_coordinates(self.conquest_input)
        self.content = self.vesta_path.read_text()
        self.parse_cellp()
        self.parse_struc()
        self.parse_vectr()
        self.assign_species()
        conquest_writer(dest=self.output_path, coords=self.conq_coords)
        print(f"Written {self.conq_coords.natoms} atoms to: {self.output_path}")

    def parse_cellp(self) -> None:
        match = re.search(
            r"CELLP\s*\n\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
            self.content,
        )
        if not match:
            raise ValueError("CELLP block not found.")
        a, b, c, _, _, _ = (float(match.group(i)) for i in range(1, 7))
        # CONQUEST only supports orthorhombic cells, assume angles are 90 degrees
        a_b = a * ANGSTROM_TO_BOHR
        b_b = b * ANGSTROM_TO_BOHR
        c_b = c * ANGSTROM_TO_BOHR
        self.conq_coords.lattice_vectors = np.array(
            [[a_b, 0.0, 0.0], [0.0, b_b, 0.0], [0.0, 0.0, c_b]]
        )

    def parse_struc(self) -> None:
        match = re.search(r"STRUC\n(.*?)\n\s*0 0 0 0 0 0 0", self.content, re.DOTALL)
        if not match:
            raise ValueError("STRUC block not found or not properly terminated.")
        lines = match.group(1).splitlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            i += 1
            if not line:
                continue
            parts = line.split()
            # VESTA prints each of these lines identically
            if len(parts) >= 7 and parts[0].lstrip("-").isdigit():
                sp_dict = self.conquest_input.species_dict.items()
                species_int = [k for k, v in sp_dict if v == str(parts[1])][0]
                # choose first element as species int, fix spin later
                self.conq_coords.atoms.append(
                    Atom(
                        species=species_int,
                        can_move=["T", "T", "T"],
                        number=int(parts[0]),
                        coords=np.array([float(parts[4]), float(parts[5]), float(parts[6])]),
                        label=str(parts[1]),
                    )
                )
                i += 1
        self.conq_coords.natoms = str(len(self.conq_coords.atoms))

    def parse_vectr(self) -> None:
        match = re.search(r"VECTR\n(.*?)\nVECTT", self.content, re.DOTALL)
        if not match:
            return
        vector_to_sign: dict[int, int] = {}
        current_vec: int | None = None
        current_sign: int = 0
        for line in match.group(1).splitlines():
            parts = line.strip().split()
            if not parts:
                continue
            if parts == ["0", "0", "0", "0", "0"]:
                current_vec = None
                current_sign = 0
                continue
            if len(parts) == 5 and "." in parts[1]:
                current_vec = int(parts[0])
                vz = float(parts[3])
                current_sign = 1 if vz > 0 else (-1 if vz < 0 else 0)
                vector_to_sign[current_vec] = current_sign
                continue
            if current_vec is not None and len(parts) == 5:
                self.atom_to_spin[int(parts[0])] = current_sign

    def build_spin_species_map(self) -> dict[tuple[str, int], int]:
        element_to_species: dict[str, list[int]] = {}
        for sp_idx, element in self.conquest_input.species_dict.items():
            element_to_species.setdefault(element, []).append(sp_idx)
        for indices in element_to_species.values():
            indices.sort()

        spin_species_map: dict[tuple[str, int], int] = {}
        for element, indices in element_to_species.items():
            if len(indices) == 1:
                spin_species_map[(element, 1)] = indices[0]
                spin_species_map[(element, -1)] = indices[0]
                spin_species_map[(element, 0)] = indices[0]
            elif len(indices) == 2:
                spin_species_map[(element, 1)] = indices[0]
                spin_species_map[(element, -1)] = indices[1]
                spin_species_map[(element, 0)] = indices[0]
            else:
                raise ValueError(
                    f"Element '{element}' maps to more than 2 species indices {indices}; "
                    "only collinear  spin is supported."
                )
        return spin_species_map

    def assign_species(self) -> None:
        spin_species_map = self.build_spin_species_map()
        for atom in self.conq_coords.atoms:
            spin = self.atom_to_spin.get(atom.number, 0)
            key = (atom.label, spin)
            if key not in spin_species_map:
                raise ValueError(
                    f"No species found for element '{atom.label}' with spin {spin}. "
                    "Check your species_dict."
                )
            atom.species = spin_species_map[key]  # set species for spin


if __name__ == "__main__":
    species = {1: "O", 2: "Bi", 3: "Co", 4: "Co", 5: "Mn", 6: "Mn"}
    conqin = conquest_input(species_dict=species)
    vesta_to_conquest(
        "tests/data/test_vesta_to_conquest.vesta",
        "tests/data/test_vesta_to_conquest.coords",
        conqin,
    )
