from io import TextIOWrapper
from multiprocessing import Value
from typing import IO, Any, Callable, Literal
from collections.abc import Iterator
import sys

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing_extensions import override
import numpy as np
from ase.units import Bohr
import conquest2a._types as c2at
from conquest2a.constants import BOHR_TO_ANGSTROM
from conquest2a.conquest import Atom, conquest_coordinates, conquest_species, processor_base
from conquest2a.writers import conquest_writer

class coords_to_conquest(processor_base):
    def __init__(self, path: str, species: conquest_species,
        destination: str,
        precision: int = 10,err_str: str | None = None, ) -> None:
        super().__init__(path, err_str)
        self._ENDBLOCK: str = r"%ENDBLOCK"

    def _resolve_species_id(self, label: str) -> int:
        if label not in self.species.allowed_element_labels:
            raise ValueError(f"'{label}' is not a valid element.")
        matches = sorted(idx for idx, lbl in self.species.species_dict.items() if lbl == label)
        if not matches:
            raise ValueError(f"Element '{label}' is not present in the supplied species map.")
        return matches[0]

    def _resolve_species_names(self, counts: list[int], elements: list[str] | None) -> list[str]:
        if elements is not None:
            return [element.strip() for element in elements]
        labels: list[str] = []
        seen: set[str] = set()
        for _, label in sorted(self.species.species_dict.items()):
            if label not in seen:
                labels.append(label)
                seen.add(label)
        if len(labels) != len(counts):
            raise RuntimeError(
                "Number of element types and number of counts are mismatched"
            )
        return labels


class cell_to_conquest(processor_base):
    def __init__(
        self,
        path: str,
        species: conquest_species,
        destination: str,
        precision: int = 10,
    ) -> None:
        super().__init__(path)
        self._ENDBLOCK: str = r"%ENDBLOCK"
        self.block_handlers: dict[str, Callable[[Iterator[str]], None]] = {
            "LATTICE_CART": self._read_lattice_cart,
            # "LATTICE_ABC": self._read_lattice_abc,
            "POSITIONS_FRAC": self._read_positions_frac,
            # "POSITIONS_ABS": self._read_positions_abs,
        }
        self.fix_ions: bool = False
        self._atom_counter: int = 0
        self.resolve_path()
        self.coords: conquest_coordinates = conquest_coordinates(species)
        self.read_cell()
        self.coords.assign_atom_labels()
        self.coords.get_cartesian_positions()
        self.coords.natoms = str(self._atom_counter)
        conquest_writer(destination, self.coords, precision=precision)

    def _read_lattice_cart(self, lines: Iterator[str]) -> None:
        units: str = "BOHR"
        vectors: list[c2at.REAL_ARRAY] = []

        for line in self._block_lines(lines):
            line = line.strip()
            if len(vectors) == 0 and line.isalpha():
                units = line
                continue
            vector = np.fromstring(line, sep=" ")

            if len(vector) != 3:
                raise RuntimeError("Expected three position components.")
            vectors.append(vector)

        if len(vectors) != 3:
            raise RuntimeError("LATTICE_CART must contain exactly three lattice vectors.")

        lattice = np.vstack(vectors)

        self.coords.lattice_vectors = lattice if units.upper() == "BOHR" else lattice / Bohr

    def _read_positions_frac(self, lines: Iterator[str]) -> None:
        units = None
        for line in self._block_lines(lines):
            atom: Atom = self._parse_atom(line)
            self.coords.atoms.append(atom)
        # Sort by atom number for convenience
        self.coords.atoms.sort(key=lambda d: d.number)

    def _skip_block(self, lines: Iterator[str]) -> None:
        for line in lines:
            if line.strip().upper().startswith(self._ENDBLOCK):
                return

    def _parse_atom(self, line: str) -> Atom:
        """
        Parse an atom positions line of the form:

        ::

            Element label xpos ypos zpos
            Element label xpos ypos zpos SPIN spinval

        A line of the form ``Element label xpos ypos zpos SPIN sx sy sz`` will raise an error as non-collinear spin is unsupported in CONQUEST.

        :param line: the atom line with or without `SPIN` flag
        :type line: ``str ``
        :return: ``Atom`` instance of the atom
        :rtype: ``Atom``
        :raises ValueError: Exits if the line does not have the correct format or atom is wrong
        :raises ValueError: Exits if non-collinear spin is detected
        """
        parts: list[str] = line.split()
        if len(parts) < 4:
            raise ValueError(f"Atom line is formatted incorrectly: {line!r}")

        label = parts[0]
        if label not in self.coords.conquest_input.allowed_element_labels:
            raise ValueError("Element label was not a valid element")

        position = np.array([float(x) for x in parts[1:4]])
        # No spin by default
        spin: c2at.REAL_ARRAY = np.array([0.0, 0.0, 0.0])
        rest: list[str] = parts[4:]
        if len(rest) > 0:
            if rest[0].upper() != "SPIN":
                raise ValueError(f"Unexpected trailing content in atom line: {line!r}")
            spin_values = rest[1:]
            if len(spin_values) == 1:
                spin = np.array([0.0, 0.0, float(spin_values[0])])
            elif len(spin_values) == 3:
                raise ValueError(f"Non-collinear spin is not supported: {line!r}")
            else:
                raise ValueError(
                    f"SPIN flag must be followed by 1 (collinear) or 3 (non-collinear) "
                    + f"values, got {len(spin_values)}: {line!r}"
                )
        move_line = ["F", "F", "F"] if self.fix_ions else ["T", "T", "T"]
        spin_species: list[int] = self.coords.conquest_input.element_to_species_dict[label]
        atom_line: Atom = Atom(
            species=spin_species[-1] if spin[2] < 0.0 else spin_species[0],
            number=self._atom_counter,
            coords=position,
            label=label,
            can_move=move_line,
            spins=spin,
        )
        self._atom_counter += 1

        return atom_line

    def _block_lines(
        self,
        lines: Iterator[str],
    ) -> Iterator[str]:
        for line in lines:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if line.upper().startswith(self._ENDBLOCK):
                return

            yield line

    def read_cell(self) -> None:
        with open(self.abs_input_path) as f:
            lines = iter(f)
            for raw_line in lines:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                upper = line.upper()
                if upper.startswith("FIX_ALL_IONS"):
                    self.fix_ions = True if upper.split()[1] == "TRUE" else False
                    continue
                if upper.startswith("%BLOCK"):
                    block = upper.split()[1]
                    handler = self.block_handlers.get(block, self._skip_block)
                    handler(lines)

class poscar_to_conquest(coords_to_conquest):
    def __init__(
        self,
        path: str,
        species: conquest_species,
        destination: str,
        precision: int = 10,
    ) -> None:
        self.species: conquest_species = species
        super().__init__(path, self.species, path)
        self._atom_counter: int = 0
        self.resolve_path()
        self.coords: conquest_coordinates = conquest_coordinates(species)
        self.read_poscar()
        self.coords.assign_atom_labels()
        self.coords.get_cartesian_positions()
        self.coords.natoms = str(self._atom_counter)
        conquest_writer(destination, self.coords, precision=precision)

    def read_poscar(self) -> None:
        with open(self.abs_input_path) as f:
            lines = [line.rstrip("\n") for line in f]
        # Skip comment line
        line_idx = 1
        scale = lines[line_idx].split()
        if len(scale) > 1:
            raise ValueError("There should only be one universal scale factor")
        scale = int(scale[0])
        line_idx += 1
        # Lattice vectors
        x = [float(v) for v in lines[line_idx].split()[:3]]
        y = [float(v) for v in lines[line_idx + 1].split()[:3]]
        z = [float(v) for v in lines[line_idx + 2].split()[:3]]
        lattice_raw = np.array([x,y,z])
        line_idx += 3

        lattice_bohr = lattice_raw * scale / Bohr
        self.coords.lattice_vectors = lattice_bohr 
        # Ion species and numbers
        elements = lines[line_idx].split()  # store elements
        line_idx += 1
        counts = [int(t) for t in lines[line_idx].split()] # store element counts
        line_idx += 1

        labels = self._resolve_species_names(counts, elements)

        selective_dynamics = lines[line_idx].strip()[:1].upper() == "S"
        if selective_dynamics:
            line_idx += 1

        mode_char = lines[line_idx].strip()[:1].upper()
        cartesian = mode_char in ("C", "D")
        line_idx += 1

        inv_lattice_t = np.linalg.inv(self.coords.lattice_vectors).T

        for label, count in zip(labels, counts):
            species_id = self._resolve_species_id(label)
            for _ in range(count):
                parts: list[str] = lines[line_idx].split()
                line_idx += 1
                coord_vals = np.array([float(v) for v in parts[:3]])
                frac =coord_vals if not cartesian else  (coord_vals  @ inv_lattice_t)
                frac = frac % 1.0
                can_move = parts[3:6] if selective_dynamics and len(parts) >= 6 else ["T", "T", "T"]
                atom = Atom(
                    species=species_id,
                    number=self._atom_counter,
                    coords=frac,
                    label=label,
                    can_move=can_move,
                )
                self.coords.atoms.append(atom)
                self._atom_counter += 1

            
def main() -> None:
    kcuf3spec = conquest_species({1: "K", 2: "Cu", 3: "Cu", 4: "F"})
    cell_to_conquest("../tests/data/files/kcuf3.cell", kcuf3spec, "../tests/data/files/kcuf3.dat")
    cell_to_conquest(
        "../tests/data/files/kcuf3_spin.cell", kcuf3spec, "../tests/data/files/kcuf3_spin.dat"
    )


if __name__ == "__main__":
    main()
