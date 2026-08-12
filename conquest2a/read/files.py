from io import TextIOWrapper
from typing import IO, Any, Callable, Literal
from collections.abc import Iterator
import sys
if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing_extensions import override
import numpy as np
import conquest2a._types as c2at
from conquest2a.constants import BOHR_TO_ANGSTROM
from conquest2a.conquest import Atom, conquest_coordinates, conquest_species, processor_base
from conquest2a.writers import conquest_writer


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
            "LATTICE_ABC": self._read_lattice_abc,
            "POSITIONS_FRAC": self._read_positions_frac,
            "POSITIONS_ABS": self._read_positions_abs,
        }
        self.fix_ions: bool = False
        self.resolve_path()
        self.coords: conquest_coordinates = conquest_coordinates(species)
        self.read_cell()
        self.coords.assign_atom_labels()
        self.coords.get_cartesian_positions()
        conquest_writer(destination, self.coords, precision=precision)

    def _read_lattice_cart(self, lines: Iterator[str]) -> None:
      units = "ANG"
      vectors: list[c2at.REAL_ARRAY] = []

      for line in self._block_lines(lines):
          line = line.strip()

          if len(vectors) == 0 and self._is_unit(line):
              units = line
              continue
          vector = np.fromstring(line, sep=" ")

          if vector.size != 3:
              raise RuntimeError(
                  "Expected three components."
              )
          vectors.append(vector)

      if len(vectors) != 3:
          raise RuntimeError(
              "LATTICE_CART must contain exactly three lattice vectors."
          )

      lattice = np.vstack(vectors)

      self.coords.lattice_vectors = self._convert_length_units(
          lattice,
          units,
      )

    def _read_positions_frac(self, lines: Iterator[str]) -> None:
        units = None
        for line in self._block_lines(lines):
          atom: Atom = self._parse_atom(line, fractional=True)
          self.coords.atoms.append(atom)

    def _skip_block(self, lines: Iterator[str]) -> None:
        for line in lines:
            if line.strip().upper().startswith(self._ENDBLOCK):
                return
    def _parse_atom(self, line) -> Atom:
      pass
#       ├── _parse_atom_line()
# ├── _parse_spin()
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
                    self.fix_ions = upper.split()[1] == "TRUE"
                    continue
                if upper.startswith("%BLOCK"):
                    block = upper.split()[1]
                    handler = self.block_handlers.get(block, self._skip_block)
                    handler(lines)
