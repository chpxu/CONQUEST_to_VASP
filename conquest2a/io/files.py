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
from conquest2a.io.writers import conquest_writer



def main() -> None:
    kcuf3spec = conquest_species({1: "K", 2: "Cu", 3: "Cu", 4: "F"})
    cell_to_conquest("../tests/data/files/kcuf3.cell", kcuf3spec, "../tests/data/files/kcuf3.dat")
    cell_to_conquest(
        "../tests/data/files/kcuf3_spin.cell", kcuf3spec, "../tests/data/files/kcuf3_spin.dat"
    )


if __name__ == "__main__":
    main()
