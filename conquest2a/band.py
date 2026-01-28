from dataclasses import dataclass, field
import re
import numpy as np
from conquest2a.conquest import block_processor
import conquest2a._types as c2at


@dataclass
class band:
    spin: int
    index: int  # band index
    kpoint: c2at.INT_ARRAY = field(default_factory=lambda: np.array([0, 0, 0]))  # kpoint index.
    energies: c2at.REAL_ARRAY = field(
        default_factory=lambda: np.empty((1, 3))
    )  # energies at that kpoint for this band


class bst_processor(block_processor):
    def __init__(self, bst_file: str) -> None:
        super().__init__()
        self.bst_file = bst_file
        self.blocks: list[c2at.REAL_ARRAY] = []
        self.bands: list[band] = []
        self.fermi_level: float = 0.0
        self.is_shifted_to_fermi: bool = True
        self.block_counter: int = 0
        self.current_spin = 0
        self.read_file(filename=self.bst_file)

    def process_headers(
        self,
        line: str,
        num_spins: int,
    ) -> None:
        if "# Spin" in line:
            self.num_spins += 1
            self.current_spin += 1
        if "# Original" in line:
            result = re.findall(self.re_float, line)
            self.fermi_level = float(result[0])
        if "# Band " in line:
            band_index = re.findall(self.re_index, line)
            self.bands.append(band(index=int(band_index[0]), spin=self.current_spin))
        if not line.startswith("# Bands shifted"):
            self.is_shifted_to_fermi = False

    def process_block(self, line: str) -> None:
        if line == "&":
            if self.current_block:
                arrayed_block = np.array(self.current_block, dtype=float)
                self.bands[self.block_counter].kpoint = arrayed_block[:, 0]
                self.bands[self.block_counter].energies = arrayed_block[:, 1]
                self.block_counter += 1
                self.current_block = []
        else:
            self.current_block.append(np.array(line.split()).astype(np.float32))
