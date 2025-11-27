from pathlib import Path
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt
import re


@dataclass
class band:
    spin: int
    index: int # band index
    kpoint: list[int] | npt.NDArray # kpoint index. 
    energies: list[float | int] | npt.NDArray # energies at that kpoint for this band

class bst_processor:
    def __init__(self, bst_file: str | Path) -> None:
        self.bst_file = bst_file
        self.blocks: list[np.ndarray | list[float | int]]= []
        self.bands: list[band] = []
        self.fermi_level: float = 0.0
        self.read_bst_file()
    def read_bst_file(self) -> None:
        with open(self.bst_file, "r", encoding="utf-8") as f:
            current_block = []
            num_spins = 0
            block_counter = 0
            for line in f:
                stripped_line = line.strip()
                if "# Spin" in stripped_line:
                    num_spins += 1
                    continue
                if "# Original" in stripped_line:
                    result = re.findall(r"\d+\.\d+", stripped_line)
                    self.fermi_level = float(result[0])
                    continue
                if "# Band " in stripped_line:
                    band_index = re.findall(r"\d+", stripped_line)
                    self.bands.append(band(index=int(band_index[0]), kpoint=[], energies=[], spin = num_spins))
                    continue
                if not stripped_line or stripped_line.startswith("# Bands shifted"):
                    continue
                if stripped_line == "&":
                    if current_block:
                        # self.blocks.append(np.array(current_block, dtype=float))
                        arrayed_block  = np.array(current_block, dtype=float)
                        self.bands[block_counter].kpoint = arrayed_block[:, 0]
                        self.bands[block_counter].energies = arrayed_block[:, 1]
                        block_counter += 1
                        current_block = []
                else:
                    current_block.append(line.split())
            if current_block:
                self.blocks.append(np.array(current_block, dtype=float))
            self.num_spins = num_spins
