from pathlib import Path
from typing import Literal
from os.path import abspath
from conquest2a import conquest
import numpy as np
import numpy.typing as npt
import re
import os
from conquest2a.conquest import block_processor
import conquest2a._types as c2at


class pdos_processor(block_processor):
    def __init__(self, conquest_rundir: str | Path, lm: Literal["lm", "l", "t"] = "t") -> None:
        # self.dos_file = dos_file
        self.blocks: list[c2at.REAL_ARRAY] = []
        # self.read_pdos_file()
        self.all_pdos_files: list[str] = []
        self.conquest_rundir = conquest_rundir
        self.lm = lm
        super().__init__()
        self.resolve_path()
        self.locate_pdos_files()

    def clear_pdos(
        self, pdos: dict[str, list[c2at.REAL_ARRAY]]
    ) -> dict[str, list[c2at.REAL_ARRAY]]:
        return {key: [] for key in pdos}

    def process_headers(self, line: str, num_spins: int) -> None:
        if "# Spin" in line:
            num_spins += 1
        if "# Original" in line:
            result = re.findall(self.re_float, line)
            self.fermi_level = float(result[0])
        if not line.startswith("# DOS shifted"):
            self.is_shifted_to_fermi = False

    def process_block(self, line: str) -> None:
        if line == "&":
            if self.current_block:
                self.blocks.append(np.array(self.current_block, dtype=float))
                self.current_block = []
        else:
            self.current_block.append(np.array(line.split()).astype(np.float64))

    def resolve_path(self) -> Path:
        abs_run_path = Path(abspath(self.conquest_rundir))
        if abs_run_path.exists():
            return abs_run_path
        else:
            raise FileNotFoundError(
                f'Conquest directory specified: "{abs_run_path}", does not exist.'
            )

    def locate_pdos_files(self) -> list[str]:
        if len(self.all_pdos_files) > 0:
            # Reset pdos files array, i.e. if changing directory
            self.all_pdos_files = []
        if self.lm == "t":
            self.all_pdos_files = ["DOS.dat"]
            return self.all_pdos_files
        abs_path = self.resolve_path()
        pdos_rgx = re.compile(rf"Atom[0-9]{{7}}DOS\_{self.lm}\.dat")
        file_list: list[str] = []
        for root, dirs, files in os.walk(abs_path, topdown=True):
            file_list = files
            break
        pdos_file_list: list[str] = []
        for filename in file_list:
            res = re.match(pdos_rgx, filename)
            if res:
                pdos_file_list.append(filename)
        pdos_file_list = sorted(pdos_file_list)
        for file in pdos_file_list:
            self.all_pdos_files.append(f"{abspath(self.conquest_rundir)}/{file}")
        return self.all_pdos_files


class pdos_l_processor(pdos_processor):
    def __init__(self, conquest_rundir: str | Path) -> None:
        self.num_spins: int = 0
        super().__init__(conquest_rundir=conquest_rundir, lm="l")
        self.l_dict: dict[str, list[c2at.REAL_ARRAY]] = {}
        # PDOS file is split into blocks separated by "&" lines
        # The first column of each block is the energy values
        # The second column is the total PDOS, i.e sum of all l, at that energy
        # Subsequent columns are the PDOS values for each l component, sorted in
        # ascending order of l
        # e.g., l = 0,  l =1,  l = 2, etc.
        # So dict will be of the form {"l": [array(spin1), array(spin2), ...],}
        self.energy_values: dict[int, c2at.REAL_ARRAY] = {}
        # self.l_map()

    def l_map(self) -> None:
        # From column 3, there are only l-contributions, and is sorted by ascending l values, and each row is just ech l-contribution at that energy
        l_dict: dict[str, list[c2at.REAL_ARRAY]] = {}
        self.l_dict = self.clear_pdos(self.l_dict)
        for idx, block in enumerate(self.blocks):
            energy = block[:, 0]
            self.energy_values[idx + 1] = energy
            pdos_values = block[:, 2:]
            num_l = pdos_values.shape[1]
            for l in range(num_l):
                if str(l) not in l_dict:
                    l_dict[str(l)] = []
                l_dict[str(l)].append(pdos_values[:, l])
        self.l_dict = l_dict


class pdos_lm_processor(pdos_processor):
    def __init__(self, conquest_rundir: str | Path) -> None:
        self.num_spins: int = 0
        super().__init__(conquest_rundir=conquest_rundir, lm="lm")
        self.lm_dict: dict[str, list[c2at.REAL_ARRAY]] = {}
        # PDOS file is split into blocks separated by "&" lines
        # The first column of each block is the energy values
        # The second column is the total PDOS, i.e sum of all l and m, at that energy
        # Subsequent columns are the PDOS values for each lm component, sorted in
        # ascending order of l and m
        # e.g., l = 0,  l =1, m = -1, 0, 1, l = 2, m = -2, -1, 0, 1, 2, etc.
        # So dict will be of the form {"l,m": [array(spin1), array(spin2), ...],}
        self.energy_values: dict[int, c2at.REAL_ARRAY] = {}
        # self.lm_map()

    def lm_map(self) -> None:
        # for every l, (2l + 1) m values in interval [-l, l]
        # We sort self.blocks to a dictionary with keys denoted by "l,m" and values as the corresponding PDOS arrays
        # We also create a map of energies for each spin
        lm_dict: dict[str, list[c2at.REAL_ARRAY]] = {}
        self.lm_dict = self.clear_pdos(self.lm_dict)
        for idx, block in enumerate(self.blocks):
            energy = block[:, 0]
            self.energy_values[idx + 1] = energy
            pdos_values = block[:, 2:]
            num_lm = pdos_values.shape[1]
            l = 0
            m_count = 0
            for i in range(num_lm):
                if m_count < (2 * l + 1):
                    m = -l + m_count
                    lm_key = f"{l},{m}"
                    if lm_key not in lm_dict:
                        lm_dict[lm_key] = []
                    lm_dict[lm_key].append(pdos_values[:, i])
                    m_count += 1
                else:
                    l += 1
                    m_count = 0
                    m = -l + m_count
                    lm_key = f"{l},{m}"
                    if lm_key not in lm_dict:
                        lm_dict[lm_key] = []
                    lm_dict[lm_key].append(pdos_values[:, i])
                    m_count += 1
        self.lm_dict = lm_dict
