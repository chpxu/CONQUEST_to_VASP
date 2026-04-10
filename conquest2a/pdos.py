from pathlib import Path
from typing import Literal, override
import os
import re
from os.path import abspath
import numpy as np
from conquest2a.conquest import block_processor
import conquest2a._types as c2at
import matplotlib.pyplot as plt

class pdos_processor(block_processor):
    def __init__(self, conquest_rundir: str | Path, lm: Literal["lm", "l", "t"] = "t") -> None:
        # self.dos_file = dos_file
        self.blocks: list[c2at.REAL_ARRAY] = []
        self.all_pdos_files: list[str] = []
        self.conquest_rundir = conquest_rundir
        self.lm = lm
        super().__init__()
        self.fermi_level: float = 0.0
        self.is_shifted_to_fermi: bool = True
        self.num_spins: int = 0
        self.resolve_path()
        self.locate_pdos_files()

    def clear_pdos(
        self, pdos: dict[str, list[c2at.REAL_ARRAY]]
    ) -> dict[str, list[c2at.REAL_ARRAY]]:
        return {key: [] for key in pdos}

    def process_headers(self, line: str) -> None:
        if "# Spin" in line:
            self.num_spins += 1
        if "# Original" in line:
            result = re.findall(self.re_float, line)
            self.fermi_level = float(result[0])
        if not line.startswith("# DOS shifted"):
            self.is_shifted_to_fermi = False

    def process_block(self, line: str) -> None:
        if line == "&":
            if self.current_block:
                self.blocks.append(np.array(self.current_block, dtype=np.float64))
                self.current_block = []
        else:
            self.current_block.append(np.array(line.split()).astype(np.float64))

    def resolve_path(self) -> Path:
        abs_run_path = Path(abspath(self.conquest_rundir))
        if abs_run_path.exists():
            return abs_run_path
        raise FileNotFoundError(f'Conquest directory specified: "{abs_run_path}", does not exist.')

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
        for _, _, files in os.walk(abs_path, topdown=True):
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
    def plot_pdos(self, *args, **kwargs) -> None:
        if self.lm == "t":
            self.energy_values: dict[int, c2at.REAL_ARRAY] = {}
            tdos: c2at.REAL_ARRAY
            ldos: c2at.REAL_ARRAY
            fig = plt.figure(figsize=())
            for idx, block in enumerate(self.blocks):
                energy = block[:, 0]
                self.energy_values[idx + 1] = energy
                tdos = block[:, 1]
                ldos = block[:, 2]
                plt.plot(energy, tdos)
                plt.savefig("DOS.png")
                plt.clf()
                plt.plot(energy, ldos)
                plt.savefig("LDOS.png")

class pdos_l_processor(pdos_processor):
    def __init__(self, conquest_rundir: str | Path) -> None:
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
    def plot_pdos(self) -> None:
        
        super().plot_pdos()


class pdos_lm_processor(pdos_processor):
    def __init__(self, conquest_rundir: str | Path) -> None:
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
        for idx, block in enumerate(self.blocks):
            energy = block[:, 0]
            self.energy_values[idx + 1] = energy
            pdos_values = block[:, 2:]
            num_lm = pdos_values.shape[1]
            l = 0
            m_count = 0
            m = 0
            for i in range(num_lm):
                if m_count >= (2 * l + 1):
                    l += 1
                    m_count = 0
                m = -l + m_count
                lm_key = f"{l},{m}"
                if lm_key not in lm_dict:
                    lm_dict[lm_key] = []
                lm_dict[lm_key].append(pdos_values[:, i])
                m_count += 1
        self.lm_dict = lm_dict
    @override
    def plot_pdos(self, orbitals: list[str]) -> None:
        self.orbital_dict = {
            "0,0": "blue",
            "1,0": "magenta",
            "1,-1": "cyan",
            "1,0": "cyan",
            "1,1": "cyan",
            "2,-2": "red",
            "2,-1": "cyan",
            "2,0": "pink",
            "2,1": "blue",
            "2,2": "black",
        }

        self.label_dict = {
            "0,0": r"$s$",
            "1,0": r"$p_z$",
            "1,-1":r"$p_y$",
            "1,1": r"$p_x$",
            "2,-2": r"$d_{xy}$",
            "2,-1": r"$d_{yz}$",
            "2,0": r"$d_{z^2}$",
            "2,1": r"$d_{xz}$",
            "2,2": r"$d_{x^2 - y^2}$",
        }
        super().plot_pdos()