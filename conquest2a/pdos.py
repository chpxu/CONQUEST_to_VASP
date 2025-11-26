from pathlib import Path
import re
import numpy as np
dos_file = "/home/chung-ping/Projects/CONQUEST_TO_VASP/tests/Atom00000002DOS_lm.dat"
# pdos_lm = np.genfromtxt(dos_file, delimiter="&")
# print(pdos_lm)

blocks = []
current_block = []

class pdos_processor:
    def __init__(self, dos_file: str | Path) -> None:
        self.dos_file = dos_file
        self.blocks: list[np.ndarray | list[float | int]]= []
        self.read_pdos_file()

    def read_pdos_file(self) -> None:
        with open(self.dos_file, "r", encoding="utf-8") as f:
            current_block = []
            num_spins = 0
            for line in f:
                stripped_line = line.strip()
                if "# Spin" in stripped_line:
                    num_spins += 1
                    continue
                if not stripped_line or stripped_line.startswith("#"):
                    continue
                if stripped_line == "&":
                    if current_block:
                        self.blocks.append(np.array(current_block, dtype=float))
                        current_block = []
                else:
                    current_block.append(line.split())
            if current_block:
                self.blocks.append(np.array(current_block, dtype=float))
            self.num_spins = num_spins
    
class pdos_l_processor(pdos_processor):
    def __init__(self, dos_file: str | Path) -> None:
        self.num_spins: int = 0
        super().__init__(dos_file=dos_file)
        # PDOS file is split into blocks separated by "&" lines
        # The first column of each block is the energy values
        # The second column is the total PDOS, i.e sum of all l, at that energy
        # Subsequent columns are the PDOS values for each l component, sorted in
        # ascending order of l 
        # e.g., l = 0,  l =1,  l = 2, etc.
        # So dict will be of the form {"l": [array(spin1), array(spin2), ...],}
        self.energy_values: dict[int, np.ndarray| list[float]] = {}
        self.l_map()
    def l_map(self) -> None:
        # From column 3, there are only l-contributions, and is sorted by ascending l values, and each row is just ech l-contribution at that energy
        l_dict = {}
        for block in self.blocks:
            energy = block[:, 0]
            pdos_values = block[:, 2:]
            num_l = pdos_values.shape[1]
            for l in range(num_l):
                if str(l) not in l_dict:
                    l_dict[str(l)] = []
                l_dict[str(l)].append(pdos_values[:, l])

class pdos_lm_processor(pdos_processor):
    def __init__(self, dos_file: str | Path) -> None:
        self.num_spins: int = 0
        super().__init__(dos_file=dos_file)
        # PDOS file is split into blocks separated by "&" lines
        # The first column of each block is the energy values
        # The second column is the total PDOS, i.e sum of all l and m, at that energy
        # Subsequent columns are the PDOS values for each lm component, sorted in
        # ascending order of l and m
        # e.g., l = 0,  l =1, m = -1, 0, 1, l = 2, m = -2, -1, 0, 1, 2, etc.
        # So dict will be of the form {"l,m": [array(spin1), array(spin2), ...],}
        self.energy_values: dict[int, np.ndarray| list[float]] = {}
        self.lm_map()
    
    def lm_map(self) -> None:
        # for every l, (2l + 1) m values in interval [-l, l]
        # We sort self.blocks to a dictionary with keys denoted by "l,m" and values as the corresponding PDOS arrays
        # We also create a map of energies for each spin
        lm_dict = {}
        for block in self.blocks:
            energy = block[:, 0]
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


pdos_lm_processor_instance = pdos_lm_processor(dos_file=dos_file)
print(pdos_lm_processor_instance.lm_dict)