from dataclasses import dataclass, field
from itertools import islice
import re
from pathlib import Path
import numpy as np
import conquest2a._types as c2at
from conquest2a.conquest import processor_base
from conquest2a.constants import HARTREE_TO_EV


@dataclass
class k_point_blocks:
    k_index: int
    weight: c2at.FLOAT
    eigenvalues: c2at.REAL_ARRAY = field(default_factory=lambda: np.zeros((10)))
    eigenvalues_with_fermi: c2at.REAL_ARRAY = field(default_factory=lambda: np.zeros((10)))
    k_vector: c2at.REAL_ARRAY = field(default_factory=lambda: np.array([0.0, 0.0, 0.0]))


class eigenvalues_processor(processor_base):
    def __init__(self, file: str) -> None:
        self.file = file
        super().__init__(path=file)
        self.resolve_path()
        self.eig_blocks: list[k_point_blocks] = []
        self.num_eigenvalues_per_block: int = 0
        self.fermi_energies: c2at.REAL_ARRAY
        self.open_file()

    def open_file(self) -> None:
        with open(self.abs_input_path, "r", encoding="utf-8") as eigfile:
            line1 = next(eigfile)
            misc = re.findall(self.re_index, line1)
            self.num_eigenvalues_per_block = int(misc[0])
            line2 = next(eigfile)
            self.fermi_energies = np.array(re.findall(self.re_float, line2), dtype=float)
            next(eigfile)
            while True:
                block = list(islice(eigfile, 0, self.num_eigenvalues_per_block + 1))
                if not block:
                    break
                header = block[0].strip().split()
                eigvals = np.array([float(eig.split()[1]) for eig in block[1:]])
                if len(header) == 5:
                    self.eig_blocks.append(
                        k_point_blocks(
                            k_index=int(header[0]),
                            k_vector=np.array(header[1:4], dtype=float),
                            weight=float(header[-1]),
                            eigenvalues=eigvals,
                            eigenvalues_with_fermi=eigvals - self.fermi_energies[0],
                        )
                    )
        eigfile.close()

    def get_vbm_cbm_eigen(
        self,
    ) -> tuple[k_point_blocks, k_point_blocks, c2at.REAL_NUMBER, c2at.REAL_NUMBER]:
        block_with_cbm = self.eig_blocks[0]
        block_with_vbm = self.eig_blocks[0]
        vbm: c2at.REAL_NUMBER = -np.inf
        cbm: c2at.REAL_NUMBER = np.inf
        for block in self.eig_blocks:
            eigenvalues = block.eigenvalues_with_fermi
            valence_eigenvalues = eigenvalues[eigenvalues <= 0.0]
            conduction_eigenvalues = eigenvalues[eigenvalues > 0.0]
            if len(valence_eigenvalues) <= 0:
                raise ValueError(f"No valence eigenvalues found for kpoint {block.k_index}")
            if len(conduction_eigenvalues) <= 0:
                raise ValueError(f"No conduction eigenvalues found for kpoint {block.k_index}")
            if len(valence_eigenvalues) > 0:
                vmax = np.max(valence_eigenvalues)
                if vmax > vbm:
                    vbm = vmax
                    block_with_vbm = block
            if len(conduction_eigenvalues) > 0:
                cmin = np.min(conduction_eigenvalues)
                if cmin < cbm:
                    cbm = cmin
                    block_with_cbm = block
        return block_with_cbm, block_with_vbm, vbm, cbm

    def get_bandgap(self) -> tuple[k_point_blocks, k_point_blocks, c2at.REAL_NUMBER]:
        temp = self.get_vbm_cbm_eigen()
        return temp[0], temp[1], (temp[3] - temp[2]) * HARTREE_TO_EV
