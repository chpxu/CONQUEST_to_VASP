"""
This file plots a band structure from a Conquest band structure output file using matplotlib.

Run in root of repo:
    python3 examples/plot_test_pdos.py
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from conquest2a.pdos import pdos_lm_processor
# from conquest2a.band import bst_processor
# Example usage of pdos_l_processor and bst_processor to plot PDOS and band structure
# Adjust the paths below to point to your actual files
data_dir = Path(__file__).parent.parent.resolve() / "tests"
dos_file = data_dir / "Atom00000002DOS_lm.dat"
pdos_lm_processor_instance = pdos_lm_processor(dos_file=dos_file)
print(pdos_lm_processor_instance.lm_dict)

fig = plt.figure(figsize=(8,6))
plt.plot(pdos_lm_processor_instance.energy_values[1], pdos_lm_processor_instance.lm_dict["0,0"][0], label = "l = 0,m = 0, spin1")
plt.plot(pdos_lm_processor_instance.energy_values[1], pdos_lm_processor_instance.lm_dict["0,0"][1], label = "l = 0,m = 0, spin2", color = "red")
plt.savefig(data_dir / "test_pdos.png", dpi=600, bbox_inches="tight")