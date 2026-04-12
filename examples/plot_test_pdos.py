"""
This file plots a band structure from a Conquest band structure output file using matplotlib.

Run in root of repo:
    python3 examples/plot_test_pdos.py
"""

from pathlib import Path
from conquest2a.pdos import pdos_lm_processor
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots  # pip3 install SciencePlots

# Load SciencePlots configuration

plt.style.use(["science", "no-latex"])
# Override some of their settings because I don't like their grid params
# I also want a solid legend
mpl.rcParams.update(
    {
        # Muted grey gridlines
        "axes.grid": True,
        "grid.color": "#CCCCCC",
        "grid.linewidth": 0.5,
        "grid.linestyle": "--",
        "grid.alpha": 0.7,
        "axes.axisbelow": True,  # grid behind data
        # Flat / frameless legend
        "legend.frameon": True,
        "legend.framealpha": 0.5,
        "legend.edgecolor": "black",
        "legend.fancybox": False,
        "legend.loc": "best",
        "legend.handlelength": 1.5,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.0,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "figure.figsize": (3.5, 2.8),
    }
)

# Some labels for colouring
color_dict = {
    "0,0": "blue",
    "1,0": "magenta",
    "1,-1": "brown",
    "1,1": "cyan",
    "2,-2": "red",
    "2,-1": "cyan",
    "2,0": "pink",
    "2,1": "blue",
    "2,2": "black",
}

orbital_dict = {
    "0,0": "blue",
    "1,0": "magenta",
    "1,-1": "red",
    "1,1": "cyan",
    "2,-2": "red",
    "2,-1": "cyan",
    "2,0": "pink",
    "2,1": "blue",
    "2,2": "black",
}

label_dict = {
    "0,0": r"$s$",
    "1,0": r"$p_z$",
    "1,-1": r"$p_y$",
    "1,1": r"$p_x$",
    "2,-2": r"$d_{xy}$",
    "2,-1": r"$d_{yz}$",
    "2,0": r"$d_{z^2}$",
    "2,1": r"$d_{xz}$",
    "2,2": r"$d_{x^2 - y^2}$",
}
# from conquest2a.band import bst_processor
# Example usage of pdos_l_processor and bst_processor to plot PDOS and band structure
# Adjust the paths below to point to your actual files
data_dir = Path(__file__).parent.parent.resolve() / "tests"
pdos_lm_processor_instance = pdos_lm_processor(conquest_rundir=data_dir)
pdos_lm_processor_instance.lm_map()
plt.plot(
    pdos_lm_processor_instance.energy_values[1],
    pdos_lm_processor_instance.lm_dict["0,0"][0],
    label="l = 0,m = 0, spin1",
)
plt.plot(
    pdos_lm_processor_instance.energy_values[1],
    pdos_lm_processor_instance.lm_dict["0,0"][1],
    label="l = 0,m = 0, spin2",
    color="red",
)
plt.savefig(data_dir / "test_pdos.png", dpi=600, bbox_inches="tight")
