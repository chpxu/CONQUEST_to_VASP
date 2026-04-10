import importlib.resources
from ase.units import Bohr, Hartree
import matplotlib.pyplot as plt
import matplotlib as mpl
import scienceplots

plt.style.use(["science", "no-latex"])

BOHR_TO_ANGSTROM = Bohr
ANGSTROM_TO_BOHR = 1 / BOHR_TO_ANGSTROM
BOHR_TO_ANGSTROM_VOLUME = (BOHR_TO_ANGSTROM) ** 3
HARTREE_TO_EV = Hartree
EV_TO_HARTREE = 1 / HARTREE_TO_EV
LIBRARY = importlib.resources.files("conquest2a")


# Some custom MPL params to use with SciencePlots
MPLGENERIC = {
    # Muted grey gridlines
    "axes.axisbelow": True,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "figure.figsize": (3.375 * 2, 2.8),
    # Flat / frameless legend
    "legend.frameon": True,
    "legend.framealpha": 0.5,
    "legend.edgecolor": "black",
    "legend.fancybox": False,
    "legend.loc": "best",
    "legend.handlelength": 1.5,
    "legend.handletextpad": 0.5,
    "legend.columnspacing": 1.0,
}

MPLGRID = {
    "axes.grid": True,
    "grid.color": "#CCCCCC",
    "grid.linewidth": 0.5,
    "grid.linestyle": "--",
    "grid.alpha": 0.7,
}

mpl.rcParams.update(MPLGENERIC)
