from pathlib import Path
from typing import Literal
import os
import re
from os.path import abspath
import numpy as np
from conquest2a.conquest import block_processor
import conquest2a._types as c2at
import matplotlib.pyplot as plt


class pdos_processor(block_processor):
    """Initialise generic PDOS processor class.
    
    CONQUEST can produce a ``DOS.dat`` containing the total DOS and the local DOS, which ``lm="t"`` will process. To process :math:`l` and :math:`lm`-resolved PDOS files, initialise an instance of :class:`pdos_l_processor` and :class:`pdos_lm_processor` respectively.

    :param conquest_rundir: String or Path to the directory containing the PDOS files generated from CONQUEST's PostProcessing tool.
    :type conquest_rundir: ``string | Path``
    :param lm: Determines the file-processing mode. Defaults to ``"t"``. 
    
    :type lm: ``Literal["lm", "l", "t"]``, optional
    """
    def __init__(self, conquest_rundir: str | Path, lm: Literal["lm", "l", "t"] = "t") -> None:
        # self.dos_file = dos_file
        self.lm = lm
        self.blocks: list[c2at.REAL_ARRAY] = []
        self.filename_regex = rf"Atom[0-9]{{7}}DOS\_{self.lm}\.dat"
        self.all_pdos_files: list[str] = []
        self.pdos_atoms: list[int] = []
        self.conquest_rundir = conquest_rundir
        super().__init__()
        self.fermi_level: float = 0.0
        self.is_shifted_to_fermi: bool = True
        self.num_spins: int = 0
        self.resolve_path()
        self.locate_pdos_files()

    def process_headers(self, line: str) -> None:
        """Processes lines in PDOS output files starting with #

        :param line: line
        :type line: ``str``
        """
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
        """Checks whether the directory containing files exists.

        :raises FileNotFoundError: Raises error if directory does not exist
        :return: Path to the directory
        :rtype: ``Path``
        """
        abs_run_path = Path(abspath(self.conquest_rundir))
        if abs_run_path.exists():
            return abs_run_path
        raise FileNotFoundError(f'Conquest directory specified: "{abs_run_path}", does not exist.')

    def locate_pdos_files(self) -> list[str]:
        """Gets the paths to all PDOS files of the right type and stores it in a list.

        :return: List of all paths to PDOS files.
        :rtype: ``list[str]``
        """
        if len(self.all_pdos_files) > 0:
            # Reset pdos files array, i.e. if changing directory
            self.all_pdos_files = []
        if self.lm == "t":
            self.all_pdos_files = ["DOS.dat"]
            return self.all_pdos_files
        abs_path = self.resolve_path()
        file_list: list[str] = []
        for _, _, files in os.walk(abs_path, topdown=True):
            file_list = files
            break
        pdos_file_list: list[str] = []
        # The search is done because a directory can contain both lm, l resolved pDOS files
        for filename in file_list:
            match = re.search(f"([0-9]{{7}})", filename)
            if match:
                self.pdos_atoms.append(int(match.group(1)))
            res = re.match(self.filename_regex, filename)
            if res:
                pdos_file_list.append(filename)
        pdos_file_list = sorted(pdos_file_list)
        for file in pdos_file_list:
            self.all_pdos_files.append(f"{abspath(self.conquest_rundir)}/{file}")
        return self.all_pdos_files

    def get_pdos(self, atom: int) -> None:
        """Reads in the PDOS data of a file corresponding to an atom.

        CONQUEST outputs pdos filenames pf the form ``AtomNNNNNNNDOS_lm.dat`` or ``AtomNNNNNNNDOS_l.dat``
        where NNNNNNN is a zero-padded atom number (as ordered in the coordinates file). This method is modified by the children :class:`pdos_l_processor` and :class:`pdos_lm_processor` for their purposes.

        :param atom: The atom to find and read in the PDOS for.
        :type atom: ``int``
        :raises ValueError: If the atom chosen does not have a PDOS file associated to it.
        """
        if atom not in self.pdos_atoms:
            raise ValueError("Chosen atom for pdos was not in the atom list")
        id = f"{atom:07d}"
        # self.all_pdos_files will be absolute oaths, use searchc
        for filename in self.all_pdos_files:
            match = re.search(rf"Atom{id}DOS_{self.lm}\.dat", filename)
            if match:
                self.read_file(filename)
                return

    def plot_pdos(self, *args, **kwargs) -> None:  # type: ignore
        """Method which plots the PDOS and LDOS inside `DOS.dat`.
        """
        if self.lm == "t":
            self.energy_values: dict[int, c2at.REAL_ARRAY] = {}
            for filename in self.all_pdos_files:
                self.read_file(filename)
            tdos: c2at.REAL_ARRAY
            ldos: c2at.REAL_ARRAY
            fig = plt.figure(figsize=(3,2))
            # print(self.blocks)
            energy = self.blocks[0][:, 0]
            self.energy_values[0] = energy
            tdos_up = self.blocks[0][:, 1]
            tdos_dn = self.blocks[1][:, 1]
            ldos_up = self.blocks[0][:, 2]
            ldos_dn = self.blocks[1][:, 2]

            plt.plot(energy, tdos_up, color="red", label="Spin up")
            plt.plot(energy, -1 * tdos_dn, color="blue", label="Spin down")
            plt.savefig("DOS.png")

            plt.plot(energy, ldos_up, color="red", label="Spin up")
            plt.plot(energy, -1 * ldos_dn, color="blue", label="Spin down")
            plt.savefig("LDOS.png")


class pdos_l_processor(pdos_processor):
    def __init__(self, conquest_rundir: str | Path) -> None:
        """Class to process and plot :math:`l`-resolved PDOS.
            
        * PDOS file is split into blocks separated by "&" lines. The first block is the spin-up and second is spin-down.
        * Column 1 records the energy in electronvolts
        * Column 2 records the sum over all :math:`l`-PDOS at that energy
        * From column 3 onwards, records specific :math:`l`-contributions, and columns are sorted by ascending :math:`l` values.
        
        :param conquest_rundir:  String or Path to the directory containing the PDOS files generated from CONQUEST's `PostProcess` binary.
        :type conquest_rundir: ``str | Path``
        """
        super().__init__(conquest_rundir=conquest_rundir, lm="l")
        self.l_dict: dict[str, list[c2at.REAL_ARRAY]] = {}
       
        # e.g., l = 0,  l =1,  l = 2, etc.
        # So dict will be of the form {"l": [array(spin1), array(spin2), ...],}
        self.energy_values: dict[int, c2at.REAL_ARRAY] = {}
        self.color_dict = {
            "0": "blue",
            "1": "cyan",
            "2": "magenta",
            "3": "red",
        }

        self.label_dict = {
            "0": r"$l = 0$",
            "1": r"$l = 1$",
            "2": r"$l = 2$",
            "3": r"$l = 3$",
        }

    def l_map(self) -> None:
        """Reads and stores the columns of an :math:`l`-resolved PDOS file.
        """
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

    def get_pdos(self, atom: int) -> None:
        """Method which clears the data, calls :func:`~get_pdos` on the ``atom``, and :func:`~l_map` in one go.

        :param atom: The atom to get the PDOS for.
        :type atom: ``int``
        :raises ValueError: If the chosen atom to plot does not have a pDOS file, the method will abort without doing anything.
        """
        if atom not in self.pdos_atoms:
            raise ValueError(f"Some chosen atoms for pdos plotting was not in the atom list")
        self.l_dict = {}
        super().get_pdos(atom)
        self.l_map()

    def plot_pdos(
        self,
        atomno: int,
        ang_mom: int,
        x1: float | None,
        x2: float | None,
        y1: float | None,
        y2: float | None,
        filename: str,
    ) -> None:
        """Plots the :math:`l` pDOS for an atom.

        :param atomno: The atom number as defined in the coordinates.
        :type atomno: ``int``
        :param ang_mom: The angular momentum value to plot.
        :type ang_mom: ``int``
        :param x1: Lower energy limit. If ``None``, defaults to the lowest energy in the data.
        :type x1: ``float | None``
        :param x2: Upper energy limit. If ``None``, defaults to the highest energy in the data.
        :type x2: ``float | None``
        :param y1: Lower y-limit
        :type y1: ``float | None``
        :param y2: Upper y-limit
        :type y2: ``float | None``
        :param filename: Filename of plot.
        :type filename: ``str``
        """
        x_label = r"$E - E_F~[\text{eV}]$" if self.is_shifted_to_fermi else r"$E~[\text{eV}]$"
        y_label = r"$\text{DOS} [\text{states/eV}]$"
        
        fig = plt.figure()

        self.get_pdos(atomno)
        key = str(ang_mom)
        plt.plot(
            self.energy_values[0],
            self.l_dict[key][0],
            label=self.label_dict[key],
            color=self.color_dict[key],
        )
        plt.plot(
            self.energy_values[1], -self.l_dict[key][1], linestyle="--", color=self.color_dict[key]
        )
        plt.legend()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.ylim(y1, y2)
        plt.xlim(x1, x2)
        plt.savefig(filename)
        plt.close()


class pdos_lm_processor(pdos_processor):
    def __init__(self, conquest_rundir: str | Path) -> None:
        """Class to process and plot :math:`lm`-resolved PDOS.
            
        * PDOS file is split into blocks separated by "&" lines. The first block is the spin-up and second is spin-down.
        * Column 1 records the energy in electronvolts
        * Column 2 records the sum over all :math:`lm`-PDOS at that energy
        * From column 3 onwards, records specific contributions, and columns are sorted by ascending :math:`l` values, and ascending :math:`m` values.
            *  E.g., (l = 0),  (l =1), m = -1, 0, 1, (l = 2), m = -2, -1, 0, 1, 2, etc.

        :param conquest_rundir:  String or Path to the directory containing the PDOS files generated from CONQUEST's `PostProcess` binary.
        :type conquest_rundir: ``str | Path``
        """
        super().__init__(conquest_rundir=conquest_rundir, lm="lm")
        self.lm_dict: dict[str, list[c2at.REAL_ARRAY]] = {}
        # PDOS file is split into blocks separated by "&" lines
        # The first column of each block is the energy values
        # The second column is the total PDOS, i.e sum of all l and m, at that energy
        # Subsequent columns are the PDOS values for each lm component, sorted in
        # ascending order of l and m
       
        # So dict will be of the form {"l,m": [array(spin1), array(spin2), ...],}
        self.energy_values: dict[int, c2at.REAL_ARRAY] = {}

        self.color_dict = {
            "0,0": "blue",
            "1,-1": "cyan",
            "1,0": "magenta",
            "1,1": "red",
            "2,-2": "red",
            "2,-1": "cyan",
            "2,0": "pink",
            "2,1": "blue",
            "2,2": "black",
        }

        self.label_dict = {
            "0,0": r"$s$",
            "1,-1": r"$p_y$",
            "1,0": r"$p_z$",
            "1,1": r"$p_x$",
            "2,-2": r"$d_{xy}$",
            "2,-1": r"$d_{yz}$",
            "2,0": r"$d_{z^2}$",
            "2,1": r"$d_{xz}$",
            "2,2": r"$d_{x^2 - y^2}$",
        }

    def lm_map(self) -> None:
        """Reads and stores the columns of an :math:`lm`-resolved PDOS file.

        * For every :math:`l`, there are :math:`2l + 1` columns of PDOS. These become the keys of a dictionary and will be accessed in the form ``["l,m"]``.
        * The values of this dictionary will be a list of NumPy arrays, ordered first by spin-up and then spin-down.
        """
       
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

    def get_pdos(self, atom: int) -> None:
        """Method which clears the data, calls :func:`~get_pdos` on the ``atom``, and :func:`~lm_map` in one go.

        :param atom: The atom to get the PDOS for.
        :type atom: ``int``
        :raises ValueError: If the chosen atom to plot does not have a pDOS file, the method will abort without doing anything.
        """
        if atom not in self.pdos_atoms:
            raise ValueError("Chosen atom for pdos was not in the atom list")
        self.lm_dict = {}
        super().get_pdos(atom)
        self.lm_map()

    def plot_pdos(
        self,
        atomnos: list[int],
        orbitals: list[str],
        x1: float | None,
        x2: float | None,
        y1: float | None,
        y2: float | None,
        filename: str,
    ) -> None:
        """Plots the :math:`l` pDOS for an atom.

        :param atomnos: The atom numbers as defined in the coordinates, to plot.
        :type atomnos: ``int``
        :param orbitals: The orbitals to plot. These are a list of the keys in the dictionary. 
        :type orbitals: ``list[int]``
        :param x1: Lower energy limit. If ``None``, defaults to the lowest energy in the data.
        :type x1: ``float | None``
        :param x2: Upper energy limit. If ``None``, defaults to the highest energy in the data.
        :type x2: ``float | None``
        :param y1: Lower y-limit
        :type y1: ``float | None``
        :param y2: Upper y-limit
        :type y2: ``float | None``
        :param filename: Filename of plot.
        :type filename: ``str``
        """
        x_label = r"$E - E_F~[\text{eV}]$" if self.is_shifted_to_fermi else r"$E~[\text{eV}]$"
        y_label = r"$\text{DOS} [\text{states/eV}]$"
        super().plot_pdos()
        if not set(atomnos).issubset(self.pdos_atoms):
            raise ValueError(f"Some chosen atoms for pdos plotting was not in the atom list")
        # energy_mask = np.ma.masked_inside(self.energy_values[1], x1, x2).mask #type: ignore
        # x_energy = self.energy_values[1][energy_mask]
        # Plot the same orbitals from each atom on the same plot
        fig = plt.figure()

        for atom in atomnos:
            self.get_pdos(atom)
            for key in orbitals:
                plt.plot(
                    self.energy_values[0],
                    self.lm_dict[key][0],
                    label=self.label_dict[key],
                    color=self.color_dict[key],
                )
                plt.plot(
                    self.energy_values[1],
                    -self.lm_dict[key][1],
                    linestyle="--",
                    color=self.color_dict[key],
                )
        plt.legend()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.ylim(y1, y2)
        plt.xlim(x1, x2)
        plt.savefig(filename)
        plt.close()
