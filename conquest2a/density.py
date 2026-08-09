from collections.abc import Sequence
from typing import Any
from pathlib import Path
import os
from os.path import abspath, basename
import re
import numpy as np
from ase.atoms import Atoms
from ase.io.cube import read_cube
from ase.units import Bohr
from scipy.ndimage import map_coordinates
import matplotlib as mpl
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable as mal
import matplotlib.pyplot as plt
import scienceplots
from conquest2a._types import INT_ARRAY, REAL_ARRAY
from conquest2a.conquest import processor_base
from conquest2a.constants import MPLGENERIC

mpl.rcParams.update(MPLGENERIC)  # type: ignore
plt.style.use(["science", "no-latex"])

_ELEMENT_COLOURS: dict[str, str] = {
    # Alkali metals
    "Li": "#cc80ff",
    "Na": "#ab5cf2",
    "K": "#8f40d4",
    "Rb": "#702eb0",
    "Cs": "#57178f",
    "Fr": "#420066",
    # Others
    "La": "#70d4ff",
    # Alkaline earth metals
    "Be": "#c2ff00",
    "Mg": "#8aff00",
    "Ca": "#3dff00",
    "Sr": "#00ff00",
    "Ba": "#00c900",
    "Ra": "#007d00",
    # 3d transition metals
    "Sc": "#e6e6e6",
    "Ti": "#bfc2c7",
    "V": "#a6a6ab",
    "Cr": "#8a99c7",
    "Mn": "#9c7ac7",
    "Fe": "#e06633",
    "Co": "#f090a0",
    "Ni": "#50d050",
    "Cu": "#c88033",
    "Zn": "#7d80b0",
    # 4d transition metals
    "Y": "#94ffff",
    "Zr": "#94e0e0",
    "Nb": "#73c2c9",
    "Mo": "#54b5b5",
    "Tc": "#3b9e9e",
    "Ru": "#248f8f",
    "Rh": "#0a7d8c",
    "Pd": "#006985",
    "Ag": "#c0c0c0",
    "Cd": "#ffd98f",
    # 5d transition metals
    "Hf": "#4dc2ff",
    "Ta": "#4da6ff",
    "W": "#2194d6",
    "Re": "#267dab",
    "Os": "#266696",
    "Ir": "#175487",
    "Pt": "#d0d0e0",
    "Au": "#ffd123",
    "Hg": "#b8b8d0",
    # Post-transition metals
    "Al": "#bfa6a6",
    "Ga": "#c28f8f",
    "In": "#a67573",
    "Sn": "#668080",
    "Tl": "#a6544d",
    "Pb": "#575961",
    "Bi": "#9e4fb5",
    # Metalloids
    "B": "#ffb5b5",
    "Si": "#f0c8a0",
    "Ge": "#668f8f",
    "As": "#bd80e3",
    "Sb": "#9e63b5",
    "Te": "#d47a00",
    "Po": "#ab5c00",
    # Halogens
    "F": "#90e050",
    "Cl": "#1ff01f",
    "Br": "#a62929",
    "I": "#940094",
    # Reactive non-metals
    "H": "#ffffff",
    "C": "#404040",
    "N": "#3050f8",
    "O": "#ff0d0d",
    "P": "#ff8000",
    "S": "#ffff30",
    "Se": "#ffa100",
    # Noble gases
    "He": "#d9ffff",
    "Ne": "#b3e3f5",
    "Ar": "#80d1e3",
    "Kr": "#5cb8d1",
    "Xe": "#429eb0",
    "Rn": "#428296",
}


class density(processor_base):
    """Process volumetric density data

    CONQUEST can output any number of ``.cube`` files, e.g.:
        - ``chden_up.cube``/``chden_dn.cube`` for a spin polarised calculation
        - a single cube file, e.g. ``chden.cube`` for unpolarised calculations
        - Band densities

    This class uses ASE's ``.cube`` file processor to get atoms and data.
    It is therefore suitable for generic ``.cube`` file processing.

    An arbitrary number of ``.cube`` files can be supplied via ``paths``. If more than
    one path is given, ``operations`` is a string describing how each subsequent
    file is combined with the running total, applied left to right, one character
    per file after the first, much like VESTA:

    - ``"+"`` add
    - ``"-"`` subtract
    - ``"x"`` multiply
    - ``"/"`` divide

    Example: for ``paths=[a, b, c, d]`` and ``operations="+/-"``, the resulting data is
    ``((a + b) / c) - d``.

    :param hkl: The :math:`hkl` slice of the crystal to plot charge densities in.
    :type hkl: :ref:`INT ARRAY <types>`
    :param offset: The :math:`hkl` direction defines a family of planes.
    Use ``offset`` to select which one (i.e. wherein the unit cell).
    :type offset: ``float``
    :param paths: Path(s) to one or more charge density (``.cube``) files. At least
        one path must be provided.
    :type paths: ``Sequence[str]``
    :param operations: How to combine each file in ``paths`` with the running total of
        those before it, as a string with one of ``"+"``, ``"-"``, ``"x"``, ``"/"`` per
        character. Must have exactly ``len(paths) - 1`` characters. Not required (and
        ignored) when only one path is given, defaults to ``None``.
    :type operations: ``str | None``, optional
    :raises ValueError: If all Miller indices are 0: cannot slice through the origin.
    :raises ValueError: If no paths are provided.
    :raises ValueError: If the length of ``operations`` does not match ``len(paths) - 1``.
    :raises ValueError: If ``operations`` contains a character other than ``"+-x/"``.
    :raises ValueError: If cube files being combined have mismatched grid shapes.
    :raises ZeroDivisionError: If a ``"/"`` operation would divide by a zero-valued voxel.
    """

    def __init__(
        self,
        hkl: tuple[int, int, int],
        offset: float,
        paths: Sequence[str],
        operations: str | None = None,
    ) -> None:
        self._VALID_OPERATION_CHARS: str = "+-x/"
        if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
            raise ValueError("Miller indices (h k l) cannot all be zero.")
        if len(paths) == 0:
            raise ValueError("At least one cube file path must be provided.")

        operations = operations if operations is not None else ""
        if len(paths) > 1 and len(operations) != len(paths) - 1:
            raise ValueError(
                f"Expected {len(paths) - 1} operation character(s) to combine "
                + f"{len(paths)} cube files, got {len(operations)!r}."
            )
        for op in operations:
            if op not in self._VALID_OPERATION_CHARS:
                raise ValueError(
                    f"Unknown operation {op!r}; expected one of {list(self._VALID_OPERATION_CHARS)}."
                )

        self.hkl: INT_ARRAY = np.array(hkl)
        self.offset: float = offset
        self.paths: list[str] = list(paths)
        self.operations: str = operations
        super().__init__(path=self.paths[0])

        #: The (data, atoms) pair loaded from each cube file in ``paths``, in order.
        self.densities: list[tuple[REAL_ARRAY, Atoms]] = [self.load_cube(p) for p in self.paths]
        self.atoms: Atoms = self.densities[0][1]
        self.cell: REAL_ARRAY = self.atoms.get_cell() / Bohr  # ASE cell information

        data: REAL_ARRAY = np.array(self.densities[0][0], copy=True)
        for (other_data, _), op in zip(self.densities[1:], self.operations):
            other_data = np.asarray(other_data)
            if other_data.shape != data.shape:
                raise ValueError(
                    "Cannot combine cube files with mismatched grid shapes: "
                    + f"{data.shape} vs {other_data.shape}."
                )
            if op == "/" and np.any(other_data == 0):
                raise ZeroDivisionError("Division by a cube file that contains zero-valued voxels.")
            if op == "+":
                data = np.add(data, other_data)
            elif op == "-":
                data = np.subtract(data, other_data)
            elif op == "x":
                data = np.multiply(data, other_data)
            else:  # op == "/"
                data = np.divide(data, other_data)
        self.data: REAL_ARRAY = data

    def load_cube(self, filename: str) -> tuple[REAL_ARRAY, Atoms]:
        self.resolve_path(filename=filename)
        with open(filename, "r", encoding="utf-8") as fh:
            cube: dict[str, Atoms] = read_cube(fh)
        fh.close()
        return cube["data"], cube["atoms"]

    def inplane_basis(
        self,
    ) -> tuple[REAL_ARRAY, REAL_ARRAY, REAL_ARRAY]:
        """Return two orthonormal Cartesian vectors :math:`(v_1, v_2)` spanning the :math:`[hkl]` plane
        and the unit plane-normal :math:`\\hat{n}`.

        1. The plane normal in Cartesian space is  :math:`n = ha^* + kb^* + lc^*`
        where :math:`a^*, b^*, c^*` are the reciprocal lattice vectors (rows of inv(cell)^T).
        2. Two fractional-space null vectors of :math:`[hkl]` are found via SVD.
        3. Those are converted to Cartesian, then Gram-Schmidt orthonormalised so the axes are perpendicular in real-space.

        """
        cell: REAL_ARRAY = self.cell
        recip: REAL_ARRAY = (np.linalg.inv(cell).T).astype(float)
        n_cart = self.hkl[0] * recip[0] + self.hkl[1] * recip[1] + self.hkl[2] * recip[2]
        n_hat = n_cart / np.linalg.norm(n_cart)
        _, _, Vt = np.linalg.svd(self.hkl.reshape(1, 3))
        u1_frac, u2_frac = Vt[1], Vt[2]
        u1 = u1_frac @ cell
        u2 = u2_frac @ cell

        v1 = u1 / np.linalg.norm(u1)
        v2 = u2 - np.dot(u2, v1) * v1
        v2 /= np.linalg.norm(v2)

        return v1, v2, n_hat

    def plane_origin(self) -> REAL_ARRAY:
        """Return a Cartesian point lying on the plane  :math:`ha + kb + lc =` ``offset``.

        ``offset`` is a dimensionless fractional intercept (0-1 spans one
        interplanar period). Pick the simplest fractional coordinate
        satisfying the plane equation.

        :returns: Origin of the slice.
        :rtype: :ref:`REAL ARRAY <types>`
        """
        n: REAL_ARRAY = self.hkl / (self.hkl @ self.hkl)  # normal direction in fractional space
        centre: REAL_ARRAY = np.array([0.5, 0.5, 0.5])
        frac: REAL_ARRAY = centre - (centre @ self.hkl) * n + self.offset * n
        return frac @ self.cell

    def inplane_range(self, v1: REAL_ARRAY, v2: REAL_ARRAY) -> tuple[float, float]:
        length_1: float = sum(abs(self.cell[i] @ v1) for i in range(3))
        length_2: float = sum(abs(self.cell[i] @ v2) for i in range(3))
        return length_1, length_2

    def extract_slice(
        self,
        n_points: int = 1000,
        interp_order: int = 5,
    ) -> tuple[Any, REAL_ARRAY, REAL_ARRAY, REAL_ARRAY, REAL_ARRAY, REAL_ARRAY]:
        """Sample the charge density on the :math:`[hkl]` plane at fractional ``offset``.

        For each point on a 2D :math:`(t_1, t_2)` grid centred on the plane origin:

        1. Compute its Cartesian position:
            :math:`p = \\text{origin} + t_1 v_1 + t_2 v_2`
        2. Convert to fractional coordinates: ``s = p @ inv(cell)``
        3. Wrap into the unit cell to enforce periodic boundaries.
        4. Map fractional -> voxel index and interpolate via ``scipy.ndimage.map_coordinates``.

        :param n_points: Number of points to sample on the slice, defaults to 1000
        :type n_points: ``int``, optional
        :param interp_order: The polynomial degree for interpolation, defaults to 5
        :type interp_order: ``int``, optional
        :returns: A tuple containing:

            - **density** (:ref:`REAL ARRAY <types>`) -- Charge density in :math:`e/a_0^3`, shape ``(n_points, n_points)``
            - **v1** (:ref:`REAL ARRAY <types>`) -- First unit vector spanning the :math:`[hkl]` plane
            - **v2** (:ref:`REAL ARRAY <types>`) -- Second unit vector spanning the :math:`[hkl]` plane
            - **t1** (:ref:`REAL ARRAY <types>`) -- Scalar grid along :math:`v_1`
            - **t2** (:ref:`REAL ARRAY <types>`) -- Scalar grid along :math:`v_2`
            - **origin** (``float``) -- Cartesian slice origin in Bohr
        """
        v1, v2, _ = self.inplane_basis()
        origin = self.plane_origin()

        length_1, length_2 = self.inplane_range(v1, v2)
        t1 = np.linspace(-length_1 / 2, length_1 / 2, n_points)
        t2 = np.linspace(-length_2 / 2, length_2 / 2, n_points)
        grid_1, grid_2 = np.meshgrid(t1, t2, indexing="ij")
        pts_cart = origin + grid_1[..., None] * v1 + grid_2[..., None] * v2
        pts_frac = (pts_cart @ np.linalg.inv(self.cell)) % 1.0
        vox = pts_frac * self.data.shape
        data_padded = np.pad(self.data, interp_order + 1, mode="wrap")

        coords = vox.reshape(-1, 3).T + interp_order + 1
        density: REAL_ARRAY = map_coordinates(
            data_padded, coords, order=interp_order, mode="nearest"
        ).reshape(n_points, n_points)

        return density, v1, v2, t1, t2, origin

    def project_atoms(
        self,
        v1: REAL_ARRAY,
        v2: REAL_ARRAY,
        n_hat: REAL_ARRAY,
        origin: REAL_ARRAY,
        thickness: float = 0.5,
    ) -> tuple[REAL_ARRAY, REAL_ARRAY, list[str]]:
        """Project atoms within ``thickness`` of the slice plane onto :math:`(v_1, v_2)` axes.

        :param v1: First unit vector spanning the :math:`[hkl]` plane
        :type v1: :ref:`REAL ARRAY <types>`
        :param v2: Second unit vector spanning the :math:`[hkl]` plane
        :type v2: :ref:`REAL ARRAY <types>`
        :param n_hat: Unit vector defining the slice
        :type n_hat: :ref:`REAL ARRAY <types>`
        :param origin: The origin of the slice
        :type origin: :ref:`REAL ARRAY <types>`
        :param thickness: The Cartesian distance perpendicular to the plane to consider atoms as lying on the slice, defaults to 0.5 Bohr.
        :type thickness: ``float``, optional
        :returns: A tuple containing:

            - **t1_proj** (:ref:`REAL ARRAY <types>`) -- Positions of label along :math:`v_1`
            - **t2_proj** (:ref:`REAL ARRAY <types>`) -- Positions of label along :math:`v_2`
            - **syms** (``list[str]``) -- List of atom labels
        """
        positions = self.atoms.get_positions() / Bohr
        symbols = self.atoms.get_chemical_symbols()
        inv_cell: REAL_ARRAY = np.linalg.inv(self.cell).astype(float)
        disp = positions - origin  # displacement from slice origin
        disp_frac: REAL_ARRAY = disp @ inv_cell
        disp_frac -= np.round(disp_frac)
        disp = disp_frac @ self.cell
        dist_norm = disp @ n_hat  # signed distance to plane

        mask = np.abs(dist_norm) < thickness
        t1_proj = disp[mask] @ v1
        t2_proj = disp[mask] @ v2
        syms: list[Any] = [symbols[i] for i in np.where(mask)[0]]
        return t1_proj, t2_proj, syms


class chden(processor_base):
    """Class which looks for charge density files in a directory.

    :param directory: Path to directory containing band density files
    :type density_instance: ``str``
    """

    def __init__(
        self,
        directory: str,
        hkl: tuple[int, int, int],
        offset: float,
        operations: str,
        charge_stub: str = "chden",
        spin: int | None = None,
    ) -> None:
        self.all_dens_files: list[str] = []
        self.filtered_dens_files: list[str] = []
        self.directory: str = directory
        self.charge_stub: str = charge_stub
        self._charge_regex: re.Pattern[str] = re.compile(rf"{charge_stub}\.cube")
        self._charge_spin_regex: re.Pattern[str] = re.compile(rf"{charge_stub}_(up|dn)\.cube")
        self.spin: int | None = spin
        super().__init__(path=directory)
        self.locate_chden_files()
        self._filter_chden_files()
        self.density: density = density(
            hkl, offset, paths=self.filtered_dens_files, operations=operations
        )

    def locate_chden_files(self) -> list[str]:
        """Gets the paths to all charge density files.

        For spin-unpolarised calculation, CONQUEST will output just ``chden.cube``.
        For spin-polarised calculation, CONQUEST will output ``chden_up.cube`` and ``chden_dn.cube``.

        This output may change if the user has set ``Process.ChargeStub`` in ``Conquest_input`` which defaults to ``chden``. Since reading ``Conquest_input`` directly is unsupported, just assume the user will pass in any stub.

        :raises FileNotFoundError: If ``self.directory`` does not exist.
        :return: List of all paths to band density files.
        :rtype: ``list[str]``
        """
        abs_path: Path = Path(abspath(self.directory))
        if not abs_path.is_dir():
            raise FileNotFoundError(f'Directory specified: "{abs_path}", does not exist.')

        all_file_list: list[str] = []
        file_list: list[str] = []
        if len(self.all_dens_files) > 0:
            # Reset files array, i.e. if changing directory
            self.all_dens_files = []
        for _, _, files in os.walk(abs_path, topdown=True):
            file_list = files
            break
        for filename in file_list:
            if re.match(self._charge_regex, filename) or re.match(
                self._charge_spin_regex, filename
            ):
                all_file_list.append(filename)

        for file in all_file_list:
            full_path = f"{abs_path}/{file}"
            self.resolve_path(filename=full_path)
            self.all_dens_files.append(full_path)
        return self.all_dens_files

    def _filter_chden_files(self) -> list[str]:
        """Method to choose the spin-unpolarised .cube file or the spin-polarised files.

        :return: List of paths to selected charge density files.
        :rtype: ``list[str]``
        """
        pattern: re.Pattern[str]
        if self.spin is None:
            pattern = self._charge_regex
        else:
            pattern = self._charge_spin_regex

        self.filtered_dens_files = [path for path in self.all_dens_files if pattern.search(path)]
        return self.filtered_dens_files


class bandden(processor_base):
    """Class which looks for band density files in a directory.

    :param directory: Path to directory containing band density files
    :type density_instance: ``str``
    """

    def __init__(
        self,
        directory: str,
        hkl: tuple[int, int, int],
        offset: float,
        operations: str,
        band: int,
        spin: int | None = None,
        kpt: int | None = None,
    ) -> None:
        self.all_dens_files: list[str] = []
        self.filtered_dens_files: list[str] = []
        self.bands: list[int] = []
        self._bykpt_regex: re.Pattern[str] = re.compile(
            r"Band([0-9]{6})den_kp([0-9]{3})S(\d)\.cube"
        )
        self._sumkpt_regex: re.Pattern[str] = re.compile(r"Band([0-9]{6})den_totS(\d)\.cube")

        self.directory: str = directory
        super().__init__(path=directory)
        self.locate_banddens_files()
        self._filter_banddens_files(band, spin=spin, kpt=kpt)
        self.density: density = density(
            hkl, offset, paths=self.filtered_dens_files, operations=operations
        )

    def locate_banddens_files(self) -> list[str]:
        """Gets the paths to all band density files, sorted by band number, :math:`k` point then spin.

        :raises FileNotFoundError: If ``self.directory`` does not exist.
        :return: List of all paths to band density files.
        :rtype: ``list[str]``
        """
        abs_path: Path = Path(abspath(self.directory))
        if not abs_path.is_dir():
            raise FileNotFoundError(f'Directory specified: "{abs_path}", does not exist.')

        all_file_list: list[str] = []
        file_list: list[str] = []
        if len(self.all_dens_files) > 0:
            # Reset files array, i.e. if changing directory
            self.all_dens_files = []
            self.bands = []
        for _, _, files in os.walk(abs_path, topdown=True):
            file_list = files
            break
        for filename in file_list:
            if re.match(self._bykpt_regex, filename) or re.match(self._sumkpt_regex, filename):
                all_file_list.append(filename)
                match: re.Match[str] | None = re.search(r"Band([0-9]{6})", filename)
                if match:
                    self.bands.append(int(match.group(1)))

        self.bands = sorted(set(self.bands))
        all_file_list = sorted(all_file_list)
        for file in all_file_list:
            full_path = f"{abs_path}/{file}"
            self.resolve_path(filename=full_path)
            self.all_dens_files.append(full_path)
        return self.all_dens_files

    def _filter_banddens_files(
        self, band: int, spin: int | None = None, kpt: int | None = None
    ) -> list[str]:
        """Method to filter band densities by band number, and optionally spin and k-point.

        :param band: Band index
        :type band: ``int``
        :param spin: Spin index. 1 is usually up and 2 is down. ``None`` fetches all spins associated with a band and k-point
        :type spin: ``int | None``
        :param kpt: k-point number. If ``None``, will look for band density files that have summed over k-points instead.
        :type kpt: ``int | None``
        :raises ValueError: If ``band`` index supplied is not found.
        :raises ValueError: If ``spin`` index supplied is less than 1.
        :raises ValueError: If ``kpt`` index supplied is less than 1.
        :return: List of paths to band density files matching the given filters.
        :rtype: ``list[str]``
        """
        if spin is not None and spin < 1:
            raise ValueError("Spin index cannot be less than 1.")
        if kpt is not None and kpt < 1:
            raise ValueError("k-point index cannot be less than 1.")
        if band not in self.bands:
            raise ValueError("Selected band(s) is not in the directory")

        # Filter by kpt or kpt sum
        pattern: re.Pattern[str] = self._bykpt_regex if kpt is not None else self._sumkpt_regex
        compiled: re.Pattern[str] = re.compile(pattern)

        category_matches: list[tuple[str, re.Match[str]]] = []
        for filepath in self.all_dens_files:
            match: re.Match[str] | None = compiled.search(basename(filepath))
            if match is None:
                continue
            if kpt is not None and int(match.group(2)) != kpt:
                continue
            category_matches.append((filepath, match))

        # band number
        band_matches: list[tuple[str, re.Match[str]]] = [
            (filepath, match) for filepath, match in category_matches if int(match.group(1)) == band
        ]

        # If spin is none, just return the remaining filtered files
        if spin is None:
            self.filtered_dens_files = [filepath for filepath, _ in band_matches]
            return self.filtered_dens_files

        # Otherwise filter for specific files with matching spin from remaining files
        spin_group: int = 3 if kpt is not None else 2
        self.filtered_dens_files = [
            filepath for filepath, match in band_matches if int(match.group(spin_group)) == spin
        ]
        return self.filtered_dens_files


class plot_densities:
    """Helper class to slice, analyse, and plot volumetric density data.

    Takes a :class:`density` instance directly, or with
    :class:`chden` or :class:`bandden` that exposes a ``.density`` attribute

    :param source: The density data to plot.
    :type source: :class:`density` | :class:`chden` | :class:`bandden`
    :param show_atoms: Whether to show atom labels on the plot, defaults to False.
    :type show_atoms: ``bool``, optional
    :param extension: File extension used for auto-generated filenames (leading
        "." optional), defaults to "png".
    :type extension: ``str``, optional
    :param cbar_label: Label for the colour bar. If not given, a sensible default
        is chosen based on the type of ``source`` -- charge density and band
        density are conventionally expressed in different units.
    :type cbar_label: ``str | None``, optional
    """

    _DEFAULT_CBAR_LABELS: dict[type, str] = {
        chden: r"$\rho$  [$e\,a_0^{-3}$]",
        bandden: r"$|\psi_{nk}|^2$  [$a_0^{-3}$]",
        density: r"$\rho$  [arb. units]",
    }

    def __init__(
        self,
        source: "density | chden | bandden",
        show_atoms: bool = False,
        extension: str = "png",
        cbar_label: str | None = None,
    ) -> None:
        self.density: density = source.density if isinstance(source, (chden, bandden)) else source
        self.show_atoms: bool = show_atoms
        self.extension: str = extension.lstrip(".")
        self.cbar_label: str = cbar_label or self._DEFAULT_CBAR_LABELS.get(
            type(source), self._DEFAULT_CBAR_LABELS[density]
        )

    def _miller_str(self, idx: int) -> str:
        if idx >= 0:
            return str(idx)
        return rf"$\overline{{{abs(idx)}}}$"

    def _vec_str(self, v: REAL_ARRAY) -> str:
        parts: list[str] = [self._miller_str(x) for x in v]
        return f"Position along [{','.join(parts)}]"

    def _default_filename(self) -> str:
        hkl = self.density.hkl
        filename = f"{hkl[0]}{hkl[1]}{hkl[2]}_{self.density.offset:.3f}"
        if self.density.operations:
            filename += f"_{self.density.operations}"
        return f"{filename}.{self.extension}"

    @staticmethod
    def _orient_for_plot(
        grid: REAL_ARRAY, t1: REAL_ARRAY, t2: REAL_ARRAY, v1: REAL_ARRAY, v2: REAL_ARRAY
    ) -> tuple[REAL_ARRAY, REAL_ARRAY, REAL_ARRAY, REAL_ARRAY, REAL_ARRAY, bool]:
        """Private method to rotate the plot so the horizontal direction is longest"""
        if (t2[-1] - t2[0]) > (t1[-1] - t1[0]):
            return grid.T, t2, t1, v2, v1, True
        return grid, t1, t2, v1, v2, False

    def _plot_atoms(
        self,
        ax: Any,
        atom_data: tuple[REAL_ARRAY, REAL_ARRAY, list[str]],
        transposed: bool,
    ) -> None:
        t1s, t2s, syms = atom_data
        for t1a, t2a, sym in zip(t1s, t2s, syms):
            if transposed:
                t1a, t2a = t2a, t1a
            color = _ELEMENT_COLOURS.get(sym, "#00ff80")
            ax.scatter(t1a, t2a, s=160, color=color, edgecolors="white", linewidths=0.8, zorder=5)
            ax.text(
                t1a,
                t2a,
                sym,
                ha="center",
                va="center",
                fontsize=5,
                color="black",
                zorder=6,
                fontweight="bold",
            )

    def plot_slice(
        self,
        density_grid: REAL_ARRAY,
        t1: REAL_ARRAY,
        t2: REAL_ARRAY,
        v1: REAL_ARRAY,
        v2: REAL_ARRAY,
        atom_data: (
            tuple[REAL_ARRAY, REAL_ARRAY, list[str]] | None
        ) = None,  # (t1s, t2s, symbols) tuple or None
        log_scale: bool = False,
        vmin: float | None = 0.0,
        vmax: float | None = None,
        output: str | None = None,
        figsize: tuple[float, float] | None = None,
        savefig_kwargs: dict[str, Any] | None = None,
        **imshow_kwargs: Any,
    ) -> None:
        """Create the plot and save it to disk.

        :param density_grid: Sliced density data, shape ``(n_points, n_points)``.
        :type density_grid: :ref:`REAL ARRAY <types>`
        :param t1: Scalar grid along :math:`v_1`
        :type t1: :ref:`REAL ARRAY <types>`
        :param t2: Scalar grid along :math:`v_2`
        :type t2:  :ref:`REAL ARRAY <types>`
        :param v1: First unit vector spanning the :math:`[hkl]` plane
        :type v1: :ref:`REAL ARRAY <types>`
        :param v2: Second unit vector spanning the :math:`[hkl]` plane
        :type v2: :ref:`REAL ARRAY <types>`
        :param atom_data: Data for atom labels, see :func:`density.project_atoms`, defaults to None
        :type atom_data: ``tuple[REAL_ARRAY, REAL_ARRAY, list[str]] | None``, optional
        :param log_scale: Whether to plot the density on a base-10 logarithmic scale -
            useful for revealing details without colour clipping, defaults to False
        :type log_scale: ``bool``, optional
        :param vmin: Minimum value to set the colour scale at (ignored if ``log_scale``), defaults to 0.0
        :type vmin: ``float | None``, optional
        :param vmax: Maximum value to set the colour scale at (ignored if ``log_scale``), defaults to None
        :type vmax: ``float | None``, optional
        :param output: Filename to save as. Will save with a useful name if not provided, defaults to None
        :type output: ``str | None``, optional
        :param figsize: Override the auto-computed figure size, defaults to None
        :type figsize: ``tuple[float, float] | None``, optional
        :param savefig_kwargs: Extra keyword arguments forwarded to
            :func:`matplotlib.pyplot.savefig` (e.g. ``dpi``, ``transparent``,
            ``bbox_inches``), defaults to None
        :type savefig_kwargs: ``dict[str, Any] | None``, optional
        :param imshow_kwargs: Any further keyword arguments (e.g. ``cmap``,
            ``interpolation``, ``alpha``, ``aspect``) are forwarded to
            :func:`matplotlib.pyplot.imshow`, overriding the defaults.
        """
        density_grid, t1, t2, v1, v2, transposed = self._orient_for_plot(
            density_grid, t1, t2, v1, v2
        )
        l1 = t1[-1] - t1[0]
        l2 = t2[-1] - t2[0]

        if figsize is None:
            fig_w = 5.0
            figsize = (fig_w, fig_w * (l2 / l1) + 0.5)
        fig, ax = plt.subplots(figsize=figsize)
        divider = mal(ax)

        imshow_args: dict[str, Any] = {
            "origin": "lower",
            "extent": (t1[0], t1[-1], t2[0], t2[-1]),
            "cmap": "viridis",
            "aspect": "equal",
            "interpolation": "lanczos",
        }
        # vmax, vmin and Norm interact separately
        if log_scale:
            imshow_args["norm"] = colors.LogNorm()
        else:
            imshow_args["vmin"] = 0.0 if vmin is None else float(vmin)
            imshow_args["vmax"] = float(np.max(density_grid)) if vmax is None else float(vmax)
        imshow_args.update(imshow_kwargs)

        im = ax.imshow(density_grid, **imshow_args)
        ax.tick_params(direction="out", which="both")
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(im, cax=cax, fraction=0.04, pad=0.1)
        cbar.set_label(
            (r"$\log_{10}$ " + self.cbar_label) if log_scale else self.cbar_label,
            fontsize=10,
        )

        if atom_data is not None:
            self._plot_atoms(ax, atom_data, transposed)

        ax.set_xlabel(f"{self._vec_str(v1)}" + r"$[a_0]$", fontsize=8)
        ax.set_ylabel(f"{self._vec_str(v2)}" + r"$[a_0]$", fontsize=8)

        fig.tight_layout()
        output = output or self._default_filename()
        plt.savefig(output, **(savefig_kwargs or {}))
        plt.close(fig)
        print(f"Saved: {output}")

    def run(
        self,
        filename: str | None = None,
        thickness: float = 0.5,
        log_scale: bool = False,
        vmin: float | None = 0.0,
        vmax: float | None = None,
        figsize: tuple[float, float] | None = None,
        savefig_kwargs: dict[str, Any] | None = None,
        **imshow_kwargs: Any,
    ) -> None:
        """Runs the full sequence of steps:

        1. Extracts a slice from the underlying :class:`density`
        2. Optionally projects nearby atoms onto the slice
        3. Plots and saves the result

        :param filename: Filename to save as. Will save with a useful name if not provided, defaults to None
        :type filename: ``str | None``, optional
        :param thickness: The Cartesian distance perpendicular to the plane to consider atoms as lying on the slice, defaults to 0.5 Bohr.
        :type thickness: ``float``, optional
        :param log_scale: Whether to plot the density on a base-10 logarithmic scale, defaults to False
        :type log_scale: ``bool``, optional
        :param vmin: Minimum value to set the colour scale at, defaults to 0.0
        :type vmin: ``float | None``, optional
        :param vmax: Maximum value to set the colour scale at, defaults to None
        :type vmax: ``float | None``, optional
        :param figsize: Override the auto-computed figure size, defaults to None
        :type figsize: ``tuple[float, float] | None``, optional
        :param savefig_kwargs: Extra keyword arguments forwarded to :func:`matplotlib.pyplot.savefig`, defaults to None
        :type savefig_kwargs: ``dict[str, Any] | None``, optional
        :param imshow_kwargs: Any further keyword arguments (e.g. ``cmap``,
            ``interpolation``) are forwarded to :func:`matplotlib.pyplot.imshow`.
        """
        density_grid, v1, v2, t1, t2, origin = self.density.extract_slice()
        atom_data: tuple[REAL_ARRAY, REAL_ARRAY, list[str]] | None = None
        if self.show_atoms:
            _, _, n_hat = self.density.inplane_basis()
            atom_data = self.density.project_atoms(v1, v2, n_hat, origin, thickness=thickness)
            print(f"  Atoms within {thickness} (Bohr) of plane: {len(atom_data[0])}")
        self.plot_slice(
            density_grid,
            t1,
            t2,
            v1,
            v2,
            atom_data,
            log_scale=log_scale,
            vmin=vmin,
            vmax=vmax,
            output=filename,
            figsize=figsize,
            savefig_kwargs=savefig_kwargs,
            **imshow_kwargs,
        )
