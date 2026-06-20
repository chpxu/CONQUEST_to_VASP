from typing import Any, Literal
import numpy as np
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

mpl.rcParams.update(MPLGENERIC)
plt.style.use(["science", "no-latex"])

_ELEMENT_COLOURS = {
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


class chden(processor_base):
    """Process charge density data

    CONQUEST supplies up to two chden .cube files:
        - chden_up/dn.cube for a spin polarised calculation
        - a single cube file, e.g. chden.cube for unpolarised calculations

    This class uses ASE's .cube file processor to get atoms and data.
    It combines chden up and dn if both arguments are supplied, to make a total charge density
    It also does up - dn to see the spin difference, if both arguments are supplied

    :param hkl: The :math:`hkl` slice of the crystal to plot charge densities in.
    :type hkl: :ref:`INT ARRAY <types>`
    :param offset: The :math:`hkl` direction defines a family of planes. Use ``offset`` to select which one (i.e. where in the unit cell).
    :type offset: ``float``
    :param ch1: Path to a charge density (``.cube``) file.
    :type ch1: ``str``
    :param ch2: Path to another charge density (``.cube``) file, defaults to ``None``.
    :type ch2: ``str | None``, optional
    :param mode: If two charge densities are supplied, whether to sum the data or find the difference, defaults to ``None``
    :type mode: ``Literal["sum", "diff"] | None``, optional
    :raises ValueError: Cannot provide a mode if only one file is specified.
    :raises ValueError: If all Miller indices are 0: cannot slice through the origin
    """

    def __init__(
        self,
        hkl: INT_ARRAY,
        offset: float,
        ch1: str,
        ch2: str | None = None,
        mode: Literal["sum", "diff"] | None = None,
    ) -> None:
        self.hkl = hkl
        self.offset = offset
        self.mode = mode
        super().__init__(path=ch1)
        self.dens1 = self.load_cube(ch1)
        self.data = self.dens1[0]
        self.dens2 = None
        self.atoms = self.dens1[1]
        self.cell = self.atoms.get_cell() / Bohr  # ASE cell information
        if self.mode is not None and ch2 is None:
            raise ValueError("Cannot have sum/diff mode if a second file is not provided!")
        if ch2 is not None:
            self.dens2 = self.load_cube(ch2)
            if self.mode == "sum":
                self.data += self.dens2[0]  # add data arrays together
            elif self.mode == "diff":
                self.data -= self.dens2[0]  # subtract data arrays together

        if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
            raise ValueError("Miller indices (h k l) cannot all be zero.")

    def load_cube(self, filename: str) -> tuple[Any, Any]:
        self.resolve_path(filename=filename)
        with open(filename, "r", encoding="utf-8") as fh:
            cube = read_cube(fh)  # type: ignore
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
        cell = self.cell
        recip = np.linalg.inv(cell).T
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

    def plane_origin(self) -> Any:
        """Return a Cartesian point lying on the plane  :math:`ha + kb + lc =` ``offset``.

        ``offset`` is a dimensionless fractional intercept (0-1 spans one
        interplanar period). Pick the simplest fractional coordinate
        satisfying the plane equation.

        :returns: Origin of the slice.
        :rtype: :ref:`REAL ARRAY <types>`
        """
        n = self.hkl / (self.hkl @ self.hkl)  # normal direction in fractional space
        centre = np.array([0.5, 0.5, 0.5])
        frac = centre - (centre @ self.hkl) * n + self.offset * n
        return frac @ self.cell

    def inplane_range(self, v1: REAL_ARRAY, v2: REAL_ARRAY) -> tuple[float, float]:
        length_1 = sum(abs(self.cell[i] @ v1) for i in range(3))
        length_2 = sum(abs(self.cell[i] @ v2) for i in range(3))
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
        density = map_coordinates(data_padded, coords, order=interp_order, mode="nearest").reshape(
            n_points, n_points
        )

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
        inv_cell = np.linalg.inv(self.cell)
        disp = positions - origin  # displacement from slice origin
        disp_frac = disp @ inv_cell
        disp_frac -= np.round(disp_frac)
        disp = disp_frac @ self.cell
        dist_norm = disp @ n_hat  # signed distance to plane

        mask = np.abs(dist_norm) < thickness
        t1_proj = disp[mask] @ v1
        t2_proj = disp[mask] @ v2
        syms = [symbols[i] for i in np.where(mask)[0]]
        return t1_proj, t2_proj, syms


class chden_plot:
    """Helper class to analyse and plot charge densities.

    :param chden_instance: :class:`chden` instance to use
    :type chden_instance: :class:`chden`
    :param show_atoms: Whether to show atom labels on the charge density plot, defaults to False
    :type show_atoms: ``bool``, optional
    :param extension: File extension, defaults to "png"
    :type extension: ``str``, optional
    """

    def __init__(
        self, chden_instance: chden, show_atoms: bool = False, extension: str = "png"
    ) -> None:
        self.chden = chden_instance
        self.show_atoms = show_atoms
        self.extension = extension

    def _miller_str(self, idx: int) -> str:
        if idx >= 0:
            return str(idx)
        return rf"$\overline{{{abs(idx):.2}}}$"

    def _vec_str(self, v: REAL_ARRAY) -> str:
        parts = [self._miller_str(x) for x in v]
        return f"Position along [{','.join(parts)}]"

    def plot_slice(
        self,
        density: REAL_ARRAY,
        t1: REAL_ARRAY,
        t2: REAL_ARRAY,
        v1: REAL_ARRAY,
        v2: REAL_ARRAY,
        atom_data: (
            tuple[REAL_ARRAY, REAL_ARRAY, list[str]] | None
        ) = None,  # (t1s, t2s, symbols) tuple or None
        cmap: str = "viridis",
        log_scale: bool = False,
        vmin: float | None = 0.0,
        vmax: float | None = None,
        interpolation: str = "lanczos",
        output: str | None = None,
    ) -> None:
        """Function to create the actual plot and figure instance.

        :param density: Charge density
        :type density: :ref:`REAL ARRAY <types>`
        :param t1: Scalar grid along :math:`v_1`
        :type t1: :ref:`REAL ARRAY <types>`
        :param t2: Scalar grid along :math:`v_2`
        :type t2:  :ref:`REAL ARRAY <types>`
        :param v1: First unit vector spanning the :math:`[hkl]` plane
        :type v1: :ref:`REAL ARRAY <types>`
        :param v2: Second unit vector spanning the :math:`[hkl]` plane
        :type v2: :ref:`REAL ARRAY <types>`
        :param atom_data: Data for atom labels, see :func:`chden.project_atoms`, defaults to None
        :type atom_data: ``tuple[REAL_ARRAY, REAL_ARRAY, list[str]] | None``, optional
        :param log_scale: Whether to plot the charge density on a base-10 logarithmic scale - useful for revealing details without colour clipping, defaults to False
        :type log_scale: ``bool``, optional
        :param vmin: Minimum value to set the colour scale at, defaults to 0.0
        :type vmin: ``float | None``, optional
        :param vmax: Maximum value to set the colour scale at, defaults to None
        :type vmax: ``float | None``, optional
        :param interpolation: Interpolation type, see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html, defaults to "lanczos"
        :type interpolation: ``str``, optional
        :param output: Filename to save as. Will save with a useful name if not provided, defaults to None
        :type output: ``str | None``, optional
        """
        l1 = t1[-1] - t1[0]
        l2 = t2[-1] - t2[0]
        transpose_label = False
        # Rotate so the longer axis is always horizontal
        if l2 > l1:
            density = density.T
            t1, t2 = t2, t1
            l1, l2 = l2, l1
            v1, v2 = v2, v1
            transpose_label = True

        fig_w = 5.0
        fig_h = fig_w * (l2 / l1) + 0.5
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        divider = mal(ax)
        extent = (t1[0], t1[-1], t2[0], t2[-1])
        vmin = float(vmin) if vmin is not None else 0.0
        vmax = float(vmax) if vmax is not None else np.max(density)
        imshow_args: dict[str, Any] = {
            "origin": "lower",
            "extent": extent,
            "cmap": cmap,
            "aspect": "equal",
            "interpolation": interpolation,
        }
        # vmax, vmin and Norm interact separately
        if log_scale:
            imshow_args["norm"] = colors.LogNorm()
        else:
            imshow_args["vmax"] = vmax
            imshow_args["vmin"] = vmin
        im = ax.imshow(density.T, **imshow_args)
        ax.tick_params(direction="out", which="both")
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(im, cax=cax, fraction=0.04, pad=0.1)
        # CONQUEST outputs cube file in e/Bohr^3 by default
        cbar.set_label(
            (r"$\log_{10}\rho$ " if log_scale else r"$\rho$  [$e a_0^{-3}$]"),
            fontsize=10,
        )

        if atom_data is not None:
            t1s, t2s, syms = atom_data
            for t1a, t2a, sym in zip(t1s, t2s, syms):
                color = _ELEMENT_COLOURS.get(sym, "#00ff80")
                if transpose_label:
                    t2a, t1a = t1a, t2a
                ax.scatter(
                    t1a, t2a, s=160, color=color, edgecolors="white", linewidths=0.8, zorder=5
                )
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

        ax.set_xlabel(f"{self._vec_str(v1)}" + r"$[a_0]$", fontsize=8)
        ax.set_ylabel(f"{self._vec_str(v2)}" + r"$[a_0]$", fontsize=8)

        fig.tight_layout()
        if output is None or output == "":
            # Create generic but useful filename
            filename = (
                f"{self.chden.hkl[0]}{self.chden.hkl[1]}{self.chden.hkl[2]}_{self.chden.offset:.3f}"
            )
            mode_string = self.chden.mode if self.chden.mode is not None else ""
            filename += f"_{mode_string}"
            plt.savefig(f"{filename}{self.extension}")
            print(f"Saved: {filename}")
        else:
            plt.savefig(output)
            print(f"Saved: {output}")

    def run(
        self,
        filename: str | None,
        vmin: float | None = 0.0,
        vmax: float | None = None,
        thickness: float = 0.5,
        log_scale: bool = False,
        cmap: str = "viridis",
        interpolation: str = "lanczos",
    ) -> None:
        """Runs the full sequence of steps:

        1. Fetches charge density files
        2. Extracts and analyses data
        3. Plots data and save

        :param filename: Filename to save as. Will save with a useful name if not provided, defaults to None
        :type filename: ``str | None``, optional
         :param vmin: Minimum value to set the colour scale at, defaults to 0.0
        :type vmin: ``float | None``, optional
        :param vmax: Maximum value to set the colour scale at, defaults to None
        :type vmax: ``float | None``, optional
        :param interpolation: Interpolation type, see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html, defaults to "lanczos"
        :type interpolation: ``str``, optional
        :param thickness: The Cartesian distance perpendicular to the plane to consider atoms as lying on the slice, defaults to 0.5 Bohr.
        :type thickness: ``float``, optional
        :param log_scale: Whether to plot the charge density on a base-10 logarithmic scale - useful for revealing details without colour clipping, defaults to False
        :type log_scale: ``bool``, optional
        :param cmap: Colour map to use, see https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.Colormap.html#matplotlib.colors.Colormap, defaults to "viridis"
        :type cmap: ``str``, optional
        """
        density, v1, v2, t1, t2, origin = self.chden.extract_slice()
        atom_data = None
        if self.show_atoms:
            _, _, n_hat = self.chden.inplane_basis()
            atom_data = self.chden.project_atoms(v1, v2, n_hat, origin, thickness=thickness)
            print(f"  Atoms within {thickness} (Bohr) of plane: {len(atom_data[0])}")
        self.plot_slice(
            density,
            t1,
            t2,
            v1,
            v2,
            atom_data,
            output=filename,
            vmax=vmax,
            vmin=vmin,
            log_scale=log_scale,
            cmap=cmap,
            interpolation=interpolation,
        )


def main() -> None:
    example_chden = chden(
        np.array([1, 0, 0]),
        0.0,
        ch1="tests/data/chden_up.cube",
        ch2="tests/data/chden_dn.cube",
        mode="sum",
    )
    filename = None
    chden_plot(example_chden, False).run(filename, log_scale=True)


if __name__ == "__main__":
    main()
