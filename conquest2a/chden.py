from ctypes import ArgumentError
import os
from pathlib import Path
from typing import Any, Literal
import numpy as np
from ase.io.cube import read_cube
from scipy.ndimage import map_coordinates
from conquest2a._types import INT_ARRAY, REAL_ARRAY, REAL_NUMBER
from conquest2a.constants import MPLGENERIC

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


class chden:
    def __init__(
        self,
        hkl: INT_ARRAY,
        offset: float,
        ch1: str,
        ch2: str | None = None,
        mode: Literal["sum", "diff"] | None = None,
    ) -> None:
        """Initialise charge density data

        CONQUEST supplies up to two chden .cube files:
            - chden_up/dn.cube for a spin polarised calculation
            - a single cube file, e.g. chden.cube
        This class uses ASE's .cube file processor to get atoms and data.
        It combines chden up and dn if both arguments are supplied, to make a total charge density
        It also does up - dn to see the spin difference, if both arguments are supplied
        Args:
            chup (str): Path to a CONQUEST chden file
            chdn (str): Path to another CONQUEST chden file
        """
        self.hkl = hkl
        self.offset = offset
        self.mode = mode
        self.dens1 = self.load_cube(ch1)
        self.data = self.dens1[0]
        self.dens2 = None
        self.atoms = self.dens1[1]
        self.cell = self.atoms.get_cell()  # ASE cell information
        if self.mode is not None and ch2 is None:
            raise ArgumentError("Cannot have sum/diff mode if a second file is not provided!")
        if ch2 is not None:
            self.dens2 = self.load_cube(ch2)
            if self.mode == "sum":
                self.data += self.dens2[0]  # add data arrays together
            elif self.mode == "diff":
                self.data -= self.dens2[0]  # subtract data arrays together

        if hkl[0] == 0 and hkl[1] == 0 and hkl[2] == 0:
            raise ValueError("Miller indices (h k l) cannot all be zero.")

    def resolve_path(self, filename: str) -> None:
        abs_coord_path: Path = Path(filename)
        if not abs_coord_path.exists():
            raise FileNotFoundError(f"{abs_coord_path} not found.")
        if abs_coord_path.is_file() and os.stat(abs_coord_path).st_size <= 0:
            raise RuntimeError(f"{abs_coord_path} was an existing file, but has no file contents.")
        self.abs_input_path = abs_coord_path

    def load_cube(self, filename: str) -> Any:
        self.resolve_path(filename=filename)
        with open(filename, "r") as fh:
            cube = read_cube(fh)
        fh.close()
        return cube["data"], cube["atoms"]

    def inplane_basis(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return two orthonormal Cartesian vectors (v1, v2) spanning the (hkl) plane
        and the unit plane-normal n_hat.

        Steps
        --------
        1. The plane normal in Cartesian space is  n = ha* + kb* + lc*
        where a*, b*, c* are the (non-2π) reciprocal basis vectors
        (rows of inv(cell)^T).
        2. Two fractional-space null vectors of [h k l] are found via SVD.
        3. Those are converted to Cartesian, then Gram-Schmidt
        orthonormalised so the axes are perpendicular in real-space.

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
        """
        Return a Cartesian point lying on the plane  ha + kb + lc = offset.

        `offset` is a dimensionless fractional intercept (0-1 spans one
        interplanar period). Pick the simplest fractional coordinate
        satisfying the plane equation, then convert to angstrom.
        """
        n = self.hkl / (self.hkl @ self.hkl)  # normal direction in fractional space
        centre = np.array([0.5, 0.5, 0.5])
        frac = centre - (centre @ self.hkl) * n + self.offset * n
        return frac @ self.cell

    def inplane_range(self, v1: np.ndarray, v2: np.ndarray) -> tuple[float, float]:
        L1 = sum(abs(self.cell[i] @ v1) for i in range(3))
        L2 = sum(abs(self.cell[i] @ v2) for i in range(3))
        return L1, L2

    def extract_slice(
        self,
        n_points: int = 1000,
        interp_order: int = 5,
    ) -> Any:
        """
        Sample the charge density on the (hkl) plane at fractional offset `offset`.

        Algorithm
        ---------
        For each point on a 2-D (t1, t2) grid centred on the plane origin:
        1. Compute its Cartesian position:
            p = origin + t1*v1 + t2*v2
        2. Convert to fractional coordinates:
            s = p @ inv(cell)
        3. Wrap into the unit cell to enforce PBC.
        4. Map fractional -> voxel index and interpolate via
            scipy.ndimage.map_coordinates.

        Returns
        -------
        density: (n_points, n_points) ndarray - charge density in e/a_0^3
        v1, v2: unit in-plane Cartesian vectors (Å)
        t1, t2: 1-D coordinate axes (Å, symmetric about 0)
        origin: Cartesian slice origin (Å)
        """
        Nx, Ny, Nz = self.data.shape
        cell_inv = np.linalg.inv(self.cell)

        v1, v2, _ = self.inplane_basis()
        origin = self.plane_origin()

        L1, L2 = self.inplane_range(v1, v2)
        t1 = np.linspace(-L1 / 2, L1 / 2, n_points)
        t2 = np.linspace(-L2 / 2, L2 / 2, n_points)
        T1, T2 = np.meshgrid(t1, t2, indexing="ij")
        pts_cart = origin + T1[..., None] * v1 + T2[..., None] * v2
        pts_frac = (pts_cart @ cell_inv) % 1.0
        vox = pts_frac * np.array([Nx, Ny, Nz])
        pad = interp_order + 1
        data_padded = np.pad(self.data, pad, mode="wrap")

        coords = vox.reshape(-1, 3).T + pad
        density = map_coordinates(data_padded, coords, order=interp_order, mode="nearest").reshape(
            n_points, n_points
        )

        return density, v1, v2, t1, t2, origin

    def project_atoms(
        self,
        v1: np.ndarray,
        v2: np.ndarray,
        n_hat: np.ndarray,
        origin: np.ndarray,
        thickness: float = 0.5,
    ) -> Any:
        """
        Project atoms within `thickness` Å of the slice plane onto (v1, v2) axes.

        Returns (t1_coords, t2_coords, symbol_list).
        """
        positions = self.atoms.get_positions()  # (N, 3) Å
        symbols = self.atoms.get_chemical_symbols()
        disp = positions - origin  # displacement from slice origin
        dist_norm = disp @ n_hat  # signed distance to plane

        mask = np.abs(dist_norm) < thickness
        t1_proj = disp[mask] @ v1
        t2_proj = disp[mask] @ v2
        syms = [symbols[i] for i in np.where(mask)[0]]
        return t1_proj, t2_proj, syms


class chden_plot:
    """
    Helper class to plot output
    """

    def __init__(
        self, chden_instance: chden, show_atoms: bool = False, format: str = "png"
    ) -> None:
        self.chden = chden_instance
        self.show_atoms = show_atoms
        self.format = format
        # Check for matplotlib and SciencePlots
        import importlib.util

        mpl_spec = importlib.util.find_spec("matplotlib")
        mpl_found = mpl_spec is not None
        sp_spec = importlib.util.find_spec("scienceplots")
        sp_found = sp_spec is not None
        if mpl_found:
            import matplotlib.pyplot as plt

            self.plt = plt
        else:
            raise ImportError("matplotlib could not be found")
        if sp_found and mpl_found:
            import scienceplots

            self.plt.style.use(["science", "no-latex"])
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1 import make_axes_locatable as mal

        self.mal = mal
        mpl.rcParams.update(MPLGENERIC)

    def miller_str(self, idx: int) -> str:
        if idx >= 0:
            return str(idx)
        else:
            return rf"$\overline{{{abs(idx):.2}}}$"

    def vec_str(self, v: REAL_ARRAY) -> str:
        parts = [self.miller_str(x) for x in v]
        return f"Position along [{','.join(parts)}]"

    def plot_slice(
        self,
        density: np.ndarray,
        t1: np.ndarray,
        t2: np.ndarray,
        v1: np.ndarray,
        v2: np.ndarray,
        atom_data: Any = None,  # (t1s, t2s, symbols) tuple or None
        cmap: str = "inferno",
        log_scale: bool = False,
        vmin: REAL_NUMBER | None = 0.0,
        vmax: REAL_NUMBER | None = None,
        interpolation: str = "lanczos",
        output: str | None = None,
    ) -> None:
        L1 = t1[-1] - t1[0]
        L2 = t2[-1] - t2[0]
        transpose_label = False
        # Rotate so the longer axis is always horizontal
        if L2 > L1:
            density = density.T
            t1, t2 = t2, t1
            L1, L2 = L2, L1
            v1, v2 = v2, v1
            transpose_label = True

        fig_w = 5.0
        fig_h = fig_w * (L2 / L1) + 0.5
        fig, ax = self.plt.subplots(figsize=(fig_w, fig_h))
        divider = self.mal(ax)
        plot_data = np.log1p(np.clip(density, 0, None)) if log_scale else density
        extent = (t1[0], t1[-1], t2[0], t2[-1])
        im = ax.imshow(
            plot_data.T,
            origin="lower",
            extent=extent,
            cmap=cmap,
            aspect="equal",
            interpolation=interpolation,
            vmin=float(vmin) if vmin is not None else 0.0,
            vmax=float(vmax) if vmax is not None else np.max(plot_data),
        )
        ax.tick_params(direction="out", which="both")
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(im, cax=cax, fraction=0.04, pad=0.1)
        cbar.set_label(
            (
                r"$\log\,(1+\rho)$   [$\mathrm{eV} a_0^{-3}$]"
                if log_scale
                else r"$\rho$  [$\mathrm{eV} a_0^{-3}$]"
            ),
            fontsize=12,
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

        ax.set_xlabel(f"{self.vec_str(v1)}" + r"$[\mathrm{\AA}]$", fontsize=8)
        ax.set_ylabel(f"{self.vec_str(v2)}" + r"$[\mathrm{\AA}]$", fontsize=8)

        fig.tight_layout()
        if output is None or output == "":
            # Create generic but useful filename
            filename = f"{self.chden.hkl[0]}{self.chden.hkl[1]}{self.chden.hkl[2]}_{self.chden.offset:.3f}" 
            filename += f"{self.chden.mode if self.chden.mode is not None else ""}"
            extension = f".{self.format}"
            self.plt.savefig(f"{filename}{extension}")
            print(f"Saved: {filename}")
        else:
            self.plt.savefig(output)
            print(f"Saved: {output}")

    def run(self, filename: str | None, vmin: REAL_NUMBER | None = 0.0, vmax: REAL_NUMBER | None = None, thickness: float = 0.5, log_scale: bool=False, cmap: str = "inferno", interpolation: str = "lanczos") -> None:
        """
        Runs the full sequence of steps including plotting
        """
        density, v1, v2, t1, t2, origin = self.chden.extract_slice()
        atom_data = None
        if self.show_atoms:
            _, _, n_hat = self.chden.inplane_basis()
            atom_data = self.chden.project_atoms(v1, v2, n_hat, origin, thickness=thickness)
            print(f"  Atoms within {thickness} (ang) of plane: {len(atom_data[0])}")
        self.plot_slice(density, t1, t2, v1, v2, atom_data, output=filename, vmax=vmax, vmin=vmin, log_scale=log_scale, cmap=cmap, interpolation=interpolation)


def main() -> None:
    example_chden = chden(np.array([1, 0, 0]), 0.0, ch1="tests/data/chden_up.cube", ch2="tests/data/chden_dn.cube", mode="sum")
    filename = None
    chden_plot(example_chden, False).run(filename, vmax=0.1)


if __name__ == "__main__":
    main()
