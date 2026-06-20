from dataclasses import dataclass, field
import re
import numpy as np
from conquest2a.conquest import block_processor
import conquest2a._types as c2at


@dataclass
class band:
    spin: int
    index: int  # band index
    kpoint: c2at.INT_ARRAY = field(default_factory=lambda: np.array([0, 0, 0]))  # kpoint index.
    energies: c2at.REAL_ARRAY = field(
        default_factory=lambda: np.empty((1, 3))
    )  # energies at that kpoint for this band


class bst_processor(block_processor):
    def __init__(self, bst_file: str) -> None:
        """Class to store data from CONQUEST's BandStructure.dat
        
        :param bst_file: Path to a file containing CONQUEST band structure data.
        :type bst_file: ``str``
        """
        super().__init__()
        self.bst_file = bst_file
        self.blocks: list[c2at.REAL_ARRAY] = []
        self.bands: list[band] = []
        self.fermi_level: float = 0.0
        self.is_shifted_to_fermi: bool = True
        self.block_counter: int = 0
        self.current_spin = 0
        self.read_file(filename=self.bst_file)
    
    def process_headers(
        self,
        line: str,
    ) -> None:
        if "# Spin" in line:
            self.num_spins += 1
            self.current_spin += 1
        if "# Original" in line:
            result = re.findall(self.re_float, line)
            self.fermi_level = float(result[0])
        if "# Band " in line:
            band_index = re.findall(self.re_index, line)
            self.bands.append(band(index=int(band_index[0]), spin=self.current_spin))
        if not line.startswith("# Bands shifted"):
            self.is_shifted_to_fermi = False

    def process_block(self, line: str) -> None:
        if line == "&":
            if self.current_block:
                arrayed_block = np.array(self.current_block, dtype=float)
                self.bands[self.block_counter].kpoint = arrayed_block[:, 0]
                self.bands[self.block_counter].energies = arrayed_block[:, 1]
                self.block_counter += 1
                self.current_block = []
        else:
            self.current_block.append(np.array(line.split()).astype(np.float64))


class bst:
    DEFAULT_SPIN_COLORS: dict[int, str] = {1: "black", 2: "red"}
 
    def __init__(
        self,
        processor,
        spin_colors: dict[int, str] | None = None,
    ) -> None:
        """Plot a bandstructure using ``BandStructure.dat``.
    
        Bands are plotted as energy (eV) against :math:`k`-point index. CONQUEST only supports collinear spin only.
    
        Energies may be plotted shifted to the Fermi level or not.

        :param processor: A :class:`bst_processor` instance whose ``bands``
            attribute is a :class:`list` of :class:`band` objects.
        :param spin_colors: Mapping from spin channel (``1`` = up, ``2`` = down)
            to a Matplotlib colour string.  Missing keys fall back to
            :attr:`DEFAULT_SPIN_COLORS`.  Pass ``None`` to use the defaults (black, red)
            unchanged.
        :type spin_colors: ``dict[int, str]`` or ``None``
        """
        self._processor = processor
        self._spin_colors = {**self.DEFAULT_SPIN_COLORS, **(spin_colors or {})}
  
    def plot(
        self,
        band_range: tuple[int, int] | None = None,
        energy_range: tuple[float, float] | None = None,
        figsize: tuple[float, float] = (8, 6),
        ylabel: str = r"Energy (eV)",
        save_path: str,
    ) -> None:
        """Generate save the band structure.
 
        :param band_range: Inclusive ``(min_index, max_index)`` range of band
            indices to plot.  Bands whose :attr:`~band.index` falls outside
            this range are excluded.  Pass ``None`` to include all bands.
        :type band_range: ``tucple[int, int]`` or ``None``
 
        :param energy_range: :math:`(E_{\\min}, E_{\\max})` energy window in eV.
            Only bands that contain at least one :math:`k` -point with an energy
            inside the closed interval :math:`[E_{\\min},\\, E_{\\max}]` are plotted, and the
            y-axis is clipped to this range.  Pass ``None`` (default) to  plot all bands.

            Note: energy filtering is done based on the energy ranges inside :class:`bst_processor`. Shifting the Fermi level is purely done for plotting purposes.
        :type energy_range: ``tuple[float, float]`` or ``None``
 
        :param figsize: ``(width, height)`` of the figure in inches, forwarded
            directly to :func:`matplotlib.pyplot.subplots`.
        :type figsize: ``tuple[float, float]``
 
        :param ylabel: Y-axis label.
        :type ylabel: ``str``
  
        :param save_path: Filesystem path at which the figure is saved via
            :meth:`~matplotlib.figure.Figure.savefig`.  Format is inferred from the file extension.
        :type save_path: ``str``
 
        :raises ValueError: If no bands remain after the *band_range* and/or
            *energy_range* filters are applied.
         """
        bands = self._filter_bands(band_range, energy_range)
        if not bands:
            raise ValueError(
                "No bands remain after filtering. "
                "Check your band_range / energy_range arguments."
            )
 
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_ylabel(ylabel)
        self._draw_bands(ax, bands)
 
        if energy_range is not None:
            ax.set_ylim(*energy_range)
        ax.legend()
        plt.savefig(save_path, bbox_inches="tight")
 
    def _filter_bands(
        self,
        band_range: tuple[int, int] | None,
        energy_range: tuple[float, float] | None,
    ) -> list:
        """Return the subset of bands that satisfy all active filters.
 
        :param band_range: Inclusive ``(min_index, max_index)`` band-index
            filter, or ``None`` to disable.
        :type band_range: tuple[int, int] or None
 
        :param energy_range: ``(E_{\\min}, E_{\\max})`` energy filter in eV,
            or ``None`` to disable.  A band passes if *any* of its energy
            values falls within the closed interval
            :math:`[E_{\\min},\\, E_{\\max}]`.
        :type energy_range: tuple[float, float] or None
 
        :returns: Filtered list of :class:`band` objects.
        :rtype: ``list[band]``
        """
        result = []
 
        for b in self._processor.bands:
            if band_range is not None:
                lo, hi = band_range
                if not (lo <= b.index <= hi):
                    continue
 
            if energy_range is not None:
                e_lo, e_hi = energy_range
                energies = np.asarray(b.energies).ravel()
                if not np.any((energies >= e_lo) & (energies <= e_hi)):
                    continue
 
            result.append(b)
 
        return result
 
    def _draw_bands(self, ax: Axes, bands: list) -> None:
        """Plot each band as energy versus sequential k-point index.
 
        Spin-up bands (``spin == 1``) and spin-down bands (``spin == 2``) are
        drawn in the colours defined by :attr:`_spin_colors`.
 
        :param ax: The :class:`~matplotlib.axes.Axes` on which to draw.
        :param bands: Pre-filtered list of :class:`band` objects to plot.
        :type bands: ``list[band]``
        """
        _spin_labels = {1: "Spin up", 2: "Spin down"}
        _seen: set[int] = set()
        for b in bands:
            energies = np.asarray(b.energies).ravel()
            k_indices = np.arange(len(energies))
            color = self._spin_colors.get(b.spin, "black")
            label = _spin_labels.get(b.spin) if b.spin not in _seen else None
            _seen.add(b.spin)

            ax.plot(k_indices, energies, color=color, linewidth=0.8)
