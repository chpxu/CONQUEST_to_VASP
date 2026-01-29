# CONQUEST2a

A [CONQUEST](https://github.com/OrderN/CONQUEST-release/) post-processing tool written in Python to do multiple, useful things:
- Convert CONQUEST coordinates format into `.vasp` and `.(ext)xyz` formats for quick and easy visualisation, e.g. in [VESTA](https://jp-minerals.org/vesta/en/).
- Create `xsf` files using `AtomCharge.dat` to visualise net spins
- Create supercells (larger cells formed of repeats of a unit cell)
- Process and sort (p)DOS files into something easy to use for plotting via matplotlib
- Process and sort `BandStructure.dat` into something easy to use for plotting via matplotlib
- Nearest-neighbour searching via periodic KDTree implementation
- Calculation of dihedral and planar angles

## Installation From 0.2.0
Usage is simple. In your `venv`, simply
```
pip3 install conquest2a numpy
```
If you are attempting to integrate this directly into your Nix devShell, you will have to manually build the package with `buildPythonPackage`. Support for this as a standalone package will come soon. The [devflake](https://github.com/chpxu/development-flake) in this repo automatically builds and adds it to the devshell environment.

## Usage

1. [Initialising your input](#initialising-your-input)
2. [Bandstructures](#bandstructures)
3. [pDOS](#pdos)
4. [kNN](#k-nearest-neighbours)
5. [Quantities](#quantities)

These steps assume you are already in the directory where `Conquest_input` and other relevant files sit. There is however, file path checking + absolute path resolution, for implementing when using in your own scripts, so relative paths _shouldn't_ be an issue.

### Initialising your input

First, import everything you might want to use:
```py
from conquest2a.conquest import * # necessary, (1)
from conquest2a.supercell import * # for supercell creation
from conquest2a.writers import * # to write output files to disk
from conquest2a.pdos import * # to process (p)DOS
from conquest2a.band import * # to process BandStructure.dat
from conquest2a.read.quantities import * # to process static output files without ASE
from conquest2a.algo.kdtree import periodic_kdtree # for nearest-neighbour searching
```
Next, get the path to your Conquest coordinates file, and instantiate `(1)` as
```py
test_input = conquest_input({1: "Bi", 2: "Mn", 3: "O"}) # replace this dict with your dict
test_coords_proc = conquest_coordinates_processor("./tests/data/test.dat", test_input)
```

Your `dict` inside `conquest_input()` will represent be the Conquest species index to element label map. Note that the `dict` integers should match the ones specified in `Conquest_input` and the coordinates file. Please ensure that the element labels represent real elements - the code will error out if it isn't.

### k-Nearest Neighbours

Traditional nearest-neighbour methods involve searching all atoms and specifying an arbitrary cutoff which is expensive for ridiculously large systems (around tens or hundreds of thousands or more atoms). 

CONQUEST2a gives each atom a number depending on their location in a Conquest coordinates file.

The algorithm used is a periodic KDTree, which automatically finds nearest neighbours using a binary tree. By specifying a number of neighbours $k$, you can then automatically get the $k$ closest neighbours, their interatomic distances, the element and their coordinates, assuming you initialised the Atoms correctly [above](#initialising-your-input). CONQUEST2a's implementation is a wrapper around [SciPy's KDTree](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html) to interface with the `Atom` class.

```py

from conquest2a.algo.nn import nearest_neighbours
conquest_map = conquest_input({1: "O", 2: "Bi", 3: "Mn", 4: "Mn", 5: "Mn", 6: "Mn"})
path = "./tests/data/test_output_input_coords.in"
coordsproc = conquest_coordinates_processor(path, conquest_map)
nn = nearest_neighbours(
    coordsproc, coordsproc.atoms[8]
)
nn.get_result(2) # Returns the interatomic distance in BOHR and the associated Atom
```

**WARNING**: to make index mapping easier, the first element is ALWAYS the atom you passed in to search around. E.g., to search for the **first** nearest neighbour, ensure the integer passed in to `get_result()` is **2**. 

### Bandstructures

First, ensure `conquest2a.band` is imported at the start of your file.

Get the path to your bandstructure file, and initialise the `bst_processor` class:

```py
from conquest2a.band import *
test_bst = bst_processor("./tests/data/test_BandStructure.dat")
```

Now, all the bands have been stored in a list of `band`, and can be accessed as  `test_bst.bands`. The default CONQUEST BandStructure outputs all k-points as the "x-axis" for all bands and the specific options `Process.BandStrucAxis` may/may not work.

You can then just iterate through the `bands` list and plot, see `examples/plot_test_band.py`.

### pDOS

In CONQUEST, there are multiple types of density of states (DOS) output:

1. total DOS (TDOS)
2. partial DOS (PDOS), either $l$ (angular-momentum) resolved or $lm$-resolved for each atom

The relevant source is located in `conquest2a/pdos.py`.

To begin, first import the `pdos_processor` classes:
```py

from conquest2a.pdos import pdos_processor, pdos_l_processor, pdos_lm_processor
```

If you care only about the total DOS and nothing else, use `pdos_processor` with `lm = "t"`, i.e.:
```py
tdos = pdos_processor(conquest_rundir="/yourpath", lm="t")
tdos.read_file(tdos.all_pdos_files[0]) # reads "DOS.dat"
```

This will search your directory, here `/yourpath`, for `DOS.dat` (if parameter `lm = "t"`, which is the default), and store the results in `self.blocks`, where each element of the list is a numpy array of pdos values, in order of the spin. In CONQUEST, spins are output in ascending order, so the first "block" is for Spin 1 etc.

Alternatively, you can get CONQUEST to output angular momentum-resolved DOS. In this case, you can either use `pdos_l_processor` or `pdos_lm_processor` depending on whether you have `AtomXXXXXXX_l.dat` or `AtomXXXXXXX_lm.dat` files respectively. Since they are also different objects, you may use both too. The valid filenames are then stored in `self.all_pdos_files` as a list of strings. To extract the data from a file, follow the same format

```py 
lmpdos = pdos_lm_processor(conquest_rundir="/yourpath") # lm is set automatically
lmpdos.read_file(lmpdos.all_pdos_files[0]) # reads, e.g. "Atom00000001_lm.dat" if that exists in your directory.
atom1 = lmpdos.blocks # NOTE: this is a SHALLOW COPY. If you do another read, this will be OVERWRITTEN
# atom1 = copy.deepcopy(lmpdos.blocks) # you may prefer to do this instead, if you need to read and store all the pdos output separately
```

Remark 1: in CONQUEST, you can choose to output pdos only for specific atoms, but nonetheless `self.all_pdos_files` will be in ascending order of atom index.

Remark 2: Reading pdos files automatically, and storing all of their data at once, is not implemented.

`pdos_lm_processor` and `pdos_l_processor` each have their own methods, `lm_map()` and `l_map()` respectively. This must be called after a read, e.g. as `lmpdos.lm_map()`. This groups each column of pdos files by their $(l,m)$ or $l$-value, thus forming a `dict` like
```py
{
  "0,0": [np.array(...), np.array(...)], # [spin 1, spin 2,....]
  "1,-1": [np.array(...), np.array(...)],
  "1,0": [np.array(...), np.array(...)],
  "1,1": [np.array(...), np.array(...)],
  # etc
}
```
where again the numpy arrays are in ascending order of spins. These dicts can be accessed as `processor_instance.lm_dict` or `processor_instance.l_dict`.

See `examples/plot_test_pdos.py` for an example of plotting the data obtained from a pDOS file.
### Quantities

As detailed on the [CONQUEST docs](https://conquest.readthedocs.io/en/latest/ase-conquest.html), you can manage it with [ASE](https://ase-lib.org/) indirectly by setting the flag `IO.WriteOutToASEFile True` in your `Conquest_input` file. Sometimes, you just forget to add flags when you need to, and then proceed to do numerous calculations without ASE, and thus `conquest2a/read/quantities.py` was created. 

By pointing to a file from a static run, this module will fetch the free energy, Harris-Foulkes energy, DFT total energy, forces on each atom (and assign them to the right `Atom` instances), max force and total stresses from near the end of the file.

First, load your species dictionary correctly, according to your coordinates file. Then,

```py
from conquest2a.conquest import *
from conquest2a.read.quantities import *
test_input = conquest_input({1: "Bi", 2: "Mn", 3: "O"})
test_coords_proc = conquest_coordinates_processor("./tests/data/test.dat", test_input)
output = read_static_file("tests/data/test_output.txt", test_coords_proc) # will do all the quantity fetching automatically

output.dft_energy
output.harris_foulkes_energy
#...
```

Forces can be accessed per atom.

### 0.1.0 and older
Usage is simple and there are **no external library dependencies** (currently). Either use `main.py` from the Releases tab, or clone the repo and copy `main.py` to your desired location.

1. Define a `dict` mapping your CONQUEST species to elements. E.g., if `Conquest_input` has a block like this:
  ```
  %block ChemicalSpeciesLabel
    1 208.9800000 Bi_SpinUp
    2 208.9800000 Bi_SpinDown
    3 16.000000 O
  %endblock
  ```
  Your `dict` will be `{1: "Bi", 2: "Bi", 3: "O"}`. Note that the `dict` integers should match the ones specified in `Conquest_input` and are completely arbitrary. please ensure however, that the element labels represent real elements (see `main.py`, `CONQUEST_INPUT` class)!

2. Create an instance of `CONQUEST_INPUT`, e.g. `CONQUEST_INPUT({1: "Bi", 2: "Bi", 3: "O"})`

3. Create an instance of `CONQUEST_COORDINATES`, feeding in the coordinate file you want to post-process and the instance of `CONQUEST_INPUT` created in Step 2:
  ```py
  conq = CONQUEST_COORDINATES(
      "./test/test.dat", CONQUEST_input=CONQUEST_INPUT({1: "Bi", 2: "Bi", 3: "O"})
  )
  ```
4. You may now call any of the writers with the path to the desination file, in this case, `test/test.ABC`:
   1. `vasp_writer("test/test.vasp", data=conq)`
   2. `xyz_writer("test/test.xyz", data=conq)`
   3. `extxyz_writer("test/test.extxyz", data=conq)`


## CONTRIBUTING

Edit `main.py` as necessary. Files may be refactored in the future. 

Dependencies: tested and developed on Python 3.12. Expected to work on Python >= 3.10.