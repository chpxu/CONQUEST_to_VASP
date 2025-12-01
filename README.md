# CONQUEST2a

A [CONQUEST](https://github.com/OrderN/CONQUEST-release/) post-processing tool written in Python to do multiple, useful things:
- Convert CONQUEST coordinates format into `.vasp` and `.(ext)xyz` formats for quick and easy visualisation, e.g. in [VESTA](https://jp-minerals.org/vesta/en/).
- Create `xsf` files using `AtomCharge.dat` to visualise net spins
- Create supercells (larger cells formed of repeats of a unit cell)
- Process and sort (p)DOS files into something easy to use for plotting via matplotlib
- Process and sort `BandStructure.dat` into something easy to use for plotting via matplotlib

## Installation From 0.2.0
Usage is simple. In your `venv`, simply
```
pip3 install conquest2a numpy
```
If you are attempting to integrate this directly into your Nix devShell, you will have to manually build the package with `buildPythonPackage`. Support for this as a standalone package will come soon.

## Usage

These steps assume you are already in the directory where `Conquest_input` and other relevant files sit. There is however, file path checking + absolute path resolution, for implementing when using in your own scripts, so relative paths _shouldn't_ be an issue.

### (partial) Density of States

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
tdos.read_pdos_file(tdos.all_pdos_files[0]) # reads "DOS.dat"
```

This will search your directory, here `/yourpath`, for `DOS.dat` (if parameter `lm = "t"`, which is the default), and store the results in `self.blocks`, where each element of the list is a numpy array of pdos values, in order of the spin. In CONQUEST, spins are output in ascending order, so the first "block" is for Spin 1 etc.

Alternatively, you can get CONQUEST to output angular momentum-resolved DOS. In this case, you can either use `pdos_l_processor` or `pdos_lm_processor` depending on whether you have `AtomXXXXXXXX_l.dat` or `AtomXXXXXXXX_lm.dat` files respectively. Since they are also different objects, you may use both too. The valid filenames are then stored in `self.all_pdos_files` as a list of strings. To extract the data from a file, follow the same format

```py 
lmpdos = pdos_lm_processor(conquest_rundir="/yourpath") # lm is set automatically
lmpdos.read_pdos_file(lmpdos.all_pdos_files[0]) # reads, e.g. "Atom00000001_lm.dat" if that exists in your directory.
atom1 = lmpdos.blocks # NOTE: this is a SHALLOW COPY. If you do another read, this will be OVERWRITTEN
# atom1 = copy.deepcopy(lmpdos.blocks) # you may prefer to do this instead, if you need to read and store all the pdos output separately
```

Remark 1: in CONQUEST, you can choose to output pdos only for specific atoms, but nonetheless `self.all_pdos_files` will be in ascending order of atom index.

Remark 2: Reading pdos files automatically, and storing all of their data at once, is not implemented.

`pdos_lm_processor` and `pdos_l_processor` each have their own methods, `lm_map()` and `l_map()` respectively, that is called whenever a class instance is defined. This automatically groups each column of pdos files by their $(l,m)$ or $l$-value, thus forming a `dict` like
```py
{
  "0,0": [np.array(...), np.array(...)], # [spin 1, spin 2,....]
  "1,-1": [np.array(...), np.array(...)],
  "1,0": [np.array(...), np.array(...)],
  "1,1": [np.array(...), np.array(...)],
  # etc
}
```
where again the numpy arrays are in ascending order of spins.

See `examples/plot_test_pdos.py` for an example of plotting the data obtained from a pDOS file.

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

Dependencies: tested and developed on Python 3.12. Expected to work on Python >= 3.8.