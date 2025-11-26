# CONQUEST2a

A [CONQUEST](https://github.com/OrderN/CONQUEST-release/) post-processing tool written in Python to do multiple, useful things:
- Convert CONQUEST coordinates format into `.vasp` and `.(ext)xyz` formats for quick and easy visualisation, e.g. in [VESTA](https://jp-minerals.org/vesta/en/).
- Create `xsf` files using `AtomCharge.dat` to visualise net spins
- Create supercells (larger cells formed of repeats of a unit cell)
- Process and sort (p)DOS files into something easy to use for plotting via matplotlib
- Process and sort `BandStructure.dat` into something easy to use for plotting via matplotlib
## Usage

### From 0.2.0
Usage is simple. In your `venv`, simply
```
pip3 install conquest2a
```
If you are attempting to integrate this directly into your Nix devShell, you will have to manually build the package with `buildPythonPackage`. Support for this as a standalone package will come soon.

You may also simply just copy the `conquest2a` directory and use the files as scripts.

### 0.1.0
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