# CONQUEST2a

A [CONQUEST](https://github.com/OrderN/CONQUEST-release/) post-processing tool written in Python to convert its coordinates format into `.vasp` and `.(ext)xyz` formats for quick and easy visualisation, e.g. in [VESTA](https://jp-minerals.org/vesta/en/).

## Usage

Usage is simple and there are **no external library dependencies**. Either use `main.py` from the Releases tab, or clone the repo and copy `main.py` to your desired location.

1. Define a `dict` mapping your CONQUEST species to elements. E.g., if `Conquest_input` has a block like this:
  ```
  %block ChemicalSpeciesLabel
    1 208.9800000 Bi_SpinUp
    2 208.9800000 Bi_SpinDown
    3 16.000000 O
  %endblock
  ```
  Your `dict` will be `{1: "Bi", 2: "O"}`. Note that the `dict` integers do not match the ones specified in `Conquest_input` and are completely arbitrary please ensure however, that the element labels represent real elements (see `main.py`, `CONQUEST_INPUT` class)!
2. Create an instance of `CONQUEST_INPUT`, e.g. `CONQUEST_INPUT({1: "Bi", 2: "O"})`
3. Create an instance of `CONQUEST_COORDINATES`, feeding in the coordinate file you want to post-process and the instance of `CONQUEST_INPUT` created in Step 2:
  ```py
  conq = CONQUEST_COORDINATES(
      "./test/test.dat", CONQUEST_input=CONQUEST_INPUT({1: "Bi", 2: "O"})
  )
  ```
4. You may now call any of the writers:
   1. `vasp_writer("test/test.vasp", data=conq)`
   2. `xyz_writer("test/test.xyz", data=conq)`
   3. `extxyz_writer("test/test.extxyz", data=conq)`


## CONTRIBUTING

Edit `main.py` as necessary. Files may be refactored in the future. 

Dependencies: tested and developed on Python 3.12. Expected to work on Python >= 3.8.