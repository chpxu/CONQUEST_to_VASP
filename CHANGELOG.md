# 0.2.0

## Features
- Refactored project structure to prepare for package release
- Split file writer classes into `src/writers.py`
- Classes dealing with processing the Conquest coordinates file is in `src/conquest.py`
- Classes making Conquest coordinates-compliant supercells are in `src/supercell.py` (NEW)
<!-- - Classes generating a Conquest_input file with sane defaults is in `src/generate_run.py` (NEW)-->
- `xsf` file-format with spin (NEW)
    - Modify the `Atom` class with a `spin` attribute that can optionally be filed
    - Create a `atom_charge` class to process `AtomCharge.dat` 
    - create an `xsf_writer` class in `src/writers.py` to write this data, including the spin
- pDOS processing (NEW!)
    - Supports reading pDOS files for $l,m$ and $l$ decomposed situations
    - Reads `Conquest_input` to search for the used coordinates file, and maps each file to its atom number.
## Misc

- Added flake8, pylint, mypy and black to the repo
- Added some initial tests in `tests/`

# 0.1.0 - Initial release

- Script for converting coordinate files to and from VASP and (ext)XYZ format.