"""This file is to used to play around with the test files, and does not form part of actual unit tests
"""
from conquest2a.conquest import *
from conquest2a.supercell import *
from conquest2a.writers import *
test_input = conquest_input({"1": "Bi", "2": "Mn", "3": "Mn", "4": "O"})
test_coords_proc = conquest_coordinates_processor("test.dat", test_input)

xsf_writer_spins("./test_with_spins.xsf", atom_charge(test_coords_proc, "test_original_AtomCharge.dat"), xsf_file="./test_orig_xsf.xsf")