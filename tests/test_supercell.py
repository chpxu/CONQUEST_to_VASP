from conquest2a.conquest import *
from conquest2a.supercell import *
from conquest2a.writers import *
import pytest

test_input = conquest_input({"1": "Bi", "2": "Mn", "3": "Mn", "4": "O"})
test_coords_proc = conquest_coordinates_processor("tests/data/test.dat", test_input)


def test_integer_params() -> None:
    single_cell: supercell = supercell(
        repeats_x=1, repeats_y=1, repeats_z=1, coords=test_coords_proc
    )
    assert single_cell.repeats_x is not int
    assert single_cell.repeats_y is not int
    assert single_cell.repeats_z is not int
    assert single_cell.coords is not None


def test_at_least_zero() -> None:
    with pytest.raises(ValueError):
        supercell(repeats_x=-1, repeats_y=1, repeats_z=1, coords=test_coords_proc)
        supercell(repeats_x=1, repeats_y=-32, repeats_z=1, coords=test_coords_proc)
        supercell(repeats_x=1, repeats_y=1, repeats_z=-1, coords=test_coords_proc)


def non_integer_repeat() -> None:
    with pytest.raises(TypeError):
        supercell(repeats_x=-0.5, repeats_y=1, repeats_z=1, coords=test_coords_proc)
        supercell(repeats_x=0, repeats_y=1.5, repeats_z=1, coords=test_coords_proc)
        supercell(repeats_x=-0.5, repeats_y=1, repeats_z=0.000001, coords=test_coords_proc)


def test_single_repeat_x() -> None:
    new_cell: supercell = supercell(repeats_x=1, repeats_y=0, repeats_z=0, coords=test_coords_proc)
    num_repeats = (new_cell.repeats_x + new_cell.repeats_y + new_cell.repeats_z) + 1
    assert int(new_cell.supercell_coords.natoms) == num_repeats * int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (
        (new_cell.repeats_x + 1) * 5.92770000
        < new_cell.supercell_coords.lattice_vectors[0][0]
        < (new_cell.repeats_x + 1) * 5.9290000
    )  # account for floating point error
    assert (
        (new_cell.repeats_y + 1) * 7.439000
        < new_cell.supercell_coords.lattice_vectors[1][1]
        < (new_cell.repeats_y + 1) * 7.44100000
    )  # account for floating point error
    assert (
        (new_cell.repeats_z + 1) * 5.371999999
        < new_cell.supercell_coords.lattice_vectors[2][2]
        < (new_cell.repeats_z + 1) * 5.372100000
    )  # account for floating point error


def test_single_repeat_y() -> None:
    new_cell: supercell = supercell(repeats_x=0, repeats_y=1, repeats_z=0, coords=test_coords_proc)
    num_repeats = (new_cell.repeats_x + new_cell.repeats_y + new_cell.repeats_z) + 1
    assert int(new_cell.supercell_coords.natoms) == num_repeats * int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (
        (new_cell.repeats_x + 1) * 5.92770000
        < new_cell.supercell_coords.lattice_vectors[0][0]
        < (new_cell.repeats_x + 1) * 5.9290000
    )  # account for floating point error
    assert (
        (new_cell.repeats_y + 1) * 7.439000
        < new_cell.supercell_coords.lattice_vectors[1][1]
        < (new_cell.repeats_y + 1) * 7.44100000
    )  # account for floating point error
    assert (
        (new_cell.repeats_z + 1) * 5.371999999
        < new_cell.supercell_coords.lattice_vectors[2][2]
        < (new_cell.repeats_z + 1) * 5.372100000
    )  # account for floating point error


def test_single_repeat_z() -> None:
    new_cell: supercell = supercell(repeats_x=0, repeats_y=0, repeats_z=1, coords=test_coords_proc)
    num_repeats = (new_cell.repeats_x + 1) * (new_cell.repeats_y + 1) * (new_cell.repeats_z + 1)
    assert int(new_cell.supercell_coords.natoms) == num_repeats * int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (
        (new_cell.repeats_x + 1) * 5.92770000
        < new_cell.supercell_coords.lattice_vectors[0][0]
        < (new_cell.repeats_x + 1) * 5.9290000
    )  # account for floating point error
    assert (
        (new_cell.repeats_y + 1) * 7.439000
        < new_cell.supercell_coords.lattice_vectors[1][1]
        < (new_cell.repeats_y + 1) * 7.44100000
    )  # account for floating point error
    assert (
        (new_cell.repeats_z + 1) * 5.371999999
        < new_cell.supercell_coords.lattice_vectors[2][2]
        < (new_cell.repeats_z + 1) * 5.372100000
    )  # account for floating point error


def test_single_repeat_all() -> None:
    new_cell: supercell = supercell(repeats_x=2, repeats_y=1, repeats_z=0, coords=test_coords_proc)
    num_repeats = (new_cell.repeats_x + 1) * (new_cell.repeats_y + 1) * (new_cell.repeats_z + 1)
    assert int(new_cell.supercell_coords.natoms) == num_repeats * int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (
        (new_cell.repeats_x + 1) * 5.927799999
        < new_cell.supercell_coords.lattice_vectors[0][0]
        < (new_cell.repeats_x + 1) * 5.9290000
    )  # account for floating point error
    assert (
        (new_cell.repeats_y + 1) * 7.439000
        < new_cell.supercell_coords.lattice_vectors[1][1]
        < (new_cell.repeats_y + 1) * 7.44100000
    )  # account for floating point error
    assert (
        (new_cell.repeats_z + 1) * 5.371999999
        < new_cell.supercell_coords.lattice_vectors[2][2]
        < (new_cell.repeats_z + 1) * 5.372100000
    )  # account for floating point error
