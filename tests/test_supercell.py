from src.conquest import *
from src.supercell import *
from src.writers import *
import pytest
test_input = CONQUEST_INPUT({"1": "Bi", "2": "Mn", "3": "Mn", "4": "O"})
test_coords_proc = CONQUEST_COORDINATES_PROCESSOR("tests/test.dat", test_input)

def test_integer_params() -> None:
    single_cell: SUPERCELL = SUPERCELL(Nx=1, Ny=1, Nz=1, coords=test_coords_proc)
    assert single_cell.Nx is not int
    assert single_cell.Ny is not int
    assert single_cell.Nz is not int
    assert single_cell.coords is not None

def test_at_least_zero() -> None:
    with pytest.raises(ValueError):
        SUPERCELL(Nx=-1, Ny=1, Nz=1, coords=test_coords_proc)
        SUPERCELL(Nx=1, Ny=-32, Nz=1, coords=test_coords_proc)
        SUPERCELL(Nx=1, Ny=1, Nz=-1, coords=test_coords_proc)

    
def non_integer_repeat() -> None:
    with pytest.raises(TypeError):
        SUPERCELL(Nx=-0.5, Ny=1, Nz=1, coords=test_coords_proc)
        SUPERCELL(Nx=0, Ny=1.5, Nz=1, coords=test_coords_proc)
        SUPERCELL(Nx=-0.5, Ny=1, Nz=0.000001, coords=test_coords_proc)

def test_single_repeat_x() -> None:
    new_cell: SUPERCELL = SUPERCELL(Nx=1, Ny=0, Nz=0, coords=test_coords_proc)
    num_repeats = (new_cell.Nx + new_cell.Ny + new_cell.Nz)  + 1
    assert int(new_cell.supercell_coords.natoms) ==  num_repeats* int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (new_cell.Nx + 1) * 5.92770000 < new_cell.supercell_coords.lattice_vectors[0][0] < (new_cell.Nx + 1) * 5.9290000 # account for floating point error
    assert (new_cell.Ny + 1) * 7.439000 < new_cell.supercell_coords.lattice_vectors[1][1] < (new_cell.Ny + 1) * 7.44100000 # account for floating point error
    assert (new_cell.Nz + 1) * 5.371999999 < new_cell.supercell_coords.lattice_vectors[2][2] < (new_cell.Nz + 1) * 5.372100000 # account for floating point error


def test_single_repeat_y() -> None:
    new_cell: SUPERCELL = SUPERCELL(Nx=0, Ny=1, Nz=0, coords=test_coords_proc)
    num_repeats = (new_cell.Nx + new_cell.Ny + new_cell.Nz)  + 1
    assert int(new_cell.supercell_coords.natoms) ==  num_repeats* int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (new_cell.Nx + 1) * 5.92770000 < new_cell.supercell_coords.lattice_vectors[0][0] < (new_cell.Nx + 1) * 5.9290000 # account for floating point error
    assert (new_cell.Ny + 1) * 7.439000 < new_cell.supercell_coords.lattice_vectors[1][1] < (new_cell.Ny + 1) * 7.44100000 # account for floating point error
    assert (new_cell.Nz + 1) * 5.371999999 < new_cell.supercell_coords.lattice_vectors[2][2] < (new_cell.Nz + 1) * 5.372100000 # account for floating point error


def test_single_repeat_z() -> None:
    new_cell: SUPERCELL = SUPERCELL(Nx=0, Ny=0, Nz=1, coords=test_coords_proc)
    num_repeats = ((new_cell.Nx + 1) * (new_cell.Ny + 1) * (new_cell.Nz + 1))
    assert int(new_cell.supercell_coords.natoms) ==  num_repeats* int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (new_cell.Nx + 1) * 5.92770000 < new_cell.supercell_coords.lattice_vectors[0][0] < (new_cell.Nx + 1) * 5.9290000 # account for floating point error
    assert (new_cell.Ny + 1) * 7.439000 < new_cell.supercell_coords.lattice_vectors[1][1] < (new_cell.Ny + 1) * 7.44100000 # account for floating point error
    assert (new_cell.Nz + 1) * 5.371999999 < new_cell.supercell_coords.lattice_vectors[2][2] < (new_cell.Nz + 1) * 5.372100000 # account for floating point error

def test_single_repeat_all() -> None:
    new_cell: SUPERCELL = SUPERCELL(Nx=2, Ny=1, Nz=0, coords=test_coords_proc)
    num_repeats = ((new_cell.Nx + 1) * (new_cell.Ny + 1) * (new_cell.Nz + 1))
    assert int(new_cell.supercell_coords.natoms) ==  num_repeats* int(new_cell.coords.natoms)
    assert len(new_cell.supercell_coords.Atoms) == num_repeats * len(new_cell.coords.Atoms)
    assert int(new_cell.supercell_coords.natoms) == len(new_cell.supercell_coords.Atoms)
    assert (new_cell.Nx + 1) * 5.927799999 < new_cell.supercell_coords.lattice_vectors[0][0] < (new_cell.Nx + 1) * 5.9290000 # account for floating point error
    assert (new_cell.Ny + 1) * 7.439000 < new_cell.supercell_coords.lattice_vectors[1][1] < (new_cell.Ny + 1) * 7.44100000 # account for floating point error
    assert (new_cell.Nz + 1) * 5.371999999 < new_cell.supercell_coords.lattice_vectors[2][2] < (new_cell.Nz + 1) * 5.372100000 # account for floating point error

