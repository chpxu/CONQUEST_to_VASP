from pathlib import Path
import copy
import numpy as np
import pytest

from conquest2a.conquest import Atom, conquest_coordinates, conquest_coordinates_processor, conquest_species
from conquest2a.cell.transform import transform_unit_cell
from conquest2a._types import REAL_ARRAY
data_dir = "tests/data/transforms/"
input_cell = f"{data_dir}/unitcell.dat"
final_cell = f"{data_dir}/expected_transformed_cell.dat"

BMO = conquest_species({1: "O", 2: "Bi", 3:"Mn", 4:"Mn"})
P = np.array([[1 / 2, 0, 1 / 2], [0, 1, 0], [-1 / 2, 0, 1 / 2]])
identity = np.eye(3)


@pytest.fixture
def species() ->conquest_species:
    return BMO


@pytest.fixture
def init_cq_proc(species: conquest_species) ->conquest_coordinates_processor:
    return conquest_coordinates_processor(input_cell, species)

@pytest.fixture
def final_cq_proc(species: conquest_species) ->conquest_coordinates_processor:
    return conquest_coordinates_processor(final_cell, species)

def volume(lattice_vectors: REAL_ARRAY) -> float:
    return abs(np.dot(lattice_vectors[0], np.cross(lattice_vectors[1], lattice_vectors[2])))


# test input_cell is not touched by transform_unit_cell

def test_same_atom_count(init_cq_proc: conquest_coordinates_processor) -> None:
    orig_natoms = init_cq_proc.coords.natoms
    num_atoms = len(init_cq_proc.coords.atoms)
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(P)
    assert int(orig_natoms) == num_atoms
    assert int(init_cq_proc.coords.natoms) == int(orig_natoms)
    assert num_atoms == len(init_cq_proc.coords.atoms)

def test_same_lattice_vect(init_cq_proc: conquest_coordinates_processor)-> None:
    orig_vects = init_cq_proc.coords.lattice_vectors
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(P)
    # Float comparison
    assert np.allclose(init_cq_proc.coords.lattice_vectors, orig_vects, atol=1e-7)

def test_same_atoms(init_cq_proc: conquest_coordinates_processor)-> None:
    """
    Test that each atom in the list of atoms is unchanged after calling tuc.transform()
    Check labels, species number and position
    """
    orig_atoms = copy.deepcopy(init_cq_proc.coords.atoms)
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(P)
    # Float comparison
    for i in range(1, len(orig_atoms)):
        assert orig_atoms[i].label == init_cq_proc.coords.atoms[i].label
        assert orig_atoms[i].number == init_cq_proc.coords.atoms[i].number
        assert orig_atoms[i].species == init_cq_proc.coords.atoms[i].species
        assert orig_atoms[i].can_move == init_cq_proc.coords.atoms[i].can_move
        assert np.allclose(init_cq_proc.coords.atoms[i].coords, orig_atoms[i].coords, atol=1e-7)

# null test: identity leaves everything unchanged


def test_identity_same_atom_count(init_cq_proc: conquest_coordinates_processor) -> None:
    num_atoms = len(init_cq_proc.coords.atoms)
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(identity)
    assert int(tuc.transformed_cell_coords.natoms) == num_atoms
    assert int(init_cq_proc.coords.natoms) == int(tuc.transformed_cell_coords.natoms)
    assert num_atoms == len(tuc.transformed_cell_coords.atoms)

def test_identity_same_lattice_vect(init_cq_proc: conquest_coordinates_processor)-> None:
    orig_vects = init_cq_proc.coords.lattice_vectors
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(identity)
    # Float comparison
    assert np.allclose(orig_vects, tuc.transformed_cell_coords.lattice_vectors, atol=1e-7)

def test_identity_same_atoms(init_cq_proc: conquest_coordinates_processor)-> None:
    """
    Test that each atom in the list of atoms is unchanged after calling tuc.transform(identity)
    Check labels, species number and position
    """
    orig_atoms = copy.deepcopy(init_cq_proc.coords.atoms)
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(identity)
    new_unit_cell = tuc.transformed_cell_coords
    # Float comparison
    for i in range(1, len(orig_atoms)):
        assert orig_atoms[i].label == new_unit_cell.atoms[i].label
        assert orig_atoms[i].number == new_unit_cell.atoms[i].number
        assert orig_atoms[i].species == new_unit_cell.atoms[i].species
        assert orig_atoms[i].can_move == new_unit_cell.atoms[i].can_move
        assert np.allclose(new_unit_cell.atoms[i].coords, orig_atoms[i].coords, atol=1e-7)

@pytest.mark.parametrize(
    "scale_matrix,multiplier",
    [
        (np.diag([2.0, 1.0, 1.0]), 2),
        (np.diag([1.0, 3.0, 1.0]), 3),
        (np.diag([1.0, 1.0, 2.0]), 2),
    ],
)
def test_supercell_expansion_multiplies_atom_count(init_cq_proc, scale_matrix, multiplier):
    original_count = len(init_cq_proc.coords.atoms)
    transformer = transform_unit_cell(init_cq_proc.coords)
    result = transformer.transform(scale_matrix)
    assert int(result.natoms) == original_count * multiplier

# Now start applying a general transform
# Under this transform, the 80 atom cell becomes a 40 atom one
def test_P_atom_count(init_cq_proc: conquest_coordinates_processor, final_cq_proc: conquest_coordinates_processor) -> None:
    num_atoms = len(init_cq_proc.coords.atoms)
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(P)
    final_atoms = len(tuc.transformed_cell_coords.atoms)

    assert final_atoms * 2 == num_atoms
    assert int(init_cq_proc.coords.natoms) == 2* int(tuc.transformed_cell_coords.natoms)
    assert len(final_cq_proc.coords.atoms) * 2 == num_atoms
    assert int(final_cq_proc.coords.natoms) == int(tuc.transformed_cell_coords.natoms)

def test_P_lattice_vect(init_cq_proc: conquest_coordinates_processor, final_cq_proc: conquest_coordinates_processor)-> None:
    final_vects = np.array([[14.3929907019,	0.0000000000,	0.0000000000],
        [0.0000000000,	13.9739181051,0.0000000000],
        [-0.3077087491	,0.0000000000,14.3897010626]])
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(P)
    assert np.allclose(final_cq_proc.coords.lattice_vectors, final_vects, atol=1e-7)
    assert np.allclose(final_cq_proc.coords.lattice_vectors, tuc.transformed_cell_coords.lattice_vectors, atol=1e-7)

def test_P_atoms(init_cq_proc: conquest_coordinates_processor, final_cq_proc: conquest_coordinates_processor)-> None:
    """
    Test that each atom in the list of atoms is the same as the expected (read in)
    """
    tuc = transform_unit_cell(init_cq_proc.coords)
    tuc.transform(P)
    new_unit_cell = tuc.transformed_cell_coords
    expected_atoms = final_cq_proc.coords.atoms
    for i in range(1, len(tuc.transformed_cell_coords.atoms)):
        assert expected_atoms[i].label == new_unit_cell.atoms[i].label
        assert expected_atoms[i].species == new_unit_cell.atoms[i].species
        assert expected_atoms[i].can_move == new_unit_cell.atoms[i].can_move
        assert np.allclose(new_unit_cell.atoms[i].coords, expected_atoms[i].coords, atol=1e-7)
        # We do not test atom number because it is carried over from original unit cell
        # There is no guarantee atom number is unchanged, but it may be useful for
        #comparisons by the user later
        # assert expected_atoms[i].number == new_unit_cell.atoms[i].number


# test errors