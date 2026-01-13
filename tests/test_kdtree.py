import numpy as np
from conquest2a.conquest import *
from conquest2a.algo.kdtree import *
import pytest

conquest_map = conquest_input({1: "O", 2: "Bi", 3: "Mn", 4: "Mn", 5: "Mn", 6: "Mn"})
path = "./tests/data/test_output_input_coords.in"
coordsproc = conquest_coordinates_processor(path, conquest_map)
tree = periodic_kdtree(
    coordsproc.Atoms, box=np.array([coordsproc.lattice_vectors[i][i] for i in range(3)])
)


def test_box_lengths() -> None:
    # Test lattice parameters are read in correctly
    # Allow for some floating point error
    assert np.allclose(tree.box, np.array([14.1965497350, 14.03059009363509, 14.30494587850491860]))


def test_fractional_points() -> None:
    # Test points are read in correctly and in order

    # First atom in coord
    assert np.isclose(0.01414703960386700 * tree.box[0], tree.points[0, 0])
    assert np.isclose(0.97857461399301104 * tree.box[1], tree.points[0, 1])
    assert np.isclose(0.01810295041438003 * tree.box[2], tree.points[0, 2])
    # last atom in coord
    assert np.isclose(0.68909265560069932 * tree.box[0], tree.points[-1, 0])
    assert np.isclose(0.78501800639685138 * tree.box[1], tree.points[-1, 1])
    assert np.isclose(0.54878158969709301 * tree.box[2], tree.points[-1, 2])


def test_num_points() -> None:
    # Ensure this test is BEFORE any calls to knn or search
    assert len(tree.idx) == 40


# def test_initial_tree() -> None:
#     # Ensure this test is BEFORE any calls to knn or search
#     print(tree.root)


nn4 = tree.knn(coordsproc.Atoms[0], k=4)
nn1 = tree.knn(coordsproc.Atoms[0], k=10)

# print(nn)
# [(Atom(species=1, coords=array([0.18786065, 0.28189486, 0.05115697]), can_move=['T', 'T', 'T'], number=39, label='O', forces=array([0., 0., 0.]), spins=array([0., 0., 0.])), np.float64(4.941346059224852)), (Atom(species=1, coords=array([0.18959664, 0.66201423, 0.04796157]), can_move=['T', 'T', 'T'], number=38, label='O', forces=array([0., 0., 0.]), spins=array([0., 0., 0.])), np.float64(5.110146385688078)), (Atom(species=4, coords=array([0.27018768, 0.72324595, 0.78566734]), can_move=['T', 'T', 'T'], number=13, label='Mn', forces=array([0., 0., 0.]), spins=array([0., 0., 0.])), np.float64(6.091109178241822))]


# Analyse the result for k = 1, k = 4
def test_search_nn1() -> None:
    # Each element of the list is a tuple of the form (atom, dist)
    # assert len(nn1) == 1
    for atom in nn1:
        assert isinstance(atom[0], Atom)

    atom = nn1[1][0]
    # assert np.allclose(atom.coords, [0.18786065, 0.28189486, 0.05115697])
    # assert atom.label == "O"
    # assert atom.species == 1
    assert nn1[1][1] == 4.941346059224852


def test_search_nn4() -> None:
    # Analyse the result for k = 4
    # Check the first element of k = 4 is the same atom as k = 1
    assert len(nn4) == 4
    for atom in nn4:
        assert isinstance(atom[0], Atom)
    atom = nn4[1][0]
    assert np.allclose(atom.coords, [0.18786065, 0.28189486, 0.05115697])
    assert atom.label == "O"
    assert atom.species == 1
    assert nn4[1][1] == 4.941346059224852
