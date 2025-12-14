from conquest2a.eigenvalue import eigenvalues_processor, k_point_blocks
import pytest
import numpy as np

data_file = "tests/data/test_eigenvalues.dat"
target_bandgap = 0.12999176310841296
proc = eigenvalues_processor(file=data_file)


# Check initialisation of default k_point_blocks
def test_default_kptsblock() -> None:
    block = k_point_blocks(k_index=1, weight=0.00003)
    assert np.array_equal(block.eigenvalues, np.array([0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    assert np.array_equal(block.eigenvalues_with_fermi, np.array([0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    assert np.array_equal(block.k_vector, np.array([0.0, 0.0, 0.0]))


def test_file_exists() -> None:
    proc.resolve_path()


def test_num_eigs() -> None:
    assert proc.num_eigenvalues_per_block == 304


def test_fermis() -> None:
    assert proc.fermi_energies[0] == -0.2237425698
    assert proc.fermi_energies[1] == -0.2237425698


def test_num_eigs_in_each_block() -> None:
    for block in proc.eig_blocks:
        assert len(block.eigenvalues) == 304
        assert len(block.eigenvalues_with_fermi) == 304


def test_num_kpts() -> None:
    # Spin up + spin down -> 2x kpts
    assert len(proc.eig_blocks) == 48


def test_bandgap() -> None:
    assert proc.get_bandgap()[-1] == target_bandgap


def test_read_first_header() -> None:
    assert proc.eig_blocks[0].k_index == 1
    assert np.array_equal(proc.eig_blocks[0].k_vector, np.array([0.06988, -0.14267, -0.22642]))
    assert proc.eig_blocks[0].weight == 0.0416666667


def test_read_second_header() -> None:
    assert proc.eig_blocks[1].k_index == 2
    assert np.array_equal(proc.eig_blocks[1].k_vector, np.array([0.06988, -0.14267, -0.07547]))
    assert proc.eig_blocks[1].weight == 0.0416666667


def test_read_spin_down_header() -> None:
    assert proc.eig_blocks[24].k_index == 1
    assert np.array_equal(proc.eig_blocks[24].k_vector, np.array([0.06988, -0.14267, -0.22642]))
    assert proc.eig_blocks[24].weight == 0.0416666667


def test_read_last_header() -> None:
    assert proc.eig_blocks[-1].k_index == 24
    assert np.array_equal(proc.eig_blocks[-1].k_vector, np.array([0.20963, 0.14267, 0.22642]))
    assert proc.eig_blocks[-1].weight == 0.0416666667
