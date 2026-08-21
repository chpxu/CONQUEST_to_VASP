from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase.units import Bohr

from conquest2a.read.files import cell_to_conquest
from conquest2a.conquest import Atom, conquest_species


@pytest.fixture
def kcuf3_species() -> conquest_species:
    """Species map with duplicated Cu ids (2, 3) for the spin up/down channels."""
    return conquest_species({1: "K", 2: "Cu", 3: "Cu", 4: "F"})


@pytest.fixture
def simple_species() -> conquest_species:
    """A species map with no duplicated (spin) ids, for non-spin tests."""
    return conquest_species(species_dict={1: "K", 2: "Cu", 3: "F"})


def write_cell(tmp_path: Path, content: str, name: str = "test.cell") -> Path:
    """Write ``content`` to ``tmp_path/name`` and return the path."""
    cell_path = tmp_path / name
    cell_path.write_text(content)
    return cell_path


LATTICE_BOHR_BLOCK = """\
%BLOCK LATTICE_CART
Bohr
  11.259286961753   0.000000000000   0.000000000000
   0.000000000000  11.259286961753   0.000000000000
   0.000000000000   0.000000000000  15.098196958861
%ENDBLOCK LATTICE_CART
"""


# data files located in tests/data/files
DATA_DIR = Path(__file__).parent / "data" / "files"
KCUF3_CELL_PATH = DATA_DIR / "kcuf3.cell"
KCUF3_SPIN_CELL_PATH = DATA_DIR / "kcuf3_spin.cell"


def test_bohr_units_no_conversion(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    cell = write_cell(
        tmp_path,
        LATTICE_BOHR_BLOCK
        + "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n",
    )
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), kcuf3_species, str(dest))

    expected = np.array(
        [
            [11.259286961753, 0.0, 0.0],
            [0.0, 11.259286961753, 0.0],
            [0.0, 0.0, 15.098196958861],
        ]
    )
    np.testing.assert_allclose(converter.coords.lattice_vectors, expected)

def test_ang_to_bohr(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = (
        "%BLOCK LATTICE_CART\n"
        +"Ang\n"
        +"5.0 0.0 0.0\n"
        +"0.0 5.0 0.0\n"
        +"0.0 0.0 5.0\n"
        +"%ENDBLOCK LATTICE_CART\n"
        +"%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), kcuf3_species, str(dest))

    expected = np.diag([5.0, 5.0, 5.0]) / Bohr
    np.testing.assert_allclose(converter.coords.lattice_vectors, expected)

def test_no_units_defaults_bohr(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    """No alphabetic units line is present, so no conversion is applied."""
    content = (
        "%BLOCK LATTICE_CART\n"
        + "10.0 0.0 0.0\n"
        + "0.0 10.0 0.0\n"
        + "0.0 0.0 10.0\n"
        + "%ENDBLOCK LATTICE_CART\n"
        + "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), kcuf3_species, str(dest))

    np.testing.assert_allclose(converter.coords.lattice_vectors, np.diag([10.0, 10.0, 10.0]))

def test_too_few_vectors_raises(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = (
        "%BLOCK LATTICE_CART\n"
        + "Bohr\n"
        + "10.0 0.0 0.0\n"
        + "0.0 10.0 0.0\n"
        + "%ENDBLOCK LATTICE_CART\n"
        + "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(RuntimeError, match="exactly three lattice vectors"):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))

def test_bad_vector_format(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = (
        "%BLOCK LATTICE_CART\n"
        + "Bohr\n"
        + "10.0 0.0\n"
        + "0.0 10.0 0.0\n"
        + "0.0 0.0 10.0\n"
        + "%ENDBLOCK LATTICE_CART\n"
        + "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(RuntimeError, match="three position components"):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))


def test_basic_atom_parsed_correctly(tmp_path: Path, simple_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
         +"K 0.1 0.2 0.3\n"
         +"%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))

    assert len(converter.coords.atoms) == 1
    atom = converter.coords.atoms[0]
    assert isinstance(atom, Atom)
    assert atom.label == "K"
    assert atom.species == 1
    assert atom.number == 0
    np.testing.assert_allclose(atom.coords, [0.1, 0.2, 0.3])
    assert atom.can_move == ["T", "T", "T"]
    np.testing.assert_allclose(atom.spins, [0.0, 0.0, 0.0])

def test_atom_numbers_sequential_and_sorted(tmp_path: Path, simple_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        + "K 0.0 0.0 0.0\n"
        + "Cu 0.5 0.5 0.5\n"
        + "F 0.25 0.25 0.25\n"
        + "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))

    numbers = [atom.number for atom in converter.coords.atoms]
    assert numbers == [0, 1, 2]
    assert converter.coords.natoms == "3"

def test_spin_up_selects_lower_species_id(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0 0.0 SPIN 2\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), kcuf3_species, str(dest))

    atom = converter.coords.atoms[0]
    assert atom.species == 2  # lower of the two Cu ids {2, 3}
    np.testing.assert_allclose(atom.spins, [0.0, 0.0, 2.0])

def test_spin_down_selects_higher_species_id(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0 0.0 SPIN -2\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), kcuf3_species, str(dest))

    atom = converter.coords.atoms[0]
    assert atom.species == 3  # higher of the two Cu ids {2, 3}
    np.testing.assert_allclose(atom.spins, [0.0, 0.0, -2.0])

def test_spin_flag_is_case_insensitive(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0 0.0 spin 2\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), kcuf3_species, str(dest))
    np.testing.assert_allclose(converter.coords.atoms[0].spins, [0.0, 0.0, 2.0])

def test_zero_spin_selects_lower_species_id(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    """spin[2] < 0.0 is the only test used to pick the down channel, so an
    explicit SPIN of exactly 0 should still resolve to the 'up' id."""
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0 0.0 SPIN 0\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), kcuf3_species, str(dest))
    assert converter.coords.atoms[0].species == 2

def test_noncollinear_spin_raises(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0 0.0 SPIN 1.0 1.0 1.0\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(ValueError, match="Non-collinear spin is not supported"):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))

def test_spin_with_two_values_raises(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0 0.0 SPIN 1.0 2.0\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(ValueError, match="1 \\(collinear\\) or 3 \\(non-collinear\\)"):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))

def test_unexpected_trailing_token_raises(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0 0.0 FOO 2\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(ValueError, match="Unexpected trailing content"):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))

def test_malformed_atom_line_raises(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Cu 0.0 0.0\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(ValueError, match="formatted incorrectly"):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))

def test_fake_element_label_raises(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "Xx 0.0 0.0 0.0\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(ValueError, match="not a valid element"):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))

def test_real_element_missing_from_species_map_raises(tmp_path: Path, kcuf3_species: conquest_species) -> None:
    """'He' is a real element, so it passes the allowed_element_labels
    check, but it was never given an id in kcuf3_species, so resolving
    its species id fails with a KeyError."""
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\n"
        "He 0.0 0.0 0.0\n"
        "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    with pytest.raises(KeyError):
        cell_to_conquest(str(cell), kcuf3_species, str(dest))


def test_default_atoms_can_move(tmp_path: Path, simple_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))
    assert converter.coords.atoms[0].can_move == ["T", "T", "T"]

def test_fix_all_ions_true_fixes_atoms(tmp_path: Path, simple_species: conquest_species) -> None:
    content = "FIX_ALL_IONS TRUE\n" + LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))
    assert converter.coords.atoms[0].can_move == ["F", "F", "F"]

def test_fix_all_ions_false_leaves_atoms_free(tmp_path: Path, simple_species: conquest_species) -> None:
    content = "FIX_ALL_IONS FALSE\n" + LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))
    assert converter.coords.atoms[0].can_move == ["T", "T", "T"]

def test_fix_all_ions_only_affects_atoms_parsed_after_it(tmp_path: Path, simple_species: conquest_species) -> None:
    """FIX_ALL_IONS is a simple flag set while scanning top-to-bottom, so
    it must appear before POSITIONS_FRAC to have any effect."""
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    ) + "FIX_ALL_IONS TRUE\n"
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))
    assert converter.coords.atoms[0].can_move == ["T", "T", "T"]

def test_unknown_block_is_skipped(tmp_path: Path, simple_species: conquest_species) -> None:
    content = (
        "%BLOCK SPECIES_MASS\n"
        + "ang\n"
        + "K 1 39.0983\n"
        + "%ENDBLOCK SPECIES_MASS\n"
        + LATTICE_BOHR_BLOCK
        + "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))
    assert len(converter.coords.atoms) == 1
    np.testing.assert_allclose(
        converter.coords.lattice_vectors, np.diag([11.259286961753, 11.259286961753, 15.098196958861])
    )

def test_comment_lines_are_ignored(tmp_path: Path, simple_species: conquest_species) -> None:
    content = (
        "# a top-level comment\n"
        + LATTICE_BOHR_BLOCK
        + "%BLOCK POSITIONS_FRAC\n"
        + "# a comment inside the positions block\n"
        + "K 0.0 0.0 0.0\n"
        + "%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    converter = cell_to_conquest(str(cell), simple_species, str(dest))
    assert len(converter.coords.atoms) == 1

def test_missing_file_raises(tmp_path: Path, simple_species: conquest_species) -> None:
    missing = tmp_path / "does_not_exist.cell"
    dest = tmp_path / "out.dat"
    with pytest.raises(FileNotFoundError):
        cell_to_conquest(str(missing), simple_species, str(dest))

def test_output_file_is_written(tmp_path: Path, simple_species: conquest_species) -> None:
    content = LATTICE_BOHR_BLOCK + (
        "%BLOCK POSITIONS_FRAC\nK 0.0 0.0 0.0\n%ENDBLOCK POSITIONS_FRAC\n"
    )
    cell = write_cell(tmp_path, content)
    dest = tmp_path / "out.dat"
    cell_to_conquest(str(cell), simple_species, str(dest))
    assert dest.exists()
    lines = dest.read_text().splitlines()
    # 3 lattice lines + natoms line + 1 atom line
    assert len(lines) == 5
    assert lines[3] == "1"

tmp_path = DATA_DIR
kcuf3: conquest_species = conquest_species({1: "K", 2: "Cu", 3: "Cu", 4: "F"})
dest = tmp_path / "kcuf3_spin.dat"
celltoconquest = cell_to_conquest(str(KCUF3_SPIN_CELL_PATH), kcuf3, str(dest))

def test_atom_count():
    assert len(celltoconquest.coords.atoms) == 20
    assert celltoconquest.coords.natoms == "20"

def test_element_counts():
    labels = [atom.label for atom in celltoconquest.coords.atoms]
    assert labels.count("K") == 4
    assert labels.count("Cu") == 4
    assert labels.count("F") == 12

def test_cu_spin_channels_split_correctly():
    cu_atoms = [atom for atom in celltoconquest.coords.atoms if atom.label == "Cu"]
    assert len(cu_atoms) == 4
    up_atoms = [atom for atom in cu_atoms if atom.species == 2]
    down_atoms = [atom for atom in cu_atoms if atom.species == 3]
    assert len(up_atoms) == 2
    assert len(down_atoms) == 2
    for atom in up_atoms:
        np.testing.assert_allclose(atom.spins, [0.0, 0.0, 2.0])
    for atom in down_atoms:
        np.testing.assert_allclose(atom.spins, [0.0, 0.0, -2.0])

def test_non_spin_atoms_have_zero_spin():
    for atom in celltoconquest.coords.atoms:
        if atom.label != "Cu":
            np.testing.assert_allclose(atom.spins, [0.0, 0.0, 0.0])

def test_lattice_vectors_unchanged_bohr_units():
    expected = np.diag([11.259286961753, 11.259286961753, 15.098196958861])
    np.testing.assert_allclose(celltoconquest.coords.lattice_vectors, expected)

def test_first_and_last_atom_coordinates():
    first = celltoconquest.coords.atoms[0]
    last = celltoconquest.coords.atoms[-1]
    assert first.label == "K"
    np.testing.assert_allclose(first.coords, [0.99999120, 0.00000044, 0.24988723])
    assert last.label == "F"
    np.testing.assert_allclose(last.coords, [0.49997694, 0.99997336, 0.74993634])

def test_cartesian_positions_populated() -> None:
    # get_cartesian_positions() is called in __init__; every atom should
    # have consistent cart_coords derived from its fractional coords.
    for atom in celltoconquest.coords.atoms:
        expected_cart = atom.coords @ celltoconquest.coords.lattice_vectors.T
        np.testing.assert_allclose(atom.cart_coords, expected_cart)

def test_output_file_line_count() -> None:
    dest_files = list(tmp_path.glob("kcuf3.dat"))
    assert len(dest_files) == 1
    lines = dest_files[0].read_text().splitlines()
    # 3 lattice lines + 1 natoms line + 20 atom lines
    assert len(lines) == 24
    assert lines[3] == "20"

def test_output_file_atom_line_format()-> None:
    lines = (tmp_path / "kcuf3.dat").read_text().splitlines()
    # First atom line corresponds to the first K atom, species 1, free to move.
    first_atom_line = lines[4]
    assert first_atom_line.split()[3] == "1"
    assert first_atom_line.split()[4:7] == ["T", "T", "T"]