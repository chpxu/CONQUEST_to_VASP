import numpy as np
import pytest

from conquest2a.conquest import conquest_species
from conquest2a.io import read_coords


def make_species():
    return conquest_species({1: "K", 2: "Cu", 3: "Cu", 4: "F"})


def write_file(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


def make_reader_unopened(tmp_path, name, content, fmt):
    path = write_file(tmp_path, name, content)
    reader = object.__new__(read_coords)
    reader.path = path
    reader.species = make_species()
    reader.format = fmt
    reader.encoding = "utf-8"
    reader.cq_units = "bohr"
    reader._atom_counter = 0
    reader.fix_ions = False
    reader.abs_input_path = path
    reader.open_file()
    return reader


def test_format_inferred_from_cell_extension(tmp_path):
    path = write_file(tmp_path, "test.cell", "")
    reader = object.__new__(read_coords)
    reader.path = path
    fmt = reader._get_format_from_path()
    assert fmt == "cell"


def test_format_inferred_from_vasp_extension(tmp_path):
    path = write_file(tmp_path, "test.vasp", "")
    reader = object.__new__(read_coords)
    reader.path = path
    fmt = reader._get_format_from_path()
    assert fmt == "vasp"


def test_format_unsupported_extension_raises(tmp_path):
    path = write_file(tmp_path, "test.xyz", "")
    reader = object.__new__(read_coords)
    reader.path = path
    with pytest.raises(ValueError):
        reader._get_format_from_path()


def test_read_percent_blocks_parses_multiple_blocks(tmp_path):
    content = (
        "%block LATTICE_CART\n"
        "ang\n"
        "1.0 0.0 0.0\n"
        "0.0 1.0 0.0\n"
        "0.0 0.0 1.0\n"
        "%ENDBLOCK\n"
        "%block POSITIONS_FRAC\n"
        "K 0.0 0.0 0.0\n"
        "F 0.5 0.5 0.5\n"
        "%endblock\n"
    )
    reader = make_reader_unopened(tmp_path, "test.cell", content, "cell")
    blocks = reader._read_percent_blocks()
    assert "lattice_cart" in blocks
    assert "positions_frac" in blocks
    assert blocks["lattice_cart"] == ["ang", "1.0 0.0 0.0", "0.0 1.0 0.0", "0.0 0.0 1.0"]
    assert blocks["positions_frac"] == ["K 0.0 0.0 0.0", "F 0.5 0.5 0.5"]


def test_read_percent_blocks_endblock_without_name(tmp_path):
    content = "%block LATTICE_ABC\n" "ang\n" "5.0 5.0 5.0\n" "90.0 90.0 90.0\n" "%endblock\n"
    reader = make_reader_unopened(tmp_path, "test.cell", content, "cell")
    blocks = reader._read_percent_blocks()
    assert "lattice_abc" in blocks
    assert len(blocks["lattice_abc"]) == 3


def test_read_percent_blocks_duplicate_name_raises(tmp_path):
    content = (
        "%block POSITIONS_FRAC\n"
        "K 0.0 0.0 0.0\n"
        "%endblock\n"
        "%block POSITIONS_FRAC\n"
        "F 0.5 0.5 0.5\n"
        "%endblock\n"
    )
    reader = make_reader_unopened(tmp_path, "test.cell", content, "cell")
    with pytest.raises(ValueError):
        reader._read_percent_blocks()


def test_read_lattice_vectors_correct_shape(tmp_path):
    reader = make_reader_unopened(tmp_path, "test.cell", "", "cell")
    lines = ["1.0 0.0 0.0", "0.0 2.0 0.0", "0.0 0.0 3.0"]
    vectors = reader._read_lattice_vectors(lines)
    expected = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    assert np.allclose(vectors, expected)


def test_read_lattice_vectors_wrong_count_raises(tmp_path):
    reader = make_reader_unopened(tmp_path, "test.cell", "", "cell")
    lines = ["1.0 0.0 0.0", "0.0 2.0 0.0"]
    with pytest.raises(ValueError):
        reader._read_lattice_vectors(lines)


def test_make_triclinic_lattice_cubic_case(tmp_path):
    reader = make_reader_unopened(tmp_path, "test.cell", "", "cell")
    half_pi = np.pi / 2
    lattice = reader._make_triclinic_lattice(2.0, 2.0, 2.0, half_pi, half_pi, half_pi)
    expected = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
    assert np.allclose(lattice, expected, atol=1e-8)


def test_read_cell_lattice_cart_and_positions_frac(tmp_path):
    content = (
        "%block LATTICE_CART\n"
        "ang\n"
        "10.0 0.0 0.0\n"
        "0.0 10.0 0.0\n"
        "0.0 0.0 10.0\n"
        "%endblock LATTICE_CART\n"
        "%block POSITIONS_FRAC\n"
        "K 0.0 0.0 0.0\n"
        "F 0.5 0.5 0.5\n"
        "%endblock POSITIONS_FRAC\n"
    )
    path = write_file(tmp_path, "kf.cell", content)
    reader = read_coords(path, species=make_species(), format="cell")
    assert len(reader.coords.atoms) == 2
    assert reader.coords.lattice_vectors.shape == (3, 3)


def test_read_cell_lattice_abc_degrees(tmp_path):
    content = (
        "%block LATTICE_ABC\n"
        "ang\n"
        "5.0 5.0 5.0\n"
        "90.0 90.0 90.0\n"
        "%endblock LATTICE_ABC\n"
        "%block POSITIONS_FRAC\n"
        "K 0.0 0.0 0.0\n"
        "%endblock POSITIONS_FRAC\n"
    )
    path = write_file(tmp_path, "cubic.cell", content)
    reader = read_coords(path, species=make_species(), format="cell")
    diag = np.diag(reader.coords.lattice_vectors)
    assert np.allclose(diag, diag[0])
    off_diag = reader.coords.lattice_vectors - np.diag(diag)
    assert np.allclose(off_diag, 0.0, atol=1e-6)


def test_read_cell_missing_positions_block_raises(tmp_path):
    content = (
        "%block LATTICE_CART\n"
        "ang\n"
        "5.0 0.0 0.0\n"
        "0.0 5.0 0.0\n"
        "0.0 0.0 5.0\n"
        "%endblock LATTICE_CART\n"
    )
    path = write_file(tmp_path, "nopositions.cell", content)
    with pytest.raises(Exception):
        read_coords(path, species=make_species(), format="cell")


def test_read_poscar_direct_coordinates(tmp_path):
    content = (
        "comment line\n"
        "1.0\n"
        "5.0 0.0 0.0\n"
        "0.0 5.0 0.0\n"
        "0.0 0.0 5.0\n"
        "K F\n"
        "1 1\n"
        "Direct\n"
        "0.0 0.0 0.0\n"
        "0.5 0.5 0.5\n"
    )
    path = write_file(tmp_path, "POSCAR", content)
    reader = read_coords(path, species=make_species(), format="vasp")
    assert len(reader.coords.atoms) == 2
    assert reader.coords.atoms[0].label == "K"
    assert reader.coords.atoms[1].label == "F"
    assert np.allclose(reader.coords.atoms[1].coords, [0.5, 0.5, 0.5])


def test_read_poscar_negative_scale_raises(tmp_path):
    content = (
        "comment line\n"
        "-1.0\n"
        "5.0 0.0 0.0\n"
        "0.0 5.0 0.0\n"
        "0.0 0.0 5.0\n"
        "K\n"
        "1\n"
        "Direct\n"
        "0.0 0.0 0.0\n"
    )
    path = write_file(tmp_path, "POSCAR", content)
    with pytest.raises(NotImplementedError):
        read_coords(path, species=make_species(), format="vasp")


def test_read_poscar_wrong_scale_count_raises(tmp_path):
    content = (
        "comment line\n"
        "1.0 1.0\n"
        "5.0 0.0 0.0\n"
        "0.0 5.0 0.0\n"
        "0.0 0.0 5.0\n"
        "K\n"
        "1\n"
        "Direct\n"
        "0.0 0.0 0.0\n"
    )
    path = write_file(tmp_path, "POSCAR", content)
    with pytest.raises(ValueError):
        read_coords(path, species=make_species(), format="vasp")


def test_read_dispatches_to_correct_reader(tmp_path, monkeypatch):
    content = "%block LATTICE_CART\nang\n1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%endblock\n%block POSITIONS_FRAC\nK 0.0 0.0 0.0\n%endblock\n"
    path = write_file(tmp_path, "dispatch.cell", content)
    called = []
    reader = object.__new__(read_coords)
    reader._FORMATS = {"cell": "read_cell", "vasp": "read_poscar"}
    reader.format = "cell"
    reader.read_cell = lambda: called.append("cell")
    reader.read_poscar = lambda: called.append("vasp")
    reader.read()
    assert called == ["cell"]


def test_read_unsupported_format_raises():
    reader = object.__new__(read_coords)
    reader._FORMATS = {"cell": "read_cell", "vasp": "read_poscar"}
    reader.format = "xyz"
    with pytest.raises(ValueError):
        reader.read()
