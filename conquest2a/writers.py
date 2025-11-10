from io import TextIOWrapper
from pathlib import Path
from typing import IO, Any
from conquest2a.conquest import conquest_coordinates_processor, conquest_coordinates, atom_charge

BOHR_TO_ANGSTROM = 0.529177249
class file_writer:
    def __init__(self, dest: Path | str, mode: str = "w", encoding: str = "utf-8", is_angstrom: bool = False) -> None:
        self.mode = mode
        self.dest_path = dest
        self.encoding = encoding
        self.file = self.open_file()
        self.B2A = BOHR_TO_ANGSTROM if not is_angstrom else 1.0
    def open_file(self) -> TextIOWrapper | IO[Any]:
        file = open(self.dest_path, mode=self.mode, encoding=self.encoding)
        return file

    def close_file(self, file: TextIOWrapper | IO[Any]) -> None:
        file.close()
    
    def write(self) -> None:
        pass

class conquest_writer(file_writer):
    def __init__(self, dest: Path | str, coords: conquest_coordinates, encoding: str = "utf-8", precision: int = 15):
        super().__init__(dest=dest, encoding=encoding)
        self.coords = coords
        self.precision = precision
        if self.precision < 1:
            raise ValueError("Cannot have less than 1 decimal of float precision.")
        self.close_file(file=self.file)
    def write(self) -> None:
        self.file.write(f"{self.coords.lattice_vectors[0][0]:.{self.precision}f} {0.0:.{self.precision}f} {0.0:.{self.precision}f}\n")
        self.file.write(f"{0.0:.{self.precision}f} {self.coords.lattice_vectors[1][1]:.{self.precision}f} {0.0:.{self.precision}f}\n")
        self.file.write(f"{0.0:.{self.precision}f} {0.0:.{self.precision}f} {self.coords.lattice_vectors[2][2]:.{self.precision}f}\n")
        self.file.write(self.coords.natoms)
        for atom in self.coords.Atoms:
            self.file.write(f"{atom.coords[0]:.{self.precision}f} {atom.coords[1]:.{self.precision}f} {atom.coords[2]:.{self.precision}f} {atom.species} {atom.can_move[0]} {atom.can_move[1]} {atom.can_move[2]}")
            self.file.write("\n")
class vasp_writer(file_writer):
    def __init__(
        self,
        dest: Path | str,
        data: conquest_coordinates_processor,
        encoding: str = "utf-8",
        *args, **kwargs
    ) -> None:
        super().__init__(dest=dest, encoding=encoding, *args, **kwargs)
        self.data = data
        self.write()
        self.close_file(file=self.file)

    def create_atoms_str(self) -> tuple[str, str]:
        num_ele = self.data.number_of_elements()
        ele_string = " ".join(list(num_ele.keys()))
        num_string = " ".join(str(x) for x in num_ele.values())
        return ele_string, num_string

    def write(self) -> None:
        ele_string, num_string = self.create_atoms_str()
        with self.file as file:
            file.write(f"{ele_string}\n")
            file.write("1.0\n")
            for lattice_vect in self.data.lattice_vectors:
                file.write(rf'  {" ".join(str(x * self.B2A) for x in lattice_vect)}')
                file.write("\n")
            file.write(f"{ele_string}\n")
            file.write(f"{num_string}\n")
            file.write("Direct\n")
            for atoms in self.data.element_map:
                for atom in self.data.element_map[atoms]:
                    file.write(rf' {" ".join(str(x) for x in atom.coords)}')
                    file.write("\n")


class xyz_writer(file_writer):
    """Class to generate .xyz for basic XYZ file format."""

    def __init__(
        self,
        dest: Path | str,
        data: conquest_coordinates_processor,
        encoding: str = "utf-8",
        comment_line: str = "comment line",
    ) -> None:
        super().__init__(dest=dest, encoding=encoding)
        self.data = data
        self.comment_line = comment_line
        self.write()
        self.close_file(file=self.file)

    def create_comment_line(self) -> str:
        return self.comment_line

    def create_atoms_str(self) -> tuple[str, str]:
        num_ele = self.data.number_of_elements()
        ele_string = " ".join(list(num_ele.keys()))
        num_string = " ".join(str(x) for x in num_ele.values())
        return ele_string, num_string

    def fractional_to_cartesian(self, vector: list[float]) -> list[float]:
        return [
            vector[0] * self.data.lattice_vectors[0][0],
            vector[1] * self.data.lattice_vectors[1][1],
            vector[0] * self.data.lattice_vectors[2][2],
        ]

    def write(self) -> None:
        """XYZ format expects cells in Cartesian coordinates. Since CONQUEST only deals with orthorhombic unit cells, we can just multiply the non-zero coordinates of the lattice vectors with the fractional coordinates"""
        # ele_string, num_string = self.create_atoms_str()
        with self.file as file:
            file.write(f"{self.data.natoms}")
            file.write(f"{self.create_comment_line()}\n")
            for atoms in self.data.element_map:
                for atom in self.data.element_map[atoms]:
                    file.write(
                        rf'{atoms} {" ".join(str(x * self.B2A) for x in self.fractional_to_cartesian(atom.coords))}'
                    )
                    file.write("\n")


class extxyz_writer(xyz_writer):
    def __init__(
        self,
        dest: Path | str,
        data: conquest_coordinates_processor,
        encoding: str = "utf-8",
        time: float = 0.0,
    ) -> None:
        self.time = time
        super().__init__(dest=dest, data=data, encoding=encoding)

    def create_comment_line(self) -> str:
        lattice: list[float] = []
        for single_vector in self.data.lattice_vectors:
            lattice.append(single_vector[0])
            lattice.append(single_vector[1])
            lattice.append(single_vector[2])
        return f'Lattice=\"{" ".join(str(x) for x in lattice)}\" Properties=species:S:1:pos:R:3 Time={str(self.time)}'


class xsf_writer(file_writer):
    def __init__(
        self,
        dest: Path | str,
        data: conquest_coordinates_processor,
        encoding: str = "utf-8",
    ) -> None:
        super().__init__(dest=dest, encoding=encoding)
        self.data = data
        self.write()
        self.close_file(file=self.file)

    def write(self) -> None:
        with self.file as file:
            file.write("CRYSTAL\n")
            file.write("PRIMVEC\n")
            for lattice_vect in self.data.lattice_vectors:
                file.write(rf' {" ".join(str(x * self.B2A) for x in lattice_vect)}')
                file.write("\n")
            file.write("PRIMCOORD\n")
            file.write(f"{self.data.natoms} 1\n")
            for atoms in self.data.element_map:
                for atom in self.data.element_map[atoms]:
                    file.write(rf' {atoms} {" ".join(str(x * self.B2A) for x in atom.coords)}')
                    file.write("\n")
    
class xsf_writer_spins(file_writer):
    def __init__(self, dest: Path | str, charges: atom_charge, xsf_file: Path | str, encoding: str = "utf-8") -> None:
        """
        XSF writer class that includes spin information from the AtomCharge.dat file, by editing an existing XSF file generated by PostProcessCQ.
        """
        self.charges = charges
        self.original_xsf_file = xsf_file
        self.encoding = encoding
        if dest == self.original_xsf_file:
            raise ValueError("Destination path for modified XSF file cannot be the same as the original XSF file path.")
        super().__init__(dest, encoding=self.encoding)
        self.process_xsf_file()
    def process_xsf_file(self) -> None:
        """
        Process an existing XSF file to include spin information.
        """
        header_data: list[str] = []
        lines: list[str] = []
        with open(self.original_xsf_file, mode="r", encoding=self.encoding) as original_xsf:
            header_data = [next(original_xsf) for _ in range(7)]
            lines = original_xsf.readlines()
        original_xsf.close()
        with open(self.dest_path, "w", encoding=self.encoding) as modified_xsf:
            atom_index = 0
            for header_line in header_data:
                modified_xsf.write(header_line)
            for line in lines:
                if line.strip() and atom_index < len(self.charges.coordinates.Atoms):
                    atom = self.charges.coordinates.Atoms[atom_index]
                    spin_info = f"{atom.spins[0]:.10f} {atom.spins[1]:.10f} {atom.spins[2]:.10f}"
                    modified_line = f"{line.strip()} {spin_info}\n"
                    modified_xsf.write(modified_line)
                    atom_index += 1
                else:
                    modified_xsf.write(line)
        modified_xsf.close()

