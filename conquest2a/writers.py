from io import TextIOWrapper
from typing import IO, Any, Literal, override
from conquest2a.conquest import Atom, conquest_coordinates, atom_charge
from conquest2a.constants import BOHR_TO_ANGSTROM


class file_writer:
    """Generic parent class to define file operations and variables.

    :param dest: _description_
    :type dest: str
    :param mode: _description_, defaults to "w"
    :type mode: str, optional
    :param encoding: _description_, defaults to "utf-8"
    :type encoding: str, optional
    :param is_angstrom: _description_, defaults to False
    :type is_angstrom: bool, optional
    """

    def __init__(
        self, dest: str, mode: str = "w", encoding: str = "utf-8", is_angstrom: bool = False
    ) -> None:
        self.mode: str = mode
        self.dest_path: str = dest.strip()
        self.encoding: str = encoding
        self.is_ang: bool = is_angstrom
        self.file: IO[Any] = self.open_file()

    def open_file(self) -> IO[Any]:
        file: IO[Any] = open(self.dest_path, mode=self.mode, encoding=self.encoding)
        return file

    def close_file(self, file: TextIOWrapper | IO[Any]) -> None:
        file.close()

    def write(self) -> None:
        pass


class conquest_writer(file_writer):
    """Class to write a CONQUEST coordinates file given a :class:`conquest_coordinates` instance.

    :param dest: Path to write the new coordinates file.
    :type dest: ``str``
    :param coords: :class:`conquest_coordinates` instance to write.
    :type coords: ``conquest_coordinates``
    :param encoding: File encoding, defaults to "utf-8"
    :type encoding: ``str``, optional
    :param precision: Float precision, defaults to 10
    :type precision: ``int``, optional
    :raises ValueError: If the float precision is not at least 1.
    """

    def __init__(
        self,
        dest: str,
        coords: conquest_coordinates,
        encoding: str = "utf-8",
        precision: int = 10,
    ) -> None:
        super().__init__(dest=dest, encoding=encoding)
        self.coords: conquest_coordinates = coords
        self.precision: int = precision
        if self.precision < 1:
            raise ValueError("Cannot have less than 1 decimal of float precision.")
        self.write()
        self.close_file(file=self.file)

    @override
    def write(self) -> None:
        prec = self.precision
        self.file.write(
            f"{self.coords.lattice_vectors[0][0]:.{prec}f} {0.0:.{prec}f} {0.0:.{prec}f}\n"
        )
        self.file.write(
            f"{0.0:.{prec}f} {self.coords.lattice_vectors[1][1]:.{prec}f} {0.0:.{prec}f}\n"
        )
        self.file.write(
            f"{0.0:.{prec}f} {0.0:.{prec}f} {self.coords.lattice_vectors[2][2]:.{prec}f}\n"
        )
        self.file.write(self.coords.natoms)
        self.file.write("\n")
        for atom in self.coords.atoms:
            move_str: str = f"{atom.can_move[0]} {atom.can_move[1]} {atom.can_move[2]}"
            atom_str: str = (
                f"{atom.coords[0]:.{prec}f} {atom.coords[1]:.{prec}f} {atom.coords[2]:.{prec}f}"
            )
            self.file.write(f"{atom_str} {atom.species} {move_str}")
            self.file.write("\n")


class vasp_writer(file_writer):
    """Class to write a VASP file given a :class:`~conquest.conquest_coordinates` instance.

    :param dest: Path to write the new coordinates file.
    :type dest: ``str``
    :param data: :class:`~conquest.conquest_coordinates` instance to write.
    :type data: ``conquest_coordinates``
    :param encoding: File encoding, defaults to "utf-8"
    :type encoding: ``str``, optional
    :param is_angstrom: Whether the data in ``conquest_coordinates`` is already in angstroms instead of Bohrs, defaults to ``False``.
    :type is_angstrom: ``bool``, optional
    """

    def __init__(
        self,
        dest: str,
        data: conquest_coordinates,
        encoding: str = "utf-8",
        is_angstrom: bool = False,
    ) -> None:
        super().__init__(dest=dest, encoding=encoding, is_angstrom=is_angstrom)
        self.data: conquest_coordinates = data
        self.write()
        self.close_file(file=self.file)

    def create_atoms_str(self) -> tuple[str, str]:
        num_ele: dict[str, int] = self.data.number_of_elements()
        ele_string: str = " ".join(list(num_ele.keys()))
        num_string: str = " ".join(str(x) for x in num_ele.values())
        return ele_string, num_string

    @override
    def write(self) -> None:
        ele_string, num_string = self.create_atoms_str()
        with self.file as file:
            file.write(f"{ele_string}\n")
            file.write("1.0\n")
            if self.is_ang:
                for lattice_vect in self.data.lattice_vectors:
                    file.write(rf'  {" ".join(str(x) for x in lattice_vect)}')
                    file.write("\n")
            else:
                for lattice_vect in self.data.lattice_vectors:
                    file.write(rf'  {" ".join(str(x * BOHR_TO_ANGSTROM) for x in lattice_vect)}')
                    file.write("\n")
            file.write(f"{ele_string}\n")
            file.write(f"{num_string}\n")
            file.write("Direct\n")
            for atoms in self.data.element_map:
                for atom in self.data.element_map[atoms]:
                    file.write(rf' {" ".join(str(x) for x in atom.coords)}')
                    file.write("\n")


class xyz_writer(file_writer):
    """Class to write a ``.xyz`` for a basic XYZ file given a :class:`~conquest.conquest_coordinates` instance.

    :param dest: Path to write the new coordinates file.
    :type dest: ``str``
    :param data: :class:`~conquest.conquest_coordinates` instance to write.
    :type data: ``conquest_coordinates``
    :param encoding: File encoding, defaults to "utf-8"
    :type encoding: ``str``, optional
    :param comment_line: What string to write as the comment line
    :type comment_line: ``str``, optional
    """

    def __init__(
        self,
        dest: str,
        data: conquest_coordinates,
        encoding: str = "utf-8",
        comment_line: str = "comment line",
    ) -> None:
        super().__init__(dest=dest, encoding=encoding)
        self.data: conquest_coordinates = data
        self.comment_line: str = comment_line
        self.write()
        self.close_file(file=self.file)

    def create_comment_line(self) -> str:
        return self.comment_line

    def create_atoms_str(self) -> tuple[str, str]:
        """Creates a line for an atom, for the XYZ file.

        :return: Components of the line to write: the element and the numbers
        :rtype: tuple[str, str]
        """
        num_ele: dict[str, int] = self.data.number_of_elements()
        ele_string: str = " ".join(list(num_ele.keys()))
        num_string: str = " ".join(str(x) for x in num_ele.values())
        return ele_string, num_string

    @override
    def write(self) -> None:
        """XYZ format expects cells in Cartesian coordinates.
        Since CONQUEST only deals with orthorhombic unit cells,
        we can just multiply the non-zero coordinates of the lattice vectors
        with the fractional coordinates"""
        # ele_string, num_string = self.create_atoms_str()
        with self.file as file:
            file.write(f"{self.data.natoms}")
            file.write(f"{self.create_comment_line()}\n")
            for atoms in self.data.element_map:
                for atom in self.data.element_map[atoms]:
                    file.write(
                        rf'{atoms} {" ".join(str(x * BOHR_TO_ANGSTROM) for x in atom.cart_coords)}'
                    )
                    file.write("\n")


class extxyz_writer(xyz_writer):
    """Class to write a ``.extxyz`` for a basic XYZ file given a :class:`~conquest.conquest_coordinates` instance.

    The main advantage of `.extxyz` is the ability to specify columns and the time, which is very useful for animations. I recommend just using CONQUEST's ability to output ``.extxyz`` files at different timesteps.

    :param dest: Path to write the new coordinates file.
    :type dest: ``str``
    :param data: :class:`~conquest.conquest_coordinates` instance to write.
    :type data: ``conquest_coordinates``
    :param encoding: File encoding, defaults to "utf-8".
    :type encoding: ``str``, optional
    :param time: The instance of time of the cell, defaults to `0.0`.
    :type time: ``float``, optional
    """

    def __init__(
        self,
        dest: str,
        data: conquest_coordinates,
        encoding: str = "utf-8",
        time: float = 0.0,
    ) -> None:

        self.time: float = time
        super().__init__(dest=dest, data=data, encoding=encoding)

    @override
    def create_comment_line(self) -> str:
        """Creates the ``.extxyz`` comment line: specifies columns, formats and time.

        :return: The file's comment line.
        :rtype: ``str``
        """
        lattice: list[float] = []
        for single_vector in self.data.lattice_vectors:
            lattice.append(single_vector[0])
            lattice.append(single_vector[1])
            lattice.append(single_vector[2])
        property_str = "Properties=species:S:1:pos:R:3"
        time_str: str = f"Time={str(self.time)}"
        lat: str = " ".join(str(x) for x in lattice)
        return f'Lattice="{lat}" {property_str} {time_str}'


class xsf_writer(file_writer):
    """Class to write `.xsf` files given a :class:`~conquest.conquest_coordinates` instance.

    XSF files support an extra 3 columns alongside the 3 columns used for position. These columns specify a vector associated with each atom.
    In CONQUEST, the relevant vectors are forces and spins.
    Note that CONQUEST only supports *collinear spin* (up and down), so the spin vector is visualised along the :math:`c`-axis of the simulation cell.

    :param dest: Path to write the new coordinates file.
    :type dest: ``str``
    :param data: :class:`~conquest.conquest_coordinates` instance to write.
    :type data: ``conquest_coordinates``
    :param encoding: File encoding, defaults to "utf-8".
    :type encoding: ``str``, optional
    :param write_extra: Whether to extract the force vector or spin vector of atoms, defaults to "spin".
    :type write_extra: Literal["spin", "force"], optional
    """

    def __init__(
        self,
        dest: str,
        data: conquest_coordinates,
        write_extra: Literal["spin", "force"] = "spin",
        encoding: str = "utf-8",
    ) -> None:
        super().__init__(dest=dest, encoding=encoding)
        self.data: conquest_coordinates = data
        self.write_extra: Literal["spin", "force"] = write_extra
        self.write()
        self.close_file(file=self.file)

    def _format_extra(self, atom: Atom) -> str:
        if self.write_extra == "spin":
            return " ".join(str(x) for x in atom.spins)
        elif self.write_extra == "force":
            return " ".join(str(x) for x in atom.forces)
        return ""

    @override
    def write(self) -> None:
        with self.file as file:
            file.write("CRYSTAL\n")
            file.write("PRIMVEC\n")
            for lattice_vect in self.data.lattice_vectors:
                file.write(rf' {" ".join(str(x * BOHR_TO_ANGSTROM) for x in lattice_vect)}')
                file.write("\n")
            file.write("PRIMCOORD\n")
            natom_line: str = f'{" ".join(self.data.natoms.split())} 1\n'
            file.write(natom_line)
            for element, atoms in self.data.element_map.items():
                for atom in atoms:
                    pos_string: str = " ".join(str(x * BOHR_TO_ANGSTROM) for x in atom.cart_coords)
                    extra: str = self._format_extra(atom)
                    file.write(f" {element} {pos_string} {extra}\n")


class xsf_writer_spins(file_writer):
    """
    XSF writer class that includes spin information from the AtomCharge.dat file,
    by editing an existing XSF file generated by PostProcessCQ.

    :param dest: Path to write the new coordinates file.
    :type dest: ``str``
    :param data: :class:`~conquest.conquest_coordinates` instance to write.
    :type data: ``conquest_coordinates``
    :param encoding: File encoding, defaults to "utf-8".
    :type encoding: ``str``, optional
    :param charges: :class:`~conquest.atom_charge` instance to read spin data from.
    :type charges: :class:`~conquest.atom_charge`
     :param xsf_file: XSF file to modify.
    :type xsf_file: ``str``
    """

    def __init__(
        self, dest: str, charges: atom_charge, xsf_file: str, encoding: str = "utf-8"
    ) -> None:
        self.charges: atom_charge = charges
        self.original_xsf_file: str = xsf_file.strip()
        self.encoding: str = encoding
        if dest == self.original_xsf_file:
            raise ValueError(
                "Destination path for modified XSF file "
                + "cannot be the same as the original XSF file path"
            )
        super().__init__(dest=dest, encoding=self.encoding)
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
                if line.strip() and atom_index < len(self.charges.coordinates.atoms):
                    atom: Atom = self.charges.coordinates.atoms[atom_index]
                    spin_info: str = (
                        f"{atom.spins[0]:.10f} {atom.spins[1]:.10f} {atom.spins[2]:.10f}"
                    )
                    modified_line: str = f"{line.strip()} {spin_info}\n"
                    modified_xsf.write(modified_line)
                    atom_index += 1
                else:
                    modified_xsf.write(line)
        modified_xsf.close()
