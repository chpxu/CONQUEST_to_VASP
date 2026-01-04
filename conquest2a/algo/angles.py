import numpy as np
from conquest2a.conquest import conquest_coordinates_processor, Atom
import conquest2a._types as c2at


class angles:
    def __init__(self, conquest_coords_processor: conquest_coordinates_processor) -> None:
        self.conquest_coords_processor = conquest_coords_processor

    def find_unit_normal_of_three_points(self, atoms: tuple[Atom, Atom, Atom]) -> c2at.REAL_ARRAY:
        d01 = atoms[0].coords - atoms[1].coords
        d02 = atoms[2].coords - atoms[1].coords
        normal_vector: c2at.REAL_ARRAY = np.cross(d01, d02)
        # Coefficients of x ,y, z normalised
        normal_unit_vector: c2at.REAL_ARRAY = normal_vector / np.linalg.norm(normal_vector)
        return normal_unit_vector

    def interior_angle(self, atoms: tuple[Atom, Atom, Atom]) -> c2at.REAL_NUMBER:
        """Find angle between 3 atoms, with the second atom at the centre"""
        d01 = atoms[0].coords - atoms[1].coords
        d02 = atoms[2].coords - atoms[1].coords
        dot = np.dot(d01, d02)
        mag01 = np.linalg.norm(d01)
        mag02 = np.linalg.norm(d02)
        return np.float64(np.arccos(dot / (mag01 * mag02)))

    def dihedral_angle(
        self, normal1: c2at.REAL_ARRAY, normal2: c2at.REAL_ARRAY
    ) -> c2at.REAL_NUMBER:
        return np.float64(np.arccos(np.dot(normal1, normal2)))
