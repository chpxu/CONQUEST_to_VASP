from collections.abc import Sequence
from typing import Any
from scipy.spatial import KDTree
from conquest2a.conquest import conquest_coordinates_processor, Atom


class nearest_neighbours:
    def __init__(
        self,
        conquest_coordinates_processor: conquest_coordinates_processor,
        atom: Atom,
    ) -> None:
        self.coords_proc = conquest_coordinates_processor
        self.kdtree = self.build_kdtree()
        self.atom_to_query = atom

    def build_kdtree(self) -> KDTree:
        return KDTree(
            self.coords_proc.cart_position_vectors,
            copy_data=True,
            boxsize=[self.coords_proc.lattice_vectors[i][i] for i in range(0, 3)],
        )

    def _knn(self, atom_query: Atom, num_neighbours: int | Sequence[int]) -> tuple[Any, Any]:
        """Perform the KDTree query on number_of_neighbours around the specific Atom.

        Args:
            query (Atom): The Atom to find nearest neighbours of
        Returns:
            list[Atom]: The list of nearest-neighbour Atoms
        """
        atom_query_cart_coords = atom_query.coords @ self.coords_proc.lattice_vectors.T
        distances, indices = self.kdtree.query(
            x=atom_query_cart_coords, k=num_neighbours, p=2, workers=-1  # type: ignore
        )  # type: ignore
        return distances, indices

    def get_result(self, num_neighbours: int | Sequence[int]) -> list[tuple[float, Atom]]:
        distances, indices = self._knn(atom_query=self.atom_to_query, num_neighbours=num_neighbours)
        # since k is int | Sequence[int], result is either squeezed or unsqueezed
        # Will handles these separately for the purposes of finding the Atom,
        # But should return the same results if the same parameters are entered
        if isinstance(num_neighbours, int) and num_neighbours == 1:
            # Single atom, 1 neighbour -> single float and single index
            return [(distances, self.coords_proc.atoms[indices])]
        # if isinstance(self.num_neighbours, int):
        # Here, k > 1, or k is a list of ints. Then expect lists
        # Since we are querying for one atom only, do not need to handle lists of lists
        all_atoms: list[Any] = []
        for pair in zip(distances, indices):
            assoc_atom = (pair[0], self.coords_proc.atoms[pair[1]])
            all_atoms.append(assoc_atom)
        return all_atoms
