from collections.abc import Sequence
from typing import Any
from scipy.spatial import KDTree
from conquest2a.conquest import conquest_coordinates_processor, Atom


class nearest_neighbours:
    def __init__(
        self,
        number_of_neighbours: int | Sequence[int],
        conquest_coordinates_processor: conquest_coordinates_processor,
        atom: Atom
    ) -> None:
        self.num_neighbours: int | Sequence[int] = number_of_neighbours
        self.coords_proc = conquest_coordinates_processor
        self.kdtree = self.build_kdtree()
        self.distances, self.indices = self._knn(atom_query=atom)
        self.result = self.get_result()
    def build_kdtree(self) -> KDTree:
        return KDTree(
            self.coords_proc.cart_position_vectors,
            copy_data=True,
            boxsize=[self.coords_proc.lattice_vectors[i][i] for i in range(0, 3)],
        )

    def _knn(self, atom_query: Atom) -> tuple[Any, Any]:
        """Perform the KDTree query on number_of_neighbours around the specific Atom.

        Args:
            query (Atom): The Atom to find nearest neighbours of
        Returns:
            list[Atom]: The list of nearest-neighbour Atoms
        """
        atom_query_cart_coords = atom_query.coords @ self.coords_proc.lattice_vectors.T
        distances, indices = self.kdtree.query(
            x=atom_query_cart_coords, k=self.num_neighbours, p=2.0, workers=-1
        ) # type: ignore
        return distances, indices
    def get_result(self) -> tuple[Any, Any]:
        # since k is int | Sequence[int], result is either squeezed or unsqueezed
        # Will handles these separately for the purposes of finding the Atom,
        # But should return the same results if the same parameters are entered
        if isinstance(self.num_neighbours, int) and self.num_neighbours == 1:
            # Single atom, 1 neighbour -> single float and single index
            return self.distances, self.coords_proc.atoms[self.indices]
        # if isinstance(self.num_neighbours, int):
        # Here, k > 1, or k is a list of ints. Then expect lists
        # Since we are querying for one atom only, do not need to handle lists of lists
        all_atoms: list[Any] = []
        for pair in zip(self.distances, self.indices):
            assoc_atom = (pair[0], self.coords_proc.atoms[pair[1]])
            all_atoms.append(assoc_atom)
        return tuple(all_atoms)
