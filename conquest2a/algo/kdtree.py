from __future__ import annotations
import heapq
import copy
from typing import Sequence, Any
import numpy as np
import conquest2a._types as c2at
from conquest2a.conquest import Atom


class kdnode:
    def __init__(self, index: c2at.INTEGER, axis: c2at.INTEGER) -> None:
        self.index = index
        self.axis = axis


class kdbranch(kdnode):
    def __init__(
        self,
        index: c2at.INTEGER,
        axis: c2at.INTEGER,
        left: kdbranch | kdnode | None,
        right: kdbranch | kdnode | None,
    ) -> None:
        super().__init__(index=index, axis=axis)
        self.left = left
        self.right = right


KD = kdbranch | kdnode | None


class periodic_kdtree:
    def __init__(self, atoms: list[Atom], box: c2at.REAL_ARRAY) -> None:
        self.atoms = atoms
        self.points = np.array([atom.coords for atom in atoms])  # fractional
        self.box = box
        self.points[:, 0] *= box[0]
        self.points[:, 1] *= box[1]
        self.points[:, 2] *= box[2]
        idx = np.arange(len(self.points))
        self.idx = idx
        self.root = self._build_tree(idx, depth=0)

    def _build_tree(self, idx: c2at.INT_ARRAY, depth: c2at.INTEGER) -> KD:
        if len(idx) == 0:
            return None
        axis = depth % 3
        if len(idx) == 1:
            return kdnode(index=idx[0], axis=axis)
        sorted_idx = idx[np.argsort(self.points[idx, axis])]
        mid = len(sorted_idx) // 2
        left_new_branch: KD = self._build_tree(sorted_idx[:mid], depth + 1)
        right_new_branch: KD = self._build_tree(sorted_idx[(mid + 1) :], depth + 1)
        # if left_new_branch is not None and right_new_branch is not None:
        return kdbranch(
            index=sorted_idx[mid], axis=axis, left=left_new_branch, right=right_new_branch
        )

    def squared_distance(self, p: c2at.REAL_ARRAY, q: c2at.REAL_ARRAY) -> c2at.REAL_NUMBER:
        d = p - q
        d -= np.rint(d / self.box) * self.box  # type: ignore
        return np.float64(np.dot(d, d))

    def add_to_heap(
        self, d_sq: c2at.REAL_NUMBER, idx: c2at.INTEGER, heap: list[Any], k: c2at.INTEGER
    ) -> None:
        if len(heap) < k:
            heapq.heappush(heap, (-d_sq, idx))
        else:
            if d_sq < -heap[0][0]:
                heapq.heapreplace(heap, (-d_sq, idx))

    def search(
        self,
        node: kdbranch | kdnode | None,
        box: c2at.REAL_ARRAY,
        heap: list[Any],
        k: c2at.INTEGER,
        query: Atom,
    ) -> None:
        if node is None:
            return
        point = self.points[node.index]
        if not np.allclose(query.coords, point):
            d_sq = self.squared_distance(query.coords, point)
            self.add_to_heap(d_sq=d_sq, idx=node.index, heap=heap, k=k)
            if not isinstance(node, kdbranch):
                return

            axis = node.axis
            diff = query.coords[axis] - point[axis]
            # map onto canonical unit cell
            diff -= np.round(diff / box[axis]) * box[axis]
            if diff < 0:
                first, second = node.left, node.right
            else:
                first, second = node.right, node.left
            self.search(first, box, heap, k, query)
            plane_dist_sq = diff * diff
            if len(heap) < k or plane_dist_sq < -heap[0][0]:
                self.search(second, box, heap, k, query)

    def knn(self, query: Atom, k: c2at.INTEGER) -> Sequence[tuple[Atom, c2at.REAL_NUMBER]]:
        heap: list[Any] = []
        if self.root is None:
            raise TypeError("Your root tree is None - unexpected!")
        if k < 1:
            raise ValueError(
                "The number of nearest neighbours to search for cannot be less than 1."
            )
        query_cart = copy.deepcopy(query)
        query_cart.coords[0] *= self.box[0]
        query_cart.coords[1] *= self.box[1]
        query_cart.coords[2] *= self.box[2]

        self.search(node=self.root, box=self.box, heap=heap, k=k, query=query_cart)
        out = []
        for neg_d_sq, idx in heap:
            d = np.sqrt(-neg_d_sq)
            out.append((self.atoms[idx], d))

        out.sort(key=lambda x: x[1])
        return out
