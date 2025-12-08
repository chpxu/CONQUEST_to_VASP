from __future__ import annotations
import numpy as np
import numpy.typing as npt
import conquest2a._types as c2at
from conquest2a.conquest import *
import heapq
from typing import Sequence, Any
import copy

class KDNode:
    def __init__(self, index: c2at.INTEGER, axis: c2at.INTEGER) -> None:
        self.index = index
        self.axis = axis


class KDBranch(KDNode):
    def __init__(
        self, index: c2at.INTEGER, axis: c2at.INTEGER, left: KDBranch, right: KDBranch
    ) -> None:
        super().__init__(index=index, axis=axis)
        self.left = left
        self.right = right


class PeriodicKDTree:
    def __init__(self, atoms: list[Atom], box: c2at.REAL_ARRAY) -> None:
        self.atoms = atoms
        self.points = np.array([atom.coords for atom in atoms]) # fractional
        self.box = box
        self.points[:,0] *= box[0]
        self.points[:,1] *= box[1]
        self.points[:,2] *= box[2]
        print(self.points)
        idx = np.arange(len(self.points))
        self.root = self._build_tree(idx, depth=0)

    def _build_tree(self, idx: np.ndarray[tuple[int]], depth: c2at.INTEGER) -> KDBranch | KDNode | None:
        if len(idx) == 0:
            return None
        axis = depth % 3
        if len(idx) == 1:
            return KDNode(index=idx[0], axis=axis)

        sorted_idx = idx[np.argsort(self.points[idx, axis])]
        mid = len(sorted_idx) // 2

        return KDBranch(
            index=sorted_idx[mid],
            axis=axis,
            left=self._build_tree(sorted_idx[:mid], depth + 1),
            right=self._build_tree(sorted_idx[(mid + 1):], depth + 1),
        )

    def squared_distance(
        self, p: np.ndarray[tuple[int]], q: np.ndarray[tuple[int]]
    ) -> c2at.FLOAT:
        d = p - q
        d -= np.rint(d / self.box) * self.box
        return np.dot(d, d)

    def add_to_heap(
        self, d_sq: c2at.FLOAT, idx: c2at.INTEGER, heap: list[Any], k: int | np.integer
    ) -> None:
        if len(heap) < k:
            heapq.heappush(heap, (-d_sq, idx))
        else:
            if d_sq < -heap[0][0]:
                heapq.heapreplace(heap, (-d_sq, idx))

    def search(
        self,
        node: KDBranch | KDNode | None,
        box: c2at.REAL_ARRAY,
        heap: list[Any],
        k: c2at.INTEGER,
        query: Atom,
    ) -> None:
        if node is None or not isinstance(node, KDBranch):
            return
        

        point = self.points[node.index]
        print("Visiting node:", node.index, "at point", point)
        d_sq = self.squared_distance(query.coords, point)
        print(d_sq)
        self.add_to_heap(d_sq=d_sq, idx=node.index, heap=heap, k=k)

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
        
    def knn(
        self, query: Atom, k: c2at.INTEGER
    ) -> Sequence[tuple[int | np.integer, float | np.floating]]:
        heap: list[Any] = []
        if self.root is None:
            raise TypeError("Your root tree is None - unexpected!")
        query_cart = copy.deepcopy(query)
        query_cart.coords[0] *= self.box[0]
        query_cart.coords[1] *= self.box[1]
        query_cart.coords[2] *= self.box[2]

        self.search(node=self.root, box=self.box, heap=heap, k=k, query=query_cart)
        out = []
        # out = [(idx, np.sqrt(d_sq)) for (d_sq, idx) in [(-h[0], h[1]) for h in heap]]
        for neg_d_sq, idx in heap:
            d = np.sqrt(-neg_d_sq)
            out.append((self.atoms[idx], d))

        out.sort(key=lambda x: x[1])
        return out