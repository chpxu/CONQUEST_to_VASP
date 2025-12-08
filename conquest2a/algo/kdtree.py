from __future__ import annotations
import numpy as np
import numpy.typing as npt
from conquest2a._types import *
from conquest2a.conquest import *
import heapq
from typing import Sequence, Any


class KDNode:
    def __init__(self, index: INTEGER, axis: INTEGER) -> None:
        self.index = index
        self.axis = axis


class KDBranch(KDNode):
    def __init__(
        self, index: INTEGER, axis: INTEGER, left: KDBranch, right: KDBranch
    ) -> None:
        super().__init__(index=index, axis=axis)
        self.left = left
        self.right = right


class PeriodicKDTree:
    def __init__(self, points: REAL_ARRAY, box: REAL_ARRAY) -> None:
        self.points = points
        self.box = box
        idx = np.arange(len(points))
        self.root = self._build_tree(idx, depth=0)

    def _build_tree(self, idx: np.ndarray[tuple[int]], depth: INTEGER) -> KDBranch | KDNode | None:
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
        self, p: np.ndarray[tuple[int]], q: np.ndarray[tuple[int]], box: npt.NDArray[np.floating]
    ) -> float | np.floating:
        d = p - q
        d -= np.rint(d / box) * box
        return np.dot(d, d)

    def add_to_heap(
        self, d_sq: float | np.floating, idx: int | np.integer, heap: list[Any], k: int | np.integer
    ) -> None:
        if len(heap) < k:
            heapq.heappush(heap, (-d_sq, idx))
        else:
            if d_sq < -heap[0][0]:
                heapq.heapreplace(heap, (-d_sq, idx))

    def search(
        self,
        node: KDBranch,
        points: npt.NDArray,
        box: npt.NDArray[np.floating],
        heap: list[Any],
        k: int | np.integer,
        query: np.ndarray,
    ) -> None:
        if node is None:
            return
        if not isinstance(node, KDBranch):
            return

        point = points[node.index]
        d_sq = self.squared_distance(query, point, box)
        self.add_to_heap(d_sq=d_sq, idx=node.index, heap=heap, k=k)

        axis = node.axis
        diff = query[axis] - point[axis]
        # map onto canonical unit cell
        diff -= np.round(diff / box[axis]) * box[axis]
        if diff < 0:
            first, second = node.left, node.right
        else:
            first, second = node.right, node.left
        self.search(first, points, box, heap, k, query)
        plane_dist_sq = diff * diff
        if len(heap) < k or plane_dist_sq < -heap[0][0]:
            self.search(second, points, box, heap, k, query)
        
    def knn(
        self, query: np.ndarray, k: int | np.integer
    ) -> Sequence[tuple[int | np.integer, float | np.floating]]:
        heap: list[Any] = []
        if self.root is None:
            raise TypeError("Your root tree is None - unexpected!")
        self.search(node=self.root, points=self.points, box=self.box, heap=heap, k=k, query=query)
        out = [(idx, np.sqrt(d_sq)) for (d_sq, idx) in [(-h[0], h[1]) for h in heap]]
        out.sort(key=lambda x: x[1])
        return out


np.random.seed(0)
points = np.random.rand(10000, 3) * 10.0
box = np.array([10.0, 10.0, 10.0])

tree = PeriodicKDTree(points=points, box=box)

q = np.array([1.0, 2.0, 3.0])
print(tree.knn(q, k=3))
