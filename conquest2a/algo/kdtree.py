
import numpy as np
import numpy.typing as npt
import heapq

class KDNode:
    def __init__(self, index, axis, left=None, right=None) -> None:
        self.index = index
        self.axis = axis
        self.left = left
        self.right = right

class PeriodicKDTree:
    def __init__(self, points: npt.NDArray[np.floating], box: npt.ArrayLike) -> None:
        self.points = np.asarray(points, dtype=float)
        self.box = np.asarray(box, dtype=float)
        idx = np.arange(len(points))
        self.root = self._build_tree(idx, depth=0)
    def _build_tree(self, idx: np.ndarray, depth: int | np.integer) -> KDNode:
        if len(idx) == 0:
            raise ValueError("Cannot have a length 0 tree index")
        axis = depth % 3
        sorted_idx = idx[np.argsort(self.points[idx, axis])]
        mid = len(sorted_idx) // 2

        return KDNode(
            index=sorted_idx[mid],
            axis=axis,
            left=self._build_tree(sorted_idx[:mid], depth+1),
            right=self._build_tree(sorted_idx[mid+1:], depth+1)
        )
    def squared_distance(self, p: np.ndarray, q: np.ndarray, box: np.ndarray) -> float | np.floating:
        d = p - q
        d -= np.rint(d / box) * box
        return np.dot(d, d)
    
    def knn(self, query: np.ndarray, k: int | np.integer) -> list:
        query = np.asarray(query, dtype=np.float32)
        box = self.box
        pts = self.points
        heap = []
        
        def add_to_heap(d_sq: float | np.floating, idx: int | np.integer):
            if len(heap) < k:
                heapq.heappush(heap, (-d_sq, idx))
            else:
                if d_sq < -heap[0][0]:
                    heapq.heapreplace(heap, (-d_sq, idx))

        def search(node) -> None:
            if node is None:
                return
            pt = pts[node.index]
            d_sq = self.squared_distance(query, pt, box)
            add_to_heap(d_sq, node.index)

            axis = node.axis
            diff = query[axis] - pt[axis]
            # map onto canonical unit cell
            diff -= np.round(diff / box[axis]) * box[axis]
            if diff < 0:
                first, second = node.left, node.right
            else:
                first, second = node.right, node.left
            search(first)
            plane_dist_sq = diff * diff
            if len(heap) < k or plane_dist_sq < -heap[0][0]:
                search(second)
        search(self.root)
        out = [(idx, np.sqrt(d_sq)) for (d_sq, idx) in [(-h[0], h[1]) for h in heap]]
        out.sort(key=lambda x: x[1])
        return out
