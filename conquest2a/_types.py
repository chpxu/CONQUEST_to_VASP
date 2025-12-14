import typing
import numpy as np
import numpy.typing as npt

E = typing.TypeVar("E", bound=np.integer)
INTEGER = int | np.integer
FLOAT = float | np.floating
REAL_NUMBER = INTEGER | FLOAT
NUMBER_ARRAY = npt.NDArray[np.number]
REAL_ARRAY = npt.NDArray[np.floating | np.integer]
INT_ARRAY = npt.NDArray[np.integer]
