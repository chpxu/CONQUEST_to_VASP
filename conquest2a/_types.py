import typing
import numpy as np
import numpy.typing as npt
from pathlib import Path

INTEGER = int | np.int64
FLOAT = float | np.float64
REAL_NUMBER = INTEGER | FLOAT | np.floating[typing.Any]
NUMBER_ARRAY = npt.NDArray[np.number]
REAL_ARRAY = npt.NDArray[np.floating[typing.Any]]
INT_ARRAY = npt.NDArray[np.int64]

FILE_PATH = str | Path
