from io import TextIOWrapper
import os
from os.path import abspath
from pathlib import Path
from typing import IO, Any, TextIO
from dataclasses import dataclass
from writers import *
from conquest import *
bohr_to_angstrom = 0.529177249