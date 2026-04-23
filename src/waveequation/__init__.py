"""Simple module for numerically solving the 2D wave equation"""

from .PDESolver2d import WaveSolver  # noqa: F401
from importlib.metadata import version

__version__ = version("waveequation")
