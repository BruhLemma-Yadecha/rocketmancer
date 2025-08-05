"""
Rocketmancer: Multi-stage rocket optimization library
"""

from .stage import Stage
from .rocket import Rocket
from .exceptions import SolverError

__version__ = "0.1.0"
__all__ = ["Stage", "Rocket", "SolverError"]
