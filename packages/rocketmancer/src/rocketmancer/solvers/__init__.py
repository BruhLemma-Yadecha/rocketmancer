from .greedy_solver import GreedySolver
from .base import BaseSolver

_SOLVERS = {
    "greedy": GreedySolver,
}


def get_solver(name: str, **kwargs) -> BaseSolver:
    try:
        return _SOLVERS[name](**kwargs)
    except KeyError:
        raise ValueError(f"Unknown solver '{name}'. Valid options: {list(_SOLVERS)}")
