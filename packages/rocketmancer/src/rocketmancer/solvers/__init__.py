from .lagrangian_solver import LagrangianSolver
from .base import BaseSolver

_SOLVERS = {
    "lagrangian": LagrangianSolver,
}


def get_solver(name: str, **kwargs) -> BaseSolver:
    try:
        return _SOLVERS[name](**kwargs)
    except KeyError:
        raise ValueError(f"Unknown solver '{name}'. Valid options: {list(_SOLVERS)}")
