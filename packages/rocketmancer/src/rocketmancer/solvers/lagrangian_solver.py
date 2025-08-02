from .base import BaseSolver, calculate_gross_mass, validate_solver_inputs
from typing import Sequence, Tuple
import jax
import jax.numpy as jnp
from scipy.optimize import minimize
import numpy as np
from scipy.constants import g as g0


def _validate_problem_feasibility(
    total_delta_v: float,
    specific_impulses: Sequence[float],
    propellant_mass_fractions: Sequence[float],
) -> None:
    """
    Pre-validate that the required delta-v is theoretically achievable.

    Args:
        total_delta_v: Required total Δv (m/s).
        specific_impulses: Sequence of per-stage Isp values (s).
        propellant_mass_fractions: Sequence of per-stage φ values.

    Raises:
        SolverError: If the required delta-v exceeds theoretical maximum.
    """
    from ..exceptions import SolverError

    # Calculate theoretical maximum delta-v for each stage
    max_possible_dv = sum(
        isp * g0 * np.log(1 / (1 - pmf))
        for isp, pmf in zip(specific_impulses, propellant_mass_fractions)
    )

    if total_delta_v > max_possible_dv:
        raise SolverError(
            f"Required delta-v ({total_delta_v:.1f} m/s) exceeds theoretical maximum "
            f"({max_possible_dv:.1f} m/s) for given stages. Consider adding stages or "
            f"improving stage performance parameters."
        )


def _validate_solution_physical(mass: float, delta_v_split: np.ndarray) -> None:
    """
    Validate that the optimization result is physically meaningful.

    Args:
        mass: Calculated total mass from the solution.
        delta_v_split: The delta-v allocation from optimization.

    Raises:
        SolverError: If the solution is not physically valid.
    """
    from ..exceptions import SolverError

    if mass <= 0:
        raise SolverError(
            "No physically valid solution exists - optimization converged to negative mass. "
            "This may indicate the required delta-v is too high for the given stage parameters."
        )

    if np.any(delta_v_split < 0):
        raise SolverError(
            "Invalid solution: negative delta-v allocation found. "
            "This suggests the optimization constraints may be incompatible."
        )


class LagrangianSolver(BaseSolver):
    """
    Finds the optimal Δv-split for a multi-stage rocket using the method of
    Lagrange Multipliers.

    This solver correctly models the cascading mass of the rocket and finds the
    delta-v distribution that minimizes the total Gross Lift-Off Weight (GLOW).
    """

    def __init__(self, **kwargs):
        """
        Initialize the Lagrangian solver.

        Args:
            **kwargs: Additional parameters for compatibility (ignored).
        """
        super().__init__(**kwargs)

    def solve(
        self,
        payload: float,
        total_delta_v: float,
        specific_impulses: Sequence[float],
        propellant_mass_fractions: Sequence[float],
    ) -> Tuple[jnp.ndarray, float]:
        """
        Solve for optimal delta-v fractions using Lagrangian optimization.

        Args:
            payload: Final payload mass (kg).
            total_delta_v: Required total Δv (m/s).
            specific_impulses: Sequence of per-stage Isp values (s).
            propellant_mass_fractions: Sequence of per-stage φ values.

        Returns:
            x: np.ndarray (N,), the fraction of total_delta_v allocated to each stage.
            m0: The minimum possible total wet-mass at liftoff (kg) under this optimal split.
        """

        validate_solver_inputs(
            payload, total_delta_v, specific_impulses, propellant_mass_fractions
        )

        # Pre-validate that the problem is theoretically solvable
        _validate_problem_feasibility(
            total_delta_v, specific_impulses, propellant_mass_fractions
        )

        specific_impulses = np.array(specific_impulses, dtype=float)
        propellant_mass_fractions = np.array(propellant_mass_fractions, dtype=float)
        N = specific_impulses.size

        # Convert to JAX arrays for automatic differentiation
        jax_specific_impulses = jnp.array(specific_impulses)
        jax_propellant_mass_fractions = jnp.array(propellant_mass_fractions)

        # Define objective function (mass to minimize)
        def objective(dv_split):
            return float(
                calculate_gross_mass(
                    jnp.array(dv_split),
                    jax_specific_impulses,
                    jax_propellant_mass_fractions,
                    payload,
                )
            )

        # Define gradient function using JAX automatic differentiation
        def objective_grad(dv_split):
            grad_fn = jax.grad(
                lambda dv: calculate_gross_mass(
                    dv, jax_specific_impulses, jax_propellant_mass_fractions, payload
                )
            )
            return np.array(grad_fn(jnp.array(dv_split)))

        # Define constraints
        def delta_v_sum_constraint(dv_split):
            """Sum of delta-v must equal total_delta_v"""
            return np.sum(dv_split) - total_delta_v

        def positive_mass_constraint(dv_split):
            """Ensure total mass is positive"""
            mass = objective(dv_split)
            return mass  # Must be > 0

        # Initial guess: even split of delta-v among stages (unbiased)
        initial_dv_split = np.full(N, total_delta_v / N)

        # Set up constraints
        constraints = [
            {"type": "eq", "fun": delta_v_sum_constraint},
            {"type": "ineq", "fun": positive_mass_constraint},
        ]

        # Relaxed bounds: allow flexible allocation while ensuring minimum delta-v per stage
        min_dv_per_stage = 100.0  # Minimum 100 m/s per stage
        max_dv_per_stage = total_delta_v - (N - 1) * min_dv_per_stage
        bounds = [(min_dv_per_stage, max_dv_per_stage) for _ in range(N)]

        # Solve using constrained optimization with gradient information
        result = minimize(
            objective,
            initial_dv_split,
            method="trust-constr",
            jac=objective_grad,  # Provide gradient information
            constraints=constraints,
            bounds=bounds,
            options={"gtol": 1e-8, "xtol": 1e-12, "disp": False, "maxiter": 1000},
        )

        # Check convergence
        if not result.success:
            from ..exceptions import SolverError

            raise SolverError(f"Lagrangian solver failed to converge: {result.message}")

        optimal_dv_split = result.x

        # Calculate the final results based on the optimal split
        m0 = calculate_gross_mass(
            optimal_dv_split, specific_impulses, propellant_mass_fractions, payload
        )

        # Validate that the solution is physically meaningful
        _validate_solution_physical(float(m0), optimal_dv_split)

        x = optimal_dv_split / total_delta_v
        x = jnp.array(x)

        return x, float(m0)
