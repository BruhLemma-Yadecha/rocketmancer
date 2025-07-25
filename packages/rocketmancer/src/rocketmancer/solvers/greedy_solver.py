from typing import Sequence, Tuple

import jax.numpy as jnp
import numpy as np
from scipy.constants import g as g0

from .base import BaseSolver, validate_solver_inputs


class GreedySolver(BaseSolver):
    """
    Closed-form greedy Δv-split solver for a multi-stage rocket.

    This solver greedily allocates delta-v to stages starting with the highest specific
    impulse (lowest linearized coefficient), subject to per-stage upper bounds.

    Parameters:
        safety_margin: Factor applied to theoretical max delta-v to avoid numerical instability (default: 0.999)
        convergence_tol: Tolerance for considering remaining allocation as zero (default: 1e-15)
    """

    def __init__(
        self, safety_margin: float = 0.999, convergence_tol: float = 1e-15, **kwargs
    ):
        """
        Initialize the greedy solver.

        Args:
            safety_margin: Factor to apply to maximum delta-v to avoid numerical instability (default: 0.999)
            convergence_tol: Tolerance for considering remaining allocation as zero (default: 1e-15)
            **kwargs: Additional parameters for compatibility (ignored)
        """
        if not (0.5 <= safety_margin <= 1.0):
            raise ValueError(
                f"safety_margin must be in [0.5, 1.0], got {safety_margin}"
            )
        if not (1e-20 <= convergence_tol <= 1e-3):
            raise ValueError(
                f"convergence_tol must be in [1e-20, 1e-3], got {convergence_tol}"
            )

        self.safety_margin = safety_margin
        self.convergence_tol = convergence_tol

    def solve(
        self,
        payload: float,
        total_delta_v: float,
        specific_impulses: Sequence[float],
        propellant_mass_fractions: Sequence[float],
    ) -> Tuple[jnp.ndarray, float]:
        """
        Solve for optimal delta-v fractions using greedy allocation.

        Args:
            payload: Mass (kg) entering the uppermost stage. Must be ≥ 0.
            total_delta_v: Required total Δv (m/s). Must be > 0.
            specific_impulses: Sequence of per-stage Isp values (s). All must be > 0.
            propellant_mass_fractions: Sequence of per-stage φᵢ = m_propᵢ/m_wetᵢ. All ∈ (0,1).

        Returns:
            x: jnp.ndarray (N,), the fraction of total_delta_v allocated to each stage (∑xᵢ=1).
            m0: Total wet-mass at liftoff (kg) under this split.

        Raises:
            ValueError: On invalid inputs or if ∑ dv_maxᵢ < total_delta_v (infeasible).
        """
        # Validate inputs using shared validation logic
        validate_solver_inputs(
            payload, total_delta_v, specific_impulses, propellant_mass_fractions
        )

        # Additional validation for greedy solver
        if payload < 0:
            raise ValueError(f"payload must be ≥ 0, got {payload}")
        if total_delta_v <= 0:
            raise ValueError(f"total_delta_v must be > 0, got {total_delta_v}")

        # Convert to numpy arrays for computation
        isps = np.array(specific_impulses, dtype=float)
        phis = np.array(propellant_mass_fractions, dtype=float)

        if isps.ndim != 1 or phis.ndim != 1 or isps.shape != phis.shape:
            raise ValueError(
                "specific_impulses and propellant_mass_fractions must be 1-D of equal length"
            )

        N = isps.size

        # Validate per-stage parameters
        if np.any(isps <= 0):
            raise ValueError(f"All Isp values must be > 0, got {isps}")
        if np.any((phis <= 0) | (phis >= 1)):
            raise ValueError(f"All φ values must be in (0,1), got {phis}")

        # Compute per-stage max Δv and upper bounds on xᵢ
        # dv_maxᵢ = Ispᵢ·g0·ln(1/(1–φᵢ))
        # Apply safety margin to avoid numerical instability at the boundary
        dv_max = (
            isps * g0 * np.log1p(-phis).copy() * -1.0 * self.safety_margin
        )  # log1p(-φ) = ln(1–φ)
        ub = dv_max / total_delta_v

        # Infeasibility check
        if ub.sum() < 1.0:
            tot = ub.sum() * total_delta_v
            raise ValueError(
                f"Infeasible problem: sum(dv_max)={tot:.3f} < total_delta_v={total_delta_v:.3f}"
            )

        # Greedy allocation on linearized coefficients
        # bᵢ = total_delta_v / (Ispᵢ·g0)
        b = total_delta_v / (isps * g0)

        # Sort indices by ascending bᵢ (i.e. highest‐Isp first)
        idx = np.argsort(b)
        x = np.zeros(N, dtype=float)
        remaining = 1.0

        for i in idx:
            take = min(float(ub[i]), remaining)
            x[i] = take
            remaining -= take
            if remaining <= self.convergence_tol:
                break

        # Compute total wet mass: m0 = payload·∏[exp(xᵢ·total_delta_v/(Ispᵢ·g0)) / (1–φᵢ)]
        multipliers = np.exp(x * (total_delta_v / (isps * g0))) / (1.0 - phis)
        m0 = payload * np.prod(multipliers)

        # Convert to JAX arrays for consistency with other solvers
        return jnp.array(x), float(m0)
