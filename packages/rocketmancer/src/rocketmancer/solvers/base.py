import jax
import jax.numpy as jnp
from typing import Protocol, Sequence, Tuple
from scipy.constants import g as g0


class BaseSolver(Protocol):
    def solve(
        self,
        payload: float,
        total_delta_v: float,
        specific_impulses: Sequence[float],
        propellant_mass_fractions: Sequence[float],
    ) -> Tuple[jnp.ndarray, float]: ...


def total_mass_objective(
    x: jnp.ndarray,
    payload: float,
    total_delta_v: float,
    specific_impulses: jnp.ndarray,
    propellant_mass_fractions: jnp.ndarray,
) -> jnp.ndarray:
    """
    Shared total mass calculator for all solvers.
    Computes wet mass via the rocket equation with propellant fractions.
    """
    stage_delta_vs = total_delta_v * x

    def stage_fn(carry, args):
        stage_payload, dv_i, isp_i, phi_i = carry, *args
        r_i = jnp.exp(dv_i / (isp_i * g0))
        stage_total_mass = stage_payload * (r_i - 1) / (1 - r_i * (1 - phi_i))
        stage_wet_mass = stage_total_mass + stage_payload
        return stage_wet_mass, None

    m_init, _ = jax.lax.scan(
        stage_fn,
        payload,
        (stage_delta_vs, specific_impulses, propellant_mass_fractions),
    )
    return m_init


def validate_solver_inputs(
    payload: float,
    total_delta_v: float,
    specific_impulses: Sequence[float],
    propellant_mass_fractions: Sequence[float],
) -> None:
    """Shared validation logic for all solvers."""
    n = len(specific_impulses)
    if n < 1 or n != len(propellant_mass_fractions):
        raise ValueError(
            "specific_impulses and propellant_mass_fractions must be same non-zero length"
        )
    if payload < 0 or total_delta_v < 0:
        raise ValueError("payload and total_delta_v must be non-negative")
