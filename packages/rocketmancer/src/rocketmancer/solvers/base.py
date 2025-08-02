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
    if payload < 0:
        raise ValueError(f"payload must be ≥ 0, got {payload}")
    if total_delta_v <= 0:
        raise ValueError(f"total_delta_v must be > 0, got {total_delta_v}")

    # Validate per-stage parameters
    if any(isp <= 0 for isp in specific_impulses):
        raise ValueError(
            f"All specific_impulse values must be > 0, got {specific_impulses}"
        )
    if any(pmf <= 0 or pmf >= 1 for pmf in propellant_mass_fractions):
        raise ValueError(
            f"All propellant_mass_fraction values must be in (0,1), got {propellant_mass_fractions}"
        )


@jax.jit
def calculate_gross_mass(
    delta_v_split: jnp.ndarray,
    specific_impulses: jnp.ndarray,
    propellant_mass_fractions: jnp.ndarray,
    payload: float,
) -> float:
    """
    Calculates the Gross Lift-Off Weight (GLOW) of a multi-stage rocket
    for a given delta-v allocation.

    This function is Just-In-Time (JIT) compiled with JAX for maximum performance
    and is fully differentiable, making it suitable for use in gradient-based
    optimization solvers.

    Args:
        delta_v_split: A 1D JAX array (N,) of delta-v values (m/s) for each stage.
        specific_impulses: A 1D JAX array (N,) of specific impulse values (s) for each stage.
        propellant_mass_fractions: A 1D JAX array (N,) of each stage's propellant mass
                                   fraction (φ = m_prop / m_wet).
        payload: The final payload mass (kg) entering the uppermost stage.

    Returns:
        The total initial takeoff mass (GLOW) of the rocket in kg.
    """

    # Define the body of the loop for a single stage.
    # This function will be scanned over all stages.
    def _scan_body(current_payload_from_above, stage_data):
        """
        Calculates the mass of a single stage and adds it to the payload.

        Args:
            current_payload_from_above: The mass carried from the stage above (the 'carry').
            stage_data: A slice containing [dv, isp, pmf] for the current stage.
        """
        stage_dv, stage_isp, stage_pmf = stage_data

        mass_ratio = jnp.exp(stage_dv / (stage_isp * g0))

        stage_wet_mass = (
            current_payload_from_above
            * (mass_ratio - 1)
            / (1 - mass_ratio * (1 - stage_pmf))
        )

        # IMPORTANT: If a stage has 0 dV, its mass is 0. This avoids 0/0 errors.
        stage_wet_mass = jnp.where(stage_dv > 1e-6, stage_wet_mass, 0.0)

        payload_for_next_stage = current_payload_from_above + stage_wet_mass

        return payload_for_next_stage, stage_wet_mass

    stacked_data = jnp.stack(
        [delta_v_split, specific_impulses, propellant_mass_fractions], axis=-1
    )
    reversed_stacked_data = jnp.flip(stacked_data, axis=0)
    final_payload, _ = jax.lax.scan(_scan_body, payload, reversed_stacked_data)

    return final_payload
