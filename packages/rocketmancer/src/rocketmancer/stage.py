"""
Stage dataclass with Pydantic validation for rocket stages
"""

from typing import Optional
from pydantic import BaseModel, Field, ConfigDict
import jax.numpy as jnp
from scipy.constants import g as g0


class Stage(BaseModel):
    """
    Represents a single stage of a multi-stage rocket.

    Attributes:
        specific_impulse: Specific impulse in seconds (≥0)
        propellant_mass_fraction: Ratio of propellant mass to total stage mass (0 ≤ φ ≤ 1)
        propellant_mass: Propellant mass (populated by optimization)
        structural_mass: Structural mass (populated by optimization)
        total_mass: Total stage mass (populated by optimization)
        delta_v: Delta-v contribution (populated by optimization)
    """

    model_config = ConfigDict(
        validate_assignment=True,
        arbitrary_types_allowed=True,
    )

    # Input parameters (validated at creation)
    specific_impulse: float = Field(
        ..., ge=0.0, description="Specific impulse in seconds"
    )
    propellant_mass_fraction: float = Field(
        ..., ge=0.0, le=1.0, description="Propellant mass fraction (0 ≤ φ ≤ 1)"
    )

    # Computed fields (populated by optimization)
    payload_mass: Optional[float] = Field(None, description="Payload mass")
    propellant_mass: Optional[float] = Field(None, description="Propellant mass")
    structural_mass: Optional[float] = Field(None, description="Structural mass")
    total_mass: Optional[float] = Field(None, description="Total stage mass")
    delta_v: Optional[float] = Field(None, description="Delta-v contribution")
    mass_ratio: Optional[float] = Field(None, description="Mass ratio")

    def is_optimized(self) -> bool:
        """Check if this stage has been optimized (computed fields populated)."""
        return self.total_mass is not None

    def reset_stage_computed_fields(self):
        """Reset all computed fields to None"""
        self.payload_mass = None
        self.propellant_mass = None
        self.structural_mass = None
        self.total_mass = None
        self.delta_v = None
        self.mass_ratio = None

    def hydrate(self, payload: float, delta_v: float) -> None:
        """Populate computed fields based on payload and delta-v."""
        self.payload_mass = payload
        self.delta_v = delta_v
        self.mass_ratio = float(jnp.exp(delta_v / (self.specific_impulse * g0)))
        self.total_mass = (
            payload
            * (self.mass_ratio - 1)
            / (1 - self.mass_ratio * (1 - self.propellant_mass_fraction))
        )
        self.propellant_mass = self.total_mass * self.propellant_mass_fraction
        self.structural_mass = self.total_mass * (1 - self.propellant_mass_fraction)
