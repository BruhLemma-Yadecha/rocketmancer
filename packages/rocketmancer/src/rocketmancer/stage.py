"""
Stage dataclass with Pydantic validation for rocket stages
"""

from typing import Optional, Tuple
from pydantic import BaseModel, Field, ConfigDict, field_validator


class Stage(BaseModel):
    """
    Represents a single stage of a multi-stage rocket.

    Attributes:
        isp: Specific impulse in seconds (≥0)
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
    isp: float = Field(..., ge=0.0, description="Specific impulse in seconds")
    propellant_mass_fraction: float = Field(
        ..., ge=0.0, le=1.0, description="Propellant mass fraction (0 ≤ φ ≤ 1)"
    )

    # Computed fields (populated by optimization)
    propellant_mass: Optional[float] = Field(None, description="Propellant mass")
    structural_mass: Optional[float] = Field(None, description="Structural mass")
    total_mass: Optional[float] = Field(None, description="Total stage mass")
    delta_v: Optional[float] = Field(None, description="Delta-v contribution")

    @field_validator("propellant_mass", "structural_mass", "total_mass", "delta_v")
    @classmethod
    def validate_computed_fields(cls, v):
        """
        Context-aware validation for computed fields.
        Only enforces non-negative constraints if values are set.
        """
        if v is not None and v < 0:
            raise ValueError("Computed field values must be non-negative")
        return v

    def json(self, **kwargs) -> str:
        """Return JSON representation of the stage"""
        return self.model_dump_json(**kwargs)

    def validate_computed_fields_consistency(self):
        """
        Validate that computed fields are consistently set (all or none).
        Raises ValueError if inconsistencies are found.
        """
        computed_fields = [
            self.propellant_mass,
            self.structural_mass,
            self.total_mass,
            self.delta_v,
        ]

        set_fields = [field for field in computed_fields if field is not None]

        # Either all should be None or all should be set
        if 0 < len(set_fields) < len(computed_fields):
            raise ValueError(
                "Computed fields must be consistently set: either all None or all populated. "
                f"Currently {len(set_fields)} of {len(computed_fields)} fields are set."
            )

    def is_optimized(self) -> bool:
        """
        Check if this stage has been optimized (computed fields populated).
        Also validates that computed fields are consistently set.
        """
        self.validate_computed_fields_consistency()

        return all(
            [
                self.propellant_mass is not None,
                self.structural_mass is not None,
                self.total_mass is not None,
                self.delta_v is not None,
            ]
        )

    def reset_stage_computed_fields(self):
        """Reset all computed fields to None"""
        self.propellant_mass = None
        self.structural_mass = None
        self.total_mass = None
        self.delta_v = None

    def compute_masses(
        self, payload: float, delta_v: float
    ) -> Tuple[float, float, float]:
        """
        Compute total, propellant, and structural mass for this stage.
        Args:
            payload: Mass being lifted by this stage in kg
            delta_v: Delta-v provided by this stage in m/s
        Returns:
            Tuple of (total_mass, propellant_mass, structural_mass)
        """
        import jax.numpy as jnp
        from scipy.constants import g as g0

        r = jnp.exp(delta_v / (self.isp * g0))
        total_mass = payload * (r - 1) / (1 - r * (1 - self.propellant_mass_fraction))
        propellant_mass = total_mass * self.propellant_mass_fraction
        structural_mass = total_mass * (1 - self.propellant_mass_fraction)
        return float(total_mass), float(propellant_mass), float(structural_mass)
