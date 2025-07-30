"""
Rocket class with optimization capabilities using JAX solver
"""

from typing import List, Tuple, Optional
from pydantic import BaseModel, Field, ConfigDict, field_validator

from .stage import Stage
from .solvers import get_solver


class Rocket(BaseModel):
    """
    Multi-stage rocket with optimization capabilities.

    Attributes:
        payload: Payload mass in kg (≥0)
        total_delta_v: Total delta-v requirement in m/s (≥0)
        stages: List of Stage instances
        total_mass: Total wet mass after optimization (populated by optimize())
    """

    model_config = ConfigDict(
        validate_assignment=True,
        arbitrary_types_allowed=True,
    )

    payload: float = Field(..., ge=0.0, description="Payload mass")
    total_delta_v: float = Field(..., ge=0.0, description="Total delta-v requirement")
    stages: List[Stage] = Field(..., min_length=1, description="List of rocket stages")
    total_mass: Optional[float] = Field(
        None, description="Total wet mass after optimization"
    )

    @field_validator("stages")
    @classmethod
    def validate_stages(cls, v):
        """Ensure all stages are valid Stage instances"""
        if not v:
            raise ValueError("At least one stage is required")
        return v

    def optimize(
        self, method: str = "lagrangian", **solver_kwargs
    ) -> Tuple[list, float]:
        """
        Optimize delta-v split to minimize total wet mass.

        Args:
            method: Solver method name ("lagrangian").
            **solver_kwargs: Additional solver hyperparameters.
        Returns:
            Tuple of (optimal_delta_v_fractions, optimal_total_mass)
        Raises:
            SolverError: If optimization fails to converge
        """
        for stage in self.stages:
            stage.reset_stage_computed_fields()

        solver = get_solver(method, **solver_kwargs)
        specific_impulses = [stage.specific_impulse for stage in self.stages]
        propellant_mass_fractions = [
            stage.propellant_mass_fraction for stage in self.stages
        ]

        optimal_delta_v_fractions, optimal_total_mass = solver.solve(
            self.payload,
            self.total_delta_v,
            specific_impulses,
            propellant_mass_fractions,
        )

        # Populate computed fields for each stage
        self._populate_stage_results(optimal_delta_v_fractions)
        self.total_mass = optimal_total_mass

        return optimal_delta_v_fractions.tolist(), optimal_total_mass

    def _populate_stage_results(self, optimal_delta_v_fractions):
        """Populate computed fields for all stages based on optimization results."""
        stage_delta_vs = self.total_delta_v * optimal_delta_v_fractions
        current_payload = self.payload

        # Process stages from top to bottom (last stage first)
        for i in reversed(range(len(self.stages))):
            stage = self.stages[i]
            stage.hydrate(current_payload, float(stage_delta_vs[i]))
            current_payload = stage.total_mass + current_payload

    def is_optimized(self) -> bool:
        """Check if the rocket has been optimized."""
        return self.total_mass is not None

    def reset_optimization(self):
        """Reset all computed fields"""
        self.total_mass = None
        for stage in self.stages:
            stage.reset_stage_computed_fields()
