"""
Rocket class with optimization capabilities using JAX solver
"""

from typing import List, Tuple, Optional
from pydantic import BaseModel, Field, ConfigDict, validator

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

    @validator("stages")
    def validate_stages(cls, v):
        """Ensure all stages have valid parameters"""
        if not v:
            raise ValueError("At least one stage is required")

        # Additional context-aware validation
        for i, stage in enumerate(v):
            if not isinstance(stage, Stage):
                raise ValueError(f"Stage {i + 1} must be a Stage instance")

        return v

    @validator("payload")
    def validate_payload(cls, v):
        """Validate payload mass"""
        if v < 0:
            raise ValueError("Payload mass must be non-negative")
        return v

    @validator("total_delta_v")
    def validate_total_delta_v(cls, v):
        """Validate total delta-v requirement"""
        if v < 0:
            raise ValueError("Total delta-v must be non-negative")
        return v

    @validator("total_mass")
    def validate_total_mass(cls, v, values):
        """
        Context-aware validation for total mass.
        Only enforces constraints if the rocket has been optimized.
        """
        if v is not None:
            if v < 0:
                raise ValueError("Total mass must be non-negative")
        return v

    def optimize(self, method: str = "greedy", **solver_kwargs) -> Tuple[list, float]:
        """
        Optimize delta-v split to minimize total wet mass.

        Args:
            method: Solver method name ("greedy").
            **solver_kwargs: Additional solver hyperparameters.
        Returns:
            Tuple of (optimal_delta_v_fractions, optimal_total_mass)
        Raises:
            SolverError: If optimization fails to converge
        """
        for stage in self.stages:
            stage.reset_stage_computed_fields()

        solver = get_solver(method, **solver_kwargs)
        specific_impulses = [stage.isp for stage in self.stages]
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
        """
        Populate computed fields for all stages based on optimization results.
        """
        stage_delta_vs = self.total_delta_v * optimal_delta_v_fractions

        # Start with the final payload and work backwards through stages
        current_payload = self.payload

        # Process stages from top to bottom (last stage first)
        for i in reversed(range(len(self.stages))):
            stage = self.stages[i]
            dv_i = float(stage_delta_vs[i])

            # Compute stage masses using the current payload
            stage_total_mass, stage_propellant_mass, stage_structural_mass = (
                stage.compute_masses(current_payload, dv_i)
            )

            # Populate computed fields
            stage.delta_v = dv_i
            stage.total_mass = stage_total_mass
            stage.propellant_mass = stage_propellant_mass
            stage.structural_mass = stage_structural_mass

            # The payload for the next lower stage is the current stage's wet mass plus current payload
            current_payload = current_payload + stage_total_mass

    def validate_optimization_consistency(self):
        """
        Validate that the optimization state is consistent between rocket and stages.
        Raises ValidationError if inconsistencies are found.
        """
        rocket_optimized = self.total_mass is not None

        for i, stage in enumerate(self.stages):
            stage_optimized = stage.is_optimized()

            if rocket_optimized != stage_optimized:
                if rocket_optimized:
                    raise ValueError(
                        f"Stage {i + 1} is not optimized, but rocket has total_mass set. "
                        "All stages must be optimized when rocket is optimized."
                    )
                else:
                    raise ValueError(
                        f"Stage {i + 1} is optimized, but rocket total_mass is not set. "
                        "Rocket must be optimized when stages have computed fields."
                    )

    def is_optimized(self) -> bool:
        """
        Check if the rocket has been optimized.
        Also validates optimization state consistency.
        """
        optimized = self.total_mass is not None and all(
            stage.is_optimized() for stage in self.stages
        )

        # Only validate consistency if we're checking optimization state
        if (
            optimized
            or self.total_mass is not None
            or any(stage.is_optimized() for stage in self.stages)
        ):
            self.validate_optimization_consistency()

        return optimized

    def json(self, **kwargs) -> str:
        """Return JSON representation of the rocket"""
        return self.model_dump_json(**kwargs)

    def reset_optimization(self):
        """Reset all computed fields"""
        self.total_mass = None
        for stage in self.stages:
            stage.reset_stage_computed_fields()
