"""
Skeleton tests for Rocket class functionality.
"""

import pytest
from rocketmancer import Rocket


class TestRocketBasics:
    """Test basic Rocket class functionality."""

    def test_rocket_creation(self, simple_two_stage_data):
        """Test that Rocket objects can be created with valid data."""
        payload, total_delta_v, stages = simple_two_stage_data
        rocket = Rocket(payload=payload, total_delta_v=total_delta_v, stages=stages)

        assert rocket.payload == payload
        assert rocket.total_delta_v == total_delta_v
        assert len(rocket.stages) == len(stages)
        assert rocket.total_mass is None  # Should be None before optimization

    def test_rocket_optimization_with_different_solvers(
        self, solver_name, simple_two_stage_data
    ):
        """Test that rocket optimization works with all solver types."""
        payload, total_delta_v, stages = simple_two_stage_data
        rocket = Rocket(payload=payload, total_delta_v=total_delta_v, stages=stages)

        # Should not raise an error
        result = rocket.optimize(method=solver_name)
        assert isinstance(result, tuple)
        assert len(result) == 2

        # Rocket should now be optimized
        assert rocket.is_optimized()
        assert rocket.total_mass is not None


class TestRocketValidation:
    """Test Rocket validation and error handling."""

    def test_rocket_requires_stages(self):
        """Test that rocket creation requires at least one stage."""
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match="List should have at least 1 item"):
            Rocket(payload=1000.0, total_delta_v=9000.0, stages=[])

    def test_rocket_reset_optimization(self, simple_two_stage_data):
        """Test that optimization can be reset."""
        payload, total_delta_v, stages = simple_two_stage_data
        rocket = Rocket(payload=payload, total_delta_v=total_delta_v, stages=stages)

        # Optimize then reset
        rocket.optimize()
        assert rocket.is_optimized()

        rocket.reset_optimization()
        assert not rocket.is_optimized()
        assert rocket.total_mass is None


class TestRocketOptimizationResults:
    """Test optimization result correctness with reference values."""

    def test_reference_mass_calculation(self):
        """Test optimization against known reference mass."""
        from rocketmancer import Stage

        # Reference case with known optimal mass
        stages = [
            Stage(specific_impulse=350.0, propellant_mass_fraction=0.85),  # Upper stage
            Stage(specific_impulse=280.0, propellant_mass_fraction=0.90),  # Lower stage
        ]

        rocket = Rocket(payload=1000.0, total_delta_v=9500.0, stages=stages)

        # Optimize and check against reference mass
        delta_v_fractions, total_mass = rocket.optimize()

        # Reference mass from external calculation
        expected_mass = 97533.4875

        # Allow for small numerical differences (within 0.1%)
        tolerance = expected_mass * 0.001
        assert abs(total_mass - expected_mass) < tolerance, (
            f"Expected mass {expected_mass:.4f} kg, got {total_mass:.4f} kg "
            f"(difference: {abs(total_mass - expected_mass):.4f} kg, {((total_mass - expected_mass) / expected_mass * 100):.2f} %)"
        )

        # Verify rocket state is consistent
        assert rocket.is_optimized()
        assert rocket.total_mass == total_mass

        # Verify delta-v fractions sum to 1
        assert abs(sum(delta_v_fractions) - 1.0) < 1e-6

        # Verify all stages are optimized
        assert all(stage.is_optimized() for stage in rocket.stages)
