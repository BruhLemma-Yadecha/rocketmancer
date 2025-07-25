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


# Placeholder for more detailed tests
class TestRocketOptimizationResults:
    """Test optimization result correctness - implement in next iteration."""

    def test_optimization_mass_conservation_placeholder(self):
        """Placeholder for mass conservation tests."""
        # TODO: Implement mass conservation validation
        pass
