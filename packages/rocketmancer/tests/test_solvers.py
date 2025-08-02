"""
Tests for the solver factory and base functionality.
"""

import pytest
from rocketmancer.solvers import get_solver


class TestSolverFactory:
    """Test the solver factory and registry."""

    def test_get_solver_returns_correct_types(self, solver_name):
        """Test that get_solver returns objects implementing BaseSolver."""
        solver = get_solver(solver_name)
        assert hasattr(solver, "solve"), f"{solver_name} solver missing solve method"
        # Note: Runtime isinstance check with Protocol requires special handling

    def test_get_solver_invalid_name_raises_error(self):
        """Test that invalid solver names raise appropriate errors."""
        with pytest.raises(ValueError, match="Unknown solver 'invalid'"):
            get_solver("invalid")

    def test_get_solver_with_kwargs(self, solver_name):
        """Test that solver-specific kwargs are passed correctly."""
        # Should not raise an error with common hyperparameters
        solver = get_solver(solver_name, maxiter=100, tol=1e-4)
        assert solver is not None


class TestSolverInterface:
    """Test that all solvers implement the required interface."""

    def test_solver_solve_method_signature(self, solver_name, simple_two_stage_data):
        """Test that solve method accepts correct arguments and returns correct types."""
        payload, total_delta_v, stages = simple_two_stage_data
        specific_impulses = [stage.specific_impulse for stage in stages]
        propellant_mass_fractions = [stage.propellant_mass_fraction for stage in stages]

        solver = get_solver(solver_name)
        result = solver.solve(
            payload, total_delta_v, specific_impulses, propellant_mass_fractions
        )

        # Check return type structure
        assert isinstance(result, tuple), "solve should return a tuple"
        assert len(result) == 2, "solve should return (delta_v_fractions, total_mass)"

        delta_v_fractions, total_mass = result
        assert hasattr(delta_v_fractions, "tolist"), (
            "delta_v_fractions should be JAX array"
        )
        assert isinstance(total_mass, float), "total_mass should be float"

    def test_solver_input_validation(self, solver_name):
        """Test that all solvers validate inputs consistently."""
        solver = get_solver(solver_name)

        # Test negative payload
        with pytest.raises(ValueError, match="payload must be ≥ 0"):
            solver.solve(-100.0, 9000.0, [350.0], [0.85])

        # Test mismatched array lengths
        with pytest.raises(ValueError, match="must be same non-zero length"):
            solver.solve(1000.0, 9000.0, [350.0, 280.0], [0.85])
