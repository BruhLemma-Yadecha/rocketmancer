"""
Test fixtures and utilities for rocketmancer library tests.
"""

import pytest
from rocketmancer import Stage


@pytest.fixture
def simple_two_stage_data():
    """Simple two-stage rocket test data for basic validation."""
    payload = 1000.0  # kg
    total_delta_v = 9500.0  # m/s
    stages = [
        Stage(isp=350.0, propellant_mass_fraction=0.85),  # Upper stage
        Stage(isp=280.0, propellant_mass_fraction=0.90),  # Lower stage
    ]
    return payload, total_delta_v, stages


@pytest.fixture
def three_stage_data():
    """Three-stage rocket test data for more complex scenarios."""
    payload = 5000.0  # kg
    total_delta_v = 12000.0  # m/s
    stages = [
        Stage(isp=420.0, propellant_mass_fraction=0.85),  # Upper stage
        Stage(isp=350.0, propellant_mass_fraction=0.88),  # Middle stage
        Stage(isp=280.0, propellant_mass_fraction=0.92),  # Lower stage
    ]
    return payload, total_delta_v, stages


@pytest.fixture(params=["greedy"])
def solver_name(request):
    """Parameterized fixture for all available solver types."""
    return request.param


# Test tolerances and constants
MASS_TOLERANCE = 1e-3  # kg
DELTA_V_TOLERANCE = 1e-6  # m/s
CONVERGENCE_TOLERANCE = 1e-4
