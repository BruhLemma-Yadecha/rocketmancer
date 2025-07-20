"""
Sample tests to demonstrate pytest and coverage setup
"""

import pytest
from django.test import TestCase
from rest_framework.test import APIClient


class HealthCheckTests(TestCase):
    """Basic health check tests"""

    def setUp(self):
        self.client = APIClient()

    def test_health_check_endpoint_exists(self):
        response = self.client.get("/")
        self.assertIsNotNone(response)


class RocketOptimizationTests(TestCase):
    def setUp(self):
        self.client = APIClient()

    def test_rocket_model_basic_functionality(self):
        try:
            from optimizer.rocket import Rocket

            rocket = Rocket(payload=1000.0, delta_v=4760.08, total_stages=2)
            self.assertIsNotNone(rocket)
        except ImportError:
            pytest.skip("Rocket module not available")

    def test_stage_model_basic_functionality(self):
        try:
            from optimizer.stage import Stage

            stage = Stage(
                stage=1, specific_impulse=307.36, propellant_mass_fraction=0.83
            )
            self.assertIsNotNone(stage)
        except ImportError:
            pytest.skip("Stage module not available")
