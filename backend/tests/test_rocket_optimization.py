"""
Tests for rocket optimization functionality
"""

from django.test import TestCase
from rocketmancer import Rocket, Stage


class TestRocketOptimization(TestCase):
    """Test rocket optimization algorithms"""

    def setUp(self):
        """Set up test fixtures"""
        self.test_payload = 50.0
        self.test_delta_v = 5000.0

        # Create test stages
        self.stage1 = Stage(isp=400.0, propellant_mass_fraction=0.9)
        self.stage2 = Stage(isp=350.0, propellant_mass_fraction=0.9)
        self.test_stages = [self.stage1, self.stage2]

    def test_rocket_creation(self):
        """Test basic rocket creation"""
        rocket = Rocket(
            payload=self.test_payload,
            total_delta_v=self.test_delta_v,
            stages=self.test_stages,
        )

        self.assertEqual(rocket.payload, self.test_payload)
        self.assertEqual(rocket.total_delta_v, self.test_delta_v)
        self.assertEqual(len(rocket.stages), 2)

    def test_stage_creation(self):
        """Test basic stage creation"""
        stage = Stage(isp=400.0, propellant_mass_fraction=0.9)

        self.assertEqual(stage.isp, 400.0)
        self.assertEqual(stage.propellant_mass_fraction, 0.9)
        # Computed fields should be None initially
        self.assertIsNone(stage.propellant_mass)
        self.assertIsNone(stage.structural_mass)
        self.assertIsNone(stage.total_mass)
        self.assertIsNone(stage.delta_v)

    def test_optimization_success(self):
        """Test successful optimization"""
        rocket = Rocket(
            payload=self.test_payload,
            total_delta_v=self.test_delta_v,
            stages=self.test_stages,
        )

        delta_v_fractions, total_mass = rocket.optimize()

        # Check that optimization completed
        self.assertIsNotNone(rocket.total_mass)
        self.assertIsNotNone(total_mass)
        self.assertGreater(total_mass, self.test_payload)
        self.assertEqual(rocket.total_mass, total_mass)

        # Check that stages were populated with optimization results
        for stage in rocket.stages:
            self.assertIsNotNone(stage.propellant_mass)
            self.assertIsNotNone(stage.structural_mass)
            self.assertIsNotNone(stage.total_mass)
            self.assertIsNotNone(stage.delta_v)
            # Assert that total_mass is >= 0 (can be 0 if stage is not used)
            assert stage.total_mass is not None
            self.assertGreaterEqual(stage.total_mass, 0.0)

    def test_optimization_with_invalid_stages(self):
        """Test optimization with no stages raises error"""
        with self.assertRaises(ValueError):
            Rocket(
                payload=self.test_payload,
                total_delta_v=self.test_delta_v,
                stages=[],  # Empty stages should raise validation error
            )

    def test_json_serialization(self):
        """Test rocket JSON serialization"""
        rocket = Rocket(
            payload=self.test_payload,
            total_delta_v=self.test_delta_v,
            stages=self.test_stages,
        )
        rocket.optimize()

        json_str = rocket.json()

        # Should return a JSON string
        self.assertIsInstance(json_str, str)

        # Should be valid JSON
        import json

        json_data = json.loads(json_str)

        self.assertIn("payload", json_data)
        self.assertIn("total_delta_v", json_data)
        self.assertIn("stages", json_data)
        self.assertIn("total_mass", json_data)
        self.assertEqual(len(json_data["stages"]), 2)


class TestStage(TestCase):
    """Test individual stage functionality"""

    def test_stage_creation(self):
        """Test stage creation"""
        stage = Stage(isp=400.0, propellant_mass_fraction=0.9)

        self.assertEqual(stage.isp, 400.0)
        self.assertEqual(stage.propellant_mass_fraction, 0.9)

    def test_stage_validation(self):
        """Test stage validation"""
        # Test negative isp
        with self.assertRaises(ValueError):
            Stage(isp=-100.0, propellant_mass_fraction=0.9)

        # Test invalid propellant mass fraction
        with self.assertRaises(ValueError):
            Stage(isp=400.0, propellant_mass_fraction=1.5)

        with self.assertRaises(ValueError):
            Stage(isp=400.0, propellant_mass_fraction=-0.1)

    def test_stage_json_serialization(self):
        """Test stage JSON serialization"""
        stage = Stage(isp=400.0, propellant_mass_fraction=0.9)

        json_str = stage.json()

        # Should return a JSON string
        self.assertIsInstance(json_str, str)

        # Should be valid JSON
        import json

        json_data = json.loads(json_str)

        expected_fields = ["isp", "propellant_mass_fraction"]
        for field in expected_fields:
            self.assertIn(field, json_data)
