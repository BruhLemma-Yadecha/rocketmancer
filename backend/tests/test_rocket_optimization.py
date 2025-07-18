"""
Tests for rocket optimization functionality
"""

import pytest
from django.test import TestCase
from optimizer.rocket import Rocket
from optimizer.stage import Stage


class TestRocketOptimization(TestCase):
    """Test rocket optimization algorithms"""

    def setUp(self):
        """Set up test fixtures"""
        self.test_payload = 50
        self.test_delta_v = 5000
        self.test_stages = 2

    def test_rocket_creation(self):
        """Test basic rocket creation"""
        rocket = Rocket(self.test_payload, self.test_delta_v, self.test_stages)
        
        self.assertEqual(rocket.payload, self.test_payload)
        self.assertEqual(rocket.delta_v, self.test_delta_v)
        self.assertEqual(rocket.total_stages, self.test_stages)
        self.assertEqual(len(rocket.stages), 0)

    def test_stage_addition(self):
        """Test adding stages to rocket"""
        rocket = Rocket(self.test_payload, self.test_delta_v, self.test_stages)
        
        rocket.add_stage(400, 0.9)
        rocket.add_stage(350, 0.9)
        
        self.assertEqual(len(rocket.stages), 2)
        self.assertEqual(rocket.stages[0].specific_impulse, 400)
        self.assertEqual(rocket.stages[1].specific_impulse, 350)

    def test_optimization_success(self):
        """Test successful optimization"""
        rocket = Rocket(self.test_payload, self.test_delta_v, self.test_stages)
        rocket.add_stage(400, 0.9)
        rocket.add_stage(350, 0.9)
        
        rocket.optimize()
        
        # Check that optimization completed
        self.assertIsNotNone(rocket.total_mass)
        self.assertGreater(rocket.total_mass, self.test_payload)
        
        # Check that stages were built
        for stage in rocket.stages:
            self.assertIsNotNone(stage.wet_mass)
            self.assertIsNotNone(stage.dry_mass)
            self.assertGreater(stage.wet_mass, 0)

    def test_optimization_with_invalid_stages(self):
        """Test optimization with no stages raises error"""
        rocket = Rocket(self.test_payload, self.test_delta_v, 0)
        
        with self.assertRaises(ValueError):
            rocket.optimize()

    def test_json_serialization(self):
        """Test rocket JSON serialization"""
        rocket = Rocket(self.test_payload, self.test_delta_v, self.test_stages)
        rocket.add_stage(400, 0.9)
        rocket.add_stage(350, 0.9)
        rocket.optimize()
        
        json_data = rocket.to_json()
        
        self.assertIn('payload', json_data)
        self.assertIn('totalDeltaV', json_data)
        self.assertIn('totalStages', json_data)
        self.assertIn('totalMass', json_data)
        self.assertIn('stages', json_data)
        self.assertEqual(len(json_data['stages']), 2)


class TestStage(TestCase):
    """Test individual stage functionality"""

    def test_stage_creation(self):
        """Test stage creation"""
        stage = Stage(0, 400, 0.9)
        
        self.assertEqual(stage.stage, 0)
        self.assertEqual(stage.specific_impulse, 400)
        self.assertEqual(stage.propellant_mass_fraction, 0.9)

    def test_stage_build(self):
        """Test stage building with payload and delta-V"""
        stage = Stage(0, 400, 0.9)
        stage.build(100, 2500)
        
        self.assertEqual(stage.payload_mass, 100)
        self.assertEqual(stage.delta_v, 2500)
        self.assertGreater(stage.wet_mass, stage.payload_mass)
        self.assertGreater(stage.wet_mass, stage.dry_mass)

    def test_stage_properties(self):
        """Test stage calculated properties"""
        stage = Stage(0, 400, 0.9)
        stage.build(100, 2500)
        
        # Test that all properties are calculated
        self.assertIsNotNone(stage.wet_mass)
        self.assertIsNotNone(stage.dry_mass)
        self.assertIsNotNone(stage.structural_mass)
        self.assertIsNotNone(stage.propellant_mass)
        self.assertIsNotNone(stage.exhaust_velocity)
        self.assertIsNotNone(stage.mass_ratio)
        
        # Test mass relationships
        self.assertAlmostEqual(
            stage.wet_mass,
            stage.dry_mass + stage.propellant_mass,
            places=5
        )

    def test_stage_json_serialization(self):
        """Test stage JSON serialization"""
        stage = Stage(0, 400, 0.9)
        stage.build(100, 2500)
        
        json_data = stage.to_json()
        
        expected_fields = [
            'stage', 'specificImpulse', 'propellantMassFraction',
            'payloadMass', 'deltaV', 'wetMass', 'dryMass',
            'structuralMass', 'propellantMass', 'exhaustVelocity', 'massRatio'
        ]
        
        for field in expected_fields:
            self.assertIn(field, json_data)