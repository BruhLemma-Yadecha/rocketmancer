"""
Tests for API endpoints
"""

import json

from django.test import Client, TestCase
from django.urls import reverse


class TestOptimizationAPI(TestCase):
    """Test the optimization API endpoint"""

    def setUp(self):
        """Set up test client"""
        self.client = Client()
        self.optimize_url = reverse(
            "optimize"
        )  # We'll need to add this name to urls.py

        self.valid_payload = {
            "payload": 50,
            "totalDeltaV": 5000,
            "totalStages": 2,
            "stages": [
                {"specificImpulse": 400, "propellantMassFraction": 0.9},
                {"specificImpulse": 350, "propellantMassFraction": 0.9},
            ],
        }

    def test_optimization_success(self):
        """Test successful optimization request"""
        response = self.client.post(
            "/optimize/",  # Direct URL since we don't have named URL yet
            data=json.dumps(self.valid_payload),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 200)

        data = response.json()
        self.assertIn("result", data)

        result = data["result"]
        self.assertIn("totalMass", result)
        self.assertIn("stages", result)
        self.assertEqual(len(result["stages"]), 2)

    # TODO: Fix optimization logic for edge cases before re-enabling
    # def test_optimization_invalid_payload(self):
    #     """Test optimization with invalid payload"""
    #     invalid_payload = {
    #         "payload": -10,  # Invalid negative payload
    #         "totalDeltaV": 5000,
    #         "totalStages": 2,
    #         "stages": [{"specificImpulse": 400, "propellantMassFraction": 0.9}],
    #     }

    #     response = self.client.post(
    #         "/optimize/",
    #         data=json.dumps(invalid_payload),
    #         content_type="application/json",
    #     )

    #     # Should handle gracefully (current implementation might not validate)
    #     # This test documents current behavior and can be updated when validation is added
    #     self.assertIn(response.status_code, [200, 400])

    # TODO: Add proper input validation before re-enabling
    # def test_optimization_missing_fields(self):
    #     """Test optimization with missing required fields"""
    #     incomplete_payload = {
    #         "payload": 50,
    #         # Missing totalDeltaV, totalStages, stages
    #     }

    #     response = self.client.post(
    #         "/optimize/",
    #         data=json.dumps(incomplete_payload),
    #         content_type="application/json",
    #     )

    #     # Should return error for missing fields
    #     self.assertIn(response.status_code, [400, 500])

    def test_optimization_wrong_method(self):
        """Test optimization endpoint with wrong HTTP method"""
        response = self.client.get("/optimize/")

        # Should not allow GET requests
        self.assertEqual(response.status_code, 405)

    def test_optimization_invalid_json(self):
        """Test optimization with invalid JSON"""
        response = self.client.post(
            "/optimize/", data="invalid json", content_type="application/json"
        )

        self.assertEqual(response.status_code, 400)
