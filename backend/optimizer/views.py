import json

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.request import Request
from rest_framework.response import Response

from rocketmancer import Rocket, SolverError, Stage


@api_view(["POST"])  # type: ignore[misc]
def optimize(request: Request) -> Response:
    """
    Optimize rocket delta-v split using the rocketmancer library.

    Expected request format:
    {
        "payload": float,
        "totalDeltaV": float,
        "stages": [
            {
                "specificImpulse": float,
                "propellantMassFraction": float
            },
            ...
        ]
    }
    """
    try:
        # Check if request has valid JSON data
        if not hasattr(request, "data") or not request.data:
            return Response(
                {"error": "No JSON data provided"}, status=status.HTTP_400_BAD_REQUEST
            )

        data = request.data

        # Create Stage instances using the rocketmancer library
        stages = []
        for stage_data in data["stages"]:
            stage = Stage(
                isp=float(stage_data["specificImpulse"]),
                propellant_mass_fraction=float(stage_data["propellantMassFraction"]),
            )
            stages.append(stage)

        # Create Rocket instance
        rocket = Rocket(
            payload=float(data["payload"]),
            total_delta_v=float(data["totalDeltaV"]),
            stages=stages,
        )

        # Optimize and return JSON - library owns serialization
        delta_v_fractions, total_mass = rocket.optimize()

        # Return the rocket data with result in camelCase format
        result_data = json.loads(rocket.json())

        # Convert snake_case to camelCase for API response
        camel_case_result = {
            "payload": result_data["payload"],
            "totalDeltaV": result_data["total_delta_v"],
            "totalMass": result_data["total_mass"],
            "stages": [],
        }

        for stage in result_data["stages"]:
            camel_stage = {
                "specificImpulse": stage["isp"],
                "propellantMassFraction": stage["propellant_mass_fraction"],
                "propellantMass": stage["propellant_mass"],
                "structuralMass": stage["structural_mass"],
                "totalMass": stage["total_mass"],
                "deltaV": stage["delta_v"],
            }
            camel_case_result["stages"].append(camel_stage)

        response_data = {"result": camel_case_result}
        return Response(response_data, status=status.HTTP_200_OK)

    except SolverError as e:
        return Response(
            {"error": f"Optimization failed: {e!s}"}, status=status.HTTP_400_BAD_REQUEST
        )
    except (KeyError, ValueError, TypeError, json.JSONDecodeError) as e:
        return Response(
            {"error": f"Invalid input data: {e!s}"}, status=status.HTTP_400_BAD_REQUEST
        )
    except Exception as e:
        # Check if it's a JSON parsing error from DRF
        if "JSON parse error" in str(e) or "Expecting value" in str(e):
            return Response(
                {"error": f"Invalid input data: {e!s}"},
                status=status.HTTP_400_BAD_REQUEST,
            )
        return Response(
            {"error": f"Unexpected error: {e!s}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
