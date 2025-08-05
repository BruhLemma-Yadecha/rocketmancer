import json
from typing import Any

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.request import Request
from rest_framework.response import Response
from rocketmancer import Rocket, SolverError, Stage


@api_view(["POST"])  # type: ignore[misc]
def optimize(request: Request) -> Response:
    """
    Optimize rocket delta-v split using the rocketmancer library.
    """
    try:
        if not hasattr(request, "data") or not request.data:
            return Response(
                {"error": "No JSON data provided"}, status=status.HTTP_400_BAD_REQUEST
            )

        data = request.data

        stages = []
        for stage_data in data["stages"]:
            stage = Stage(
                specific_impulse=float(stage_data["specificImpulse"]),
                propellant_mass_fraction=float(stage_data["propellantMassFraction"]),
            )
            stages.append(stage)

        rocket = Rocket(
            payload=float(data["payload"]),
            total_delta_v=float(data["totalDeltaV"]),
            stages=stages,
        )

        delta_v_fractions, total_mass = rocket.optimize()
        print(delta_v_fractions, total_mass)

        # Return the rocket data with result in camelCase format
        result_data = json.loads(rocket.json())

        def switch_case(s: str) -> str:
            """Convert snake_case to camelCase."""
            parts = s.split("_")
            return parts[0] + "".join(part.capitalize() for part in parts[1:])

        def json_switch_case(data: Any) -> Any:
            """Recursively convert keys in a dictionary from snake_case to camelCase."""
            if isinstance(data, dict):
                return {
                    switch_case(key): json_switch_case(value)
                    for key, value in data.items()
                }
            elif isinstance(data, list):
                return [json_switch_case(item) for item in data]
            else:
                return data

        response_data = {"result": json_switch_case(result_data)}
        return Response(response_data, status=status.HTTP_200_OK)

    except SolverError as e:
        print(f"SolverError: {e!r}")
        return Response(
            {"error": f"Optimization failed: {e!s}"}, status=status.HTTP_400_BAD_REQUEST
        )
    except (KeyError, ValueError, TypeError, json.JSONDecodeError) as e:
        print(f"Input Data Error: {type(e).__name__}: {e!r}")
        return Response(
            {"error": f"Invalid input data: {e!s}"}, status=status.HTTP_400_BAD_REQUEST
        )
    except Exception as e:
        print(f"General Exception: {type(e).__name__}: {e!r}")
        if "JSON parse error" in str(e) or "Expecting value" in str(e):
            print(f"DRF JSON parse error: {e!r}")
            return Response(
                {"error": f"Invalid input data: {e!s}"},
                status=status.HTTP_400_BAD_REQUEST,
            )
        return Response(
            {"error": f"Unexpected error: {e!s}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )
