from rest_framework.decorators import api_view
from rest_framework.response import Response

from .rocket import Rocket


@api_view(["POST"])
def optimize(request):
    data = request.data
    rocket0 = Rocket(
        payload=data["payload"],
        delta_v=data["totalDeltaV"],
        total_stages=data["totalStages"],
    )
    for stage in data["stages"]:
        print(stage)
        rocket0.add_stage(
            specific_impulse=float(stage["specificImpulse"]),
            propellant_mass_fraction=float(stage["propellantMassFraction"]),
        )
    rocket0.optimize()
    return Response({"result": rocket0.to_json()})
