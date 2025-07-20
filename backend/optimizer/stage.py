from math import exp
from typing import Any

from scipy.constants import g


class Stage:
    def __init__(
        self, stage: int, specific_impulse: float, propellant_mass_fraction: float
    ) -> None:
        self.stage = stage
        self.specific_impulse = specific_impulse
        self.propellant_mass_fraction = propellant_mass_fraction

    def build(self, payload_mass: float, delta_v: float) -> None:
        self.payload_mass = payload_mass
        self.delta_v = delta_v
        self.stage_mass = (self.payload_mass * (1 - self.mass_ratio)) / (
            self.mass_ratio * (1 - self.propellant_mass_fraction) - 1
        )

    @property
    def wet_mass(self) -> float:
        return self.stage_mass + self.payload_mass

    @property
    def dry_mass(self) -> float:
        return self.stage_mass * (1 - self.propellant_mass_fraction) + self.payload_mass

    @property
    def structural_mass(self) -> float:
        return self.stage_mass * (1 - self.propellant_mass_fraction)

    @property
    def propellant_mass(self) -> float:
        return self.stage_mass * self.propellant_mass_fraction

    @property
    def exhaust_velocity(self) -> float:
        return float(self.specific_impulse * g)

    @property
    def mass_ratio(self) -> float:
        return exp(self.delta_v / (self.exhaust_velocity))

    def to_json(self) -> dict[str, Any]:
        return {
            "stage": self.stage,
            "specificImpulse": self.specific_impulse,
            "propellantMassFraction": self.propellant_mass_fraction,
            "payloadMass": self.payload_mass,
            "deltaV": self.delta_v,
            "wetMass": self.wet_mass,
            "dryMass": self.dry_mass,
            "structuralMass": self.structural_mass,
            "propellantMass": self.propellant_mass,
            "exhaustVelocity": self.exhaust_velocity,
            "massRatio": self.mass_ratio,
        }
