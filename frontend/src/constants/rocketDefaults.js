/**
 * Default values and constants for rocket configurations
 */

// Default rocket configuration
export const DEFAULT_ROCKET_CONFIG = {
  payload: 50,
  totalDeltaV: 5000,
  totalStages: 2,
  stages: [
    {
      specificImpulse: 400,
      propellantMassFraction: 0.9
    },
    {
      specificImpulse: 350,
      propellantMassFraction: 0.9
    }
  ]
};

// Common specific impulse values for different propellant types
export const COMMON_ISP_VALUES = {
  SOLID: { min: 200, max: 300, typical: 250 },
  HYPERGOLIC: { min: 280, max: 340, typical: 310 },
  KEROSENE_LOX: { min: 300, max: 370, typical: 340 },
  HYDROGEN_LOX: { min: 380, max: 465, typical: 450 },
  METHANE_LOX: { min: 330, max: 380, typical: 355 }
};

// Typical mass fraction ranges
export const TYPICAL_MASS_FRACTIONS = {
  SOLID: { min: 0.85, max: 0.95, typical: 0.90 },
  LIQUID: { min: 0.80, max: 0.95, typical: 0.88 },
  PRESSURE_FED: { min: 0.75, max: 0.90, typical: 0.82 },
  PUMP_FED: { min: 0.85, max: 0.95, typical: 0.90 }
};

// Mission profiles with typical delta-V requirements
export const MISSION_PROFILES = {
  LEO: { deltaV: 9400, description: 'Low Earth Orbit' },
  GTO: { deltaV: 12500, description: 'Geostationary Transfer Orbit' },
  LUNAR: { deltaV: 15000, description: 'Lunar Transfer' },
  MARS: { deltaV: 16000, description: 'Mars Transfer' },
  ESCAPE: { deltaV: 17000, description: 'Solar System Escape' }
};

// Validation limits
export const VALIDATION_LIMITS = {
  PAYLOAD: { min: 0.1, max: 1000000 },
  DELTA_V: { min: 100, max: 20000 },
  STAGES: { min: 1, max: 10 },
  ISP: { min: 50, max: 500 },
  MASS_FRACTION: { min: 0.01, max: 0.99 }
};