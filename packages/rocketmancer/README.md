# Rocketmancer

Multi-stage rocket optimization library for delta-v allocation.

## Overview

Rocketmancer is a Python library for optimizing multi-stage rocket configurations. It uses a greedy allocation algorithm to distribute delta-v across rocket stages to minimize total vehicle mass while meeting mission requirements.

## Features

- **Multi-stage optimization**: Optimize delta-v allocation across any number of rocket stages
- **Greedy solver**: Efficient closed-form solver with configurable parameters
- **Robust validation**: Comprehensive input validation and error handling
- **Type safety**: Full type hints and runtime validation with Pydantic
- **High performance**: Built on JAX and NumPy for efficient numerical computation

## Installation

```bash
pip install rocketmancer
```

## Quick Start

```python
from rocketmancer import Rocket, Stage

# Define rocket stages (from top to bottom)
stages = [
    Stage(isp=350.0, propellant_mass_fraction=0.85),  # Upper stage
    Stage(isp=280.0, propellant_mass_fraction=0.90),  # Lower stage
]

# Create rocket with mission requirements
rocket = Rocket(
    payload=1000.0,        # kg
    total_delta_v=9500.0,  # m/s
    stages=stages
)

# Optimize delta-v allocation
result = rocket.optimize()
print(f"Optimized total mass: {rocket.total_mass:.1f} kg")

# Access detailed stage information
for i, stage in enumerate(rocket.stages):
    print(f"Stage {i+1}: {stage.delta_v:.1f} m/s, {stage.total_mass:.1f} kg")
```

## API Reference

### Core Classes

#### `Stage`

Represents a single rocket stage with propulsion characteristics.

**Parameters:**

- `isp` (float): Specific impulse in seconds
- `propellant_mass_fraction` (float): Ratio of propellant mass to total stage mass (0 < φ < 1)

#### `Rocket`

Represents a complete multi-stage rocket system.

**Parameters:**

- `payload` (float): Payload mass in kg
- `total_delta_v` (float): Total required delta-v in m/s
- `stages` (List[Stage]): List of stages from top to bottom

### Solver Configuration

The greedy solver supports configurable parameters for fine-tuning optimization behavior:

```python
# Use custom solver parameters
result = rocket.optimize(
    method='greedy',
    safety_margin=0.995,     # Numerical stability factor (default: 0.999)
    convergence_tol=1e-12    # Allocation precision (default: 1e-15)
)
```

### Error Handling

```python
from rocketmancer import SolverError

try:
    result = rocket.optimize()
except SolverError as e:
    print(f"Optimization failed: {e}")
except ValueError as e:
    print(f"Invalid input: {e}")
```

## Requirements

- Python ≥ 3.11
- JAX ≥ 0.7.0
- NumPy ≥ 2.3.1
- Pydantic ≥ 2.0.0
- SciPy ≥ 1.11

## License

MIT License - see LICENSE file for details.
