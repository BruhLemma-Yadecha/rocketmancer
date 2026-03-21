<div align="center">

<img src="assets/icon.png" alt="rocketmancer" width="128" height="128">

# rocketmancer

![License](https://img.shields.io/badge/license-MIT-green?style=flat)
![CI](https://img.shields.io/github/actions/workflow/status/BruhLemma-Yadecha/rocketmancer/ci.yml?branch=main&style=flat&label=CI)
![React](https://img.shields.io/badge/React-18+-61DAFB?style=flat&logo=react&logoColor=white)
![Vite](https://img.shields.io/badge/Vite-7-646CFF?style=flat&logo=vite&logoColor=white)

*A multi-stage rocket delta-v optimizer.*

[**Live Demo**](https://rocketmancer.thaumics.org)

<img src="assets/screenshot.png" alt="rocketmancer screenshot" width="800">

</div>

## Background

rocketmancer is a multi-stage rocket delta-v optimizer based on the [original rocketmancer](https://github.com/BruhLemma-Yadecha/rocketmancer-cpp) and its predecessor, [multistage](https://github.com/BruhLemma-Yadecha/multistage). It takes trip and vehicle parameters and finds the delta-v split across stages that minimizes total vehicle mass.

The primary goal was to fix the usability issues of those legacy projects by providing a proper web interface. The solver itself was also reworked: rather than using a general-purpose constrained optimizer, the problem is solved analytically via Lagrange multipliers, reducing the full optimization to a simple bisection on a single scalar.

## How It Works

The backbone is the [Tsiolkovsky Rocket Equation](https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation). For a multi-stage rocket, total mass is a product of per-stage mass amplification factors, each depending only on that stage's delta-v allocation. Taking the log decomposes the objective into a sum of independent terms, and the Lagrange optimality condition gives a closed-form expression for each stage's mass ratio given a multiplier λ. Bisecting on λ until the delta-v allocations sum to the target solves the problem exactly.

When a minimum stage contribution is set, the solver uses an active-set method: stages that fall below the floor are pinned to the minimum, and the remaining delta-v is re-optimized across unpinned stages. Dual variable checks ensure pinned stages are released if they become cost-effective, guaranteeing the constrained optimum.

## Quick Start

```bash
git clone https://github.com/BruhLemma-Yadecha/rocketmancer.git
cd rocketmancer
npm install
npm run dev
```

## References

The solver's formulation draws from the following:

- J.W. Cornelisse, H.F.R. Schoyer & K.F. Wakker, *Rocket Propulsion and Spaceflight Dynamics* (1979) — mass ratio decomposition, structural coefficient definitions, and Lagrange multiplier formulation for optimal staging
- [NASA Glenn Research Center — Mass Ratios](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-ratios/) — reference for stage mass ratio, structural coefficient, and payload ratio relationships
- J. Nocedal & S.J. Wright, *Numerical Optimization*, 2nd ed. (2006) — active-set methods for inequality-constrained convex problems (ch. 16), KKT conditions and dual variable interpretation

## License

[MIT](LICENSE)
