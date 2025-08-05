<div align="center">

<img src="assets/icon.png" alt="Rocketmancer Icon" width="128" height="128">

# rocketmancer-web

![Version](https://img.shields.io/badge/version-0.1.1-blue?style=flat)
![License](https://img.shields.io/badge/license-MIT-green?style=flat)
![Status](https://img.shields.io/badge/status-beta-orange?style=flat)
![CI](https://img.shields.io/github/actions/workflow/status/BruhLemma-Yadecha/rocketmancer-web/pr-checks.yml?branch=main&style=flat&label=CI)

![Python](https://img.shields.io/badge/Python-3.13+-306998?style=flat&logo=python&logoColor=white)
![React](https://img.shields.io/badge/React-18+-61DAFB?style=flat&logo=react&logoColor=white)
![Django](https://img.shields.io/badge/Django-5.1+-092E20?style=flat&logo=django&logoColor=white)
![JAX](https://img.shields.io/badge/JAX-0.7+-FF6F00?style=flat&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-supported-2496ED?style=flat&logo=docker&logoColor=white)
![uv](https://img.shields.io/badge/uv-package%20manager-306998?style=flat&logo=python&logoColor=white)

*A multi-stage delta-v optimizer for rockets, built with modern web technologies.*

</div>

## Background
rocketmancer-web is a multi-stage rocket delta-V optimizer based on the [original rocketmancer](https://github.com/BruhLemma-Yadecha/rocketmancer) and its own predecessor, [multistage](https://github.com/BruhLemma-Yadecha/multistage). It takes a set of trip and vehicle parameters to make the initial mass of the rocket as small as possible.

With rocketmancer-web, I wanted to fix the primary issue of both of those legacy projects: the lack of an easy-to-use, friendly user experience. By using JAX and modern optimization techniques, I've also further optimized the algorithm that generates the final configuration, making this version best-in-class.

## Concept

The backbone of rocketmancer is the Tsiolkovsky Rocket Equation (see the [Wikipedia article](https://en.wikipedia.org/wiki/Tsiolkovsky_rocket_equation)). The parameters that are often at hand when planning a mission and designing a vehicle are total delta-v, stage specific impulse values, and achievable ranges for the propellant mass fraction. The formula is trivial for a single stage vehicle and vehicles with stages that have the same mass fraction and specific impulse, but requires finding the optimal split when each stage performs differently. This is the primary optimization task rocketmancer does.

## Quick Start

### Using Docker (Recommended)

1. **Clone and setup**

   ```bash
   git clone <repository-url>
   cd rocketmancer-web
   cp .env.example .env
   ```

2. **Start the application**

   ```bash
   make dev
   # or
   docker-compose up --build
   ```

3. **Access the application**
   - Frontend: <http://localhost>
   - Backend API: <http://localhost:8000/api/v1/>
   - Admin Panel: <http://localhost:8000/admin>

## Development

Use `make help` to see all available development commands. The project uses uv for Python package management and Docker for containerization.

## Project Structure

```text
rocketmancer-web/
├── backend/                # Django REST API
├── frontend/               # React SPA  
├── packages/rocketmancer/  # Core optimization library
├── nginx/                  # Reverse proxy config
├── docs/                   # Documentation
└── Makefile               # Development commands
```

## License

This project is licensed under the [MIT License](https://opensource.org/license/mit) - see the [LICENSE](LICENSE) file for details.
