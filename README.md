# Cloud Modeling System (CMS)

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Status](https://img.shields.io/badge/Status-Phase%201%3A%20Foundation-yellow)](DEVELOPMENT.md)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

## Overview
The **Cloud Modeling System (CMS)** is a high-performance, modular Python framework designed for simulating atmospheric thermodynamics, cloud microphysics, and weather modification (seeding) scenarios.

It is built on a **modular Eulerian-Lagrangian** hybrid architecture, prioritizing:
*   **Separation of Concerns:** Dynamics, Microphysics, and Diffusion are decoupled.
*   **Vectorization:** Core numerical operations utilize `numpy` for performance.
*   **Configurability:** Physics constants and grid parameters are strictly isolated.

## Project Structure
```text
cms/
├── core/           # Numerical backbone (Grid, Stencils, RK-Integrators)
├── dynamics/       # Navier-Stokes solvers (Advection, Pressure, Buoyancy)
├── microphysics/   # Phase changes (Kessler, Morrison, Seeding interaction)
├── diffusion/      # Turbulence modeling (LES, Smagorinsky)
└── utils/          # Data I/O (NetCDF), Logging, Diagnostics
```

## Documentation
*   **[IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md)**: The **Source of Truth** for all physics equations, constants, and algorithms.
*   **[DEVELOPMENT.md](DEVELOPMENT.md)**: Current roadmap, active phase, and implementation status.

## Getting Started

### Prerequisites
*   Python 3.10 or higher
*   Virtual environment (recommended)

### Installation
1.  Clone the repository.
2.  Create and activate a virtual environment:
    ```bash
    python -m venv venv
    source venv/bin/activate
    ```
3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Running Tests
To verify the installation and currently implemented modules:
```bash
python -m unittest discover tests -v
```

## Contributing
Please refer to `DEVELOPMENT.md` for the current tasks and coding standards.
*   **Branching:** Use `feature/<topic>` branches.
*   **Commits:** Use descriptive, atomic commits.
*   **Validation:** All changes must pass unit tests and adhere to the physics guide.

## License
MIT License
