# Development Guide: Cloud Modeling System (CMS)

## Project Overview
This project implements a numerical model for cloud dynamics, microphysics, and seeding effects, based on the specifications in `IMPLEMENTATION_GUIDE.md`. The system is designed to be modular, extensible, and performant.

## Architecture

The project follows a modular architecture to separate concerns between dynamics, microphysics, and dispersion.

```
cloud-modeling-system/
├── cms/                      # Main package source
│   ├── __init__.py
│   ├── config.py             # Configuration parameters (Grid, Physics constants)
│   ├── core/                 # Core numerical infrastructure
│   │   ├── __init__.py
│   │   ├── grid.py           # Grid generation and management
│   │   └── time_integration.py # Runge-Kutta schemes
│   ├── dynamics/             # Fluid dynamics
│   │   ├── __init__.py
│   │   ├── advection.py      # WENO schemes
│   │   └── navier_stokes.py  # Momentum equations
│   ├── microphysics/         # Cloud microphysics
│   │   ├── __init__.py
│   │   ├── warm.py           # Warm rain processes
│   │   ├── ice.py            # Ice phase processes
│   │   └── seeding.py        # AgI seeding interactions
│   ├── diffusion/            # Turbulence and Dispersion
│   │   ├── __init__.py
│   │   └── turbulence.py     # Smagorinsky & Eddy diffusivity
│   └── utils/                # Utilities
│       ├── __init__.py
│       ├── io.py             # NetCDF I/O
│       └── visualization.py  # Basic plotting tools
├── tests/                    # Unit and integration tests
├── main.py                   # Entry point for simulations
├── requirements.txt          # Python dependencies
├── DEVELOPMENT.md            # This file
└── IMPLEMENTATION_GUIDE.md   # Physics and Math reference
```

## Development Phases

### Phase 1: Foundation (Current)
- **Goal**: Establish the grid, basic time stepping, and simple advection.
- **Tasks**:
    1.  Implement `cms/config.py` with physical constants.
    2.  Create the 3D grid class in `cms/core/grid.py`.
    3.  Implement a basic Runge-Kutta integrator in `cms/core/time_integration.py`.
    4.  Implement a simple advection test in `tests/test_advection.py`.

### Phase 2: Dynamics
- **Goal**: Implement the compressible Navier-Stokes solver.
- **Tasks**:
    1.  Add momentum and continuity equations in `cms/dynamics/`.
    2.  Implement boundary conditions (sponge layer, no-slip).

### Phase 3: Basic Microphysics
- **Goal**: Add warm rain physics (Kessler scheme).
- **Tasks**:
    1.  Implement `cms/microphysics/warm.py`.
    2.  Coupling microphysics with dynamics.

### Phase 4: Advanced Physics & Seeding
- **Goal**: Full double-moment microphysics and AgI seeding.
- **Tasks**:
    1.  Implement ice phase in `cms/microphysics/ice.py`.
    2.  Add seeding dispersion in `cms/diffusion/`.

## Setup & Workflow

### Prerequisites
- Python 3.10+
- `numpy`, `scipy`, `matplotlib` (for visualization), `netCDF4` (for I/O)

### Installation
```bash
pip install -r requirements.txt
```

### Running Tests
```bash
python -m unittest discover tests
```

## Coding Standards
- **Style**: PEP 8.
- **Documentation**: Docstrings for all public functions/classes.
- **Type Hinting**: Use Python type hints strictly.
- **Testing**: Every module must have a corresponding unit test.

## Git Workflow
- Main branch: `main`
- Feature branches: `feature/<name>`
- Commit messages: Descriptive and imperative (e.g., "Add WENO5 advection scheme").
