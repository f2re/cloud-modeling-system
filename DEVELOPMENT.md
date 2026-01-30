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

### Phase 1: Foundation (Completed)
- **Goal**: Establish the grid, basic time stepping, and simple advection.
- **Tasks**:
    1.  [x] Implement `cms/config.py` with physical constants.
    2.  [x] Create the 3D grid class in `cms/core/grid.py`.
    3.  [x] Implement a basic Runge-Kutta integrator in `cms/core/time_integration.py`.
    4.  [x] Implement a simple advection test in `tests/test_advection.py`.

### Phase 2: Dynamics (Completed)
- **Goal**: Implement the compressible Navier-Stokes solver and WENO-5 advection.
- **Tasks**:
    1.  [x] Implement 5th-order WENO advection scheme in `cms/dynamics/advection.py`.
    2.  [x] Add momentum and continuity equations in `cms/dynamics/navier_stokes.py`.
    3.  [x] Implement boundary conditions (sponge layer, no-slip).

### Phase 3: Basic Microphysics (Completed)
- **Goal**: Add warm rain physics (Kessler scheme).
- **Tasks**:
    1.  [x] Implement `cms/microphysics/warm.py`.
    2.  [x] Coupling microphysics with dynamics.

### Phase 4: Advanced Physics & Seeding (Completed)
- **Goal**: Full double-moment microphysics and AgI seeding.
- **Tasks**:
    1.  [x] Implement ice phase in `cms/microphysics/ice.py`.
    2.  [x] Add seeding dispersion in `cms/diffusion/`.

### Phase 5: Advanced Features (Completed)
- **Goal**: Aerosol-cloud interactions and LES turbulence model refinement.
- **Tasks**:
    1.  [x] Aerosol-cloud interactions (CCN activation).
    2.  [x] Refine LES turbulence model (Smagorinsky implemented in Phase 4).

### Phase 6: Validation (Completed)
- **Goal**: Radar forward operator and multi-sensor comparison.
- **Tasks**:
    1.  [x] Implement Radar forward operator in `cms/utils/diagnostics.py`.
    2.  [x] Create multi-sensor validation script in `cms/utils/validation.py`.

### Phase 7: Optimization & Deployment (In Progress)
- **Goal**: Performance and production readiness.
- **Tasks**:
    1.  [x] Profiling and Numba optimization for core stencils (Advection).
    2.  [x] Containerization (Dockerfile).
    3.  [ ] Comprehensive documentation and tutorials.

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
