# Cloud Modeling System (CMS) - Context & Instructions

## Project Mission
The Cloud Modeling System (CMS) is a high-performance Python framework for simulating atmospheric thermodynamics, cloud microphysics, and weather modification (seeding) scenarios. The goal is to provide a scientifically rigorous, modular, and extensible tool for researchers.

## Architecture & Design Philosophy

The system is built on a **modular Eulerian-Lagrangian** hybrid architecture (where applicable) or purely Eulerian for the baseline.
- **Separation of Concerns:** Dynamics (fluid flow) are decoupled from Microphysics (phase changes), allowing different schemes to be swapped easily.
- **NumPy-First:** All core numerical operations must be vectorized. Avoid explicit Python loops over grid points.
- **Configurability:** Physics constants and grid parameters are isolated in configuration files, ensuring reproducibility.

```text
cms/
├── core/           # Numerical backbone (Grid, Stencils, RK-Integrators)
├── dynamics/       # Navier-Stokes solvers (Advection, Pressure, Buoyancy)
├── microphysics/   # Phase changes (Kessler, Morrison, Seeding interaction)
├── diffusion/      # Turbulence modeling (LES, Smagorinsky)
└── utils/          # Data I/O (NetCDF), Logging, Diagnostics
```

## Development Standards & Best Practices

### 1. Scientific Coding Standards
*   **Variable Naming:** Use descriptive names that reflect physical meaning and units.
    *   *Bad:* `t`, `p`, `val`
    *   *Good:* `time_s`, `pressure_pa`, `mixing_ratio_kg_kg`
*   **Math-to-Code Mapping:** In complex numerical functions, explicitly cite the equation number from `IMPLEMENTATION_GUIDE.md` or the source paper in the docstring/comments.
    ```python
    def compute_terminal_velocity(diameter_m):
        """
        Calculates terminal velocity using the power law.
        Ref: IMPLEMENTATION_GUIDE.md Eq 3.4
        """
        ...
    ```
*   **Type Hinting:** Strictly use `typing` (or `numpy.typing`) for all function signatures to clarify expected tensor shapes and data types.

### 2. Performance Guidelines
*   **Vectorization:** Use `numpy` array operations for all field updates.
    *   *Avoid:* `for i in range(nx): for j in range(ny): ...`
    *   *Use:* `field[1:-1, 1:-1] += dt * tendency[...]`
*   **Memory Efficiency:** Avoid unnecessary copying of large 3D arrays. Use in-place operations (`+=`, `*=`) where safe.
*   **Profiling:** Use `cProfile` or `line_profiler` before attempting optimization. Don't guess where the bottleneck is.

### 3. Verification & Validation (V&V)
*   **Verification (Code Correctness):**
    *   **Unit Tests:** Test individual math functions against analytical solutions (e.g., advecting a square wave).
    *   **Conservation Checks:** The model **must** conserve mass and energy (within machine precision/truncation error). Implement automated checks that sum total water/energy in the domain at each step.
*   **Validation (Physical Realism):**
    *   Compare model output against idealized test cases (e.g., "warm bubble" rise) before running complex real-world scenarios.
    *   Sanity checks: Ensure no negative masses or absolute temperatures below 0K.

### 4. Implementation Workflow
1.  **Plan:** Break down a physical process into inputs, outputs, and governing equations.
2.  **Prototype:** Implement the math in a standalone script or notebook to verify the algorithm.
3.  **Integrate:** Move the verified code into the `cms` package structure.
4.  **Test:** Add a unit test in `tests/` that covers the new functionality.
5.  **Document:** Update docstrings and `DEVELOPMENT.md` if architectural changes occurred.

## Operational Commands

### Environment
Ensure your virtual environment is active.
```bash
source venv/bin/activate  # or equivalent
```

### Installation
```bash
pip install -r requirements.txt
```

### Testing
Run the full suite with verbosity to see which physics modules are passing.
```bash
python -m unittest discover tests -v
```

### Static Analysis
Use standard tools to maintain code quality.
```bash
# Check for PEP 8 compliance and potential errors
flake8 cms/ tests/
# Check static types
mypy cms/
```
## Comments
*   **All comments must be on Russian** and describe variables and functions logics

## Key Documentation Map
*   `IMPLEMENTATION_GUIDE.md`: **The Source of Truth** for physics equations and constants.
*   `DEVELOPMENT.md`: Current implementation status and roadmap.
*   `tests/`: Examples of how to use the modules in isolation.

## Git Strategy
*   **Main Branch (`main`):** Always stable, passing all tests.
*   **Feature Branches (`feature/topic`):** Used for developing new physics modules (e.g., `feature/ice-microphysics`).
*   **Commits:** Atomic commits grouping logical changes. Start messages with the module changed (e.g., `dynamics: fix boundary condition bug`).