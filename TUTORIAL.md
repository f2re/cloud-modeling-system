# Tutorial: Working with Cloud Modeling System (CMS)

This guide walks you through setting up, running simulations, and validating results using the CMS framework.

## 1. Environment Setup

Ensure you have Python 3.10+ and the required dependencies installed:

```bash
# Clone the repository
git clone https://github.com/f2re/cloud-modeling-system.git
cd cloud-modeling-system

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

## 2. Running a Basic Simulation

The main entry point is `main.py` (ensure it exists or create a simple one). A basic simulation script looks like this:

```python
import numpy as np
from cms.model import CMSModel
from cms.config import GridConfig, PhysicsConfig

# 1. Setup Configuration
g_config = GridConfig(nx=64, ny=64, nz=40, dx=1000, dy=1000, dz=250)
p_config = PhysicsConfig()

# 2. Initialize Model
model = CMSModel(g_config, p_config)

# 3. Set Initial Conditions (e.g., a warm bubble)
# model.theta[32, 32, 5] += 2.0 

# 4. Simulation Loop
dt = 2.0  # seconds
for step in range(100):
    model.step(dt)
    if step % 10 == 0:
        print(f"Step {step} completed.")

# 5. Access Results
print(f"Max vertical velocity: {np.max(model.w)} m/s")
```

## 3. Using Diagnostics (Phase 6 Features)

You can now calculate radar reflectivity and seeding metrics:

```python
from cms.utils.diagnostics import compute_radar_reflectivity, compute_seeding_metrics

# Calculate Radar Reflectivity (dBZ)
dbz = compute_radar_reflectivity(
    model.rho, 
    model.qr, model.nr, 
    model.qi, model.ni,
    config=model.physics
)

print(f"Max Reflectivity: {np.max(dbz):.2f} dBZ")

# Calculate Seeding Metrics (if comparing two runs)
# metrics = compute_seeding_metrics(precip_seeded, precip_control, ...)
```

## 4. Multi-Sensor Validation

To validate your model against observations (e.g., from a field campaign like SNOWIE):

```python
from cms.utils.validation import run_validation_suite

# Prepare model output dictionary
model_data = {
    'rho': model.rho, 'qc': model.qc, 'qr': model.qr, 
    'qi': model.qi, 'qs': model.qs, 'qg': model.qg,
    'nc': model.nc, 'nr': model.nr, 'ni': model.ni,
    'ns': model.ns, 'ng': model.ng
}

# Provide observations (as numpy arrays matching the grid or interpolated)
obs_data = {
    'Z': obs_radar_dbz,
    'LWP': obs_radiometer_lwp
}

results = run_validation_suite(model_data, obs_data, dz=g_config.dz)

for sensor, stats in results.items():
    print(f"--- {sensor} ---")
    print(f"RMSE: {stats['rmse']:.2f}")
    print(f"Bias: {stats['bias']:.2f}")
    print(f"Correlation: {stats['correlation']:.2f}")
```

## 5. Optimization & Performance (Phase 7)

CMS uses **Numba** to accelerate numerical kernels. 

### How to ensure Numba is used:
1. Ensure `numba` is installed: `pip install numba`
2. The `WENO5` advection scheme will automatically detect Numba and use the optimized kernels in `cms/dynamics/advection_numba.py`.
3. You can verify this by checking `model.dynamics.weno.use_numba`.

### Running in Docker:
```bash
docker build -t cms-model .
docker run cms-model
```

## 6. Development Workflow with Gemini CLI

This project is optimized for AI-assisted development.

- Use `@cms-orchestrator` to plan your next steps.
- Use `@cms-coder` to implement new physics modules.
- Use `@cms-physicist` to verify your implementation against `IMPLEMENTATION_GUIDE.md`.

Example command:
`gemini "Add a new ice nucleation parameterization based on Guide section 3.5.2"`
