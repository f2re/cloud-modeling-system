import numpy as np
import time
from cms.model import CMSModel
from cms.config import GridConfig, PhysicsConfig, ComputeConfig
from cms.utils.diagnostics import compute_radar_reflectivity
from cms.utils.io import IOManager

def run_simulation():
    print("=== CMS: Starting Warm Bubble Simulation ===")
    
    # 1. Configuration
    # Smaller grid for a quick demonstration
    g_config = GridConfig(nx=40, ny=40, nz=40, dx=100, dy=100, dz=100)
    p_config = PhysicsConfig()
    c_config = ComputeConfig(use_gpu=True)
    
    # 2. Initialize Model
    model = CMSModel(g_config, p_config, c_config)
    io_manager = IOManager("output")
    
    # 3. Initial Condition: Warm Bubble
    # Add a +2K perturbation in the center to trigger convection
    mid_x, mid_y = g_config.nx // 2, g_config.ny // 2
    model.theta[mid_x-2:mid_x+3, mid_y-2:mid_y+3, 5:10] += 2.0
    # Initialize background droplets to avoid singularity
    model.nc[:] = 1e8
    
    print(f"Grid size: {g_config.nx}x{g_config.ny}x{g_config.nz}")
    print(f"Compute Mode: {'GPU (CUDA)' if model.dynamics.weno.use_gpu else ('CPU (Numba)' if model.dynamics.weno.use_numba else 'CPU (Pure Python)')}")
    
    # 4. Integration Loop
    # Reduced dt to ensure stability (CFL condition)
    # dx=100m, max_vel ~10m/s -> dt < 100/10 = 10s. 
    # But acoustic waves ~340m/s -> dt < 100/340 ~ 0.3s
    dt = 0.2  
    total_steps = 250
    
    print(f"\nRunning {total_steps} steps...")
    start_wall = time.time()
    
    for step in range(1, total_steps + 1):
        model.step(dt)
        
        if step % 50 == 0:  # Print every 10 simulation seconds
            max_w = np.max(model.w)
            max_qc = np.max(model.qc) * 1000  # g/kg
            
            # Calculate diagnostics on the fly
            dbz = compute_radar_reflectivity(
                model.rho, model.qr, model.nr, 
                model.qi, model.ni, config=model.physics
            )
            max_z = np.max(dbz)
            
            print(f"Step {step:03d} | Max W: {max_w:6.3f} m/s | Max Qc: {max_qc:6.3f} g/kg | Max Z: {max_z:5.1f} dBZ")
            
            # Save state
            io_manager.save_state(
                f"cms_output_{step:04d}.nc",
                model.grid, model.time,
                model.u, model.v, model.w,
                model.rho, model.theta,
                model.qc, model.qr, model.nc, # nc added
                model.qi, model.qs, model.qg,
                model.ni, model.ns, model.ng,
                model.c_agi, model.n_ccn, model.n_inp_nat
            )

    end_wall = time.time()
    
    print("\n=== Simulation Complete ===")
    print(f"Total Wall Time: {end_wall - start_wall:.2f}s")
    print(f"Final Max W: {np.max(model.w):.3f} m/s")

    # 5. Visualization
    from cms.utils.visualization import Visualizer
    print("\nGenerating animations...")
    viz = Visualizer("output")
    if viz.files:
        viz.create_animation("w", y_index=mid_y, output_filename="output/w_animation.gif")
        viz.create_animation("qc", y_index=mid_y, output_filename="output/qc_animation.gif")
        viz.create_animation("theta", y_index=mid_y, output_filename="output/theta_animation.gif")

if __name__ == "__main__":
    run_simulation()
