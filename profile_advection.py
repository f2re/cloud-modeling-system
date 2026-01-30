import time
import cProfile
import pstats
import numpy as np
from cms.core.grid import Grid
from cms.config import GridConfig
from cms.dynamics.advection import WENO5

def benchmark():
    # Setup a reasonably sized grid
    config = GridConfig(nx=64, ny=64, nz=64)
    grid = Grid(config)
    weno = WENO5(grid)

    # Initialize fields
    q = np.random.rand(*grid.shape)
    u = np.ones(grid.shape) * 10.0
    v = np.ones(grid.shape) * 5.0
    w = np.ones(grid.shape) * 2.0

    print(f"Benchmarking WENO5 Advection on {grid.shape} grid...")
    
    # Warmup for JIT
    print("Warming up JIT...")
    _ = weno.advect(q, u, v, w)
    
    start_time = time.time()
    # Run 10 iterations
    for _ in range(10):
        _ = weno.advect(q, u, v, w)
    end_time = time.time()
    
    print(f"Total time for 10 steps: {end_time - start_time:.4f}s")
    print(f"Time per step: {(end_time - start_time)/10:.4f}s")

if __name__ == "__main__":
    # Use cProfile for detailed stats
    profiler = cProfile.Profile()
    profiler.enable()
    benchmark()
    profiler.disable()
    
    stats = pstats.Stats(profiler).sort_stats('cumtime')
    stats.print_stats(20)
