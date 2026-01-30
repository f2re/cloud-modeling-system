#!/usr/bin/env python3
"""
Diagnostic script to check available compute backends for CMS.
Usage: python check_compute_backend.py
"""

import sys
import platform

def check_system_info():
    print("=" * 60)
    print("CMS Compute Backend Diagnostic")
    print("=" * 60)
    print(f"\nPython: {sys.version}")
    print(f"Platform: {platform.platform()}")
    print(f"Architecture: {platform.machine()}")
    print()

def check_numpy():
    try:
        import numpy as np
        print(f"✓ NumPy: {np.__version__}")
        return True
    except ImportError:
        print("✗ NumPy: NOT INSTALLED")
        return False

def check_numba():
    try:
        import numba
        print(f"✓ Numba: {numba.__version__}")
        print(f"  - JIT compilation: Available")
        print(f"  - Threading layer: {numba.config.THREADING_LAYER}")
        
        # Test compilation
        @numba.njit
        def test_func(x):
            return x * 2
        
        result = test_func(5)
        print(f"  - Test compilation: OK (result={result})")
        return True
    except ImportError:
        print("✗ Numba: NOT INSTALLED")
        print("  Install with: pip install numba>=0.58.0")
        return False
    except Exception as e:
        print(f"⚠ Numba: Installed but error during test: {e}")
        return False

def check_cuda():
    try:
        from numba import cuda
        print(f"\n✓ Numba CUDA: Available")
        
        if cuda.is_available():
            print(f"  - CUDA devices detected: {len(cuda.gpus)}")
            for i, gpu in enumerate(cuda.gpus):
                print(f"    [{i}] {gpu.name.decode()}")
                print(f"        Compute Capability: {gpu.compute_capability}")
                ctx = gpu.create_context()
                mem_info = cuda.current_context().get_memory_info()
                print(f"        Free/Total Memory: {mem_info[0]/(1024**3):.2f} GB / {mem_info[1]/(1024**3):.2f} GB")
                ctx.pop()
            return True
        else:
            print("✗ CUDA devices: None detected")
            print("  Check: nvidia-smi")
            return False
    except ImportError:
        print("\n✗ Numba CUDA: Module not available")
        return False
    except Exception as e:
        print(f"\n⚠ Numba CUDA: Error - {e}")
        return False

def check_cupy():
    try:
        import cupy as cp
        print(f"\n✓ CuPy: {cp.__version__}")
        
        if cp.cuda.is_available():
            device = cp.cuda.Device()
            print(f"  - Device: {device.compute_capability}")
            mem_info = cp.cuda.runtime.memGetInfo()
            print(f"  - Free/Total Memory: {mem_info[0]/(1024**3):.2f} GB / {mem_info[1]/(1024**3):.2f} GB")
            
            # Test computation
            a = cp.arange(10)
            b = a * 2
            print(f"  - Test computation: OK")
            return True
        else:
            print("✗ CuPy: No CUDA devices")
            return False
    except ImportError:
        print("\n✗ CuPy: NOT INSTALLED")
        print("  Install with: pip install cupy-cuda12x  # or cupy-cuda11x")
        return False
    except Exception as e:
        print(f"\n⚠ CuPy: Error - {e}")
        return False

def check_cms_modules():
    print("\n" + "=" * 60)
    print("CMS Module Check")
    print("=" * 60)
    
    try:
        from cms.dynamics.advection_numba import weno5_reconstruct_x
        print("✓ cms.dynamics.advection_numba: Available")
    except ImportError as e:
        print(f"✗ cms.dynamics.advection_numba: Import failed - {e}")
    
    try:
        from cms.dynamics.advection_cuda import weno5_reconstruct_x_kernel
        print("✓ cms.dynamics.advection_cuda: Available")
    except ImportError as e:
        print(f"✗ cms.dynamics.advection_cuda: Import failed - {e}")

def recommend_config():
    print("\n" + "=" * 60)
    print("Recommended Configuration")
    print("=" * 60)
    
    has_numba = check_numba()
    has_cuda = False
    
    if has_numba:
        has_cuda = check_cuda()
        has_cupy = check_cupy()
        has_cuda = has_cuda and has_cupy
    
    check_cms_modules()
    
    print("\n" + "=" * 60)
    print("RECOMMENDATION")
    print("=" * 60)
    
    if has_cuda:
        print("\n✅ Use GPU (CUDA) mode")
        print("\nIn main.py, set:")
        print("  c_config = ComputeConfig(use_gpu=True)")
        print("\nExpected speedup: 50-100× vs pure NumPy")
    elif has_numba:
        print("\n✅ Use CPU (Numba) mode")
        print("\nIn main.py, set:")
        print("  c_config = ComputeConfig(use_gpu=False)")
        print("\nExpected speedup: 10-20× vs pure NumPy")
        print("\nTo enable GPU acceleration:")
        print("  1. Install CUDA Toolkit: See GPU_SETUP.md")
        print("  2. Install CuPy: pip install cupy-cuda12x")
    else:
        print("\n⚠ Falling back to pure NumPy (SLOW)")
        print("\nTo improve performance:")
        print("  1. Install Numba: pip install numba>=0.58.0")
        print("  2. Rerun this script to verify")
    
    print()

if __name__ == "__main__":
    check_system_info()
    has_numpy = check_numpy()
    
    if not has_numpy:
        print("\n❌ CRITICAL: NumPy is required but not installed!")
        print("Install with: pip install numpy")
        sys.exit(1)
    
    recommend_config()
