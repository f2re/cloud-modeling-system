import netCDF4 as nc
import numpy as np
import os
from datetime import datetime
from cms.core.grid import Grid

class IOManager:
    """
    Handles NetCDF input/output for the CMS model.
    """
    def __init__(self, output_dir: str = "output"):
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def save_state(self, filename: str, grid: Grid, time_s: float, 
                   u: np.ndarray, v: np.ndarray, w: np.ndarray, 
                   rho: np.ndarray, theta: np.ndarray, 
                   qc: np.ndarray, qr: np.ndarray,
                   qi: np.ndarray, qs: np.ndarray, qg: np.ndarray,
                   ni: np.ndarray, ns: np.ndarray, ng: np.ndarray,
                   c_agi: np.ndarray):
        """
        Writes the current model state to a NetCDF file.
        """
        filepath = os.path.join(self.output_dir, filename)
        
        # Open file in write mode
        with nc.Dataset(filepath, 'w', format='NETCDF4') as ds:
            ds.description = "Cloud Modeling System Output"
            ds.history = f"Created {datetime.now().isoformat()}"
            ds.time = time_s
            
            # Define dimensions (using inner domain size)
            ds.createDimension('x', grid.nx)
            ds.createDimension('y', grid.ny)
            ds.createDimension('z', grid.nz)
            
            # Create variables
            x_var = ds.createVariable('x', 'f4', ('x',))
            y_var = ds.createVariable('y', 'f4', ('y',))
            z_var = ds.createVariable('z', 'f4', ('z',))
            
            x_var[:] = grid.x
            y_var[:] = grid.y
            z_var[:] = grid.z
            
            def save_field(name, data, units):
                var = ds.createVariable(name, 'f4', ('x', 'y', 'z'))
                var[:] = grid.get_inner(data)
                var.units = units
            
            save_field('u', u, 'm/s')
            save_field('v', v, 'm/s')
            save_field('w', w, 'm/s')
            save_field('rho', rho, 'kg/m3')
            save_field('theta', theta, 'K')
            
            # Warm
            save_field('qc', qc, 'kg/kg')
            save_field('qr', qr, 'kg/kg')
            
            # Ice
            save_field('qi', qi, 'kg/kg')
            save_field('qs', qs, 'kg/kg')
            save_field('qg', qg, 'kg/kg')
            save_field('ni', ni, 'm-3')
            save_field('ns', ns, 'm-3')
            save_field('ng', ng, 'm-3')
            
            # Seeding
            save_field('c_agi', c_agi, 'kg/kg')
            
            # Global attributes
            ds.dx = grid.dx
            ds.dy = grid.dy
            ds.dz = grid.dz

        return filepath
