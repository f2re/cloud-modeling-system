import os
import glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from typing import List, Optional

class Visualizer:
    """
    Visualization tool for CMS model output.
    """
    def __init__(self, input_dir: str = "output"):
        self.input_dir = input_dir
        self.files = sorted(glob.glob(os.path.join(input_dir, "cms_output_*.nc")))
        if not self.files:
            print(f"No CMS output files found in {input_dir}")
        else:
            print(f"Found {len(self.files)} output files.")

    def load_data(self, file_index: int):
        """Loads data from a specific file index."""
        if file_index >= len(self.files):
            raise IndexError("File index out of range")
        
        filepath = self.files[file_index]
        ds = nc.Dataset(filepath)
        return ds

    def plot_slice_xz(self, field_name: str, y_index: int, file_index: int, 
                      ax: Optional[plt.Axes] = None, vmin=None, vmax=None, cmap='viridis'):
        """
        Plots a vertical slice (X-Z) of the given field at a specific Y index.
        """
        ds = self.load_data(file_index)
        
        # Get coordinates
        x = ds.variables['x'][:]
        z = ds.variables['z'][:]
        
        # Get field (dimensions: x, y, z) -> slice -> (x, z)
        # Note: NetCDF saved as (x, y, z) based on io.py, but usually we want (z, x) for plotting with pcolormesh
        # Let's check shape.
        field = ds.variables[field_name][:, y_index, :] # Shape (nx, nz)
        
        # Transpose for plotting: (nz, nx) so Z is vertical axis
        field_t = field.T
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))
            
        X, Z = np.meshgrid(x, z)
        
        mesh = ax.pcolormesh(X, Z, field_t, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Z (m)')
        ax.set_title(f"{field_name} at Y={y_index*ds.dy:.0f}m (t={ds.time:.1f}s)")
        
        ds.close()
        return mesh

    def plot_3d_field(self, field_name: str, file_index: int, threshold: float = 1e-6, 
                      ax: Optional[Axes3D] = None, vmin=None, vmax=None, cmap='viridis', 
                      alpha: float = 0.5, s: float = 1.0):
        """
        Plots a 3D scatter plot of the given field where values exceed a threshold.
        Визуализирует 3D-скаттер-график заданного поля, где значения превышают порог.
        """
        ds = self.load_data(file_index)
        
        x = ds.variables['x'][:]
        y = ds.variables['y'][:]
        z = ds.variables['z'][:]
        field = ds.variables[field_name][:] # (nx, ny, nz)
        
        # Find indices where field values exceed the threshold
        # Находим индексы, где значения поля превышают порог
        # Indices are (x_idx, y_idx, z_idx)
        x_indices, y_indices, z_indices = np.where(field > threshold)
        
        # Map indices back to actual coordinates
        # Сопоставляем индексы с фактическими координатами
        x_coords = x[x_indices]
        y_coords = y[y_indices]
        z_coords = z[z_indices]
        
        # Get corresponding field values for coloring
        # Получаем соответствующие значения поля для раскраски
        field_values = field[x_indices, y_indices, z_indices]

        if ax is None:
            raise ValueError("An Axes3D object must be provided to plot_3d_field for animation compatibility.")

        scatter = ax.scatter(x_coords, y_coords, z_coords, 
                             c=field_values, cmap=cmap, vmin=vmin, vmax=vmax, 
                             alpha=alpha, s=s, marker='o') # 's' is size of markers
        
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title(f"3D {field_name} (threshold > {threshold:.1e}) at t={ds.time:.1f}s")
        
        ds.close()
        return scatter

    def create_animation(self, field_name: str, y_index: int, output_filename: str = "animation.gif", fps: int = 10):
        """
        Creates an animation of the X-Z slice for the given field.
        """
        if not self.files:
            return

        print(f"Creating animation for {field_name}...")
        
        # Setup figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Determine global min/max for colorbar stability
        # Reading all files might be slow, so we estimate or take first/last
        # ideally we scan, but for speed let's just dynamic or fixed.
        # Let's use dynamic for now, or scan quickly.
        print("Scanning for min/max...")
        vmax = -1e9
        vmin = 1e9
        
        # Quick scan of first, middle, last
        indices_to_scan = [0, len(self.files)//2, len(self.files)-1]
        for idx in indices_to_scan:
             ds = self.load_data(idx)
             val = ds.variables[field_name][:, y_index, :]
             vmax = max(vmax, np.max(val))
             vmin = min(vmin, np.min(val))
             ds.close()
        
        # Refine for specific fields
        if field_name == 'w':
            limit = max(abs(vmin), abs(vmax))
            vmin, vmax = -limit, limit
            cmap = 'RdBu_r'
        elif field_name in ['qc', 'qr', 'qi']:
            vmin = 0
            cmap = 'Blues'
        else:
            cmap = 'viridis'

        def update(frame_idx):
            ax.clear()
            self.plot_slice_xz(field_name, y_index, frame_idx, ax=ax, vmin=vmin, vmax=vmax, cmap=cmap)
            
        ani = animation.FuncAnimation(fig, update, frames=len(self.files), interval=1000/fps)
        
        # Save
        if output_filename.endswith('.gif'):
            ani.save(output_filename, writer='pillow', fps=fps)
        else:
            # Requires ffmpeg usually
            try:
                ani.save(output_filename, writer='ffmpeg', fps=fps)
            except:
                print("FFmpeg not found, falling back to GIF")
                ani.save(output_filename.replace('.mp4', '.gif'), writer='pillow', fps=fps)
                
        print(f"Saved animation to {output_filename}")
        plt.close(fig)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="CMS Visualization Tool")
    parser.add_argument("--dir", type=str, default="output", help="Input directory")
    parser.add_argument("--field", type=str, default="w", help="Field to visualize (w, theta, qc, etc.)")
    parser.add_argument("--y", type=int, default=-1, help="Y-index for slice (default: center)")
    
    args = parser.parse_args()
    
    viz = Visualizer(args.dir)
    
    if viz.files:
        # Determine center if not specified
        ds = viz.load_data(0)
        ny = ds.dimensions['y'].size
        y_idx = args.y if args.y >= 0 else ny // 2
        ds.close()
        
        viz.create_animation(args.field, y_idx, output_filename=f"{args.field}_animation.gif")
