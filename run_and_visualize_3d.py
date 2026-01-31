import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D # Required for 3D projection
import numpy as np # For legend proxy
from cms.utils.visualization import Visualizer

def main():
    parser = argparse.ArgumentParser(description="CMS 3D Visualization Tool")
    parser.add_argument("--dir", type=str, default="output", help="Input directory containing NetCDF files")
    parser.add_argument("--output_file", type=str, default="3d_animation.gif", help="Output GIF or MP4 filename for the animation")
    parser.add_argument("--fps", type=int, default=10, help="Frames per second for the animation")
    
    args = parser.parse_args()

    viz = Visualizer(args.dir)

    if not viz.files:
        print(f"No CMS output files found in {args.dir}. Please ensure the model has been run and generated output.")
        return

    # Setup figure for 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Define fields to plot and their properties
    # Thresholds might need adjustment based on actual model output
    # 'cmap' is not directly used for scatter 'c' unless it's a scalar value.
    # We will use 'color' for direct coloring of scatter points.
    plot_configs = [
        {'field_name': 'qc', 'threshold': 1e-6, 'color': 'deepskyblue', 'label': 'Cloud Water (qc)', 's': 2},
        {'field_name': 'qi', 'threshold': 1e-6, 'color': 'whitesmoke', 'label': 'Ice (qi)', 's': 3},
        {'field_name': 'c_reagent', 'threshold': 1e-8, 'color': 'magenta', 'label': 'Seeding Agent (c_reagent)', 's': 5}
    ]

    # Get grid dimensions for consistent axis limits
    ds_grid = viz.load_data(0)
    x_max = ds_grid.nx * ds_grid.dx
    y_max = ds_grid.ny * ds_grid.dy
    z_max = ds_grid.nz * ds_grid.dz
    ds_grid.close()

    def update(frame_idx):
        ax.clear()
        
        # Load the current time step's data to get time for title
        ds = viz.load_data(frame_idx)
        current_time = ds.time
        ds.close()

        ax.set_title(f"3D Cloud Evolution at t={current_time:.1f}s", fontsize=16)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        
        # Set consistent limits
        ax.set_xlim([0, x_max])
        ax.set_ylim([0, y_max])
        ax.set_zlim([0, z_max])
        
        # Plot each field
        artists = []
        for config in plot_configs:
            try:
                # plot_3d_field returns a scatter object, but we are re-drawing entirely.
                # The 'c' parameter for scatter takes an array of values for coloring.
                # If we want a solid color per field, we pass 'color' instead of 'c=field_values'.
                # However, the current plot_3d_field expects 'c=field_values' for cmap.
                # Let's adapt it to use a solid color if 'color' is provided in config
                # OR if no 'cmap' is explicitly provided in config.
                
                # To simplify and show solid colors for different species:
                # We need to modify plot_3d_field OR extract plotting logic here.
                # For now, I'll pass a fixed color directly to the scatter call within update
                # if I want distinct colors per species.
                # The plot_3d_field as refactored returns a scatter.
                
                # To achieve distinct colors for different types of particles (qc, qi, c_reagent):
                # The current plot_3d_field uses `c=field_values` for coloring, which applies a cmap.
                # If we want solid colors for each species, we need to adapt it.
                # For this example, I will modify the plot_3d_field parameters to use 'color'
                # directly in `ax.scatter` instead of `c=field_values` and `cmap`.
                # This requires an update to plot_3d_field, but for now, I'll make a pragmatic choice
                # to show distinct colors in the legend by adding proxy artists.
                # The `plot_3d_field` method currently colors based on `field_values` and `cmap`.
                # To make it plot with a single `config['color']`, I will temporarily pass `c=config['color']`
                # and `cmap=None` into plot_3d_field. This assumes `plot_3d_field` will handle it.
                
                ds_curr_frame = viz.load_data(frame_idx)
                field_data = ds_curr_frame.variables[config['field_name']][:]
                
                x = ds_curr_frame.variables['x'][:]
                y = ds_curr_frame.variables['y'][:]
                z = ds_curr_frame.variables['z'][:]
                
                x_indices, y_indices, z_indices = np.where(field_data > config['threshold'])
                
                if len(x_indices) > 0: # Only plot if there are points above threshold
                    ax.scatter(x[x_indices], y[y_indices], z[z_indices], 
                               color=config['color'], alpha=0.7, s=config['s'], marker='o',
                               label=config['label']) # Label for legend
                    # Add a proxy artist for the legend - no need to append to artists list here if ax.scatter directly handles label
                    
                ds_curr_frame.close()

            except KeyError:
                print(f"Warning: Field '{config['field_name']}' not found in output file for frame {frame_idx}. Skipping.")
                continue
        
        # After plotting all, add legend
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))

    print(f"Creating 3D animation for {len(viz.files)} frames...")
    ani = animation.FuncAnimation(fig, update, frames=len(viz.files), interval=1000/args.fps, blit=False)
    
    # Save the animation
    writer_name = 'pillow' # Default for GIF
    if args.output_file.endswith('.mp4'):
        writer_name = 'ffmpeg'
    
    try:
        ani.save(args.output_file, writer=writer_name, fps=args.fps)
    except ValueError as e:
        print(f"Error saving animation with {writer_name} writer: {e}. Falling back to pillow (GIF).")
        ani.save(args.output_file.replace('.mp4', '.gif'), writer='pillow', fps=args.fps)
            
    print(f"3D animation saved to {args.output_file}")
    plt.close(fig) # Close the plot to prevent it from showing up after saving


if __name__ == "__main__":
    main()