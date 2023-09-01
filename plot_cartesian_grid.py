"""
Visualization of CMG GEM results on a 3D Cartesian grid using the SR3 reader and pyvista.

Copyright 2023 Nikolai Andrianov, nia@geus.dk

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from sr3_reader import *
import numpy as np
import pandas as pd
import pyvista as pv

# Read CMG results file
infile = 'cartesian_grid_example.sr3'
sr3 = read_SR3(infile)

# Get the timeseries data from the wells
wells_ts = get_wells_timeseries(sr3)

# Get the indices of completed cells
wells_comp = get_wells_completions(sr3)

# Get the pressures and saturations
(sp_ind, sp) = get_spatial_properties(sr3, ['PRES', 'SG', 'SW'])
t_sp = sr3.times['Days'].iloc[sp_ind]

pressure = sp['PRES'][0]
sw = sp['SW'][1]

# The dimensions of grid nodes (+1 to the cells' dimensions)
dim = np.array(sr3.grid.cart_dims) + 1

# Create the pyvista grid representation
grid = pv.ExplicitStructuredGrid(dim, sr3.grid.cells.corners)
grid = grid.compute_connectivity()
grid = grid.compute_connections()
grid = grid.compute_cell_sizes(length=False, area=False, volume=True)

# Scale the y- and z-axis for better visibility
grid.scale([1, 3, 10], inplace=True)

# Hide inactive cells
grid.hide_cells(sr3.grid.cells.inactive, inplace=True)

# Initialize the pyvista plotter
pl = pv.Plotter()

# Identify the completed cells as
for well, comp in wells_comp.items():
    points = []
    top_comp = None
    for cell in comp:
        # Get the global id of the cell
        cell_id = grid.cell_id(cell)
        # Get the coordinates of the cell center
        cv = grid.cell[cell_id].points
        center = np.mean(cv, axis=0)
        points.append(center)
        # Find the coordinates of the top completion
        if top_comp is None:
            top_comp = center
        else:
            if top_comp[2] < center[2]:
                top_comp = center

    # Set the wellhead 20% away from the top of the reservoir
    # Z-coordinate of the reservoir top (z are negative by convention)
    zf = max(grid.points[:, 2])
    dist = 1.2 * (zf - top_comp[2])
    wh = np.array([top_comp[0], top_comp[1], top_comp[2] + dist])

    # Define the vertical section of the well between the top completion and the wellhead
    vert_well = np.array([top_comp, wh])
    poly = pv.PolyData()
    poly.points = vert_well
    well_ind = np.arange(0, len(vert_well), dtype=np.int_)
    well_ind = np.insert(well_ind, 0, len(vert_well))
    poly.lines = well_ind
    tube = poly.tube(radius=20)
    pl.add_mesh(tube, show_edges=False)
    # Add the well labels
    pl.add_point_labels([wh], [well])

# Identify the variables for visualization
grid.cell_data['pressure'] = pressure
grid.cell_data['sw'] = sw

# Plot pressure
pl.add_mesh(grid, scalars='pressure', opacity=0.5, scalar_bar_args=dict(color='black'))

# Show the cells
pl.show_grid(color='black')  

# Add the axes
pl.add_axes(color='black')

# Change the background color
pl.background_color = 'white'

# Actual plotting
pl.show()

