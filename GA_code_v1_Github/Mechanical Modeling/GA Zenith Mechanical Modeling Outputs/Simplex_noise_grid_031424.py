# -*- coding: utf-8 -*-
"""
Generates a Simplex noise pattern, and saves as a .mat

Also plots the results for inspection/post analysis.

@author: 10024674
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
from opensimplex import OpenSimplex
from grid import grid_config

''' PREP '''

# Create grid
grid_size = grid_config # Number of divisions in the grid.
R = 249 # Dish Radius (Aperture) in millimeter
x = np.linspace(-1*R, R, grid_size)
z = np.linspace(-1*R, R, grid_size)
xx, zz = np.meshgrid(x, z)

# Aperture Masking
dR = np.sqrt(xx**2 + zz**2)
mask = dR <= R # So anything less than R is OK.

nPV = 0.1 # Peak to Vally Scale for Noise, units in mm

''' Simplex Noise '''

tmp = OpenSimplex(seed=1) # Initialize OpenSimplex with a seed
noise_scale = 7e-3  # Spatial scale of the noise. Adjust as desired
nullx, nullz = np.meshgrid(x, z)
ny = np.zeros_like(nullx)

# Generate noise to each grid point
for i in range(grid_size):
    for j in range(grid_size):
        ny[i][j] = tmp.noise2( x[i] * noise_scale, z[j] * noise_scale )
ny = ny * nPV # Apply PV for this noise

# Apply the mask to the interpolated data, fill mask with zeros

Z_fill = np.ma.masked_where(mask, ny*0)
masked_n1 = np.ma.masked_where(~mask, ny)

Noise_Grid = np.ma.array(Z_fill.filled(1)*masked_n1.filled(1))

''' .MAT SAVE '''

n_mat = {"nz0":Noise_Grid}
savemat("simplexNoise.mat", n_mat)

# Flatten the masked arrays to calc statitics
flat_noise1 = masked_n1.compressed()

''' STATISTICS '''

# Calculate RMS (Ensure Correct Unit System)
waves = 550 # nm
allowed_rms = 1e-6*(waves/6.0) #Convert to mm for percentage calc!


PV_n1 = np.max(flat_noise1) - np.min(flat_noise1)
rms_n1 = np.sqrt(np.mean(flat_noise1 ** 2))
perc_n1= 100 * rms_n1/allowed_rms

print("")
print("Total RMS [nm]:", rms_n1*1e6)
print("Total RMS [waves]:", rms_n1*1e6/waves)
print("")
print("PV Total [nm]:", PV_n1*1e6)
print("PV Total [waves]:", PV_n1*1e6/waves)
print("")

# Function to take a line out along X or Z axis (For nice center spacing)
def take_line_out(axis, index, xx, zz, mapped_displacement):
    if axis == 'x':
        # Take a line out along the X axis at a given Z index
        x_line = xx[index, :]
        values = mapped_displacement[index, :]
    elif axis == 'z':
        # Take a line out along the Z axis at a given X index
        z_line = zz[:, index]
        values = mapped_displacement[:, index]
    else:
        raise ValueError("Axis must be 'x' or 'z'")
    return (x_line if axis == 'x' else z_line), values

# Select the index for the line out (For nice center spacing)
x_index = grid_size-int(grid_size/2)
z_index = grid_size-int(grid_size/2)

# Take line outs across mapped grid
x_line, x_values = take_line_out('x', z_index, xx, zz, masked_n1)
z_line, z_values = take_line_out('z', x_index, xx, zz, masked_n1)


''' PLOTTING '''

# Plotting the line outs
plt.figure(figsize=(5, 5))

# Plot for X line out
plt.plot(x_line, x_values, linewidth=3, color='red', label='Line Out Along X Axis')

# Plot for Z line out
plt.plot(z_line, z_values, linewidth=3, color='blue', label='Line Out Along Z Axis')
plt.title('Line Out Along Dish')
plt.xlabel('Distance (mm)')
plt.ylabel('Nodel Total Disp. Value')

plt.legend()

# Contour Plot of Interpolated Data
plt.figure(figsize=(6, 5))
plt.contourf(xx, zz, Noise_Grid, levels=10, cmap='rainbow')
plt.colorbar(label='Total Displacement (mm)')

# Add horizontal line at the Z position of the X-axis line out
plt.axhline(2*R*(x_index/grid_size) - R, color='red', linewidth=2)

# Add vertical line at the X position of the Z-axis line out
plt.axvline(2*R*(z_index/grid_size) - R, color='blue', linewidth=2)

plt.title('Simplex Noise')
plt.grid('on')
plt.xlabel('X (mm)')
plt.ylabel('Z (mm)')
plt.show()
