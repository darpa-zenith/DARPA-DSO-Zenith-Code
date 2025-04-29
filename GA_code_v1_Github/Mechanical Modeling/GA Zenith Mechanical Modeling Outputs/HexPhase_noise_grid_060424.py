# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 13:42:16 2024

@author: 10024674
"""

import numpy as np
import matplotlib.pyplot as plt
from opensimplex import OpenSimplex
from grid import grid_config
from scipy.io import savemat
import scipy.io
import scipy.interpolate
from skimage.restoration import unwrap_phase

plt.close('all')

''' PREP '''
waves = 550e-6              # mm
grid_size = grid_config     # Number of divisions in the grid. 256x256 is a typical pixel density.
R = 249                     # Dish Radius (Aperture) in millimeters

# Create a sample grid
x = np.linspace(-1 * R, R, grid_size)
z = np.linspace(-1 * R, R, grid_size)
xx, zz = np.meshgrid(x, z)  # MRS Configuration
N_D = 15                    # Sub-aps per diameter Increase this value to reduce gaps. 
                            # N_D x N_D is ~the number of iterations, and the driving number for 'i' in the loops below, so be careful.

''' Simplex Noise '''

tmp = OpenSimplex(seed=1)   # Initialize OpenSimplex with a seed
noise_scale = 10e-3         # Spatial scale of the noise. Adjust as desired
nullx, nullz = np.meshgrid(x, z)
ny = np.zeros_like(nullx)   # empty zero grid

# Generate noise to each grid point
for i in range(grid_size):
    for j in range(grid_size):
        ny[i][j] = tmp.noise2( x[i] * noise_scale, z[j] * noise_scale )
ny = ny * 5 * waves        # Apply PV for this noise

''' MASKING '''

# Aperture Masking
dR = np.sqrt(xx**2 + zz**2)
mask  = dR <= R             # So anything less than R is OK.
umask = dR <= R             # Alternative masking option

# Apply the mask to the image, zero if outside mask
masked_ny = np.where(mask, ny, 0)

k = 2.0*np.pi/waves
phase = np.angle(np.exp(1j*k*(masked_ny))) #With 2X for reflection
unwrapped_noise = unwrap_phase(phase, wrap_around=(False, False))
unwrapped_noise = np.where(umask, unwrapped_noise, 0)

# Conversions for data handling
x_flat = xx.flatten()
z_flat = zz.flatten()
masked_flat = unwrapped_noise.flatten()

''' HEXBIN PLOT '''

plt.figure(figsize=(5, 5))

hb = plt.hexbin(x_flat, -z_flat, C=masked_flat, gridsize=N_D+1, 
                reduce_C_function=np.mean, cmap='gray', 
                extent=(-(R), (R), -(R), (R) ), mincnt=1) 
                # (R - 2*R/N_D) attempts to take the mean of only the sub-aps within the total aperture. 

hex_coords = hb.get_offsets()
hex_values = hb.get_array()
hex_values = hex_values.astype(float)

''' FINAL PLOTS '''

plt.figure(figsize=(15, 3))

# Plot
plt.subplot(1, 3, 1)
plt.imshow(phase, cmap='viridis', extent=(-R, R, -R, R))
plt.title('Underlaying Noise Function - Wrapped Noise')
plt.colorbar()

plt.subplot(1, 3, 2)
plt.imshow(unwrapped_noise, cmap='viridis', extent=(-R, R, -R, R))
plt.title('Underlaying Noise Function - Unwrapped Noise')
plt.colorbar()

mapped_grid = np.zeros((grid_size, grid_size), dtype=float)

''' HEXBIN --> ARRAY '''

def hexGrab(x, z, hex_x, hex_z, r):
    geo_x = (z - hex_z) * np.sqrt(3)/3 - (x - hex_x) / 3 # Hexagon driven geometry
    geo_z = (x - hex_x) * 2/3 # for 90 hex rotation x,z are swapped
    dx = abs(geo_x)
    dz = abs(geo_z)
    dy = abs(-geo_x - geo_z)
    return max(dx, dy, dz) < r

plt.subplot(1, 3, 3)

# This complex function maps values from hexbin, back into a pixel/grid space. This is unreliable.

hex_radius = (R/N_D)/(3**.5)    # Approximate radius for each MRS element. Play around.

print("hex_coords = " + str(hex_coords))
for (i, (hex_x, hex_z)) in enumerate(hex_coords): # Searces row/col and assigns the hexbin value accociated.
    hex_z = hex_z*-1    # Invert to match coords.
    print("Grid Loop..." + str(i))
    for row in range(grid_size):
        print("Row Loop..." + str(i))
        for col in range(grid_size):
            print("Col Loop..." + str(i))
            grid_x = -R + 2 * R * col / (grid_size - 1)
            grid_z = -R + 2 * R * row / (grid_size - 1)

            if hexGrab(grid_x, grid_z, hex_x, hex_z, hex_radius): # Here we grab the hex value.
                print("if..." + str(i))
                mapped_grid[row, col] = hex_values[i]

print("Done!")

plt.title('Piston Modes - Unwrapped Phase per Sub-Ap')
plt.imshow(mapped_grid, cmap='magma', extent=(-R, R, -R, R))
plt.colorbar()
plt.show()

''' .MAT SAVE '''

n_mat = {"pz1" : mapped_grid}
savemat("pistonError.mat", n_mat)