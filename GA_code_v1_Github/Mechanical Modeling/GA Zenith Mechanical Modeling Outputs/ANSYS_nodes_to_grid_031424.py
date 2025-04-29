# -*- coding: utf-8 -*-
"""
Imports .txt file formated from ANSYS export, converts to meshgrid, then saves as .mat

Also plots the results for inspection/post analysis.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.io import savemat
from grid import grid_config

''' LOAD FEA RESULTS '''

# File path to .txt 
file_name = 'IF_OUTER_CC_ON'  # Example file path
file_type = '.txt'

file_path = file_name + file_type

# Read the nodal displacement data from the text file (last column)
# Ensure the file and units are in the same format!
data = np.loadtxt(file_path, skiprows=1, usecols=(1, 2, 3, 4))

# Extract Nodal Data
points = data[:, [0,2]] # Initial Node Positions
nodal_disp = data[:, 3]

# Create mesh grid
grid_size = grid_config # Number of divisions in the grid read from grid.py file
R = 249 # Dish Radius (Aperture) in millimeter, with 1mm subtracted to prevent edge effect. Again, has to mach unit system in .txt
x = np.linspace(-1*R, R, grid_size)
z = np.linspace(-1*R, R, grid_size)
xx, zz = np.meshgrid(x, z)

# Perform cubic spline interpolation over the regular grid
mapped = griddata(points, nodal_disp, (xx, zz), method='cubic') # Choose interpolation method.

# Aperture Masking for the results.
dR = np.sqrt(xx**2 + zz**2)
mask = dR <= R # So anything less than R is OK.

# Replace NaN values with zero (or another chosen value)
filled_mapped_disp = mapped.copy()
filled_mapped_disp[np.isnan(filled_mapped_disp)] = 0
FEA = filled_mapped_disp

# Apply the mask to the interpolated data. Do after incorporating any noise functions.
masked_FEA = np.ma.masked_where(~mask, FEA)
flat_FEA = masked_FEA.compressed() # Flatten the masked arrays to calc statitics


''' .MAT SAVE '''

n_mat = {str(file_name):FEA} # Change label as desired
savemat(str(file_name)+".mat", n_mat)

''' STATISTICS '''

# Calculate RMS (Ensure Correct Unit System!)
waves = 550 # nm
allowed_rms = 1e-6*(waves/6.0) #Convert to mm for percentage calc!

PV_FEA = np.max(flat_FEA) - np.min(flat_FEA)
rms_FEA = np.sqrt(np.mean(flat_FEA ** 2))
perc_FEA = 100 * rms_FEA/allowed_rms

print("")
print("Error Budget Components [%]: FEA = " + str(perc_FEA))
print("")
print("Total RMS [nm]:", rms_FEA*1e6)
print("Total RMS [waves]:", rms_FEA*1e6/waves)
print("")
print("PV Total [nm]:", PV_FEA*1e6)
print("PV Total [waves]:", PV_FEA*1e6/waves)
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
x_line, x_values = take_line_out('x', z_index, xx, zz, masked_FEA)
z_line, z_values = take_line_out('z', x_index, xx, zz, masked_FEA)

''' PLOTTING '''

# Plotting the line outs
plt.figure(figsize=(5, 5))

# Plot for X Z line out results
plt.plot(x_line, x_values, linewidth=3, color='red', label='Line Out Along X Axis')
plt.plot(z_line, z_values, linewidth=3, color='blue', label='Line Out Along Z Axis')

plt.title('Line Out Along Dish')
plt.xlabel('Distance (mm)')
plt.ylabel('Nodel Total Disp. Value')

plt.legend()

# Contour Plot of Interpolated Data
plt.figure(figsize=(6, 5))
plt.contourf(xx, zz, masked_FEA, levels=10, cmap='rainbow')
plt.colorbar(label='Total Displacement (mm)')

# Add horizontal line at X and Z position
plt.axhline(2*R*(x_index/grid_size) - R, color='red', linewidth=2)
plt.axvline(2*R*(z_index/grid_size) - R, color='blue', linewidth=2)

plt.title('Interpolated Nodal Results')
plt.grid('on')
plt.xlabel('X (mm)')
plt.ylabel('Z (mm)')
plt.show()
