# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 14:24:37 2023

For use on preping grid-sag data .DAT file

.DAT must be formatted as:

    nx ny delx dely unitflag xdec ydec
    z dz/dx dz/dy d2z/dxdy nodata
    .
    .
    .

@author: 10024674
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
R = 249 # Dish Radius (Aperture) in millimeter. Again, has to mach unit system in .txt
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


''' SAVE AS .DAT FOR ZEMAX GRID SAG '''

# Convert X,Z units to Zemax's X,Y Millimeters
zemax_x = x
zemax_y = z
zemax_z = masked_FEA

# Calculate the first-order derivatives
dz_dx, dz_dy = np.gradient(zemax_z, axis=(1, 0))
nx, ny = zemax_z.shape  # Grid dimensions

delx = (zemax_x[-1] - zemax_x[0]) / (nx - 1)  # Spacing in X dimension (was X)
dely = (zemax_y[-1] - zemax_y[0]) / (ny - 1)  # Spacing in Y dimension (was Z)

# Normalize the derivatives by the grid spacing
dz_dx /= delx
dz_dy /= dely

# Calculate the mixed derivative
d2z_dxdy = np.gradient(dz_dx, axis=0) / dely

unitflag = 0  # Set according to your units (e.g., 0 for millimeters)
xdec, ydec = 0, 0  # Decentering parameters
xprc, yprc = 6, 6  # Decimal precision for grid sag.

# Replace NaNs in derivative arrays, with threhold to avoid neg-zero terms
threshold = 1e-9
dz_dx[np.abs(dz_dx) < threshold] = 0
dz_dy[np.abs(dz_dy) < threshold] = 0
d2z_dxdy[np.abs(d2z_dxdy) < threshold] = 0

# Prepare the header
header = f"{nx} {ny} {delx:.{xprc}f} {dely:.{yprc}f} {unitflag} {xdec} {ydec}\n"

# Write to file
dat_path = file_name + "_" + str(grid_size) +"X" + str(grid_size)+ "_GRID_SAG.dat"
with open(dat_path, "w") as file:
    # Write the header
    file.write(header)
    
    # Write the data lines
    for i in range(nx):
        for j in range(ny):
            
            disp = zemax_z[i, j]
            grad_x = 0 if np.isnan(dz_dx[i, j]) else dz_dx[i, j]
            grad_y = 0 if np.isnan(dz_dy[i, j]) else dz_dy[i, j]
            grad_xy = 0 if np.isnan(d2z_dxdy[i, j]) else d2z_dxdy[i, j]
            
            line = f"{disp:.{yprc}f} {grad_x:.{yprc}f} {grad_y:.{yprc}f} {grad_xy:.{yprc}f} 0\n"
            file.write(line)

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
plt.contourf(xx, zz, masked_FEA, levels=10, cmap='rainbow')
plt.colorbar(label='Total Displacement (mm)')

# Add horizontal line at the Z position of the X-axis line out
plt.axhline(2*R*(x_index/grid_size) - R, color='red', linewidth=2)

# Add vertical line at the X position of the Z-axis line out
plt.axvline(2*R*(z_index/grid_size) - R, color='blue', linewidth=2)

plt.title('Interpolated Nodal Results')
plt.grid('on')
plt.xlabel('X (m)')
plt.ylabel('Z (m)')
plt.show()

# .DAT terms plotting, for insepction
plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.imshow(zemax_z, origin='lower', extent=[0, nx, 0, ny], aspect='auto')
plt.colorbar(label='z (mm)')

plt.subplot(2, 2, 2)
plt.imshow(dz_dx, origin='lower', extent=[0, nx, 0, ny], aspect='auto')
plt.colorbar(label='dz_dx (mm)')

plt.subplot(2, 2, 3)
plt.imshow(dz_dy, origin='lower', extent=[0, nx, 0, ny], aspect='auto')
plt.colorbar(label='dz_dy (mm)')

plt.subplot(2, 2, 4)
plt.imshow(d2z_dxdy, origin='lower', extent=[0, nx, 0, ny], aspect='auto')
plt.colorbar(label='d2z_dxdy (mm)')
