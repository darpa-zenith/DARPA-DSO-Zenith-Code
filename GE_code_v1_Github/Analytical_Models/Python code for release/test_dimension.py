
import numpy as np
import matplotlib.pyplot as plt

grid = 2  # 0: rectangular, 1: triangular, 2: hexagonal
circulair = True
radius = 50/2
miroir_dimension = radius*2/np.sqrt(1-(1/2)**2)  # Set your desired miroir size


def gaussian(x, y, x0, y0, s=1):
    return np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * s**2))

def positions_array(distance, n):
    positions = np.linspace(-distance / 2 , distance / 2, n, endpoint=True)
    return positions

def remove_values_outside_radius(x, y, radius):
    if circulair == True:
        distances = np.sqrt(x**2 + y**2)
        mask = distances > radius
        x[mask] = np.nan
        y[mask] = np.nan
    return x, y

def distribution_grid(positions):
    # base is rectangular
    distance_gaus = positions[1] - positions[0]
    distribution = np.meshgrid(positions, positions)

    if grid == 1:  # triangular
        odd_columns_indices = np.arange(1, len(positions), 2)
        distribution[1][:, odd_columns_indices] = distribution[1][:, odd_columns_indices] + distance_gaus / 2
    elif grid == 2:  # hexagonal
        odd_columns_indices = np.arange(1, len(positions), 2)
        distribution[1][:, odd_columns_indices] = distribution[1][:, odd_columns_indices] + distance_gaus / 2
        distribution[0] = distribution[0] * np.sqrt(1 - (1 / 2)**2)
    
    return [distribution, distance_gaus]

def fig_gaussian_point(positions):  # display geometry of Gaussian distribution
    distribution = distribution_grid(positions)
    # distributio = remove_values_outside_radius(distribution[0][0], distribution[0][1], radius)
    # print(distributio[0].size - np.count_nonzero(np.isnan(distributio[0])) )
    # print(distribution[0].size - np.count_nonzero(np.isnan(distribution[0])) )
    plt.scatter(distribution[0][0], distribution[0][1])
    if grid == 0:
        plt.title(f"Visualisation of the centered actuator on a rectangular grid with \n minimal distance between actuators {distribution[1]:.3f} ")
    if grid == 1:
        plt.title(f"Visualisation of the centered actuator on a triangular grid with \n minimal distance between actuators {distribution[1]:.3f} ")
    if grid == 2:
        plt.title(f"Visualisation of the centered actuator on a hexagonal grid with \n distance between actuators {distribution[1]:.3f} ")
    plt.xlabel(r"Miroir position [cm]")
    plt.ylabel(r"Miroir position [cm]")
    circle = plt.Circle((0, 0), radius, color='r', fill=False, linestyle='dashed')
    plt.gca().add_patch(circle)
    plt.show()

def fig_gaus(mirror,gaus_positions, n):  # display geometry of Gaussian distribution
    positions_x, positions_y = np.meshgrid(mirror, mirror)
    func = np.zeros_like(positions_x)

    distribution = distribution_grid(gaus_positions)
    dist_gaussian = distribution[1]
    # sigma = np.sqrt( - (dist_gaussian)**2 / (4* np.log(1/2))  )
    sigma = np.sqrt( - dist_gaussian**2 / (2* np.log(0.2))  )
    
    for i in range(len(distribution[0][0].flatten())):
        func += gaussian(positions_x, positions_y, distribution[0][0].flatten()[i], distribution[0][1].flatten()[i], sigma)

    min_plot = np.minimum(positions_x.min(),positions_y.min())
    max_plot = np.maximum(positions_x.max(),positions_y.max())
    
    plt.imshow(func, extent=[min_plot, max_plot, min_plot, max_plot])
    plt.colorbar()
    if grid == 0:
        plt.title(f"Response of {n:.0f} Gaussians with distance with  \n  rectangular grid of {dist_gaussian:.2f} cm of side ")
    if grid == 1:
        plt.title(f"Response of {n:.0f} Gaussians with distance with  \n  triangular grid of {dist_gaussian:.2f} cm of side ")
    if grid == 2:
        plt.title(f"Response of {n:.0f} Gaussians with distance with  \n  hexagonal grid of {dist_gaussian:.2f} cm of side ")
  
    # plt.title(f"Response of {n:.0f} Gaussians with distance with  \n  rectangular grid of {dist_gaussian:.2f} $\sigma$ of side")
    plt.xlabel(r"Miroir position [cm]")
    plt.ylabel(r"Miroir position [cm]")
    plt.show()

def response(mirror,gaus_positions, n, s1=1, s2=1):
    mean = []
    erms = []
    dd = []
    
    positions_x, positions_y = np.meshgrid(mirror, mirror)
    func = np.zeros_like(positions_x)
    
    distribution = distribution_grid(gaus_positions)
    dist_gaussian = distribution[1]
    # sigma = np.sqrt( - dist_gaussian**2 / (4* np.log(1/2))  )
    sigma = np.sqrt( - dist_gaussian**2 / (2* np.log(0.2))  )

    for i in range(len(distribution[0][0].flatten())):
        func += gaussian(positions_x, positions_y, distribution[0][0].flatten()[i], distribution[0][1].flatten()[i], sigma)
        
    m = np.mean(func)
    mean.append(m)
    erms.append(np.sqrt(np.mean(np.square(np.subtract(func, m)))))
    dd.append(dist_gaussian)

    return [dd, mean, erms]

# main
n_gaussian = 19
positions = positions_array(miroir_dimension, n_gaussian)
mirror_array = positions_array(miroir_dimension, 100)

fig_gaussian_point(positions)
fig_gaus(mirror_array,positions,n_gaussian)
data = response(mirror_array,positions, n_gaussian)

for frac in np.arange(0.35,1.05,0.1):
    mirror_array_frac = positions_array(miroir_dimension*frac, 501)
    # fig_gaus(mirror_array_frac,positions,n_gaussian)
    data = response(mirror_array_frac,positions, n_gaussian)
    print(data[0],data[1],data[2])
    plt.scatter(frac,data[2],linewidth=2)
    
plt.title(f"RMSE of response of {n_gaussian:.0f} Gaussians using a fraction of the mirror")
plt.xlabel(r" Fraction of the mirror used [-]")
plt.ylabel("RMSE of normalized intensity [-]")   
plt.show()



