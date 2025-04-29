import numpy as np
import matplotlib.pyplot as plt

def rectangular_grid(rows, cols, spacing):
    positions = []
    for i in range(rows):
        for j in range(cols):
            positions.append([i * spacing, j * spacing])
    return np.array(positions)

def triangular_grid(rows, cols, spacing):
    positions = []
    for i in range(rows):
        for j in range(cols):
            x = i * spacing
            y = j * spacing if i % 2 == 0 else j * spacing + spacing / 2
            positions.append([x, y])
    return np.array(positions)

def hexagonal_grid(rows, cols, spacing):
    positions = []
    for i in range(rows):
        for j in range(cols):
            x = i * spacing * 1.5
            y = j * spacing * np.sqrt(3)
            if i % 2 == 1:
                y += spacing * np.sqrt(3) / 2
            positions.append([x, y])
    return np.array(positions)

# Example usage:
rectangular_positions = rectangular_grid(rows=5, cols=5, spacing=1.0)
triangular_positions = triangular_grid(rows=5, cols=5, spacing=0.787)
hexagonal_positions = hexagonal_grid(rows=5, cols=5, spacing=0.5)

print("Rectangular Grid Positions:")
# print(rectangular_positions)
print(rectangular_positions[:,0])
plt.scatter(rectangular_positions[:,0],rectangular_positions[:,1])
plt.xlim(0,4)
plt.ylim(0,4)
plt.show()

print("\nTriangular Grid Positions:")
# print(triangular_positions)
plt.scatter(triangular_positions[:,0],triangular_positions[:,1])
plt.xlim(0,4)
plt.ylim(0,4)
plt.show()

# # print("\nHexagonal Grid Positions:")
# # print(hexagonal_positions)
# plt.scatter(hexagonal_positions[:,0],hexagonal_positions[:,1])
# plt.xlim(0,4)
# plt.ylim(0,4)
# plt.show()


def gaussian(x,y, x0, y0, s=1):
    return np.exp(-( (x-x0)**2 + (y-y0)**2 ) / (2 * s**2))

# def fig_gaus(dist_x, dist_y, n):
#     cmap = plt.get_cmap('coolwarm')
#     j = 0
#     for i, row in enumerate(dist_x):
#         for nn, col in enumerate(dist_y):
#             color = cmap(j / len(dist_x))
#             plt.plot(row, gaussian(row, 0), linewidth=1, color="k")
#             plt.plot(row, gaussian(row, col), linewidth=1, color=color)
#             j += 1

#     plt.title(f"{n:.0f} Gaussian responses")
#     plt.xlabel(r"Distance between Gaussian [$\sigma$]")
#     plt.ylabel("Normalized Intensity [-]")
#     plt.show()

#     for i, row in enumerate(dist_x):
#         array = np.zeros_like(row)
#         for nn, col in enumerate(dist_y):
#             array = array + gaussian(row, col)
#         plt.plot(row, array, linewidth=1, label=f"D = {np.max(row)*2/(n-1):.2f} $\sigma$")
    
#     plt.title(f"Sum of response of {n:.0f} Gaussians")
#     plt.xlabel(r"Distance between Gaussian [$\sigma$]")
#     plt.ylabel("Normalized Intensity [-]")
#     plt.legend()
#     plt.show()

# def response(pos_x, pos_y, distance_x, distance_y, s1=1, s2=1):
#     func = np.zeros_like(distance_x)
#     for px, py in zip(pos_x, pos_y):
#         func += gaussian(distance_x, px, s1) * gaussian(distance_y, py, s2)
#     mean = np.mean(func)
#     erms = np.sqrt(np.mean(np.square(np.subtract(func, mean))))
#     return mean, erms

def distance_array(xmin, xmax, step, n):
    dist = []
    if xmin == xmax:
        dist.append(np.linspace(-xmin*(n-1)/2, xmin*(n-1)/2, 1000, endpoint=True))
    else:
        for i in np.arange(xmin, xmax, step):
            dist.append(np.linspace(-i*(n-1)/2, i*(n-1)/2, 1000, endpoint=True))
    return dist

def pos_gaussian(dist, n):
    return np.linspace(dist[0], dist[-1], n, endpoint=True)

# def response_over_array(dist_x, dist_y, n):
#     mean = []
#     erms = []
#     dd = []
#     for dx, dy in zip(dist_x, dist_y):
#         pos_gaus_x = pos_gaussian(dx, n)
#         pos_gaus_y = pos_gaussian(dy, n)
#         res = response(pos_gaus_x, pos_gaus_y, dx, dy, s1=1, s2=1)
#         mean.append(res[0])
#         erms.append(res[1])
#         dd.append(np.max(dx)*2/(n-1))
#     return [dd, mean, erms]

# main
n_gaussian = 11

d_min = np.sqrt(2*np.log(2))
d_max = 2.5

# display Gaussian
dist_x = distance_array(d_min, d_max, 0.15, n_gaussian)
dist_y = distance_array(d_min, d_max, 0.15, n_gaussian)
# fig_gaus(dist_x, dist_y, n_gaussian)

# # computation
# dist_x = distance_array(d_min, d_max, 0.1, n_gaussian)
# dist_y = distance_array(d_min, d_max, 0.1, n_gaussian)
# data = response_over_array(dist_x, dist_y, n_gaussian)
# d_opt = data[0][np.argmin(data[2])]
# dist_opt_x = distance_array(d_opt, d_opt, 0.1, n_gaussian)
# dist_opt_y = distance_array(d_opt, d_opt, 0.1, n_gaussian)
# fig_gaus(dist_opt_x, dist_opt_y, n_gaussian)
# print(r"Minimal RMS error when distance =", d_opt, r"$\sigma$")
# print(np.min(data[2]))