
## 1D investigation  
## this code has been a very quick endavour, not beautiful
## it has function for the complete mirror and fonction for ROI (smaller domain)

import numpy as np
import matplotlib.pyplot as plt


def gaussian(r,r0,s=1): # definition
    return np.exp(-(r-r0)**2/(2*s**2))

def fig_gaus(dist,n):   # display geometry of Gaussian distribution
    cmap = plt.get_cmap('coolwarm')
    j = 0
    for i in dist:
        color = cmap(j / len(dist))
        for nn in np.linspace(i[0],i[-1],n,endpoint=True):
            plt.plot(i, gaussian(i,0), linewidth=1, color="k")
            plt.plot(i, gaussian(i,nn), linewidth=1, color=color)
        j += 1
    plt.title(f"{n:.0f} Gaussian responses")
    plt.xlabel(r"Distance between Gaussian [$\sigma$]")
    plt.ylabel("Normalized Intensity [-]")
    plt.show()
    
    for i in dist:
        array = np.zeros_like(i)
        for nn in np.linspace(i[0],i[-1],n,endpoint=True):
            array = array + gaussian(i,nn)
        plt.plot(i, array, linewidth=1, label= f"D = {np.max(i)*2/(n-1):.2f} $\sigma$")
    plt.title(f"Sum of response of {n:.0f} Gaussians")
    plt.xlabel(r"Distance between Gaussian [$\sigma$]")
    plt.ylabel("Normalized Intensity [-]")
    plt.legend()
    plt.show()
    
def fig_gaus_ROI(dist,dist_ROI,n): # display geometry of Gaussian distribution over smaller domain
    cmap = plt.get_cmap('coolwarm')
    j = 0
    x_lim=0
    for d in range(len(dist)):
        color = cmap(j / len(dist))
        for nn in np.linspace(dist[d][0],dist[d][-1],n,endpoint=True):
            # print(nn)
            plt.plot(dist[d], gaussian(dist[d],0), linewidth=1, color="k")
            plt.plot(dist[d], gaussian(dist[d],nn), linewidth=1, color=color)
            if x_lim <= dist_ROI[d][-1]:
                x_lim = dist_ROI[d][-1]
        j += 1
    plt.xlim(-x_lim,x_lim)
    plt.title(f"{n:.0f} Gaussian responses")
    plt.xlabel(r"Distance between Gaussian [$\sigma$]")
    plt.ylabel("Normalized Intensity [-]")
    # plt.legend()
    plt.show()
    
    for d in range(len(dist)):
        array = np.zeros_like(dist[d])
        for nn in np.linspace(dist[d][0],dist[d][-1],n,endpoint=True):
            array = array + gaussian(dist[d],nn)
        plt.plot(dist[d], array, linewidth=1, label= f"D = {np.max(dist[d])*2/(n-1):.2f} $\sigma$")
    
    # plt.xlim(-x_lim,x_lim)
    plt.title(f"Sum of response of {n:.0f} Gaussians")
    plt.xlabel(r"Distance between Gaussian [$\sigma$]")
    plt.ylabel("Normalized Intensity [-]")
    plt.legend()
    plt.show()

def response(pos,distance,s1=1,s2=1): # give mean and RMSE of the distribution
    func = np.zeros_like(distance)
    for p in pos:
        func += gaussian(distance,p,s1)
    mean = np.mean(func)
    # erms = np.sqrt( np.trapz((func-mean)**2,distance)/np.trapz(1+0*distance,distance) ) 
    erms = np.sqrt(np.mean(np.square(np.subtract(func,mean)))) # faster, same result
    return mean, erms

def distance_array(xmin,xmax,step, n): # choosing a smaller domain for evaluation of performance
    dist = []
    if xmin == xmax:
        dist.append(np.linspace(-xmin*(n-1)/2, xmin*(n-1)/2, 1000,endpoint=True))
    else:
        for i in np.arange(xmin,xmax,step):
            dist.append(np.linspace(-i*(n-1)/2, i*(n-1)/2, 1000,endpoint=True))
    return dist

def distance_array_ROI(xmin,xmax,step, n, frac): # choosing a smaller domain for evaluation of performance
    dist = []
    if xmin == xmax:
        dist.append(np.linspace(-frac*xmin*(n-1)/2, frac*xmin*(n-1)/2, 1000,endpoint=True))
    else:
        for i in np.arange(xmin,xmax,step):
            dist.append(np.linspace(-frac*i*(n-1)/2, frac*i*(n-1)/2, 1000,endpoint=True))
    return dist

def pos_gaussian(dist, n):  # gives the center position of the gaussian
    return np.linspace(dist[0],dist[-1],n,endpoint=True)

def response_over_array(dist,n): #call functions 
    mean = []
    erms = []
    dd = []
    for d in dist:
        pos_gaus = pos_gaussian(d, n)
        res = response(pos_gaus, d, s1=1,s2=1)
        mean.append(res[0])
        erms.append(res[1])
        dd.append(np.max(d)*2/(n-1))
    return [dd,mean,erms]

def response_over_array_ROI(dist,dist_ROI,n): #to compute over smaller domain
    mean = []
    erms = []
    dd = []
    for d in range(len(dist)):
        pos_gaus = pos_gaussian(dist[d], n)
        res = response(pos_gaus, dist_ROI[d], s1=1,s2=1)
        mean.append(res[0])
        erms.append(res[1])
        dd.append(np.max(dist[d])*2/(n-1))
    return [dd,mean,erms]




# main

n_gaussian = 7  #number of Gaussians (impair/odd)
 
d_min = np.sqrt(2*np.log(2)) # assuming diameter actuator = FWHM Gaussian hard limit
d_max = 2.5                  # maximum distance between adjacent Gaussian


# display Gaussian
distx = distance_array(d_min,d_max,0.15,n_gaussian)
fig_gaus(distx,n_gaussian)

#computation
dist = distance_array(d_min,d_max,0.1,n_gaussian)
data = response_over_array(dist,n_gaussian)

# plt.plot(data[0],data[2])
# plt.title(f"RMSE of response of {n_gaussian:.0f} Gaussians")
# plt.xlabel(r"Distance between Gaussians [$\sigma$]")
# plt.ylabel("RMS error [-]")
# plt.show()

#find best filling for mirror
fraction = []
value = []
for frac in np.arange(0.5,1.01,0.1):
    dist_ROI = distance_array_ROI(d_min,d_max,0.1,n_gaussian,frac)
    data_ROI = response_over_array_ROI(dist,dist_ROI,n_gaussian)
    
    plt.plot(data[0],data_ROI[2], label=f"{frac:.2f}")
    
    fraction.append(frac)
    value.append(data[0][np.argmin(data_ROI[2])])  # position of minimum of the curve
    
plt.title(f"RMSE of response of {n_gaussian:.0f} Gaussians using a fraction of the total mirror")
plt.xlabel(r"Distance between Gaussians [$\sigma$]")
plt.ylabel("RMSE of normalized intensity [-]")    
plt.legend()
plt.show()

# plt.plot(data[0],data[1])
    # plt.title(f"Mean value of response of {n_gaussian:.0f} Gaussians")
    # plt.xlabel(r"Distance between Gaussians [$\sigma$]")
    # plt.ylabel("mean [-]")
    # plt.show()
    
    # plt.plot(data[0],data_ROI[1])
    # plt.title(f"Mean value of response of {n_gaussian:.0f} Gaussians flat ROI")
    # plt.xlabel(r"Distance between Gaussians [$\sigma$]")
    # plt.ylabel("mean [-]")
    # plt.show()


plt.plot(value,fraction)
plt.title(f" Fraction of the mirror to be used when {n_gaussian:.0f} Gaussians with minimising the RMSE ")
plt.xlabel(r" Distance between Gaussians [$\sigma$]")
plt.ylabel(" Fraction used of the mirror [-]")    
# plt.legend()
plt.show()


d_opt = data[0][np.argmin(data[2])]
dist_opt = distance_array(d_opt,d_opt,0.1,n_gaussian)
fig_gaus(dist_opt,n_gaussian)
print(r"minimal RMS error when distance =",d_opt,r"$\sigma $")
print(np.min(data[2]))

# frac = 0.5

# d_opt_ROI = data[0][np.argmin(data_ROI[2])]
# dist_opt_ROI = distance_array_ROI(d_opt_ROI,d_opt_ROI,0.01,n_gaussian, frac)
# fig_gaus_ROI(dist_opt_ROI,dist_opt_ROI,n_gaussian)
# print(r"minimal RMS error when distance =",d_opt_ROI,r"$\sigma $ flat ROI")


# fig_gaus_ROI(dist,dist_ROI,n_gaussian)



