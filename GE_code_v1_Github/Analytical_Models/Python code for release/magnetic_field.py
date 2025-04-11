
# modified from https://netdenizen.com/emagnet/offaxis/off_axis_loop.html

from scipy.special import ellipk, ellipe
from numpy import pi, sqrt, linspace, NaN
from pylab import plot, xlabel, ylabel, legend, show
import numpy as np
import matplotlib.pyplot as plt

uo = 4E-7*pi     # Permeability constant - units of H/m
ur = 2.7         # relative Permeability constant - units of H/m
rho = 1210          # ferrofluid density - scaled
Bo = lambda i, a, u=uo: i*u/2./a    # Central field = f(current, loop radius, perm. constant)
al = lambda r, a: r/a               # Alpha = f(radius of measurement point, radius of loop)
be = lambda x, a: x/a               # Beta = f(axial distance to meas. point, radius of loop)
ga = lambda x, r: x/r               # Gamma = f(axial distance, radius to meas. point)
Q = lambda r, x, a: (1 + al(r,a))**2 + be(x,a)**2   # Q = f(radius, distance to meas. point, loop radius)
k = lambda r, x, a: sqrt(4*al(r,a)/Q(r,x,a))       # k = f(radius, distance to meas. point, loop radius)
K = lambda k: ellipk(k**2.0)          # Elliptic integral, first kind, as a function of k
E = lambda k: ellipe(k**2.0)          # Elliptic integral, second kind, as a function of k

# On-Axis field = f(current and radius of loop, x of measurement point)
def Baxial(i, a, x, u=uo):
    if a == 0:
        if x == 0:
            return NaN
        else:
            return 0.0
    else:
        return (u*i*a**2)/2.0/(a**2 + x**2)**(1.5)
# Axial field component = f(current and radius of loop, r and x of meas. point)
def Bz(i, a, x, r):
    if r == 0:
        if x == 0:
            return Bo(i,a)         # central field
        else:
            return Baxial(i,a,x)   # axial field
    else:                          # axial component, any location
        return Bo(i,a)*\
            (E(k(r,x,a))*((1.0-al(r,a)**2-be(x,a)**2)/(Q(r,x,a)-4*al(r,a))) + K(k(r,x,a)))\
            /pi/sqrt(Q(r,x,a))     
# Radial field component = f(current and radius of loop, r and x of meas. point)
def Br(i, a, x, r):
    if r == 0:
        return 0                   # no radial component on axis!
    else:                          # radial component, any location other than axis.
        return Bo(i,a)*ga(x,r)*\
            (E(k(r,x,a))*((1.0+al(r,a)**2+be(x,a)**2)/(Q(r,x,a)-4*al(r,a))) - K(k(r,x,a)))\
            /pi/sqrt(Q(r,x,a))
            
axiallimit = 0.2 # meters from center
radiallimit = 0.2 # maximum radius to investigate
curveqty = 101
X = linspace(0 ,axiallimit, curveqty)
R = linspace(-radiallimit, radiallimit, curveqty)


def deviation(current, dim_coil, z,r): #Eq 3.2
    return (ur-1)*(Bz(current,dim_coil,z,r)**2+ Br(current,dim_coil,z,r)**2*ur)/(uo * ur* 9.81* rho)

def deviation2(current, dim_coil, z,r): #Eq 3.2 with tilted (fixed for the tilted case - not general)
    fieldz2 =  (Bz(current,dim_coil,z,r)**2 + (Br(current,dim_coil,z,r) * -0.002/0.035)**2)/(1+(-0.002/0.035)**2)
    fieldr2 =  (Bz(current,dim_coil,z,r) * (1 - 1/np.sqrt(1+(-0.002/0.035)**2))  )**2 + (Br(current,dim_coil,z,r) * (1 + (-0.002/0.035)/np.sqrt(1+(-0.002/0.035)**2))  )**2
    return (ur-1)*(  fieldz2 + fieldr2*ur)/(uo * ur* 9.81* rho)    

plot(R, [deviation(1,0.035/2,0.03,np.abs(r)) for r in R], label="Flat")
# plot(R, [deviation(1,0.035/2,0.03-0.002*r/0.035,np.abs(r)) for r in R], label="Naive Tilted " , alpha = 0.8, linestyle= '--')
plot(R, [deviation2(1,0.035/2,0.03-0.002*r/0.035,np.abs(r)) for r in R], label="Tilted"  , alpha = 0.8, linestyle= '-.')
xlabel("Position from the center of the coil [m]")
ylabel(r" $\Omega(x,y)$ [m]")
plt.title(r" Comparison of the expected ferrofluid surface deviation ($\rho=1210$ kg/m$^3$, $\mu_r=2.7$) from the "+"\n" +"  magnetic field of a single coil (35mm diameter, 1A current, surface distance 30 mm at center) \n for a flat surface and a tilted surface of 2 mm over the coil ")
legend()
show()

# see tilted surface
plot(R, [0.03-0.001*r/(0.035/2) for r in R])
plt.vlines(-0.035/2,0,0.4)
plt.vlines(0.035/2,0,0.4)
plt.ylim(0.027,0.033)
plt.show()

























