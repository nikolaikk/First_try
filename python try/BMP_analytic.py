import numpy as np
import matplotlib.pyplot as plt 

global zR 
ratio = 1
Lambda = 600e-9   # meter
wo=0.0001

k = 2* np.pi/Lambda
x = np.linspace(-0.001,0.001,20) * ratio
z = np.linspace(0,0.3,len(x)) * ratio

zR = np.pi*wo**2/Lambda


def width(z, Lambda=Lambda, wo=wo):
#    print (wo)
#    print (Lambda)
    return wo*np.sqrt(1+(z/zR)**2)
    
def R(z, Lambda=Lambda, wo=wo):
    return np.multiply(z,(1+(zR/z)**2))
    
def psi(z):
    return np.arctan(z/zR)

# =======================================================================


#E = np.zeros([len(x),len(z)])
#for i in range(0,len(z)):
#    E[:,i] = np.multiply(np.multiply(wo/width(z[i]),np.exp(-x**2/width(z[i])**2)),np.exp(-1j*(k*z[i] + 0.5*k*x**2/R(z[i])-psi(z[i]))))  #.reshape(1,len(z))
#
#plt.imshow(np.multiply(np.conjugate(E),E))
#plt.show()


# =======================================================================

E = np.zeros([len(x),len(z)])
for i in range(0,len(z)):
    E[:,i] = np.multiply(wo/width(z[i]),np.exp(-x**2/width(z[i])**2))

plt.imshow(E**2)
plt.show()


#       PLOT 3D
from mpl_toolkits.mplot3d import Axes3D
Nz = len(z)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
u = np.linspace(0,  Nz, Nz)
x = np.outer(np.ones(np.size(u)), u)
y = x.T
z = E

# Plot the surface
ax.plot_surface(x, y, z, cmap = 'RdBu_r',shade = 'flat',linewidth = 0)
ax.set_axis_off


plt.show()