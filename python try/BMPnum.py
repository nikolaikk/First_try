import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time


t = time.time()

zo = 0                    # micrometer
zend = 600                # micrometer
z_mesh = 20                # micrometer
ratio = 1
xo = -400                 # micrometer
xend = -xo                # micrometer
x_mesh = z_mesh/ratio     # micrometer
Lambda = 10               # micrometer
wo = 25                   # micrometer
no = 1

## Start Calculation

part_of_program = 0;
part_of_program = 1;

## Convergence

x = np.arange(xo,xend,x_mesh)
z = np.arange(zo,zend,z_mesh)
Nz = len(z)
Nx = len(x)

epsilon = 8.854187817e-12
c = 299792458

ko = 2*np.pi*no/Lambda
n = 1
k = 2*np.pi*n*np.ones(len(x))/Lambda
f = np.zeros([Nx,Nz])
f[:,0] = np.exp(-(x/wo)**2)
f[0,0] = 0
f[-1,0] = 0

## Definition of L matrix

L = np.zeros([Nx,Nx])
for i in range(Nx):
    L[i,i]=-2+(ko**2-k[i]**2)*x_mesh**2    # in free space^2 - in material^2
    print(L[i,i])

for i in range(1,Nx):
    L[i,i-1]=1;
    L[i-1,i]=1;
    
L = (1/x_mesh**2)*L

## Psi evolution

E = np.eye(Nx)
const_to_simplify = (np.linalg.pinv(E-x_mesh*L/(4*1j*ko))*(E+x_mesh*L/(4*1j*ko)))
#const_to_simplify(find(const_to_simplify<1e-14))=0

for i in range(Nz-1):
    f[:,i+1]= np.matmul(const_to_simplify,f[:,i])

    f[0,i+1] = 0
    f[-1,i+1] = 0

    
periodic = np.array([np.exp(-1j*ko*z),]*Nx)
f=np.multiply(f,periodic);

I = epsilon*no*c*0.5*np.abs(f)**2

#       PLOT 3D
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
#    ax = fig.add_subplot(111,projection='3d')
Z, X = np.meshgrid(z, x)


plt.contourf(Z, X, np.real(f),20, cmap='jet')
plt.show()


###============================================================================
##       PLOT 3D
###============================================================================
#import matplotlib.gridspec as gridspec
#from mpl_toolkits.mplot3d import Axes3D
#
#gs = gridspec.GridSpec(1,5)
#fig1 = plt.subplot(gs[0,0:2])
#plt.imshow(np.real(I))
#plt.xlabel('z (m)')
#plt.ylabel('x (mm)')
#fig1.set_xticklabels(   fig1.get_xticks()*delta_z   )
#fig1.set_yticklabels(   fig1.get_yticks()*delta_x*1000   )
#
#fig2 = plt.subplot(gs[0,2:6],projection='3d')
#u = np.linspace(0,  Nz, Nz)
#x = np.outer(np.ones(np.size(u)), u)
#y = x.T
#z = np.multiply(np.conjugate(f),f)
#fig2.plot_surface(x, y, z, cmap = 'RdBu_r')
#
#plt.show()