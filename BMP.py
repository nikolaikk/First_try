import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

t = time.time()

zo = 0  # meter
zend = 0.3  # meter
Nz = 200
delta_z = (zend-zo)/Nz
#z=zo:delta_z:zend
z = np.linspace(zo,zend,Nz)

Lambda = 500e-9
no=1
n=1
k=2*np.pi*n/Lambda
ko=2*np.pi*no/Lambda

xo=-0.001
xend=-xo
delta_x=(xend-xo)/(Nz-1)
#x=np.arange(xo,xend,delta_x)
x = np.linspace(xo,xend,Nz)
wo=0.0001
f=np.zeros([Nz,Nz])
f[:,0] = np.exp(-(x/wo)**2)

 ## Define L matrix

L = np.zeros([Nz,Nz]);
for i in range(0,Nz):
     L[i,i]=-2+(ko**2-k**2)*delta_x**2 # in free space^2 - in material^2

for i in range (1,Nz):
     L[i,i-1]=1
     L[i-1,i]=1

L = (1/delta_x**2)*L


E = np.eye(Nz)

# #   Psi evolution
import scipy.io
import sys
mat = scipy.io.loadmat('const_to_simplify.mat')
const_to_simplify = mat['const_to_simplify']

#const_to_simplify = np.dot((E+(delta_z*L/(4*1j*ko))), np.linalg.inv (E-(delta_z*L/(4*1j*ko))))
#const_to_simplify[np.where(const_to_simplify<1e-15)] = 0

for i in range(0, Nz-1):
#for i in range(0, 2):
    f[:,i+1] = np.float128(np.dot(const_to_simplify,  f[:,i]))

elapsed = time.time()-t
print('Elapsed time : ',elapsed)

#f = scipy.io.loadmat('f')['f']
I = np.multiply(np.conjugate(f),f)


##============================================================================
#       PLOT 3D
##============================================================================
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D

gs = gridspec.GridSpec(1,5)
fig1 = plt.subplot(gs[0,0:2])
plt.imshow(np.real(I))
plt.xlabel('z (m)')
plt.ylabel('x (mm)')
fig1.set_xticklabels(   fig1.get_xticks()*delta_z   )
fig1.set_yticklabels(   fig1.get_yticks()*delta_x*1000   )

fig2 = plt.subplot(gs[0,2:6],projection='3d')
u = np.linspace(0,  Nz, Nz)
x = np.outer(np.ones(np.size(u)), u)
y = x.T
z = np.multiply(np.conjugate(f),f)
fig2.plot_surface(x, y, z, cmap = 'RdBu_r')

plt.show()