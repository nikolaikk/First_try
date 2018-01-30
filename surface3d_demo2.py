'''
========================
3D surface (solid color)
========================

Demonstrates a very basic plot of a 3D surface using a solid color.
'''
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

# Make data
u = np.linspace(0,  np.pi, 100)
v = np.linspace(0, 0.5*np.pi, 100)
x = 10 * np.outer(np.cos(u), np.sin(v))
y = 10 * np.outer(np.sin(u), np.sin(v))
z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface

fig.add_subplot(111, projection='3d').plot_surface(x, y, z, color='b')

plt.show()




#=========================================================================

#import numpy as np
##from mpl_toolkits.mplot3d import Axes3Ds
#import matplotlib.pyplot as plt
#
#
#def fun(x, y):
#  return x**2 + y
#
#fig = plt.figure()S
#ax = fig.add_subplot(111, projection='3d')
#x = y = np.arange(-3.0, 3.0, 0.05)
#X, Y = np.meshgrid(x, y)
#zs1 = fun(x,y)
#zs = np.array([fun(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
#Z = zs.reshape(X.shape)
#
#ax.plot_surface(X, Y, Z)
#
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')
#
#plt.show()

#=========================================================================

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
                       
# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
plt.show()