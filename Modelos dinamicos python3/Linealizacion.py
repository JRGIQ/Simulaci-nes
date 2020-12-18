# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 11:59:13 2019

@author:Jheison Rene Gutierrez gomez, JRIQ.
"""

import numpy as np
import matplotlib.pyplot as pt
import sympy as sp

sp.init_printing()
# define symbols
x,u = sp.symbols(['x','u'])
# define equation
dxdt = -x**2 + sp.log(u)

print(sp.diff(dxdt,x))
print(sp.diff(dxdt,u))

# numeric solution with Python
import numpy as np
from scipy.misc import derivative
u = 16.0
x = 2.0
def pd_x(x):
    dxdt = -x**2 + np.sqrt(u)
    return dxdt
def pd_u(u):
    dxdt = -x**2 + np.sqrt(u)
    return dxdt

print('Approximate Partial Derivatives')
print(derivative(pd_x,x,dx=1e-4))
print(derivative(pd_u,u,dx=1e-4))

print('Exact Partial Derivatives')
print(-2.0*x) # exact d(f(x,u))/dx
print(0.5 / np.sqrt(u)) # exact d(f(x,u))/du




############ PROGRAMA PARA GRAFICAR 3D MODELO LINEAL Y NO LINEAL
#
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import numpy as np
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
## Make data.
#X = np.arange(0, 4, 0.25)
#U = np.arange(0, 20, 0.25)
#X, U = np.meshgrid(X, U)
#DXDT = -X**2 + np.sqrt(U)
#LIN = -4.0 * (X-2.0) + 1.0/8.0 * (U-16.0)
#
## Plot the surface.
#surf = ax.plot_wireframe(X, U, LIN)
#surf = ax.plot_surface(X, U, DXDT, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
#
## Customize the z axis.
#ax.set_zlim(-10.0, 5.0)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
## Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
## Add labels
#plt.xlabel('x')
#plt.ylabel('u')
#
#plt.show()