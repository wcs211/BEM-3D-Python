# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 21:35:10 2015

@author: biofluids4
"""
import numpy as np
import matplotlib.pyplot as plt
from input_parameters_3D import PARAMETERS as IP
from mpl_toolkits.mplot3d import Axes3D   


N_CHORD = IP['N_CHORD']
N_SPAN = IP['N_SPAN']
C = IP['C']
SPAN = IP['SPAN']     
D = 0.2*C

# Defining mid-chord line
MC = 0.5 * C * np.ones(N_SPAN)

# Defining Chord length function
C = C * np.ones(N_SPAN)
TE = MC + 0.5 * C
LE = MC - 0.5 * C
# Initialize bottom and top x and z coordinates

xb = np.zeros((N_CHORD, N_SPAN))
xt = np.zeros((N_CHORD, N_SPAN))
x_c = np.zeros((N_CHORD, N_SPAN))
xb1 = np.zeros((N_SPAN))
xb2 = np.zeros((N_SPAN))
zb = np.zeros((N_CHORD, N_SPAN))
zt = np.zeros((N_CHORD, N_SPAN))
x  = np.zeros((2 * N_CHORD - 1, N_SPAN))
y  = np.zeros((2 * N_CHORD - 1, N_SPAN))
z  = np.zeros((2 * N_CHORD - 1, N_SPAN))
x_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
y_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
z_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))


# Stepping through each spanwise position to calculate the positions of
# the fluke neutral plane at the given time step.
for i in range(N_SPAN):   
           
    xb[:,i] = np.linspace(np.pi,0.,N_CHORD)
    xt[:,i] = np.linspace(0.,np.pi,N_CHORD)
    
    # Slopes and intersects for the line segments
    m = -D/2/(C[i] - D/2)
    b = D/2 + D**2/4/(C[i] - D/2)
    
            # Tear drop shape equation.
    x_c[:,i] = 0.5 * (1 - np.cos(xb[:,i]))
    xb[:,i] = x_c[:,i] * C[i]
    xb1 = xb[xb[:,i] <= D/2]
    xb2 = xb[xb[:,i] > D/2]
    
    zb2 = -m * xb2[:,i] - b
    zb1 = -np.sqrt((D/2)**2 - (xb1[:,i] - D/2)**2)
    zb[:,i] = np.hstack((zb2, zb1))
    
            # Tear drop shape equation.
    x_c[:,i] = 0.5 * (1 - np.cos(xt[:,i]))
    xt[:,i] = x_c[:,i] * C[i]
    xt1 = xt[xt[:,i] <= D/2]
    xt2= xt[xt[:,i] > D/2]
    
    zt1 = np.sqrt((D/2)**2 - (xt1[:,i] - D/2)**2)
    zt2 = m * xt2[:,i] + b
    zt[:,i] = np.hstack((zt1, zt2))
            
    zb[0,i] = 0.
    zt[0,i] = 0.
    zb[-1,i] = 0.
    
    # Merge top and bottom surfaces together
    x[:,i] = np.hstack((xb[:,i] , xt[1:,i]))
    z[:,i] = np.hstack((zb[:,i] , zt[1:,i]))
    
for i in range(2*N_CHORD-1):
    y[i,:] = np.linspace(0.,SPAN,N_SPAN)

for i in range(N_SPAN-2):
    x_mid[:,i] = ((x[1:,i] + x[:-1,i])/2)
    y_mid[:,i] = ((y[1:,i] + y[:-1,i])/2)
    z_mid[:,i] = ((z[1:,i] + z[:-1,i])/2)

# Plotting the tear-drop shape in 3D
figure = plt.figure(1)
plt.plot(xt,zt,'bo',xb,zb,'go')
plt.gca().set_aspect('equal')

figure = plt.figure(2)
ax = figure.gca(projection='3d')
ax.set_aspect('equal')
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
#ax.set_xlim(-1.2*SPAN, 1.2*SPAN)
#ax.set_ylim(-1.2*SPAN, 1.2*SPAN)
#ax.set_zlim(-1.2*SPAN, 1.2*SPAN)