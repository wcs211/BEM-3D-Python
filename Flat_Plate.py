
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:38:14 2015

@author: biofluids4

"""

import numpy as np
import matplotlib.pyplot as plt
from input_parameters_3D import PARAMETERS as IP
from mpl_toolkits.mplot3d import Axes3D 


# Geometry parameters
N_CHORD = IP['N_CHORD']
N_SPAN = IP['N_SPAN']
C = IP['C']
SPAN = IP['SPAN']     
D = 0.2*C

# Defining mid-chord line
MC = 0.5 * C * np.ones(N_SPAN)

# Defining Chord length function
C = C * np.ones(N_SPAN)

# Defining Leading and Trailing Edge Line
TE = MC + 0.5 * C
LE = MC - 0.5 * C

# Initialize bottom and top x and z coordinates
xb = np.zeros((N_CHORD, N_SPAN))
xt = np.zeros((N_CHORD, N_SPAN))
zb = np.zeros((N_CHORD, N_SPAN))
zt = np.zeros((N_CHORD, N_SPAN))
x  = np.zeros((2 * N_CHORD - 1, N_SPAN))
y  = np.zeros((2 * N_CHORD - 1, N_SPAN))
z  = np.zeros((2 * N_CHORD - 1, N_SPAN))
x_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
y_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
z_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
oLE = np.zeros(N_SPAN)
oTE = np.zeros(N_SPAN)

for i in range(N_SPAN):
    # Create bottom x-corrdinates
    start = 0.
    stop  = np.pi
    step  = np.copy(N_CHORD)
    theta = np.linspace(start,stop,step)
    xb[:,i]    = 0.5 * (C[i] * np.cos(theta) + C[i])
    
    # Create top x-corrdinates
    start = np.pi
    stop  = 0.
    step  = np.copy(N_CHORD)
    theta = np.linspace(start,stop,step)
    xt[:,i]    = 0.5 * (C[i] * np.cos(theta) + C[i])
    
    # Create bottom and top z-coordinates
    zb[:,i] = -0.5 * D * np.ones(N_CHORD)
    zt[:,i] =  0.5 * D * np.ones(N_CHORD)

    # LE and TE circle origins (z is assumed to be zero)
    oLE[i] = 0.5 * D
    oTE[i] = C[i] - 0.5 * D

    # Calculate new theta positions for points on the rounded ends
    count1 = np.shape(xb[xb[:,i]>=oTE[i],i])[0]
    count2 = np.shape(xb[xb[:,i]<=oLE[i],i])[0]
    count3 = np.shape(xt[xt[:,i]<=oLE[i],i])[0]
    count4 = np.shape(xt[xt[:,i]>=oTE[i],i])[0]

    # Determine angular new positions along the rounded ends 
    thetab = np.linspace(0,np.pi,count1+count2).T
    thetat = np.linspace(np.pi,2*np.pi,count3+count4).T

    # Calculate transformed leading and trailing edge points
    x1 = oTE[i] + 0.5 * D * np.cos(thetab[0:count1])
    z1 = 0.     - 0.5 * D * np.sin(thetab[0:count1])
    x2 = oLE[i] + 0.5 * D * np.cos(thetab[-1:-1-count1+1:-1])
    z2 = 0.     - 0.5 * D * np.sin(thetab[-1:-1-count1+1:-1])
    x3 = oLE[i] + 0.5 * D * np.cos(thetat[0:count3])
    z3 = 0.     - 0.5 * D * np.sin(thetat[0:count3])
    x4 = oTE[i] + 0.5 * D * np.cos(thetat[-1:-1-count3+1:-1])
    z4 = 0.     - 0.5 * D * np.sin(thetat[-1:-1-count3+1:-1])

    # Replace x transformed points
    xb[:count1,i]           = x1
    xb[-1:-1-count2+1:-1,i] = x2
    xt[:count3,i]           = x3
    xt[-1:-1-count4+1:-1,i] = x4

    # Replace z transformed points
    zb[:count1,i]           = z1
    zb[-1:-1-count2+1:-1,i] = z2
    zt[:count3,i]           = z3
    zt[-1:-1-count4+1:-1,i] = z4

    # Make sure LE and TE points are correctly enforced (no round-off error)
    zb[0,i]  = 0.
    zt[0,i]  = 0.
    zb[-1,i] = 0.

    # Merge top and bottom surfaces together
    x[:,i] = np.hstack((xb[:,i] , xb[-2::-1,i]))
    z[:,i] = np.hstack((zb[:,i] , zt[-2::-1,i]))

    # Shift points to correct LE position
    x[:,i] += LE[i]

# Create top y-corrdinates   
for i in range(2 * N_CHORD - 1):
    y[i,:] = np.linspace(0., SPAN, N_SPAN)
 
# Define panel mid-points   
for i in range(N_SPAN - 2):
    x_mid[:,i] = 0.25 * (x[1:,i] + x[:-1,i] + x[1:,i+1] + x[:-1,i+1])
    y_mid[:,i] = 0.25 * (y[1:,i] + y[:-1,i] + y[1:,i+1] + y[:-1,i+1])
    z_mid[:,i] = 0.25 * (z[1:,i] + z[:-1,i] + z[1:,i+1] + z[:-1,i+1])   

# Plot the profile points to check    
figure = plt.figure(1)
fig = plt.plot(xt, zt, 'bo', xb,zb, 'go')
plt.gca().set_aspect('equal')

# Plot the 3D Surface
fig2 = plt.figure(2)
ax = fig2.add_subplot(111, projection='3d')
ax.set_aspect('equal')
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
#ax.set_xlim(-1.2*SPAN, 1.2*SPAN)
#ax.set_ylim(-1.2*SPAN, 1.2*SPAN)
#ax.set_zlim(-1.2*SPAN, 1.2*SPAN)