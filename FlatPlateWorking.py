# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:25:42 2015

@author: biofluids4
"""
import numpy as np
import matplotlib.pyplot as plt
from TDinput import PARAMETERS as IP
from mpl_toolkits.mplot3d import Axes3D 

#Flat plate geometry
# Stepping through each spanwise position to calculate the positions of the
# fin neutral plane at the given time step.

#Parameters of the geometry
Nx = IP['Nx']
Ny = IP['Ny']
h = IP['h']
w = IP['w']    
D = 0.2*h
bw = np.copy(w)

# Rectangular shape body
c_r = np.copy(h)     #root chord 
c_ms = np.copy(h)    #mid-chord length of body
c_t = np.copy(h)     #tip-chord length of body


start = 0
stop  = np.pi
step  = np.copy(Ny)
yfluke = np.linspace(start,2*bw,step)

# Defining mid-chord line, Chord length function, leading edge and trailing edge line
# Mid-chord line definition
        
MC = (h/2*np.ones((Ny,1)))
     
# Second order approximation of the chord distribution along the span length.
        
CFluke = (2/bw**2*(c_r - 2*c_ms + c_t)*(yfluke.T - yfluke[0])**2+1/bw*(-(3*c_r - 4*c_ms + c_t))*(yfluke.T - yfluke[0])+ c_r)

TE = MC[:,0] + CFluke/2
LE = MC[:,0] - CFluke/2
C = np.copy(CFluke)

theta = np.linspace(start,stop,step)
xb = (C*np.cos(theta).T + C)/(2.)

#print xb
start = np.pi
stop  = 0
step  = np.copy(Ny)
#step  = np.copy(Nx)
theta = np.linspace(start,stop,step)
xt = (C*np.cos(theta).T + C)/(2.)

#print xt
zb = -0.5*D*np.ones(Ny)
zt =  0.5*D*np.ones(Ny)

# Circle origins ( z is assumed to be zero)
oF = 0.5 * D
oB = C[0] - 0.5 * D

# Calculate new theta positions for points on the rounded ends
count1 = np.shape(xb[xb>=oB])[0]
count2 = np.shape(xb[xb<=oF])[0]
count3 = np.shape(xt[xt<=oF])[0]
count4 = np.shape(xt[xt>=oB])[0]

thetab = np.linspace(0,np.pi,count1+count2).T
thetat = np.linspace(np.pi,2*np.pi,count3+count4).T

# Calculate transform leading and trailing edge points
x1 = oB + 0.5 * D * np.cos(thetab[0:count1])
z1 = 0  - 0.5 * D * np.sin(thetab[0:count1])
x2 = oF + 0.5 * D * np.cos(thetab[-1:-1-count1+1:-1])
z2 = 0  - 0.5 * D * np.sin(thetab[-1:-1-count1+1:-1])
x3 = oF + 0.5 * D * np.cos(thetat[0:count3])
z3 = 0  - 0.5 * D * np.sin(thetat[0:count3])
x4 = oB + 0.5 * D * np.cos(thetat[-1:-1-count3+1:-1])
z4 = 0  - 0.5 * D * np.sin(thetat[-1:-1-count3+1:-1])

# Replace x and z transformed points
xb[:count1] = x1
xb[-1:-1-count2+1:-1] = x2
xt[:count3] = x3
xt[-1:-1-count4+1:-1] = x4

zb[:count1] = z1
zb[-1:-1-count2+1:-1] = z2
zt[:count3] = z3
zt[-1:-1-count4+1:-1] = z4

zb[0] = 0
zt[0] = 0
zb[-1] = 0

# Merge top and bottom surfaces together
xb = np.hstack((xb , xb[-2::-1]))
zb = np.hstack((zb , zt[-2::-1]))

xb_col = ((xb[1:] + xb[:-1])/2)
zb_col = ((zb[1:] + zb[:-1])/2)

#thick = np.zeros((Ny,Nx))
#xf = np.zeros((Ny,Nx))
thick = np.zeros((Nx,Ny))
xf = np.zeros((Nx,Ny))

for i in xrange (Ny):
        
# Creating Ny points along the chord for each spanwise position for 
# for total of Nx*Ny points. The points are spaced more closely at the 
# leading edge and trailing edge.  
        
#    xf[:,i] = LE[i]+ C[i]*xt
    xf[i,:] = LE[i]+ C[i]*xt
#    thick[:,i] = zt
    thick[i,:] = zt
#    thick[:,-1] = np.minimum(thick[:,-1])


# Plot the figures
figure = plt.figure(1)
fig = plt.plot(xt, zt, 'bo', xb,zb, 'go')
plt.gca().set_aspect('equal')

#yff = np.ones((Ny,Nx))
yff = np.ones((Nx,Ny))

for i in xrange (Ny):
#    yff[:,i] = yfluke
    yff[i,:] = yfluke
    
fig2 = plt.figure(2)
ax = fig2.add_subplot(111, projection='3d')
surf = ax.plot_surface(xf, yff, thick, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)

