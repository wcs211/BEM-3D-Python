# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 21:35:10 2015

@author: biofluids4
"""
import numpy as np
import matplotlib.pyplot as plt
from TDinput import PARAMETERS as IP
from mpl_toolkits.mplot3d import Axes3D   

 
    #Tear-drop geometry

Nx = IP['Nx']
Ny = IP['Ny']
h = IP['h']
w = IP['w']
       
D = 0.2*h

        # Rectangular shape body
c_r = np.copy(h)     #root chord 
c_ms = np.copy(h)    #mid-chord length of body
c_t = np.copy(h)     #tip-chord length of body

#        # Trapezoidal shape body
#c_r = np.copy(h)     #root chord 
#c_ms = np.copy(h)    #mid-chord length of body
#c_t = np.copy(h)     #tip-chord length of body

bw = np.copy(w)
#tmax_f = np.copy(D)
        
        
        # Defining number of x positions 
Nfluke = np.copy(Ny)
yfluke = 2*bw*np.linspace(0,1,Nfluke)
yfluke.shape = (Nfluke,1)
        
        
        # Defining mid-chord line, Chord length function, leading edge and trailing edge line
        # Mid-chord line definition
        
MC = (h/2*np.ones((Ny,1)))
        
        # Second order approximation of the chord distribution along the span length.
        
CFluke = (2/bw**2*(c_r - 2*c_ms + c_t)*(yfluke - yfluke[0])**2+1/bw*(-(3*c_r - 4*c_ms + c_t))*(yfluke - yfluke[0])+ c_r)

TE = MC + CFluke/2
LE = MC - CFluke/2
C = np.copy(CFluke)

        # Stepping through each spanwise position to calculate the positions of
        # the fluke neutral plane at the given time step.
                
xb = np.linspace(np.pi,0.,Nx)
xt = np.linspace(0.,np.pi,Nx)

        # Slopes and intersects for the line segments
m = -D/2/(C[0] - D/2)
b = D/2 + D**2/4/(C[0] - D/2)

        # Tear drop shape equation.
x_c = 0.5 * (1 - np.cos(xb))
xb = x_c * C[0]
xb1 = xb[xb <= D/2]
xb2 = xb[xb > D/2]

zb2 = -m * xb2 - b
zb1 = -np.sqrt((D/2)**2 - (xb1 - D/2)**2)
zb = np.hstack((zb2, zb1))

        # Tear drop shape equation.
x_c = 0.5 * (1 - np.cos(xt))
xt = x_c * C[0]
xt1 = xt[xt <= D/2]
xt2 = xt[xt > D/2]

zt1 = np.sqrt((D/2)**2 - (xt1 - D/2)**2)
zt2 = m * xt2 + b
zt = np.hstack((zt1, zt2))
        
zb[0] = 0
zt[0] = 0
zb[-1] = 0

        # Merge top and bottom surfaces together
x = np.hstack((xb , xt[1:]))
z = np.hstack((zb , zt[1:]))

x_col = ((x[1:] + x[:-1])/2)
z_col = ((z[1:] + z[:-1])/2)

#tmax = tmax_f * np.ones((len(yfluke),1), dtype=np.int)
thick = np.zeros((Nx,Ny))
xf = np.zeros((Nx,Ny))

for i in xrange (0,Ny):
        
        # Creating Ny points along the chord for each spanwise position for 
        # for total of Nx*Ny points. The points are spaced more closely at the 
        # leading edge and trailing edge.  
    xf[:,i] = LE[i]+ C[i]*x_c
    thick[:,i] = zt
#    thick[:,-1] = np.minimum(thick[:,-1])
    
         # Plotting the tear-drop shape in 3D
figure = plt.figure(1)
plt.plot(xt,zt,'bo',xb,zb,'go')
plt.gca().set_aspect('equal')

yff = np.ones((Nx,Ny), dtype=np.int)
figure = plt.figure(2)
for i in xrange (0,Nx):

    yff[:,i] = yff[:,i]*yfluke[i]
    ax = figure.gca(projection='3d')
    surf = ax.plot_surface(xf[:i], yff[:i],thick[:i], rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
    plt.show()
    ax.set_zlim(-0.05, 0.05)
    ax.set_xlim(0, 0.2)
    ax.set_ylim(0, 1.2)
#    Axes3D.plot_surface(xf[:,i], yff[:,i],thick[:,i])
#   plot3(xf(:,i),yff(:,i),thick(:,i)) 
