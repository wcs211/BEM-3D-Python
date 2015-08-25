# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:38:14 2015

@author: biofluids4

"""
import numpy as np
import matplotlib.pyplot as plt
from input_parameters_3D import PARAMETERS as IP
from mpl_toolkits.mplot3d import Axes3D   

    
#N = GeoVDVParameters.N
#S = GeoVDVParameters.S
#C = GeoVDVParameters.C
#K = GeoVDVParameters.K
#EPSILON = GeoVDVParameters.EPSILON
    
N_CHORD = IP['N_CHORD']
N_SPAN = IP['N_SPAN']
C = IP['C']
SPAN = IP['SPAN']
EPSILON = IP['EPSILON']
K = IP['K']

# Defining Chord length function
C = C * np.ones(N_SPAN)
# Initialize bottom and top x and z coordinates and parameters
A = np.zeros((N_CHORD, N_SPAN))   
R1 = np.zeros((N_CHORD, N_SPAN))
R2 = np.zeros((N_CHORD, N_SPAN))
THETA1 = np.zeros((N_CHORD, N_SPAN))
THETA2 = np.zeros((N_CHORD, N_SPAN))
x_body = np.zeros((N_CHORD,N_SPAN))
x = np.zeros((2*N_CHORD-1, N_SPAN))
y = np.zeros((2*N_CHORD-1, N_SPAN))
z = np.zeros((2*N_CHORD-1, N_SPAN))
z_top = np.zeros((N_CHORD, N_SPAN))
z_bot = np.zeros((N_CHORD, N_SPAN))
x_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
y_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
z_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
# Stepping through each spanwise position to calculate the positions of
# the fluke neutral plane at the given time step.
for i in range(N_SPAN):   
           
    A[:,i] = C[i]*((1+EPSILON)**(K-1))*(2**(-K))
    THETA = np.linspace(0,np.pi,N_SPAN)
    R1[:,i] = np.sqrt((A[:,i]*np.cos(THETA)-A[:,i])**2+(A[:,i]**2)*np.sin(THETA)**2)
    R2[:,i] = np.sqrt((A[:,i]*np.cos(THETA)-EPSILON*A[:,i])**2+(A[:,i]**2)*np.sin(THETA)**2)
    
    THETA1[:,i] = np.arctan2((A[:,i]*np.sin(THETA)) , (A[:,i]*np.cos(THETA)-A[:,i]))
    THETA2[:,i] = np.arctan2(A[:,i]*np.sin(THETA) ,(A[:,i]*np.cos(THETA)-EPSILON*A[:,i]))
    
    x_body[:,i] = ((R1[:,i]**K)/(R2[:,i]**(K-1)))*(np.cos(K*THETA1[:,i])*np.cos((K-1)*THETA2[:,i]) + np.sin(K*THETA1[:,i])*np.sin((K-1)*THETA2[:,i]))
    z_top[:,i] = ((R1[:,i]**K)/(R2[:,i]**(K-1)))*(np.sin(K*THETA1[:,i])*np.cos((K-1)*THETA2[:,i]) - np.cos(K*THETA1[:,i])*np.sin((K-1)*THETA2[:,i]))
    z_bot[:,i] = -((R1[:,i]**K)/(R2[:,i]**(K-1)))*(np.sin(K*THETA1[:,i])*np.cos((K-1)*THETA2[:,i]) - np.cos(K*THETA1[:,i])*np.sin((K-1)*THETA2[:,i]))
    
    x_body[:,i] = x_body[:,i]-x_body[-1,i] # Carrying the leading edge to the origin
    x_body[0,i] = C[i]
    
    z_top[0,i] = 0.
    z_bot[0,i] = 0.
    z_bot[-1,i] = 0.
    
    # Merge top and bottom surfaces together
    x[:,i] = np.hstack((x_body[:,i] , x_body[-2::-1,i]))
    z[:,i] = np.hstack((z_bot[:,i] , z_top[-2::-1,i]))
    
for i in range(2*N_CHORD-1):
    y[i,:] = np.linspace(0.,SPAN,N_SPAN)

for i in range(N_SPAN-2):
    x_mid[:,i] = ((x[1:,i] + x[:-1,i])/2)
    y_mid[:,i] = ((y[1:,i] + y[:-1,i])/2)
    z_mid[:,i] = ((z[1:,i] + z[:-1,i])/2)
    
# Plotting the tear-drop shape in 3D
figure = plt.figure(1)
plt.plot(x_body,z_top,'bo',x_body,z_bot,'go')
plt.gca().set_aspect('equal')

figure = plt.figure(2)
ax = figure.gca(projection='3d')
ax.set_aspect('equal')
surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='hot', linewidth=0, antialiased=False)
#ax.set_xlim(-1.2*SPAN, 1.2*SPAN)
#ax.set_ylim(-1.2*SPAN, 1.2*SPAN)
#ax.set_zlim(-1.2*SPAN, 1.2*SPAN)