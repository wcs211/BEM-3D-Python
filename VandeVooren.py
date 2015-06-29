# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:38:14 2015

@author: biofluids4

"""


N = GeoVDVParameters.N
S = GeoVDVParameters.S
C = GeoVDVParameters.C
K = GeoVDVParameters.K
EPSILON = GeoVDVParameters.EPSILON

A = C*((1+EPSILON)**(K-1))*(2**(-K))

THETA = np.linspace(0,np.pi,N/2+1)

R1 = np.sqrt((A*np.cos(THETA)-A)**2+(A**2)*np.sin(THETA)**2)
R2 = np.sqrt((A*np.cos(THETA)-EPSILON*A)**2+(A**2)*np.sin(THETA)**2)

THETA1 = np.arctan2((A*np.sin(THETA)) , (A*np.cos(THETA)-A))
THETA2 = np.arctan2(A*np.sin(THETA) ,(A*np.cos(THETA)-EPSILON*A))

x = ((R1**K)/(R2**(K-1)))*(np.cos(K*THETA1)*np.cos((K-1)*THETA2) + np.sin(K*THETA1)*np.sin((K-1)*THETA2))
z_top = ((R1**K)/(R2**(K-1)))*(np.sin(K*THETA1)*np.cos((K-1)*THETA2) - np.cos(K*THETA1)*np.sin((K-1)*THETA2))
z_bot = -((R1**K)/(R2**(K-1)))*(np.sin(K*THETA1)*np.cos((K-1)*THETA2) - np.cos(K*THETA1)*np.sin((K-1)*THETA2))

x = x-x[-1] # Carrying the leading edge to the origin
x[0] = C

z_top[0] = 0
z_bot[0] = 0
z_bot[-1] = 0

# Merge top and bottom surfaces together
x = np.hstack((x , x[-2::-1]))
z = np.hstack((z_bot , z_top[-2::-1]))

x_col = ((x[1:] + x[:-1])/2)
z_col = ((z[1:] + z[:-1])/2)

BodyFrameCoordinates = PC.BodyBFC(x, z, x_col, z_col)

return Body(N, S, BodyFrameCoordinates, MotionParameters)