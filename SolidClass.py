#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np

class solid(object):
    'Toolkit for Finite Element structural analysis'
    def __init__(self, Body, N_ELEMENTS_S,tmax):
        """
        Iniitalizes object related variables needed for other class methods.
        
        Args:
            Body (object): A body object created from the swimmer class.
            N_ELEMENTS (int): Number of elements to discretize the solid body.
            t_max (float): Maximum thickness of the solid body.
        """
        x_mid = (Body.BF.x[:-1]+Body.BF.x[1:])/2
        
        self.Nelements = N_ELEMENTS_S        
        self.Nnodes = self.Nelements + 1
        self.xp_0 = np.copy(Body.BF.x)
        self.zp_0 = np.copy(Body.BF.z)
        self.pivotPoint = (0.5*tmax) / (max(self.xp_0) - min(self.xp_0))
        self.nodes = np.zeros((self.Nnodes,3))
        self.nodesNew = np.zeros((self.Nnodes,3))
        self.nodes_0 = np.zeros((self.Nnodes,3))
        self.tempNodes = np.zeros((self.Nnodes,3))
        self.tBeam = np.zeros((self.Nelements,1))
        self.ttemp = np.zeros((self.Nelements,1))
        self.tBeamStruct = np.zeros((self.Nelements,1))
        self.meanline_p0 = Body.BF.x / (np.max(Body.BF.x) - np.min(Body.BF.x))
        self.meanline_c0 = x_mid / (np.max(Body.BF.x) - np.min(Body.BF.x))
        self.fixedCounter = 0
        self.beamCounter = 0
        self.elemAsso = np.zeros(Body.BF.x.shape, int)
        
    def initThinPlate(self, tmax, c, SW_CNST_THK_BM, tConst, FLEX_RATIO):
        """
        This function initializes the element nodal positions.
        
        Args:
            tmax (float): NumPy array of solid thicknesses for each element.
            c (float): Undeformed/initial chord length.
            tConst (float):constant thickness position (0 <= tConst <= 1)
            FLEX_RATIO (float): Percent of the body to remain rigid as measured
                from the leading edge.
        """
        for i in xrange(self.Nelements):
            if self.nodes[i,0] <= 0.5*tmax:
                self.tBeam[i,0] = np.copy(tmax)
                self.tBeamStruct[i,0] = 0.1*np.copy(self.tBeam[i,0])
                self.ttemp[i,0] = 0.1*np.copy(self.tBeam[i,0])
                self.beamCounter += 1
                self.fixedCounter += 1
            elif self.nodes[i,0] >= c-0.5*tmax:
                self.tBeam[i,0] = 2*np.sqrt((0.5*tmax)**2 -(self.nodes[i,0]-(c-0.5*tmax))**2 )
                self.tBeamStruct[i,0] = 0.1*np.copy(self.tBeam[i,0])
                self.ttemp[i,0] = np.copy(self.tBeam[i,0])
            else:
                self.ttemp[i,0] = 0.1*np.copy(tmax)
                self.tBeam[i,0] = np.copy(tmax)
                self.tBeamStruct[i,0] = 0.1*np.copy(self.tBeam[i,0])
                if (SW_CNST_THK_BM == 1 and  self.nodes[i,2] >= tConst):
                    self.tBeamStruct[i,0] = np.copy(self.tBeamStruct[i-1,0])
                else:
                    self.tBeamStruct[i,0] = 0.1*np.copy(self.tBeam[i,0])
                if (self.nodes[i,2] <= FLEX_RATIO):
                    self.fixedCounter += 1
                    
    def initTearDrop(self, tmax, c, SW_CNST_THK_BM, tConst, FLEX_RATIO):
        """
        This function initializes the element nodal positions.
        
        Args:
            tmax (float): NumPy array of solid thicknesses for each element.
            c (float): Undeformed/initial chord length.
            constThickBeam (bool): Flag argument for constant thickness properties 
            tConst (float):constant thickness position (0 <= tConst <= 1)
            FLEX_RATIO (float): Percent of the body to remain rigid as measured
                from the leading edge.
        """
        for i in xrange(self.Nelements):
            if self.nodes[i,0] <= 0.5*tmax:
                self.tBeam[i,0] = np.copy(tmax)
                self.tBeamStruct[i,0] = np.copy(self.tBeam[i,0])
                self.ttemp[i,0] = np.copy(self.tBeam[i,0])
                self.beamCounter += 1
                self.fixedCounter += 1
            else:
                self.ttemp[i,0] = -tmax / (c - 0.5*tmax) * self.nodes[i,0] + tmax / (c - 0.5*tmax) * c
                self.tBeam[i-1,0] = 0.5 * (self.ttemp[i,0] + self.ttemp[i-1,0])
                self.tBeamStruct[i-1,0] = np.copy(self.tBeam[i-1,0])
                if i == self.Nelements-1:
                    self.tBeam[i,0] = 0.5 * self.ttemp[i,0]
                    self.tBeamStruct[i,0] = np.copy(self.tBeam[i,0])
                    if (SW_CNST_THK_BM == 1 and self.nodes[i,2] >= tConst):
                        self.tBeamStruct[i,0] = np.copy(self.tBeamStruct[i-1,0])
                if (SW_CNST_THK_BM == 1 and  self.nodes[i,2] >= tConst):
                    self.tBeamStruct[i,0] = np.copy(self.tBeamStruct[i-1,0])
                else:
                    self.tBeamStruct[i-1,0] = np.copy(self.tBeam[i-1,0])
                if (self.nodes[i,2] <= FLEX_RATIO):
                    self.fixedCounter += 1

    def initMesh(self):
        """
        Initializes the finite element mesh based on the object's __init__
        values. This is only valid for the undeformed structure at time t = 0.
        """
        self.nodes[:,0] = np.arange(min(self.xp_0),max(self.xp_0)+(max(self.xp_0)-min(self.xp_0))/self.Nelements,\
                                (max(self.xp_0)-min(self.xp_0))/self.Nelements)
        self.nodes[:,1] = np.zeros((self.Nnodes,1)).T
        self.nodes[:,2] = self.nodes[:,0] / (max(self.nodes[:,0])-min(self.nodes[:,1]))
        
        self.nodes_0 = np.copy(self.nodes)      
        self.nodesNew = np.copy(self.nodes) 
        
        for i in xrange(self.xp_0.shape[0]):
            for j in xrange(self.Nelements):
                if (self.xp_0[i] <= self.nodes_0[j+1,0] and self.xp_0[i] >= self.nodes_0[j,0]):
                    self.elemAsso[i] = np.copy(j)
    
    #TODO: This function exists in FSIClass.py. Remove it in a future release?    
    def rotatePts(x0, y0, theta):
        """
        Rotates a pair of points a specified angle.
        
        Args:
            x0 (float): A NumPy array of x-coordinates to be rotated.
            y0 (float): A NumPy array of z-coordinates to be rotated.
            theta (float): Angle (in radians) to rotate coordinates.
            
        Returns:
            x (float): Rotated x-coordinate NumPy array theta radians from x0.
            y (float): Rotated y-coordinate NumPy array theta radians from y0.
        """
        x = x0 * np.cos(theta) - y0 * np.sin(theta)
        y = x0 * np.sin(theta) + y0 * np.cos(theta)
        return (x, y)