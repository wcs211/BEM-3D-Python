# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 09:51:20 2015

@author: biofluids4
"""
import numpy as np

class SwimmerParameters(object):
    """A collection of parameters related to a single swimmer.
    
    Attributes:
        CE: Constant that determines the length of Edge panels.
        DELTA_CORE: Constant used in wake_rollup to avoid singularities near
            wake panels.
        SW_KUTTA: Switch for Kutta condition (explicit or unsteady).
    """
    def __init__(self, CE, DELTA_CORE,DELTA_CORE_F,SW_KUTTA):
        self.CE = CE
        self.DELTA_CORE = DELTA_CORE
        self.DELTA_CORE_F = DELTA_CORE_F
        self.SW_KUTTA = SW_KUTTA
        
class GeoVDVParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N_CHORD:     Number of chordwise panels.
        N_SPAN:      Number of spanwise panels.
        SPAN:        Half span of the rectangular body.
        C:           Chord length of the rectangular body.
        S:           Collocation point shifting coefficient.
        K, EPSILON:  Van de Vooren airfoil parameters.
    """
    def __init__(self,N_CHORD,N_SPAN,C,SPAN,S,K,EPSILON):
        self.SW_GEOMETRY = 'VDV'
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.C = C
        self.SPAN = SPAN
        self.S = S
        self.K = K
        self.EPSILON = EPSILON
        
class GeoFPParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N_CHORD:     Number of chordwise panels.
        N_SPAN:      Number of spanwise panels.
        SPAN:        Half span of the rectangular body.
        C:           Chord length of the rectangular body.
        S:           Collocation point shifting coefficient.
        K, EPSILON:  Van de Vooren airfoil parameters.
        D: Diameter of leading and trailing edge circles.
    """
    def __init__(self,N_CHORD,N_SPAN,C,SPAN,S,D):
        self.SW_GEOMETRY = 'FP'
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.C = C
        self.SPAN = SPAN
        self.S = S
        self.D = D
        
class GeoTDParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N_CHORD:     Number of chordwise panels.
        N_SPAN:      Number of spanwise panels.
        SPAN:        Half span of the rectangular body.
        C:           Chord length of the rectangular body.
        S:           Collocation point shifting coefficient.
        K, EPSILON:  Van de Vooren airfoil parameters.
        D: Diameter of leading and trailing edge circles.
    """
    def __init__(self,N_CHORD,N_SPAN,C,SPAN,S,D):
        self.SW_GEOMETRY = 'TD'
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.C = C
        self.SPAN = SPAN
        self.S = S
        self.D = D
    
class MotionParameters(object):
    """A collection of parameters related to a swimmer's path of motion.
   ###################################################################################
    Attributes:
    
        THETA_MAX: Maximum pitching angle of the body.
        HEAVE_MAX: Maximum heaving amplitude.
        T: Time.
        F: Frequency of the body's pitching motion,Hz.
        PHI: Phase of flapping time signal.
    """
    def __init__(self,X0,Y0,Z0,V0,HEAVE_MAX,THETA_MAX,F,PHI):
        self.X0 = X0
        self.Y0 = Y0
        self.Z0 = Z0
        self.V0 = V0
        self.HEAVE_MAX = HEAVE_MAX
        
        self.THETA_MAX = THETA_MAX
        self.F = F 
        self.PHI = PHI
    
class BodyBFC(object):
    """A collection of body-frame coordinates for a Body object.
    
    Attributes:
        x, z: Body-frame panel endpoint positions.
        x_col, z_col: Body-frame panel midpoint positions (unshifted).
    """
    def __init__(self, x, z, x_col, z_col):
        self.x = x
        self.z = z
        self.x_col = x_col
        self.z_col = z_col
    
class BodyAFC(object):
    """A collection of absolute-frame coordinates for a Body object.
    
    Attributes:
        x, z: Panel endpoint positions.
        x_col, z_col: Panel collocation point positions (shifted midpoints).
        x_mid, z_mid: Panel midpoint positions.
        x_neut, z_neut: Point arrays along the body's neutral axis.
        x_le, z_le: Position of the body's leading edge point.
    """
    def __init__(self, N):
        self.x = np.empty(N+1)
        self.z = np.empty(N+1)
        self.x_col = np.empty(N)
        self.z_col = np.empty(N)
        self.x_mid = np.zeros((3,N))
        self.z_mid = np.zeros((3,N))
        self.x_neut = np.empty(N)
        self.z_neut = np.empty(N)
        self.x_le = 0.
        self.z_le = 0.