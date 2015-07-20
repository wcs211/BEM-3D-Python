# -*- coding: utf-8 -*-

import numpy as np

class SwimmerParameters(object):
    """A collection of parameters related to a single swimmer.
    
    Attributes:
        CE: Constant that determines the length of Edge panels.
        DELTA_CORE: Constant used in wake_rollup to avoid singularities near
            wake panels.
        SW_KUTTA: Switch for Kutta condition (explicit or unsteady).
    """
    def __init__(self, CE, DELTA_CORE, SW_KUTTA):
        self.CE = CE
        self.DELTA_CORE = DELTA_CORE
        self.SW_KUTTA = SW_KUTTA

class GeoVDVParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N: Number of body panels.
        S: Collocation point shifting coefficient.
        C: Body chord length.
        K, EPSILON: Van de Vooren airfoil parameters.
    """
    def __init__(self, N, S, C, K, EPSILON):
        self.SW_GEOMETRY = 'VDV'
        self.N = N
        self.S = S
        self.C = C
        self.K = K
        self.EPSILON = EPSILON
        
class GeoFPParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N: Number of body panels.
        S: Collocation point shifting coefficient.
        C: Body chord length.
        D: Diameter of leading and trailing edge circles.
    """
    def __init__(self, N, S, C, D):
        self.SW_GEOMETRY = 'FP'
        self.N = N
        self.S = S
        self.C = C
        self.D = D
        
class GeoTDParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N: Number of body panels.
        S: Collocation point shifting coefficient.
        C: Body chord length.
        D: Diameter of leading edge circle.
    """
    def __init__(self, N, S, C, D):
        self.SW_GEOMETRY = 'TD'
        self.N = N
        self.S = S
        self.C = C
        self.D = D
    
class MotionParameters(object):
    """A collection of parameters related to a swimmer's path of motion.
    
    Attributes:
        X0, Z0: Initial position of the leading edge (absolute frame).
        V0: Free-stream velocity.
        THETA_MAX: Maximum pitching angle of the body.
        F: Frequency of the body's pitching motion.
        PHI: Phase offset of the body's pitching motion.
    """
    def __init__(self, X0, Z0, V0, THETA_MAX, F, PHI):
        self.X0 = X0
        self.Z0 = Z0
        self.V0 = V0
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