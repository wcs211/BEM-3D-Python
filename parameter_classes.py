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
        N_SPAN: Number of spanwise body panels.
        N_CHORD: Number of chordwise body panels.
        S: Collocation point shifting coefficient.
        C: Body chord length.
        SPAN: Span length of body.
        K, EPSILON: Van de Vooren airfoil parameters.
    """
    def __init__(self, N_CHORD, N_SPAN, S, C, SPAN, K, EPSILON):
        self.SW_GEOMETRY = 'VDV'
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.S = S
        self.C = C
        self.SPAN = SPAN
        self.K = K
        self.EPSILON = EPSILON
        
class GeoFPParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N_SPAN: Number of spanwise body panels.
        N_CHORD: Number of chordwise body panels.
        S: Collocation point shifting coefficient.
        C: Body chord length.
        SPAN: Span length of body.
        D: Diameter of leading and trailing edge circles.
    """
    def __init__(self, N_CHORD, N_SPAN, S, C, SPAN, D):
        self.SW_GEOMETRY = 'FP'
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.S = S
        self.C = C
        self.SPAN = SPAN
        self.D = D
        
class GeoTDParameters(object):
    """A collection of parameters related to a body geometry.
    
    Attributes:
        SW_GEOMETRY: String describing the geometry type (used as a switch).
        N_SPAN: Number of spanwise body panels.
        N_CHORD: Number of chordwise body panels.
        S: Collocation point shifting coefficient.
        C: Body chord length.
        SPAN: Span length of body.
        D: Diameter of leading and trailing edge circles.
    """
    def __init__(self, N_CHORD, N_SPAN, S, C, SPAN, D):
        self.SW_GEOMETRY = 'TD'
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.S = S
        self.C = C
        self.SPAN = SPAN
        self.D = D

#==============================================================================
# # Re_organize according to Emre's code.
#==============================================================================
class MotionParameters(object):
    """A collection of parameters related to a swimmer's path of motion.
    
    Attributes:
        X0, Y0, Z0: Initial position of the leading edge (absolute frame).
        V0: Free-stream velocity.
        THETA_MAX: Maximum pitching angle of the body.
        HEAVE_MAX: Maximum heaving amplitude.
        F: Frequency of the body's pitching motion.
        PHI: Phase offset of the body's pitching motion.
    """
    def __init__(self, X0, Y0, Z0, V0, THETA_MAX, HEAVE_MAX, F, PHI):
        self.X0 = X0
        self.Y0 = Y0
        self.Z0 = Z0
        self.V0 = V0
        self.THETA_MAX = THETA_MAX
        self.HEAVE_MAX = HEAVE_MAX
        self.F = F
        self.PHI = PHI
    
class BodyBFC(object):
    """A collection of body-frame coordinates for a Body object.
    
    Attributes:
        x, y, z: Body-frame panel endpoint positions.
        x_mid, y_mid, z_mid: Body-frame panel midpoint positions (unshifted).
    """
    def __init__(self, x, y, z, x_mid, y_mid, z_mid):
        self.x = x
        self.y = y
        self.z = z
        self.x_mid = x_mid
        self.y_mid = y_mid
        self.z_mid = z_mid
        
class BodyAFC(object):
    """A collection of absolute-frame coordinates for a Body object.
    
    Attributes:
        N: Total number of body panels.
        x, y, z: Panel endpoint positions.
        x_col, y_col, z_col: Panel collocation point positions (shifted midpoints).
        x_mid, y_mid, z_mid: Panel midpoint positions.
        x_neut, y_neut, z_neut: Point arrays along the body's neutral axis.
        x_le, y_le, z_le: Position of the body's leading edge point.
    """
    def __init__(self, N_CHORD, N_SPAN):
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.N = 4. * N_CHORD * N_SPAN
        self.x = np.empty((2 * N_CHORD - 1, N_SPAN))
        self.y = np.empty((2 * N_CHORD - 1, N_SPAN))
        self.z = np.empty((2 * N_CHORD - 1, N_SPAN))
        self.x_col = np.empty((2 * N_CHORD - 2, N_SPAN - 1))
        self.y_col = np.empty((2 * N_CHORD - 2, N_SPAN - 1))
        self.z_col = np.empty((2 * N_CHORD - 2, N_SPAN - 1))
        self.x_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1, 3))
        self.y_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1, 3))
        self.z_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1, 3))
        self.x_neut = np.empty((2 * N_CHORD - 1, N_SPAN))
        self.y_neut = np.empty((2 * N_CHORD - 1, N_SPAN))
        self.z_neut = np.empty((2 * N_CHORD - 1, N_SPAN))
        self.x_le = 0.
        self.y_le = 0.
        self.z_le = 0.
