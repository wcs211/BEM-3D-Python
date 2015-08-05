# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 22:39:54 2015

@author: biofluids4
"""

import numpy as np

# Constants that determine other parameters and don't actually need lookup

P = PARAMETERS = {
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Data I/O                                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  'SW_SAVE_DATA':       False
, 'SAVE_EVERY':         1
, 'OUTPUT_DIR':         '/home/wcs211/BEM-2D-Python/data'
, 'START_FROM':         'zeroTime'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Geometry Definition                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_GEOMETRY':    'FP'
, 'N_CHORD':        100                 # Number of chordwise panels.
, 'N_SPAN':         50                  # Number of spanwise panels.
, 'C':              1.0                 # Chord  length of rectangular body.
, 'SPAN':           0.21
, 'C_B':            2.0                 # Body chord length.
, 'K':              2.-(12.4/180)       # Van de Vooren airfoil parameter.
, 'EPSILON':        0.075               # Van de Vooren airfoil parameter.
, 'T_MAX':          0.002               # Maximum thickness.
, 'S':              0.1                 # Collocation point shifting coefficient.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

, 'V0':             -1.0                # Free-stream velocity, m/s.
, 'THETA_MAX':      5.73*np.pi/180      # Maximum pitch amplitude.
, 'ALPHA_MAX':      -0*np.pi/180        # Angle of attack (in radians).
, 'F':              2.5                 # Frequency of the body's pitching motion, Hz.
, 'PHI':            0                   # Phase of flapping time signal.
, 'CE':             0.4                 # Constant that determines the length of Edge panels (TE Factor).
, 'DELTA_CORE':     0.25                # Desingularization radius of wake panels in root chord lengths.
, 'DELTA_CORE_F':   0.02                # Fence distance from the body surface in body chord lengths.
, 'RHO':            998.2               # Fluid density, kg/m^3.
, 'NU':             1.004*10**-6        # Fluid kinematic viscosity, m^2/s.
, 'HEAVE_MAX':      0.31                # Heave amplitude, m.
, 'SCw':            0
, 'SW_KUTTA':       True

#** Virtual body properties **
# This defines the properties of a virtual body to simulate the added
# drag due to a fish body without explicitly defining the geometry and
# meshing the body with elements.  This is a first-order approximation
# of the presence of a body.

, 'AR_b':           3.5                 # Ellipsoidal body aspect ratio.
, 'S_wp':           10                  # Wetted area of the body+propulsor compared to planform area of the propulsor.
, 'Cd_bod':         0.01                # Virtual body drag coefficient. 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-step and Misc. Parameters                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

, 'N_STEP':         401                 # Number of time steps.
, 'N_CYCLE':        3                   # Number of flapping cycles.
, 'N_LUMP':         2                   # Number of cycles of the wake panels that are not lumped into a single element.
, 'DSTEP':          10**-5              # Time step of displacement. 
, 'TSTEP':          10**-5              # Time step of translation.
, 'VERBOSITY':      20

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Structural Solver Constants                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'INT_METHOD':         'HHT'
, 'M_TYPE':             'consistent'
, 'ALPHA':              0.02
, 'BETA':               0.25*(1+0.02)**2
, 'GAMMA':              0.5+0.02
, 'N_ELEMENTS_S':       100
, 'MATERIAL':           'Polyethylene'
, 'E':                  3.8e9
, 'RHO_S':              935
, 'FRAC_DELT':          1.0
, 'FLEX_RATIO':         0.1
, 'T_CONST':            0.95

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FSI Coupling Constants                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_FSI':             True
, 'N_OUTERCORR_MAX':    1500
, 'OUTER_CORR_TOL':     1e-7
, 'FIXED_PT_RELAX':     1e-8
, 'COUPLING_SCHEME':    'Aitken'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Solver Switches                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_FREE_SWIM':       False
, 'SW_VISC_DRAG':       False
, 'SW_INTERP_MTD':      True
, 'SW_CNST_THK_BM':     True
, 'SW_RAMP':            True
, 'BOD':                False

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-stepParameters                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

P['DEL_T']   = 1. / P['F'] / P['N_STEP']            # Time step.
P['COUNTER'] = P['N_CYC'] * P['N_STEP'] + 1

# Constants dependent on declared parameters
P['RE'] = -P['V0']*P['C']/P['NU']                   # Reynolds number based on the body chord length.
P['U0'] = -P['V0']*np.cos(P['ALPHA_MAX'])           # U component of free-stream velocity.
P['W0'] = P['V0']*np.sin(P['ALPHA_MAX'])            # W component of free-stream velocity.
P['RF'] = (2*np.pi*P['F']*P['C'])/abs(P['U0'])      # Reduced frequency based on root chord.
P['ST'] = P['F']*2*P['H_MAX']/abs(P['U0'])          # Strouhal number based on the fin tip peak-to-peak amplitude.  
P['V_CORE'] = P['C_B']*P['DELTA_CORE']              # Vortex core radius.
P['CUT_OFF'] = 1**-10*P['C_B']                      # Cutoff distance for influence coeff. calc: body panels.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Body Motion Parameters                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['T']           = [P['DEL_T'] * i for i in xrange(P['COUNTER'])]
RAMP             = [0.5*np.tanh(0.25*(P['T'][i]-4))+0.5 for i in xrange(P['COUNTER'])]
RAMP_P           = [0.5*np.tanh(0.25*((P['T'][i] + P['TSTEP'])-4))+0.5 for i in xrange(P['COUNTER'])]
RAMP_M           = [0.5*np.tanh(0.25*((P['T'][i] - P['TSTEP'])-4))+0.5 for i in xrange(P['COUNTER'])]

P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI'])  * RAMP[i] for i in xrange(P['COUNTER'])]
P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) * RAMP_M[i] for i in xrange(P['COUNTER'])]
P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) * RAMP_P[i] for i in xrange(P['COUNTER'])]
H_DOT            = [2 * np.pi * P['HEAVE_MAX'] * P['F'] * np.cos(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
H_DOT_PLUS       = [2 * np.pi * P['HEAVE_MAX'] * P['F'] * np.cos(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
H_DOT_MINUS      = [2 * np.pi * P['HEAVE_MAX'] * P['F'] * np.cos(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
P['THETA']       = [RAMP[i] * np.arctan(H_DOT[i] / P['V0']) for i in xrange(P['COUNTER'])]
P['THETA_MINUS'] = [RAMP_M[i] * np.arctan(H_DOT_MINUS[i] / P['V0']) for i in xrange(P['COUNTER'])]
P['THETA_PLUS']  = [RAMP_P[i] * np.arctan(H_DOT_PLUS[i] / P['V0']) for i in xrange(P['COUNTER'])]

#P['THETA']       = [np.tanh(P['T'][i])*5./180.*np.pi for i in xrange(P['COUNTER'])]
#P['THETA_MINUS'] = [np.tanh(P['T'][i])*5./180.*np.pi for i in xrange(P['COUNTER'])]
#P['THETA_PLUS']  = [np.tanh(P['T'][i])*5./180.*np.pi for i in xrange(P['COUNTER'])]
#P['HEAVE']       = [0. for i in xrange(P['COUNTER'])]
#P['HEAVE_MINUS'] = [0. for i in xrange(P['COUNTER'])]
#P['HEAVE_PLUS']  = [0. for i in xrange(P['COUNTER'])]

#P['THETA']       = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
#P['THETA_MINUS'] = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
#P['THETA_PLUS']  = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
#P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i]) * np.tanh(3 * P['T'][i])  * RAMP[i] for i in xrange(P['COUNTER'])]
#P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP'])) * np.tanh(3 * (P['T'][i] - P['TSTEP'])) * RAMP[i] for i in xrange(P['COUNTER'])]
#P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP'])) * np.tanh(3 * (P['T'][i] + P['TSTEP'])) * RAMP[i] for i in xrange(P['COUNTER'])]