import numpy as np

# Constants that determine other parameters and don't actually need lookup
MU = 0.001003

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
, 'SW_GEOMETRY':        'FP'
, 'N_BODY':             100
, 'C':                  0.1
, 'K':                  2.-(12.4/180)
, 'EPSILON':            0.075
, 'T_MAX':              0.0011
, 'CE':                 0.4
, 'S':                  0.01

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-step and Misc. Parameters                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_STEP':             100
, 'N_CYC':              10
, 'DSTEP':              10**-5
, 'TSTEP':              10**-5
, 'VERBOSITY':          1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'V0':                 -0.05
, 'THETA_MAX':          0.
, 'HEAVE_MAX':          0.001755
, 'F':                  0.7104278595
, 'PHI':                0.5*np.pi
, 'RHO':                998.2
, 'SW_KUTTA':           True

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Skin Friction Solver Constants                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



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
}


##### The Following parameters are based on perviously declared variables #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Multi-Swimmer Parameters                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['N_SWIMMERS']  = 1
P['X_START']     = [i * 0.0 for i in xrange(P['N_SWIMMERS'])]
P['Z_START']     = [i * 0.4 for i in xrange(P['N_SWIMMERS'])]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-stepParameters                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['DEL_T']   = 1. / P['F'] / P['N_STEP']
P['COUNTER'] = P['N_CYC'] * P['N_STEP'] + 1

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

# Constants dependent on declared parameters
P['DELTA_CORE']  = (0.005*P['THETA_MAX']+0.09)*P['C']
P['RE']          = P['RHO']*-P['V0']*P['C']/MU
