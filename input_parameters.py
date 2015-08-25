import numpy as np

# Constants that determine other parameters and don't actually need lookup
MU = 0.001003

P = PARAMETERS = {
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Data I/O                                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  'SW_SAVE_DATA':       False
, 'SW_SV_L_CYCLE':      True
, 'SW_SV_FORCES':       True
, 'SAVE_EVERY':         1
, 'OUTPUT_DIR':         '/home/wcs211/BEM-3D-Python/data'
, 'START_FROM':         'zeroTime'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Geometry Definition                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'SW_GEOMETRY':        'FP'
, 'N_CHORD':            61                  # Number of chordwise panels.
, 'N_SPAN':             61                  # Number of spanwise panels.
, 'C':                  1.0                 # Chord  length of rectangular body.
, 'SPAN':               0.21                # Span length of the body
, 'C_B':                2.0                 # Body chord length.
, 'K':                  2.-(12.4/180)       # Van de Vooren airfoil parameter.
, 'EPSILON':            0.075               # Van de Vooren airfoil parameter.
, 'T_MAX':              0.002               # Maximum thickness.
, 'S':                  0.1                 # Collocation point shifting coefficient.
, 'AR_b':               3.5                 # Ellipsoidal body aspect ratio.
, 'S_wp':               10                  # Wetted area of the body+propulsor compared to planform area of the propulsor.
, 'Cd_bod':             0.01                # Virtual body drag coefficient. 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Time-step and Misc. Parameters                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_STEP':             100                 # Number of time-steps per cycle
, 'N_CYC':              10                  # Number of cycles to simulate
, 'N_LUMP':             2                   # Number of cycles of the wake panels that are not lumped into a single element.
, 'DSTEP':              10**-5              # Displacement incremental time-step
, 'TSTEP':              10**-5              # Translation incremental time-step
, 'VERBOSITY':          1                   # Frequency of information is printed to terminal

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'V0':             -1.0                    # Free-stream velocity, m/s.
, 'THETA_MAX':      5.73*np.pi/180          # Maximum pitch amplitude.
, 'ALPHA_MAX':      -0*np.pi/180            # Angle of attack (in radians).
, 'F':              2.5                     # Frequency of the body's pitching motion, Hz.
, 'PHI':            0                       # Phase of flapping time signal.
, 'CE':             0.4                     # Constant that determines the length of Edge panels (TE Factor).
, 'DELTA_CORE':     0.25                    # Desingularization radius of wake panels in root chord lengths.
, 'DELTA_CORE_F':   0.02                    # Fence distance from the body surface in body chord lengths.
, 'RHO':            998.2                   # Fluid density, kg/m^3.
, 'NU':             1.004*10**-6            # Fluid kinematic viscosity, m^2/s.
, 'HEAVE_MAX':      0.31                    # Heave amplitude, m.
, 'SCw':            0
, 'SW_KUTTA':       True

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
# Time-step Parameters                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['DEL_T']   = 1. / P['F'] / P['N_STEP']
P['COUNTER'] = P['N_CYC'] * P['N_STEP'] + 1

# Constants dependent on declared parameters
P['N'] = 4*P['N_CHORD']*P['N_SPAN']                 # Total number of body panels.
P['RE'] = -P['V0']*P['C']/P['NU']                   # Reynolds number based on the body chord length.
P['U0'] = -P['V0']*np.cos(P['ALPHA_MAX'])           # U component of free-stream velocity.
P['W0'] = P['V0']*np.sin(P['ALPHA_MAX'])            # W component of free-stream velocity.
P['RF'] = (2*np.pi*P['F']*P['C'])/abs(P['U0'])      # Reduced frequency based on root chord.
P['ST'] = P['F']*2*P['HEAVE_MAX']/abs(P['U0'])      # Strouhal number based on the fin tip peak-to-peak amplitude.  
P['V_CORE'] = P['C_B']*P['DELTA_CORE']              # Vortex core radius.
P['CUT_OFF'] = 1**-10*P['C_B']                      # Cutoff distance for influence coeff. calc: body panels.
P['DELTA_CORE']  = (0.005*P['THETA_MAX']+0.09)*P['C']

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Body Motion Parameters                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
P['T']           = [P['DEL_T'] * i for i in xrange(P['COUNTER'])]
P['THETA']       = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i] + P['PHI']) for i in xrange(P['COUNTER'])]
P['THETA_MINUS'] = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
P['THETA_PLUS']  = [P['THETA_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP']) + P['PHI']) for i in xrange(P['COUNTER'])]
P['HEAVE']       = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * P['T'][i]) * np.tanh(3 * P['T'][i]) for i in xrange(P['COUNTER'])]
P['HEAVE_MINUS'] = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] - P['TSTEP'])) * np.tanh(3 * (P['T'][i] - P['TSTEP'])) for i in xrange(P['COUNTER'])]
P['HEAVE_PLUS']  = [P['HEAVE_MAX'] * np.sin(2 * np.pi * P['F'] * (P['T'][i] + P['TSTEP'])) * np.tanh(3 * (P['T'][i] + P['TSTEP'])) for i in xrange(P['COUNTER'])]