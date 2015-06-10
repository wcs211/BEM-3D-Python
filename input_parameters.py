import numpy as np

# Constants that determine other parameters and don't actually need lookup
RF = np.pi # Reduced frequency
MU = 0.001003

P = PARAMETERS = {

  'COUNTER':        401
, 'DEL_T':          np.pi*0.01/RF
#, 'DEL_T':          1/0.7104278595/200
, 'DSTEP':          10**-5
, 'TSTEP':          10**-5
, 'VERBOSITY':      20


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fluid Body Constants                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
, 'N_BODY':         100
, 'C':              1.0
, 'K':              2.-(12.4/180)
, 'EPSILON':        0.075
, 'V0':             -1.0
, 'THETA_MAX':      5.73*np.pi/180
#, 'THETA_MAX':      5.73*np.pi/180
, 'F':              RF/(2*np.pi)
#, 'F':              0.7104278595
, 'PHI':            0
, 'T_MAX':          0.002

, 'CE':             0.4
, 'S':              0.1
, 'RHO':            998.2

, 'SW_GEOMETRY':    'FP'
, 'SW_KUTTA':       0

}

# Constants dependent on declared parameters
P['DELTA_CORE'] = (0.005*P['THETA_MAX']+0.09)*P['C']
P['RE'] = P['RHO']*-P['V0']*P['C']/MU