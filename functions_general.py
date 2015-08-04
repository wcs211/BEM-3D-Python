import os
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy import arange, array, exp
#from swimmer_class import Swimmer

    # x,z components of each panel's tangential and normal vectors
def panel_vectors(x,z):
    lpanel = np.sqrt((x[1:]-x[:-1])**2 + (z[1:]-z[:-1])**2)
    tx = (x[1:]-x[:-1])/lpanel
    tz = (z[1:]-z[:-1])/lpanel
    nx = -tz
    nz = tx
    return (tx,tz,nx,nz,lpanel)

    # x,z components of each midpoint's/collocation point's tangential and normal vectors
def point_vectors(xdp,xdm,zdp,zdm):
    tx = (xdp-xdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    tz = (zdp-zdm)/np.sqrt((xdp-xdm)**2 + (zdp-zdm)**2)
    nx = -tz
    nz = tx
    return(tx,tz,nx,nz)

def archive(array, axis=0):
    """
    Shifts array values along an axis (row-wise by default).

    Used for arrays that keep past values for differencing with respect to time.

    Args:
        array: The array that will be manipulated.
        axis: The axis to shift values along (0==row-wise, 1==column-wise).
    """
    if len(np.shape(array)) == 1:
        array[1:] = array[:-1]
    elif len(np.shape(array)) == 2:
        if axis == 0:
            array[1:,:] = array[:-1,:]
        else:
            array[:,1:] = array[:,:-1]

# Velocity and velocity potential equations are defined in panel coordinates so a transformation should be done
# Each row of xp1/xp2/zp is a target, and each column is an influence
# NI is N influences, NT is N targets
# xi/zi is x/z of influences, xt/zt is x/z of target points
def transformation(xt,zt,xi,zi):
# Returns xp1, xp2, zp
# Others: NT, NI, tx, tz, nx, nz, dummy, x_tile, z_tile, tx_tile, tz_tile

    NT = np.size(xt)
    NI = np.size(xi)-1

    (tx,tz,nx,nz) = panel_vectors(xi,zi)[:-1]

    # Intermediary variables to reduce number of tile/repeat operations
    # From normalvectors: tx==nz, tz==-nx
    x_tile = np.repeat(xt[:,np.newaxis],NI,1) - np.repeat(xi[:-1,np.newaxis].T,NT,0)
    z_tile = np.repeat(zt[:,np.newaxis],NI,1) - np.repeat(zi[:-1,np.newaxis].T,NT,0)
    tx_tile = np.repeat(tx[:,np.newaxis].T,NT,0)
    tz_tile = np.repeat(tz[:,np.newaxis].T,NT,0)

    # Transforming left side collocation points from global to local coordinates
    xp1 = x_tile*tx_tile + z_tile*tz_tile
    zp = x_tile*(-tz_tile) + z_tile*tx_tile

    # Transforming right side panel points into local coordinate system
    dummy = (xi[1:]-xi[:-1])*tx + (zi[1:]-zi[:-1])*tz
    xp2 = xp1 - np.repeat(dummy[:,np.newaxis].T,NT,0)

    return(xp1,xp2,zp)

def absoluteToBody(Body, Solid, THETA, HEAVE):
    """Transforms absolute reference frame to body reference frame"""
    Body.BF.x = ((Body.AF.x - Body.AF.x_le) * np.cos(-1*THETA) - (Body.AF.z - Body.AF.z_le) * np.sin(-1*THETA))
    Body.BF.z = ((Body.AF.z - Body.AF.z_le) * np.cos(-1*THETA) + (Body.AF.x - Body.AF.x_le) * np.sin(-1*THETA))
    Body.BF.x_col = ((Body.BF.x[1:] + Body.BF.x[:-1])/2)
    Body.BF.z_col = ((Body.BF.z[1:] + Body.BF.z[:-1])/2)

    Solid.nodesNew[:,0] = (Solid.nodes[:,0] - Body.AF.x_le) * np.cos(-1*THETA) - (Solid.nodes[:,1] - Body.AF.z_le) * np.sin(-1*THETA)
    Solid.nodesNew[:,1] = (Solid.nodes[:,1] - Body.AF.z_le) * np.cos(-1*THETA) + (Solid.nodes[:,0] - Body.AF.x_le) * np.sin(-1*THETA)

def ramp(t, slope, startTime):
    """
    This function can generate a ramp signal based on the following inputs:

    Args:
        t: array of time samples
        slope: slope of the ramp signal
        startTime: location where the ramp turns on
    """
    # Get the number of samples in the output signal
    N = t.size

    # Initialize the ramp signal
    r = np.zeros(N)

    # Find the index where the ramp turns on
    if (np.median(np.diff(t)) > 0):
        startInd = np.min((t>=startTime).nonzero())
        popInd =np.arange(startInd,N)
    elif (np.median(np.diff(t)) < 0):
        # Time-reversed ramp
        startTime = -1. * startTime
        startInd = np.max((t>=startTime).nonzero())
        popInd = np.arange(startInd)
        slope = -1. * slope

    # For indicies greater than the start time, compute the
    # proper signal value based on slope
    r[popInd] = slope * (t[popInd] + startTime) - 2 * startTime * slope

    return (r)

def geom_setup(P, PC, Swimmer, solid=None, FSI=None, PyFEA=None):
    SwiP     = [None for x in xrange(P['N_SWIMMERS'])]
    GeoP     = [None for x in xrange(P['N_SWIMMERS'])]
    MotP     = [None for x in xrange(P['N_SWIMMERS'])]
    Swimmers = [None for x in xrange(P['N_SWIMMERS'])]
    SolidP   = [None for x in xrange(P['N_SWIMMERS'])]
    FSIP     = [None for x in xrange(P['N_SWIMMERS'])]
    PyFEAP   = [None for x in xrange(P['N_SWIMMERS'])]

    for i in xrange(P['N_SWIMMERS']):
        SwiP[i] = PC.SwimmerParameters(P['CE'], P['DELTA_CORE'], P['SW_KUTTA'])
        if (P['SW_GEOMETRY'] == 'FP'):
            GeoP[i] = PC.GeoFPParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
        elif (P['SW_GEOMETRY'] == 'TD'):
            GeoP[i] = PC.GeoTDParameters(P['N_BODY'], P['S'], P['C'], P['T_MAX'])
        elif (P['SW_GEOMETRY'] == 'VDV'):
            GeoP[i] = PC.GeoVDVParameters(P['N_BODY'], P['S'], P['C'], P['K'], P['EPSILON'])
        else:
            print 'ERROR! Invalid geometry type.'

        MotP[i] = PC.MotionParameters(P['X_START'][i], P['Z_START'][i], P['V0'], P['THETA_MAX'], P['F'], P['PHI'])

        Swimmers[i] = Swimmer(SwiP[i], GeoP[i], MotP[i], P['COUNTER']-1)

        if (P['SW_FSI'] == True):
            SolidP[i] = solid(Swimmers[i].Body, P['N_ELEMENTS_S'], P['T_MAX'])
            FSIP[i] = FSI(Swimmers[i].Body, SolidP[i])
            PyFEAP[i] = PyFEA(SolidP[i], P['FRAC_DELT'], P['DEL_T'], P['E'], P['RHO_S'])

            SolidP[i].initMesh()
            if (P['SW_GEOMETRY'] == 'FP'):
                SolidP[i].initThinPlate(P['T_MAX'],P['C'],P['SW_CNST_THK_BM'],P['T_CONST'],P['FLEX_RATIO'])
            elif (P['SW_GEOMETRY'] == 'TD'):
                SolidP[i].initTearDrop(P['T_MAX'],P['C'],P['SWITCH_CNST_THK_BM'],P['T_CONST'],P['FLEX_RATIO'])
            else:
                print 'ERROR! Invalid geometry type.'

    return (SwiP, GeoP, MotP, Swimmers, SolidP, FSIP, PyFEAP)

def simulation_startup(P, DIO, PC, Swimmer, solid=None, FSI=None, PyFEA=None):
    if (os.path.exists(P['OUTPUT_DIR']) == False or os.listdir(P['OUTPUT_DIR']) == []):
        P['START_FROM'] = 'zeroTime'

    if (P['START_FROM'] == 'latestTime'):
        startTime = 0.
        for file in os.listdir(''.join((P['OUTPUT_DIR'], '/'))):
            startTime = max(float(file), startTime)

        (sP, i, FLOWTIME, SwiP, GeoP, MotP, Swimmers, SolidP, FSIP, PyFEAP) = DIO.read_data(''.join((P['OUTPUT_DIR'], '/', '%.8f' % startTime)))
        if not (sP['DEL_T'] == P['DEL_T']) and (sP['N_SWIMMERS'] == P['N_SWIMMERS']) and (sP['N_BODY'] == P['N_BODY']):
            print 'ERROR! Inconsistent input parameters with starting data file.'

        if (Swimmers[0].Wake.x.shape[0] < P['COUNTER']):
            for Swim in Swimmers:
                Swim.Wake.x.resize(P['COUNTER'])
                Swim.Wake.z.resize(P['COUNTER'])
                Swim.Wake.mu.resize(P['COUNTER']-1)
                Swim.Wake.gamma.resize(P['COUNTER'])

        START_COUNTER = i + 1
        COUNTER = P['COUNTER']

    elif (P['START_FROM'] == 'firstTime'):
        startTime = sys.float_info.max
        for file in os.listdir(''.join((P['OUTPUT_DIR'], '/'))):
            startTime = max(float(file), startTime)

        (sP, i, FLOWTIME, SwiP, GeoP, MotP, Swimmers, SolidP, FSIP, PyFEAP) = DIO.read_data(''.join((P['OUTPUT_DIR'], '/', '%.8f' % startTime)))
        if not (sP['DEL_T'] == P['DEL_T']) and (sP['N_SWIMMERS'] == P['N_SWIMMERS']) and (sP['N_BODY'] == P['N_BODY']):
            print 'ERROR! Inconsistent input parameters with starting data file.'

        if (Swimmers[0].Wake.x.shape[0] < P['COUNTER']):
            for Swim in Swimmers:
                Swim.Wake.x.resize(P['COUNTER'])
                Swim.Wake.z.resize(P['COUNTER'])
                Swim.Wake.mu.resize(P['COUNTER']-1)
                Swim.Wake.gamma.resize(P['COUNTER'])

        START_COUNTER = i + 1
        COUNTER = P['COUNTER']

    elif (P['START_FROM'] == 'zeroTime'):
        startTime = '0.00000000'
        (SwiP, GeoP, MotP, Swimmers, SolidP, FSIP, PyFEAP) = geom_setup(P, PC, Swimmer, solid, FSI, PyFEA)

        START_COUNTER = 0
        COUNTER = P['COUNTER']

    else:
        print 'ERROR! Invalid START_FROM. Valid values are:'
        print '    latestTime'
        print '    firstTime'
        print '    zeroTime'

    return (START_COUNTER, COUNTER, SwiP, GeoP, MotP, Swimmers, SolidP, FSIP, PyFEAP)

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike