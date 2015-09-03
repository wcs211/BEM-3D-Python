import os
import sys
import numpy as np

#from swimmer_class import Swimmer

# Ns is the numper of panels in spanwise direction
# Nc is the number of panels in chordwise direction
# TODO: panel_vectors now need a second tangential vector and area of panel
def panel_vectors(x,y,z):
    
    Nc = x.shape[0]
    Ns = x.shape[1]    
    
    # Calculating the normal unit vector
    # diagonal vectors associated to each panel
    diag1 = (x[1:,:-1]-x[:-1,1:],y[:-1,1:]-y[1:,:-1],z[1:,:-1]-z[:-1,1:])
    diag2 = (x[:-1,:-1]-x[1:,1:],y[1:,1:]-y[:-1,:-1],z[:-1,:-1]-z[1:,1:])
#    diag1 = (x[:-1,1:]-x[1:,:-1],y[:-1,1:]-y[1:,:-1],z[:-1,1:]-z[1:,:-1])
#    diag2 = (x[1:,1:]-x[:-1,:-1],y[1:,1:]-y[:-1,:-1],z[1:,1:]-z[:-1,:-1])
    
    diag1 = np.asarray(diag1)
    diag2 = np.asarray(diag2)
    
    vn_x = np.zeros((Nc-1,Ns-1))
    vn_y = np.zeros((Nc-1,Ns-1))
    vn_z = np.zeros((Nc-1,Ns-1))
    
    for i in xrange(Ns-1):
        for j in xrange(Nc-1):
        
            v_n = np.cross(diag1[:,j,i],diag2[:,j,i])
            #full lengths of each normal vector component
            vn_x[j,i] = v_n[0]
            vn_y[j,i] = v_n[1]
            vn_z[j,i] = v_n[2]

    nx = vn_x/np.sqrt(vn_x**2 + vn_y**2 + vn_z**2)
    ny = vn_y/np.sqrt(vn_x**2 + vn_y**2 + vn_z**2)
    nz = vn_z/np.sqrt(vn_x**2 + vn_y**2 + vn_z**2)

    # Calculating the tangential unit vector in the span direction
    xms = 0.5*(x[:-1,:] + x[1:,:])
    yms = 0.5*(y[:-1,:] + y[1:,:])
    zms = 0.5*(z[:-1,:] + z[1:,:])
    
    lps = np.sqrt((xms[:,1:] - xms[:,:-1])**2 + (yms[:,1:] - yms[:,:-1])**2 + (zms[:,1:] - zms[:,:-1])**2)
    txs = (xms[:,1:] - xms[:,:-1]) / lps
    tys = (yms[:,1:] - yms[:,:-1]) / lps
    tzs = (zms[:,1:] - zms[:,:-1]) / lps
    
     # Calculating the tangential unit vector in the chord direction
    xmc = 0.5*(x[:,:-1] + x[:,1:])
    ymc = 0.5*(y[:,:-1] + y[:,1:])
    zmc = 0.5*(z[:,:-1] + z[:,1:])
    
    lpc = np.sqrt((xmc[1:,:] - xmc[:-1,:])**2 + (ymc[1:,:] - ymc[:-1,:])**2 + (zmc[1:,:] - zmc[:-1,:])**2)
    txc = (xmc[1:,:] - xmc[:-1,:]) / lpc
    tyc = (ymc[1:,:] - ymc[:-1,:]) / lpc
    tzc = (zmc[1:,:] - zmc[:-1,:]) / lpc

    return (nx, ny, nz, txs, tys, tzs, txc, tyc, tzc, lps, lpc)

    # x,z components of each midpoint's/collocation point's tangential and normal vectors
# TODO: point_vectors now need a second tangential vector
def point_vectors(Nc, Ns, v1,v2,v3,v4):
    vy = np.asarray(v1) - np.asarray(v2)
    vx = np.asarray(v3) - np.asarray(v4)        
    
    vn_x = np.zeros((Nc,Ns))
    vn_y = np.zeros((Nc,Ns))
    vn_z = np.zeros((Nc,Ns))

    for i in xrange(Ns):
        for j in xrange(Nc):       
            v_n = np.cross(vx[:,j,i],vy[:,j,i])
            #full lengths of each normal vector component
            vn_x[j,i] = v_n[0]
            vn_y[j,i] = v_n[1]
            vn_z[j,i] = v_n[2]
     
    vn_x_unit = vn_x/np.sqrt(vn_x**2 + vn_y**2 + vn_z**2)
    vn_y_unit = vn_y/np.sqrt(vn_x**2 + vn_y**2 + vn_z**2)
    vn_z_unit = vn_z/np.sqrt(vn_x**2 + vn_y**2 + vn_z**2)    

    return (vn_x_unit,vn_y_unit,vn_z_unit)        


def archive(array, axis=3):
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
    elif len(np.shape(array)) == 3:
        if axis == 0:
            array[1:,:,:] = array[:-1,:,:]
        elif axis == 1:
            array[:,1:,:] = array[:,:-1,:]
        else:
            array[:,:,1:] = array[:,:,:-1]
        

# Velocity and velocity potential equations are defined in panel coordinates so a transformation should be done
# Each row of xp1/xp2/zp is a target, and each column is an influence
# NI is N influences, NT is N targets
# xi/zi is x/z of influences, xt/zt is x/z of target points
def transformation(xt,yt,zt,xi,yi,zi):
    """
    This code transforms the N points (Xb, Yb, Zb)' from a global coordinate 
    system into K coordinate systems defined by the unit vectors 
    vc, vt, and vn and their origins cpts.
    
    Xb = R^(1 X N)
    Yb = R^(1 X N)
    Zb = R^(1 X N)
    
    cpts = R^(3 X K)
    vc = R^(3 X K)
    vt = R^(3 X K)
    vn = R^(3 X K)
    
    The output Xp is Xp = R^(3*K X N) matrix.  The transformed coordinate
    points are assembled in the following way:
    
               Xp = [x11 x12 ... x1N]
                    [y11 y12 ... y1N]
                    [z11 z12 ... z1N]
                    [x21 x22 ... x2N]
                    [y21 y22 ... y2N]
                    [z21 z22 ... z2N]
                    [x21 x22 ... x2N]
                    [ .   .  ...  . ]
                    [ .   .  ...  . ]
                    [ .   .  ...  . ]
                    [xK1 xK2 ... xKN]
                    [yK1 yK2 ... yKN]
                    [zK1 zK2 ... zKN]
    """
    
    (nx, ny, nz, txs, tys, tzs, txc, tyc, tzc) = panel_vectors(xi, yi, zi)[:-2]

    Xb = np.reshape(xi,(1, xi.size))
    Yb = np.reshape(yi,(1, yi.size))
    Zb = np.reshape(zi,(1, zi.size))

    K = cpts.shape[1]
    N = Xb.shape[1]
     
    i = reshape(repmat(1:3*K,3,1),9*K,1);
    j = reshape(repmat(reshape((1:3*K)',3,K),3,1),9*K,1);
    v = [vc;vt;vn];
    
    v = reshape(v,9*K,1);
    
    Qt = sparse(i,j,v,3*K,3*K);
    
    Xp = Qt*(repmat([Xb;Yb;Zb],K,1) - repmat(reshape(cpts,3*K,1),1,N));

    return(xp1,xp2,zp)

def absoluteToBody(Body, Solid, THETA, HEAVE):
    """Transforms absolute reference frame to body reference frame"""
    Body.BF.x = ((Body.AF.x - Body.AF.x_le) * np.cos(-1*THETA) - (Body.AF.z - Body.AF.z_le) * np.sin(-1*THETA))
    Body.BF.z = ((Body.AF.z - Body.AF.z_le) * np.cos(-1*THETA) + (Body.AF.x - Body.AF.x_le) * np.sin(-1*THETA))
    Body.BF.x_col = ((Body.BF.x[1:] + Body.BF.x[:-1])/2)
    Body.BF.z_col = ((Body.BF.z[1:] + Body.BF.z[:-1])/2)

    Solid.nodesNew[:,0] = (Solid.nodes[:,0] - Body.AF.x_le) * np.cos(-1*THETA) - (Solid.nodes[:,1] - Body.AF.z_le) * np.sin(-1*THETA)
    Solid.nodesNew[:,1] = (Solid.nodes[:,1] - Body.AF.z_le) * np.cos(-1*THETA) + (Solid.nodes[:,0] - Body.AF.x_le) * np.sin(-1*THETA)

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
#            N_CHORD, N_SPAN, S, C, SPAN, D
            GeoP[i] = PC.GeoFPParameters(P['N_CHORD'], P['N_SPAN'], P['S'], P['C'], P['SPAN'], P['T_MAX'])
        elif (P['SW_GEOMETRY'] == 'TD'):
            GeoP[i] = PC.GeoTDParameters(P['N_CHORD'], P['N_SPAN'], P['S'], P['C'], P['SPAN'], P['T_MAX'])
        elif (P['SW_GEOMETRY'] == 'VDV'):
            GeoP[i] = PC.GeoVDVParameters(P['N_CHORD'], P['N_SPAN'], P['S'], P['C'], P['SPAN'], P['K'], P['EPSILON'])
        else:
            print 'ERROR! Invalid geometry type.'
#        X0, Y0, Z0, V0, THETA_MAX, HEAVE_MAX, F, PHI
        MotP[i] = PC.MotionParameters(P['X_START'][i], P['Y_START'][i], P['Z_START'][i], P['V0'], P['THETA_MAX'], P['HEAVE_MAX'], P['F'], P['PHI'])

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