# -*- coding: utf-8 -*-
"""Module for the Body, Edge, and Wake classes."""

import numpy as np
from functions_general import point_vectors, panel_vectors
import parameter_classes as PC

class Edge(object):
    """An edge doublet panel located where separation occurs on a body.

    Attributes:
        N: Number of edge panels (just one).
        CE: Constant that determines the length of the edge panel.
        x, z: X- and Z-coordinates of the edge panel endpoints.
        mu: Doublet strength of the edge panel.
        gamma: Circulation at the edge panel endpoints.
    """
    def __init__(self, CE):
        """Inits Edge with all necessary parameters."""
        self.N = 1
        self.mu = np.zeros(self.N)
        self.gamma = np.zeros(self.N+1)
        self.CE = CE
        self.x = np.zeros(self.N+1)
        self.z = np.zeros(self.N+1)

class Wake(object):
    """A chain of wake doublet panels.

    Attributes:
        N: Number of wake panels.
        x, z: X- and Z-coordinates of the wake panel endpoints.
        mu: Doublet strengths of the wake panels.
        gamma: Circulations at the wake panel endpoints.
    """
    def __init__(self, COUNTER):
        """Inits Wake with all necessary parameters."""
        self.N = COUNTER-1
        self.x = np.zeros(self.N+1)
        self.z = np.zeros(self.N+1)
        self.mu = np.zeros(self.N)
        self.gamma = np.zeros(self.N+1)

class Body(object):
    """An arrangement of source/doublet panels in the shape of a swimming body.

    Attributes:
        N: Number of body panels.
        S: Parameter for shifting collocation points into the body.
        BF: A collection of various body-frame coordinates.
        AF: A collection of various absolute-frame coordinates.
        MP: A collection of parameters describing the motion of the swimmer.
        V0: The free-stream velocity (included in MP as well).
        vx, vz: X- and Z-components of body-frame surface velocities.
        sigma: Source strengths of the body panels.
        phi_s: Matrix of body source panel influences on the body.
        phi_db: Matrix of body doublet panel influences on the body.
        phi_dw: Matrix of edge and wake panel influences on the body.
        mu: Doublet strengths of the body panels.
        gamma: Circulations at the body panel endpoints.
        p: Surface pressures of the body panels.
        cp: Surface pressure coefficients of the body panels.
        mu_past: mu arrays from previous time steps for backwards differencing.
    """
    def __init__(self, N, S, BodyFrameCoordinates, MotionParameters):
        """Inits Body with all necessary parameters."""
        self.N = N
        self.S = S

        # Body-frame panel coordinates:
        # x, z, x_col, z_col
        self.BF = BodyFrameCoordinates
        # Initialize absolute-frame panel coordinates:
        # x, z, x_col, z_col, x_mid, z_mid, x_neut, z_neut
        self.AF = PC.BodyAFC(N)
        # Prescribed motion
        self.MP = MotionParameters
        self.V0 = MotionParameters.V0

        self.vx = np.zeros(N)
        self.vz = np.zeros(N)

        self.sigma = np.zeros(N)
        self.phi_s = np.zeros((N,N))
        self.phi_db = np.zeros((N,N))
        self.phi_dw = np.zeros(N)
        self.mu = np.zeros(N)
        self.gamma = np.zeros(N+1)

        self.p = np.zeros(N)
        self.cp = np.zeros(N)
        self.mu_past = np.zeros((2,N))

#   TODO: Add 3D Geometry creation here for different profile shapes
    @classmethod
    def from_van_de_vooren(cls, GeoVDVParameters, MotionParameters):
        """Creates a Body object based on a Van de Vooren airfoil geometry.

        MotionParameters are unused here, just getting passed through for the
        creation of the Body.

        Args:
            GeoVDVParameters: A collection of parameters for constructing a
                Van de Vooren geometry. (N, S, C, K, EPSILON)
            MotionParameters: Motion parameters of the swimmer.

        Returns:
            A Body object with the Van de Vooren airfoil geometry.
        """
        N = GeoVDVParameters.N
        S = GeoVDVParameters.S
        C = GeoVDVParameters.C
        K = GeoVDVParameters.K
        EPSILON = GeoVDVParameters.EPSILON

        A = C*((1+EPSILON)**(K-1))*(2**(-K))

        THETA = np.linspace(0,np.pi,N/2+1)

        R1 = np.sqrt((A*np.cos(THETA)-A)**2+(A**2)*np.sin(THETA)**2)
        R2 = np.sqrt((A*np.cos(THETA)-EPSILON*A)**2+(A**2)*np.sin(THETA)**2)

        THETA1 = np.arctan2((A*np.sin(THETA)) , (A*np.cos(THETA)-A))
        THETA2 = np.arctan2(A*np.sin(THETA) ,(A*np.cos(THETA)-EPSILON*A))

        x = ((R1**K)/(R2**(K-1)))*(np.cos(K*THETA1)*np.cos((K-1)*THETA2) + np.sin(K*THETA1)*np.sin((K-1)*THETA2))
        z_top = ((R1**K)/(R2**(K-1)))*(np.sin(K*THETA1)*np.cos((K-1)*THETA2) - np.cos(K*THETA1)*np.sin((K-1)*THETA2))
        z_bot = -((R1**K)/(R2**(K-1)))*(np.sin(K*THETA1)*np.cos((K-1)*THETA2) - np.cos(K*THETA1)*np.sin((K-1)*THETA2))

        x = x-x[-1] # Carrying the leading edge to the origin
        x[0] = C

        z_top[0] = 0
        z_bot[0] = 0
        z_bot[-1] = 0

        # Merge top and bottom surfaces together
        x = np.hstack((x , x[-2::-1]))
        z = np.hstack((z_bot , z_top[-2::-1]))

        x_col = ((x[1:] + x[:-1])/2)
        z_col = ((z[1:] + z[:-1])/2)

        BodyFrameCoordinates = PC.BodyBFC(x, z, x_col, z_col)

        return Body(N, S, BodyFrameCoordinates, MotionParameters)

    @classmethod
    #Flat plate geometry
    def flat_plate(cls, GeoFPParameters, MotionParameters):
        N = GeoFPParameters.N
        S = GeoFPParameters.S
        C = GeoFPParameters.C
        D = GeoFPParameters.D

        # Stepping through each spanwise position to calculate the positions of the
        # fin neutral plane at the given time step.
        start = 0
        stop  = np.pi
        step  = (N+2)/2
        theta = np.linspace(start,stop,step)
        xb = (C*np.cos(theta).T + C)/(2.)
#        print xb

        start = np.pi
        stop  = 0
        step  = (N+2)/2
        theta = np.linspace(start,stop,step)
        xt = (C*np.cos(theta).T + C)/(2.)
#        print xt
#        xb = np.linspace(c,0,(N+2)/2).T
#        xt = np.linspace(0,c,(N+2)/2).T
        zb = -0.5*D*np.ones((N+2)/2)
        zt =  0.5*D*np.ones((N+2)/2)

        # Circle origins ( z is assumed to be zero)
        oF = 0.5 * D
        oB = C - 0.5 * D

        # Calculate new theta positions for points on the rounded ends
        count1 = np.shape(xb[xb>=oB])[0]
        count2 = np.shape(xb[xb<=oF])[0]
        count3 = np.shape(xt[xt<=oF])[0]
        count4 = np.shape(xt[xt>=oB])[0]

        thetab = np.linspace(0,np.pi,count1+count2).T
        thetat = np.linspace(np.pi,2*np.pi,count3+count4).T

        # Calculate transform leading and trailing edge points
        x1 = oB + 0.5 * D * np.cos(thetab[0:count1])
        z1 = 0  - 0.5 * D * np.sin(thetab[0:count1])
        x2 = oF + 0.5 * D * np.cos(thetab[-1:-1-count1+1:-1])
        z2 = 0  - 0.5 * D * np.sin(thetab[-1:-1-count1+1:-1])
        x3 = oF + 0.5 * D * np.cos(thetat[0:count3])
        z3 = 0  - 0.5 * D * np.sin(thetat[0:count3])
        x4 = oB + 0.5 * D * np.cos(thetat[-1:-1-count3+1:-1])
        z4 = 0  - 0.5 * D * np.sin(thetat[-1:-1-count3+1:-1])

        # Replace x and z transformed points
        xb[:count1] = x1
        xb[-1:-1-count2+1:-1] = x2
        xt[:count3] = x3
        xt[-1:-1-count4+1:-1] = x4

        zb[:count1] = z1
        zb[-1:-1-count2+1:-1] = z2
        zt[:count3] = z3
        zt[-1:-1-count4+1:-1] = z4

        zb[0] = 0
        zt[0] = 0
        zb[-1] = 0

        # Merge top and bottom surfaces together
        xb = np.hstack((xb , xb[-2::-1]))
        zb = np.hstack((zb , zt[-2::-1]))

        xb_col = ((xb[1:] + xb[:-1])/2)
        zb_col = ((zb[1:] + zb[:-1])/2)

        BodyFrameCoordinates = PC.BodyBFC(xb, zb, xb_col, zb_col)

        return Body(N, S, BodyFrameCoordinates, MotionParameters)

    @classmethod
    #Tear-drop geometry
    def tear_drop(cls, GeoTDParameters, MotionParameters):
        N = GeoTDParameters.N
        S = GeoTDParameters.S
        C = GeoTDParameters.C
        D = GeoTDParameters.D

        # Stepping through each spanwise position to calculate the positions of
        # the fin neutral plane at the given time step.
        xb = np.linspace(np.pi, 0., (N+2.)/2)
        xt = np.linspace(0., np.pi, (N+2.)/2)

        # Slopes and intersects for the line segments
        m = -D/2/(C - D/2)
        b = D/2 + D**2/4/(C - D/2)

        # Tear drop shape equation.
        x_c = 0.5 * (1 - np.cos(xb))
        xb = x_c * C
        xb1 = xb[xb <= D/2]
        xb2 = xb[xb > D/2]

        zb2 = -m * xb2 - b
        zb1 = -np.sqrt((D/2)**2 - (xb1 - D/2)**2)
        zb = np.hstack((zb2, zb1))

        # Tear drop shape equation.
        x_c = 0.5 * (1 - np.cos(xt))
        xt = x_c * C
        xt1 = xt[xt <= D/2]
        xt2 = xt[xt > D/2]

        zt1 = np.sqrt((D/2)**2 - (xt1 - D/2)**2)
        zt2 = m * xt2 + b
        zt = np.hstack((zt1, zt2))

        zb[0] = 0
        zt[0] = 0
        zb[-1] = 0

        # Merge top and bottom surfaces together
        x = np.hstack((xb , xt[1:]))
        z = np.hstack((zb , zt[1:]))

        x_col = ((x[1:] + x[:-1])/2)
        z_col = ((z[1:] + z[:-1])/2)

        BodyFrameCoordinates = PC.BodyBFC(x, z, x_col, z_col)

        return Body(N, S, BodyFrameCoordinates, MotionParameters)

# TODO: Change neutral_axis to neutral_plane
    def neutral_axis(self, x, T, DSTEP=0, TSTEP=0):
        """Finds a body's neutral axis for a given time.

        The neutral axis is the axis which coincides with the chord line and
        divides the symmetric airfoil into two. CURRENTLY PITCHING MOTION ONLY.

        The axis that it finds is in an absolute frame of reference.

        Args:
            x: An array of body-frame x-coordinates to use as reference points.
            DSTEP, TSTEP: Small incremental distance/time offsets
                (intended for differencing).
            T: Time of the current step.
            X0, Z0: Initial position of the leading edge (absolute frame).
            THETA_MAX: Maximum pitching angle of the body.
            F: Frequency of the body's pitching motion.
            PHI: Phase offset of the body's pitching motion.
            V0: Free-stream velocity.

        Returns:
            x_neut and z_neut: X- and Z-coordinates of the neutral axis points.
        """
        X0 = self.MP.X0
        Z0 = self.MP.Z0
        THETA_MAX = self.MP.THETA_MAX
        F = self.MP.F
        PHI = self.MP.PHI
        V0 = self.MP.V0

        x_neut = X0 + (x+DSTEP)*np.cos(THETA_MAX*np.sin(2*np.pi*F*(T+TSTEP) + PHI)) + V0*T
        z_neut = Z0 + (x+DSTEP)*np.sin(THETA_MAX*np.sin(2*np.pi*F*(T+TSTEP) + PHI))

        return(x_neut, z_neut)

#   TODO: Update panel_positions to account for 3D objects
    def panel_positions(self, DSTEP, T):
        """Updates all the absolute-frame coordinates of the body.

        Args:
            DSTEP: Small incremental distance to pass into neutral_axis().
            T: Time of current step.
            bfx, bfz: Body-frame x- and z-coordinates.
            bfz_col: Body-frame collocation z-coordinates (unshifted)
            V0: Free-stream velocity.
        """
        bfx = self.BF.x
        bfz = self.BF.z
        bfz_col = self.BF.z_col
        V0 = self.V0 # Used only for x_le

        (x_neut, z_neut) = self.neutral_axis(bfx, T)

        # Infinitesimal differences on the neutral axis to calculate the tangential and normal vectors
        (xdp_s, zdp_s) = self.neutral_axis(bfx, T, DSTEP)
        (xdm_s, zdm_s) = self.neutral_axis(bfx, T, -DSTEP)

        # Absolute-frame panel endpoint positions for time t
        afx = x_neut + point_vectors(xdp_s, xdm_s, zdp_s, zdm_s)[2]*bfz
        afz = z_neut + point_vectors(xdp_s, xdm_s, zdp_s, zdm_s)[3]*bfz

        # Absolute-frame panel midpoint positions
        x_mid = (afx[:-1]+afx[1:])/2
        z_mid = (afz[:-1]+afz[1:])/2

        # Collocation points are the points where impermeable boundary condition is forced
        # They should be shifted inside or outside of the boundary depending on the dirichlet or neumann condition
        # Shifting surface collocation points some percent of the height from the neutral axis
        # Normal vectors point outward but positive S is inward, so the shift must be subtracted from the panel midpoints
        afx_col = x_mid - self.S*panel_vectors(afx, afz)[2]*np.absolute(bfz_col)
        afz_col = z_mid - self.S*panel_vectors(afx, afz)[3]*np.absolute(bfz_col)

        self.AF.x = afx
        self.AF.z = afz
        self.AF.x_col = afx_col
        self.AF.z_col = afz_col
        self.AF.x_mid[0,:] = x_mid
        self.AF.z_mid[0,:] = z_mid
        self.AF.x_neut = x_neut
        self.AF.z_neut = z_neut
        # Location of leading edge (assuming pitching motion only)
        self.AF.x_le = V0*T
        self.AF.z_le = 0.

    def surface_kinematics(self, DSTEP, TSTEP, DEL_T, T, i):
        """Calculates the body-frame surface velocities of body panels.

        Args:
            DSTEP, TSTEP: Incremental distance/time passed into neutral_axis().
            DEL_T: Time step length.
            T: Time of current step.
            i: Time step number.
            x_col, z_col: Unshifted body-frame collocation point coordinates.
        """
        if i == 0:

            x_col = self.BF.x_col
            z_col = self.BF.z_col

            # Panel midpoint velocity calculations
            # Calculating the surface positions at tplus(tp) and tminus(tm)
            (xtpneut, ztpneut) = self.neutral_axis(x_col, T, 0, TSTEP)
            (xtpdp, ztpdp) = self.neutral_axis(x_col, T, DSTEP, TSTEP)
            (xtpdm, ztpdm) = self.neutral_axis(x_col, T, -DSTEP, TSTEP)
            (xtmneut, ztmneut) = self.neutral_axis(x_col, T, 0, -TSTEP)
            (xtmdp, ztmdp) = self.neutral_axis(x_col, T, DSTEP, -TSTEP)
            (xtmdm, ztmdm) = self.neutral_axis(x_col, T, -DSTEP, -TSTEP)

            # Displaced airfoil's panel midpoints for times tplus(tp) and tminus(tm)
            xctp = xtpneut + point_vectors(xtpdp, xtpdm, ztpdp, ztpdm)[2]*z_col
            xctm = xtmneut + point_vectors(xtmdp, xtmdm, ztmdp, ztmdm)[2]*z_col

            zctp = ztpneut + point_vectors(xtpdp, xtpdm, ztpdp, ztpdm)[3]*z_col
            zctm = ztmneut + point_vectors(xtmdp, xtmdm, ztmdp, ztmdm)[3]*z_col

            # Velocity calculations on the surface panel midpoints
            self.vx = (xctp - xctm)/(2*TSTEP)
            self.vz = (zctp - zctm)/(2*TSTEP)

        elif i == 1:
            # First-order backwards differencing of body collocation point positions
            self.vx = (self.AF.x_mid[0,:]-self.AF.x_mid[1,:])/DEL_T - self.V0
            self.vz = (self.AF.z_mid[0,:]-self.AF.z_mid[1,:])/DEL_T

        else:
            # Second-order backwards differencing of body collocation point positions
            self.vx = (3*self.AF.x_mid[0,:]-4*self.AF.x_mid[1,:]+self.AF.x_mid[2,:])/(2*DEL_T) - self.V0
            self.vz = (3*self.AF.z_mid[0,:]-4*self.AF.z_mid[1,:]+self.AF.z_mid[2,:])/(2*DEL_T)

        (nx,nz) = panel_vectors(self.AF.x,self.AF.z)[2:4]
        self.sigma = nx*(self.V0 + self.vx) + nz*self.vz

    def pressure(self, RHO, DEL_T, i):
        """Calculates the pressure distribution along the body's surface.

        Args:
            RHO: Fluid density.
            DEL_T: Time step length.
            i: Time step number.
        """
        (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x,self.AF.z)

        # Tangential panel velocity dmu/dl, first-order differencing
        dmu_dl = np.empty(self.N)
        dmu_dl[0] = (self.mu[0]-self.mu[1]) / (lpanel[0]/2 + lpanel[1]/2)
        dmu_dl[1:-1] = (self.mu[:-2]-self.mu[2:]) / (lpanel[:-2]/2 + lpanel[1:-1] + lpanel[2:]/2)
        dmu_dl[-1] = (self.mu[-2]-self.mu[-1]) / (lpanel[-2]/2 + lpanel[-1]/2)

        # Potential change dmu/dt, second-order differencing after first time step
        if i == 1:
            dmu_dt = (self.mu - self.mu_past[0,:])/DEL_T
        else:
            dmu_dt = (3*self.mu - 4*self.mu_past[0,:] + self.mu_past[1,:])/(2*DEL_T)

        # Unsteady pressure calculation (from Matlab code)
        qpx_tot = dmu_dl*tx + self.sigma*nx
        qpz_tot = dmu_dl*tz + self.sigma*nz

        self.p = -RHO*(qpx_tot**2 + qpz_tot**2)/2. + RHO*dmu_dt + RHO*(qpx_tot*(self.V0+self.vx) + qpz_tot*self.vz)
        self.cp = self.p / (0.5*RHO*self.V0**2)

#    def force(self, i):
#        """Calculates drag and lift forces acting on the body.
#
#        Args:
#            i: Time step number.
#        """
#        (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x, self.AF.z)
#
#        self.drag[i-1] = np.dot(self.p[i-1,:]*lpanel, np.reshape(tx,(self.N,1)))\
#                      + np.dot(self.p[i-1,:]*lpanel, np.reshape(-tz,(self.N,1)))
#
#        self.lift[i-1] = np.dot(self.p[i-1,:]*lpanel, np.reshape(-nz,(self.N,1)))\
#                      + np.dot(self.p[i-1,:]*lpanel, np.reshape(nx,(self.N,1)))