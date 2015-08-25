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
        self.CE = CE
        self.x = np.zeros(self.N+1)
        self.z = np.zeros(self.N+1)
        self.mu = np.zeros(self.N)
        self.gamma = np.zeros(self.N+1)

class Wake(object):
    """A chain of wake doublet panels.

    Attributes:
        N: Number of wake panels.
        x, z: X- and Z-coordinates of the wake panel endpoints.
        mu: Doublet strengths of the wake panels.
        gamma: Circulations at the wake panel endpoints.
    """
    def __init__(self, N):
        """Inits Wake with all necessary parameters."""
        self.N = N
        self.x = np.zeros(N+1)
        self.z = np.zeros(N+1)
        self.mu = np.zeros(N)
        self.gamma = np.zeros(N+1)

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
        
        self.Cf = 0.
        self.Cl = 0.
        self.Ct = 0.
        self.Cpow = 0.

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
        N_CHORD = GeoVDVParameters.N_CHORD
        N_SPAN = GeoVDVParameters.N_SPAN
        S = GeoVDVParameters.S
        C = GeoVDVParameters.C
        SPAN = GeoVDVParameters.SPAN
        K = GeoVDVParameters.K
        EPSILON = GeoVDVParameters.EPSILON

        # Defining Chord length function
        C = C * np.ones(N_SPAN)
        
        # Initialize bottom and top x and z coordinates and parameters
        A = np.zeros((N_CHORD, N_SPAN))   
        R1 = np.zeros((N_CHORD, N_SPAN))
        R2 = np.zeros((N_CHORD, N_SPAN))
        THETA1 = np.zeros((N_CHORD, N_SPAN))
        THETA2 = np.zeros((N_CHORD, N_SPAN))
        x_body = np.zeros((N_CHORD,N_SPAN))
        x = np.zeros((2*N_CHORD-1, N_SPAN))
        y = np.zeros((2*N_CHORD-1, N_SPAN))
        z = np.zeros((2*N_CHORD-1, N_SPAN))
        z_top = np.zeros((N_CHORD, N_SPAN))
        z_bot = np.zeros((N_CHORD, N_SPAN))
        x_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        y_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        z_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        
        # Stepping through each spanwise position to calculate the positions of
        # the fluke neutral plane at the given time step.
        for i in range(N_SPAN):   
                   
            A[:,i] = C[i]*((1+EPSILON)**(K-1))*(2**(-K))
            THETA = np.linspace(0,np.pi,N_SPAN)
            R1[:,i] = np.sqrt((A[:,i]*np.cos(THETA)-A[:,i])**2+(A[:,i]**2)*np.sin(THETA)**2)
            R2[:,i] = np.sqrt((A[:,i]*np.cos(THETA)-EPSILON*A[:,i])**2+(A[:,i]**2)*np.sin(THETA)**2)
            
            THETA1[:,i] = np.arctan2((A[:,i]*np.sin(THETA)) , (A[:,i]*np.cos(THETA)-A[:,i]))
            THETA2[:,i] = np.arctan2(A[:,i]*np.sin(THETA) ,(A[:,i]*np.cos(THETA)-EPSILON*A[:,i]))
            
            x_body[:,i] = ((R1[:,i]**K)/(R2[:,i]**(K-1)))*(np.cos(K*THETA1[:,i])*np.cos((K-1)*THETA2[:,i]) + np.sin(K*THETA1[:,i])*np.sin((K-1)*THETA2[:,i]))
            z_top[:,i] = ((R1[:,i]**K)/(R2[:,i]**(K-1)))*(np.sin(K*THETA1[:,i])*np.cos((K-1)*THETA2[:,i]) - np.cos(K*THETA1[:,i])*np.sin((K-1)*THETA2[:,i]))
            z_bot[:,i] = -((R1[:,i]**K)/(R2[:,i]**(K-1)))*(np.sin(K*THETA1[:,i])*np.cos((K-1)*THETA2[:,i]) - np.cos(K*THETA1[:,i])*np.sin((K-1)*THETA2[:,i]))
            
            x_body[:,i] = x_body[:,i]-x_body[-1,i] # Carrying the leading edge to the origin
            x_body[0,i] = C[i]
            
            z_top[0,i] = 0.
            z_bot[0,i] = 0.
            z_bot[-1,i] = 0.
            
            # Merge top and bottom surfaces together
            x[:,i] = np.hstack((x_body[:,i] , x_body[-2::-1,i]))
            z[:,i] = np.hstack((z_bot[:,i] , z_top[-2::-1,i]))
            
        for i in range(2*N_CHORD-1):
            y[i,:] = np.linspace(0.,SPAN,N_SPAN)
        
        for i in range(N_SPAN-2):
            x_mid[:,i] = ((x[1:,i] + x[:-1,i])/2)
            y_mid[:,i] = ((y[1:,i] + y[:-1,i])/2)
            z_mid[:,i] = ((z[1:,i] + z[:-1,i])/2)

        BodyFrameCoordinates = PC.BodyBFC(x, y, z, x_mid, y_mid, z_mid)

        return Body(N_CHORD, N_SPAN, S, BodyFrameCoordinates, MotionParameters)

    @classmethod
    #Flat plate geometry
    def flat_plate(cls, GeoFPParameters, MotionParameters):
        N_CHORD = GeoFPParameters.N_CHORD
        N_SPAN = GeoFPParameters.N_SPAN
        S = GeoFPParameters.S
        C = GeoFPParameters.C
        SPAN = GeoFPParameters.SPAN
        D = GeoFPParameters.D
        
        # Defining mid-chord line
        MC = 0.5 * C * np.ones(N_SPAN)
        # Defining Chord length function
        C = C * np.ones(N_SPAN)
        # Defining Leading and Trailing Edge Line
        TE = MC + 0.5 * C
        LE = MC - 0.5 * C
        
        # Initialize bottom and top x and z coordinates
        xb = np.zeros((N_CHORD, N_SPAN))
        xt = np.zeros((N_CHORD, N_SPAN))
        zb = np.zeros((N_CHORD, N_SPAN))
        zt = np.zeros((N_CHORD, N_SPAN))
        x  = np.zeros((2 * N_CHORD - 1, N_SPAN))
        y  = np.zeros((2 * N_CHORD - 1, N_SPAN))
        z  = np.zeros((2 * N_CHORD - 1, N_SPAN))
        x_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        y_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        z_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        oLE = np.zeros(N_SPAN)
        oTE = np.zeros(N_SPAN)
        
        for i in range(N_SPAN):
            # Create bottom x-corrdinates
            start = 0.
            stop  = np.pi
            step  = np.copy(N_CHORD)
            theta = np.linspace(start,stop,step)
            xb[:,i]    = 0.5 * (C[i] * np.cos(theta) + C[i])
            
            # Create top x-corrdinates
            start = np.pi
            stop  = 0.
            step  = np.copy(N_CHORD)
            theta = np.linspace(start,stop,step)
            xt[:,i]    = 0.5 * (C[i] * np.cos(theta) + C[i])
            
            # Create bottom and top z-coordinates
            zb[:,i] = -0.5 * D * np.ones(N_CHORD)
            zt[:,i] =  0.5 * D * np.ones(N_CHORD)
        
            # LE and TE circle origins (z is assumed to be zero)
            oLE[i] = 0.5 * D
            oTE[i] = C[i] - 0.5 * D
        
            # Calculate new theta positions for points on the rounded ends
            count1 = np.shape(xb[xb[:,i]>=oTE[i],i])[0]
            count2 = np.shape(xb[xb[:,i]<=oLE[i],i])[0]
            count3 = np.shape(xt[xt[:,i]<=oLE[i],i])[0]
            count4 = np.shape(xt[xt[:,i]>=oTE[i],i])[0]
        
            # Determine angular new positions along the rounded ends 
            thetab = np.linspace(0,np.pi,count1+count2).T
            thetat = np.linspace(np.pi,2*np.pi,count3+count4).T
        
            # Calculate transformed leading and trailing edge points
            x1 = oTE[i] + 0.5 * D * np.cos(thetab[0:count1])
            z1 = 0.     - 0.5 * D * np.sin(thetab[0:count1])
            x2 = oLE[i] + 0.5 * D * np.cos(thetab[-1:-1-count1+1:-1])
            z2 = 0.     - 0.5 * D * np.sin(thetab[-1:-1-count1+1:-1])
            x3 = oLE[i] + 0.5 * D * np.cos(thetat[0:count3])
            z3 = 0.     - 0.5 * D * np.sin(thetat[0:count3])
            x4 = oTE[i] + 0.5 * D * np.cos(thetat[-1:-1-count3+1:-1])
            z4 = 0.     - 0.5 * D * np.sin(thetat[-1:-1-count3+1:-1])
        
            # Replace x transformed points
            xb[:count1,i]           = x1
            xb[-1:-1-count2+1:-1,i] = x2
            xt[:count3,i]           = x3
            xt[-1:-1-count4+1:-1,i] = x4
        
            # Replace z transformed points
            zb[:count1,i]           = z1
            zb[-1:-1-count2+1:-1,i] = z2
            zt[:count3,i]           = z3
            zt[-1:-1-count4+1:-1,i] = z4
        
            # Make sure LE and TE points are correctly enforced (no round-off error)
            zb[0,i]  = 0.
            zt[0,i]  = 0.
            zb[-1,i] = 0.
        
            # Merge top and bottom surfaces together
            x[:,i] = np.hstack((xb[:,i] , xb[-2::-1,i]))
            z[:,i] = np.hstack((zb[:,i] , zt[-2::-1,i]))
        
            # Shift points to correct LE position
            x[:,i] += LE[i]
        
        # Create top y-corrdinates   
        for i in range(2 * N_CHORD - 1):
            y[i,:] = np.linspace(0., SPAN, N_SPAN)
         
        # Define panel mid-points   
        for i in range(N_SPAN - 2):
            x_mid[:,i] = 0.25 * (x[1:,i] + x[:-1,i] + x[1:,i+1] + x[:-1,i+1])
            y_mid[:,i] = 0.25 * (y[1:,i] + y[:-1,i] + y[1:,i+1] + y[:-1,i+1])
            z_mid[:,i] = 0.25 * (z[1:,i] + z[:-1,i] + z[1:,i+1] + z[:-1,i+1])   

        BodyFrameCoordinates = PC.BodyBFC(xb, zb, x_mid, y_mid, z_mid)

        return Body(N_CHORD, N_SPAN, S, BodyFrameCoordinates, MotionParameters)

    @classmethod
    #Tear-drop geometry
    def tear_drop(cls, GeoTDParameters, MotionParameters):
        N_CHORD = GeoTDParameters.N_CHORD
        N_SPAN = GeoTDParameters.N_SPAN
        S = GeoTDParameters.S
        C = GeoTDParameters.C
        SPAN = GeoTDParameters.SPAN
        D = GeoTDParameters.D

        # Defining Chord length function
        C = C * np.ones(N_SPAN)

        # Initialize bottom and top x and z coordinates
        
        xb = np.zeros((N_CHORD, N_SPAN))
        xt = np.zeros((N_CHORD, N_SPAN))
        x_c = np.zeros((N_CHORD, N_SPAN))
        xb1 = np.zeros((N_SPAN))
        xb2 = np.zeros((N_SPAN))
        zb = np.zeros((N_CHORD, N_SPAN))
        zt = np.zeros((N_CHORD, N_SPAN))
        x  = np.zeros((2 * N_CHORD - 1, N_SPAN))
        y  = np.zeros((2 * N_CHORD - 1, N_SPAN))
        z  = np.zeros((2 * N_CHORD - 1, N_SPAN))
        x_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        y_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        z_mid = np.zeros((2 * N_CHORD - 2, N_SPAN - 1))
        
        
        # Stepping through each spanwise position to calculate the positions of
        # the fluke neutral plane at the given time step.
        for i in range(N_SPAN):   
                   
            xb[:,i] = np.linspace(np.pi,0.,N_CHORD)
            xt[:,i] = np.linspace(0.,np.pi,N_CHORD)
            
            # Slopes and intersects for the line segments
            m = -D/2/(C[i] - D/2)
            b = D/2 + D**2/4/(C[i] - D/2)
            
                    # Tear drop shape equation.
            x_c[:,i] = 0.5 * (1 - np.cos(xb[:,i]))
            xb[:,i] = x_c[:,i] * C[i]
            xb1 = xb[xb[:,i] <= D/2]
            xb2 = xb[xb[:,i] > D/2]
            
            zb2 = -m * xb2[:,i] - b
            zb1 = -np.sqrt((D/2)**2 - (xb1[:,i] - D/2)**2)
            zb[:,i] = np.hstack((zb2, zb1))
            
                    # Tear drop shape equation.
            x_c[:,i] = 0.5 * (1 - np.cos(xt[:,i]))
            xt[:,i] = x_c[:,i] * C[i]
            xt1 = xt[xt[:,i] <= D/2]
            xt2= xt[xt[:,i] > D/2]
            
            zt1 = np.sqrt((D/2)**2 - (xt1[:,i] - D/2)**2)
            zt2 = m * xt2[:,i] + b
            zt[:,i] = np.hstack((zt1, zt2))
                    
            zb[0,i] = 0.
            zt[0,i] = 0.
            zb[-1,i] = 0.
            
            # Merge top and bottom surfaces together
            x[:,i] = np.hstack((xb[:,i] , xt[1:,i]))
            z[:,i] = np.hstack((zb[:,i] , zt[1:,i]))
            
        for i in range(2*N_CHORD-1):
            y[i,:] = np.linspace(0.,SPAN,N_SPAN)
        
        for i in range(N_SPAN-2):
            x_mid[:,i] = ((x[1:,i] + x[:-1,i])/2)
            y_mid[:,i] = ((y[1:,i] + y[:-1,i])/2)
            z_mid[:,i] = ((z[1:,i] + z[:-1,i])/2)

        BodyFrameCoordinates = PC.BodyBFC(x, y, z, x_mid, y_mid, z_mid)

        return Body(N_CHORD, N_SPAN, S, BodyFrameCoordinates, MotionParameters)

# TODO: Change neutral_axis to neutral_plane
    def neutral_axis(self, x, y, T, DSTEP=0, TSTEP=0):
        """Finds a body's neutral axis for a given time.

        The neutral axis is the axis which coincides with the chord line and
        divides the symmetric airfoil into two.

        The axis that it finds is in an absolute frame of reference.

        Args:
            x: An array of body-frame x-coordinates to use as reference points.
            y: An array of body-frame y-coordinates to use as reference points.
            DSTEP, TSTEP: Small incremental distance/time offsets
                (intended for differencing).
            T: Time of the current step.
            X0, Y0, Z0: Initial position of the leading edge (absolute frame).
            THETA_MAX: Maximum pitching angle of the body.
            F: Frequency of the body's pitching motion.
            PHI: Phase offset of the body's pitching motion.
            V0: Free-stream velocity.

        Returns:
            x_neut and z_neut: X- and Z-coordinates of the neutral axis points.
        """
        X0 = self.MP.X0
        Y0 = self.MP.Y0
        Z0 = self.MP.Z0
        V0 = self.MP.V0

        x_neut = X0 + (x+DSTEP)*np.cos(THETA_MAX*np.sin(2*np.pi*F*(T+TSTEP) + PHI)) + V0*T
        y_neut = Y0 + (y +DSTEP)
        z_neut = Z0 + (x+DSTEP)*np.sin(THETA_MAX*np.sin(2*np.pi*F*(T+TSTEP) + PHI))
        
        return(x_neut, y_neut, z_neut)

#   TODO: Update panel_positions to account for 3D objects
    def panel_positions(self, DSTEP, T):
        """Updates all the absolute-frame coordinates of the body.

        Args:
            DSTEP: Small incremental distance to pass into neutral_axis().
            T: Time of current step.
            THETA: Current pitching angle.
        """
        bfx = self.BF.x
        bfz = self.BF.z
        bfz_col = self.BF.z_col
        V0 = self.V0 # Used only for x_le

        (x_neut, z_neut) = self.neutral_axis(bfx, T, THETA, HEAVE)

        # Infinitesimal differences on the neutral axis to calculate the tangential and normal vectors
        (xdp_s, zdp_s) = self.neutral_axis(bfx, T, THETA, HEAVE, DSTEP)
        (xdm_s, zdm_s) = self.neutral_axis(bfx, T, THETA, HEAVE, -DSTEP)

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
        # Location of leading edge (currently pitching motion only)
        self.AF.x_le = V0*T
        self.AF.z_le = HEAVE
        
    def fsi_panel_positions(self, FSI, T, THETA, HEAVE):
        self.AF.x = self.AF.x + (FSI.fluidNodeDispl[:,0] - FSI.fluidNodeDisplOld[:,0])
        self.AF.z = self.AF.z + (FSI.fluidNodeDispl[:,1] - FSI.fluidNodeDisplOld[:,1])                 

        self.AF.x_mid[0,:] = (self.AF.x[:-1] + self.AF.x[1:])/2
        self.AF.z_mid[0,:] = (self.AF.z[:-1] + self.AF.z[1:])/2

        self.BF.x = (self.AF.x - self.AF.x_le) * np.cos(-1*THETA) - (self.AF.z - self.AF.z_le) * np.sin(-1*THETA)
        self.BF.z = (self.AF.z - self.AF.z_le) * np.cos(-1*THETA) + (self.AF.x - self.AF.x_le) * np.sin(-1*THETA)
        self.BF.x_col = ((self.BF.x[1:] + self.BF.x[:-1])/2)
        self.BF.z_col = ((self.BF.z[1:] + self.BF.z[:-1])/2)

        (self.AF.x_neut, self.AF.z_neut) = self.neutral_axis(self.BF.x, T, THETA, HEAVE)

        self.AF.x_col = self.AF.x_mid[0,:] - self.S*panel_vectors(self.AF.x, self.AF.z)[2]*np.absolute(self.BF.z_col)
        self.AF.z_col = self.AF.z_mid[0,:] - self.S*panel_vectors(self.AF.x, self.AF.z)[3]*np.absolute(self.BF.z_col)

    def surface_kinematics(self, DSTEP, TSTEP, THETA_MINUS, THETA_PLUS, HEAVE_MINUS, HEAVE_PLUS, DEL_T, T, i):
        """Calculates the body-frame surface velocities of body panels.

        Also finds the body panel source strengths based on these surface
        velocities.

        Args:
            DSTEP, TSTEP: Incremental distance/time passed into neutral_axis().
            DEL_T: Time step length.
            T: Time of current step.
            i: Time step number.
            THETA_MINUS: Pitching angle minus a small time difference (TSTEP)
            THETA_PLUS: Pitching angle plus a small time difference (TSTEP)
        """
        if i == 0:

            x_col = self.BF.x_col
            z_col = self.BF.z_col

            # Panel midpoint velocity calculations
            # Calculating the surface positions at tplus(tp) and tminus(tm)
            (xtpneut, ztpneut) = self.neutral_axis(x_col, T, THETA_PLUS, HEAVE_PLUS, 0)
            (xtpdp, ztpdp) = self.neutral_axis(x_col, T, THETA_PLUS, HEAVE_PLUS, DSTEP)
            (xtpdm, ztpdm) = self.neutral_axis(x_col, T, THETA_PLUS, HEAVE_PLUS, -DSTEP)
            (xtmneut, ztmneut) = self.neutral_axis(x_col, T, THETA_MINUS, HEAVE_MINUS, 0)
            (xtmdp, ztmdp) = self.neutral_axis(x_col, T, THETA_MINUS, HEAVE_MINUS, DSTEP)
            (xtmdm, ztmdm) = self.neutral_axis(x_col, T, THETA_MINUS, HEAVE_MINUS, -DSTEP)

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

        # Body source strengths with normal vector pointing outward (overall sigma pointing outward)
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
        if i == 0:
            dmu_dt = self.mu / DEL_T
        if i == 1:
            dmu_dt = (self.mu - self.mu_past[0,:])/DEL_T
        else:
            dmu_dt = (3*self.mu - 4*self.mu_past[0,:] + self.mu_past[1,:])/(2*DEL_T)

        # Unsteady pressure calculation (from Matlab code)
        qpx_tot = dmu_dl*tx + self.sigma*nx
        qpz_tot = dmu_dl*tz + self.sigma*nz

        self.p = -RHO*(qpx_tot**2 + qpz_tot**2)/2. + RHO*dmu_dt + RHO*(qpx_tot*(self.V0+self.vx) + qpz_tot*self.vz)
        self.cp = self.p / (0.5*RHO*self.V0**2)

    def force(self, THETA, RHO, V0, C, B, i):
        """Calculates drag and lift forces acting on the body.

        Args:
            i: Time step number.
        """
        
        (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x, self.AF.z)

        delFx = -self.p * lpanel * B * nx
        delFz = -self.p * lpanel * B * nz
        delF = np.array([delFx, delFz])
#        delF = -np.multiply(np.kron(self.p, np.array([1, 1])), np.array([nx, nz]).T)
        delP = np.sum(-delF * np.array([self.vx.T, self.vz.T]), 1)

        force = np.sum(delF,1)
        lift = force[1] * np.cos(THETA) - force[0] * np.sin(THETA)
        thrust = -(force[1] * np.sin(THETA) + force[0] * np.cos(THETA))
        power = np.sum(delP, 0)
        
        self.Cf = np.sqrt(force[0]**2 + force[1]**2) / (0.5 * RHO * V0**2 * C *B)
        self.Cl = lift /(0.5 * RHO * V0**2 * C *B)
        self.Ct = thrust / (0.5 * RHO * V0**2 * C *B)
        self.Cpow = power /  (0.5 * RHO * V0**3 * C *B)
        
#        Body.drag[i-1] = np.dot(self.p[i-1,:]*lpanel, np.reshape(tx,(self.N,1)))\
#                      + np.dot(self.p[i-1,:]*lpanel, np.reshape(-tz,(self.N,1)))
#
#        self.lift[i-1] = np.dot(self.p[i-1,:]*lpanel, np.reshape(-nz,(self.N,1)))\
#                      + np.dot(self.p[i-1,:]*lpanel, np.reshape(nx,(self.N,1)))
