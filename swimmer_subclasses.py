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
    def __init__(self, CE, N_SPAN):
        """Inits Edge with all necessary parameters."""
        self.N = 1
        self.CE = CE
        self.x = np.zeros((self.N+1,N_SPAN))
        self.y = np.zeros((self.N+1,N_SPAN))
        self.z = np.zeros((self.N+1,N_SPAN))
        self.mu = np.zeros((self.N, N_SPAN))
        self.gamma = np.zeros((self.N+1,N_SPAN))

class Wake(object):
    """A chain of wake doublet panels.

    Attributes:
        N: Number of wake panels.
        x, z: X- and Z-coordinates of the wake panel endpoints.
        mu: Doublet strengths of the wake panels.
        gamma: Circulations at the wake panel endpoints.
    """
    def __init__(self, N, N_SPAN):
        """Inits Wake with all necessary parameters."""
        self.N = N
        self.x = np.zeros((N+1, N_SPAN))
        self.y = np.zeros((N+1, N_SPAN))
        self.z = np.zeros((N+1, N_SPAN))
        self.mu = np.zeros((N, N_SPAN))
        self.gamma = np.zeros((N+1, N_SPAN))

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
#    N_CHORD, N_SPAN, S, BodyFrameCoordinates, MotionParameters
    def __init__(self, N_CHORD, N_SPAN, S, BodyFrameCoordinates, MotionParameters):
        """Inits Body with all necessary parameters."""
        self.N_CHORD = N_CHORD
        self.N_SPAN = N_SPAN
        self.N = 4. * N_CHORD * N_SPAN
        self.S = S

        # Body-frame panel coordinates:
        # x, z, x_col, z_col
        self.BF = BodyFrameCoordinates
        # Initialize absolute-frame panel coordinates:
        # x, z, x_col, z_col, x_mid, z_mid, x_neut, z_neut
        self.AF = PC.BodyAFC(self.N)
        # Prescribed motion
        self.MP = MotionParameters
        self.V0 = MotionParameters.V0

        self.vx = np.zeros(self.N)
        self.vz = np.zeros(self.N)

        self.sigma = np.zeros(self.N)
        self.phi_s = np.zeros((self.N,self.N))
        self.phi_db = np.zeros((self.N,self.N))
        self.phi_dw = np.zeros(self.N)
        self.mu = np.zeros(self.N)
        self.gamma = np.zeros(self.N+1)

        self.p = np.zeros(self.N)
        self.cp = np.zeros(self.N)
        self.mu_past = np.zeros((2,self.N))
        
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

        BodyFrameCoordinates = PC.BodyBFC(x, y, z, x_mid, y_mid, z_mid)

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

    def neutral_plane(self, x, y, T, THETA, HEAVE, DSTEP1=0, DSTEP2=0, DSTEP3=0):
        """Finds a body's neutral plane for a given time.

        The neutral plane is the plane which coincides with the z=0 plane and
        divides the symmetric body into two. CURRENTLY PITCHING MOTION ONLY.

        The plane that it finds is in an absolute frame of reference.

        Args:
            x: A 2D array of fin-frame x-coordinates to use as reference points.
            Each column is a chord cross section and each row is a span cross section
            DSTEP: Small incremental distance offset (intended for differencing).
            T: Time of the current step.
            THETA: Current pitching angle.

        Returns:
            x_neut, y_neut and z_neut: 
            X-Y and Z-coordinates of the neutral plane points.
        """
        X0 = self.MP.X0
        Y0 = self.MP.Y0
        Z0 = self.MP.Z0
        V0 = self.MP.V0
        ##Neutral plane positions(x_neut y_neut z_neut) are going to be 
        ##set according to the type of the motion.
        
        x_neut = X0 + (x+DSTEP1)*np.cos(THETA) + V0*T
        y_neut = Y0 + y+DSTEP2
        z_neut = Z0 + (x+DSTEP3)*np.sin(THETA) + HEAVE

        return(x_neut, y_neut, z_neut)
        
    def panel_positions(self, DSTEP1, DSTEP2, DSTEP3, T, THETA, HEAVE):
        """Updates all the absolute-frame coordinates of the body.

        Args:
            DSTEP: Small incremental distance to pass into neutral_plane().
            T: Time of current step.
            THETA: Current pitching angle.

        """
#        Nc = int(0.5*(self.BF.x.shape[0]-1))
        Nc = self.BF.x.shape[0]
        Ns = self.BF.x.shape[1]
        bfx = self.BF.x
        bfy = self.BF.y                 #No y information coming from the geometry class yet
        bfz = self.BF.z
        bfz_col = self.BF.z_mid
        V0 = self.V0 # Used only for x_le

        (x_neut, y_neut, z_neut) = self.neutral_plane(bfx, bfy, T, THETA, HEAVE)

        # Infinitesimal differences on the neutral axis to calculate the tangential and normal vectors
        v1 = (xdp_y, ydp_y, zdp_y) = self.neutral_plane(bfx, bfy, T, THETA, HEAVE, 0, DSTEP2, DSTEP3)
        v2 = (xdm_y, ydm_y, zdm_y) = self.neutral_plane(bfx, bfy, T, THETA, HEAVE, 0, -DSTEP2, -DSTEP3)

        v3 = (xdp_x, ydp_x, zdp_x) = self.neutral_plane(bfx, bfy, T, THETA, HEAVE, DSTEP1, 0, 0)
        v4 = (xdm_x, ydm_x, zdm_x) = self.neutral_plane(bfx, bfy, T, THETA, HEAVE, -DSTEP1, 0, 0)
        
        # Absolute-frame panel endpoint positions for time t
        afx = x_neut + point_vectors(Nc, Ns, v1,v2,v3,v4)[0]*bfz
        afy = y_neut + point_vectors(Nc, Ns, v1,v2,v3,v4)[1]*bfz
        afz = z_neut + point_vectors(Nc, Ns, v1,v2,v3,v4)[2]*bfz

        # Absolute-frame panel midpoint positions
        x_mid = (afx[1:,:-1]+afx[:-1,:-1])/2
        y_mid = (afy[1:,:-1]+afy[:-1,:-1])/2
        z_mid = (afz[1:,:-1]+afz[:-1,:-1])/2

        # Collocation points are the points where impermeable boundary condition is forced
        # They should be shifted inside or outside of the boundary depending on the dirichlet or neumann condition
        # Shifting surface collocation points some percent of the height from the neutral axis
        # Normal vectors point outward but positive S is inward, so the shift must be subtracted from the panel midpoints
        afx_col = x_mid - self.S*panel_vectors(bfx, bfy, bfz, Ns, Nc)[0]*np.absolute(bfz_col)
        afy_col = y_mid - self.S*panel_vectors(bfy, bfy, bfz, Ns, Nc)[1]*np.absolute(bfz_col)
        afz_col = z_mid - self.S*panel_vectors(bfz, bfy, bfz, Ns, Nc)[2]*np.absolute(bfz_col)

        self.AF.x = afx
        self.AF.y = afy
        self.AF.z = afz
        self.AF.x_col = afx_col
        self.AF.y_col = afy_col
        self.AF.z_col = afz_col
        self.AF.x_mid = x_mid
        self.AF.y_mid = y_mid
        self.AF.z_mid = z_mid
        self.AF.x_neut = x_neut
        self.AF.y_neut = y_neut
        self.AF.z_neut = z_neut
        # Location of leading edge (currently pitching motion only)
        self.AF.x_le = V0*T
        self.AF.y_le = 0
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

    def surface_kinematics(self, DSTEP1, DSTEP2, DSTEP3, TSTEP, THETA_MINUS, THETA_PLUS, HEAVE_MINUS, HEAVE_PLUS, DEL_T, T, i):
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
            Nc = self.BF.x.shape[0]-1
            Ns = self.BF.x.shape[1]-1
            x_col = self.BF.x_mid
            y_col = self.BF.y_mid
            z_col = self.BF.z_mid

            # Panel midpoint velocity calculations
            # Calculating the surface positions at tplus(tp) and tminus(tm)
            (xtpneut, ytpneut, ztpneut) = self.neutral_plane(x_col, y_col, T, THETA_PLUS, HEAVE_PLUS, 0, 0, 0)
            
            (xtpdp_y, ytpdp_y, ztpdp_y) = self.neutral_plane(x_col, y_col, T, THETA_PLUS, HEAVE_PLUS, 0, DSTEP2, DSTEP3)
            (xtpdp_x, ytpdp_x, ztpdp_x) = self.neutral_plane(x_col, y_col, T, THETA_PLUS, HEAVE_PLUS, DSTEP1, 0, 0)
            
            (xtpdm_y, ytpdm_y, ztpdm_y) = self.neutral_plane(x_col, y_col, T, THETA_PLUS, HEAVE_PLUS, 0, -DSTEP2, -DSTEP3)
            (xtpdm_x, ytpdm_x, ztpdm_x) = self.neutral_plane(x_col, y_col, T, THETA_PLUS, HEAVE_PLUS, -DSTEP1, 0, 0)
            
            (xtmneut, ytmneut, ztmneut) = self.neutral_plane(x_col, y_col, T, THETA_MINUS, HEAVE_MINUS, 0, 0, 0)
            
            (xtmdp_y, ytmdp_y, ztmdp_y) = self.neutral_plane(x_col, y_col, T, THETA_MINUS, HEAVE_MINUS, 0, DSTEP2, DSTEP3)
            (xtmdp_x, ytmdp_x, ztmdp_x) = self.neutral_plane(x_col, y_col, T, THETA_MINUS, HEAVE_MINUS, DSTEP1, 0, 0)
            
            (xtmdm_y, ytmdm_y, ztmdm_y) = self.neutral_plane(x_col, y_col, T, THETA_MINUS, HEAVE_MINUS, 0, -DSTEP2, -DSTEP3)
            (xtmdm_x, ytmdm_x, ztmdm_x) = self.neutral_plane(x_col, y_col, T, THETA_MINUS, HEAVE_MINUS, -DSTEP1, 0, 0)

            # Displaced airfoil's panel midpoints for times tplus(tp) and tminus(tm)
            xctp = xtpneut + point_vectors(Nc, Ns, (xtpdp_y, ytpdp_y, ztpdp_y), (xtpdm_y, ytpdm_y, ztpdm_y), (xtpdp_x, ytpdp_x, ztpdp_x), (xtpdm_x, ytpdm_x, ztpdm_x))[0]*z_col
            xctm = xtmneut + point_vectors(Nc, Ns, (xtmdp_y, ytmdp_y, ztmdp_y), (xtmdm_y, ytmdm_y, ztmdm_y), (xtmdp_x, ytmdp_x, ztmdp_x), (xtmdm_x, ytmdm_x, ztmdm_x))[0]*z_col

            yctp = ytpneut + point_vectors(Nc, Ns, (xtpdp_y, ytpdp_y, ztpdp_y), (xtpdm_y, ytpdm_y, ztpdm_y), (xtpdp_x, ytpdp_x, ztpdp_x), (xtpdm_x, ytpdm_x, ztpdm_x))[1]*z_col
            yctm = ytmneut + point_vectors(Nc, Ns, (xtmdp_y, ytmdp_y, ztmdp_y), (xtmdm_y, ytmdm_y, ztmdm_y), (xtmdp_x, ytmdp_x, ztmdp_x), (xtmdm_x, ytmdm_x, ztmdm_x))[1]*z_col
            
            zctp = ztpneut + point_vectors(Nc, Ns, (xtpdp_y, ytpdp_y, ztpdp_y), (xtpdm_y, ytpdm_y, ztpdm_y), (xtpdp_x, ytpdp_x, ztpdp_x), (xtpdm_x, ytpdm_x, ztpdm_x))[2]*z_col
            zctm = ztmneut + point_vectors(Nc, Ns, (xtmdp_y, ytmdp_y, ztmdp_y), (xtmdm_y, ytmdm_y, ztmdm_y), (xtmdp_x, ytmdp_x, ztmdp_x), (xtmdm_x, ytmdm_x, ztmdm_x))[2]*z_col

            # Velocity calculations on the surface panel midpoints
            self.vx = (xctp - xctm)/(2*TSTEP)
            self.vy = (yctp - yctm)/(2*TSTEP)
            self.vz = (zctp - zctm)/(2*TSTEP)

        elif i == 1:
            # First-order backwards differencing of body collocation point positions
            self.vx = (self.AF.x_mid[0,:]-self.AF.x_mid[1,:])/DEL_T - self.V0
            self.vy = (self.AF.y_mid[0,:]-self.AF.y_mid[1,:])/DEL_T
            self.vz = (self.AF.z_mid[0,:]-self.AF.z_mid[1,:])/DEL_T

        else:
            # Second-order backwards differencing of body collocation point positions
            self.vx = (3*self.AF.x_mid[0,:]-4*self.AF.x_mid[1,:]+self.AF.x_mid[2,:])/(2*DEL_T) - self.V0
            self.vy = (3*self.AF.y_mid[0,:]-4*self.AF.y_mid[1,:]+self.AF.y_mid[2,:])/(2*DEL_T)
            self.vz = (3*self.AF.z_mid[0,:]-4*self.AF.z_mid[1,:]+self.AF.z_mid[2,:])/(2*DEL_T)

#        # Body source strengths with normal vector pointing outward (overall sigma pointing outward)
#        (nx,nz) = panel_vectors(self.AF.x,self.AF.z)[2:4]
#        self.sigma = nx*(self.V0 + self.vx) + nz*self.vz

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
        """Calculates drag, lift, and thrust forces acting on the body. 
        Calculates input power.

        Args:
            i: Time step number.
        """
        
        (tx,tz,nx,nz,lpanel) = panel_vectors(self.AF.x, self.AF.z)

        delFx = -self.p * lpanel * B * nx
        delFz = -self.p * lpanel * B * nz
        delF = np.array([delFx, delFz])
        delP = np.sum(-delF * np.array([self.vx.T, self.vz.T]), 1)

        force = np.sum(delF,1)
        lift = force[1] * np.cos(THETA) - force[0] * np.sin(THETA)
        thrust = -(force[1] * np.sin(THETA) + force[0] * np.cos(THETA))
        power = np.sum(delP, 0)
        
        self.Cf = np.sqrt(force[0]**2 + force[1]**2) / (0.5 * RHO * V0**2 * C *B)
        self.Cl = lift /(0.5 * RHO * V0**2 * C *B)
        self.Ct = thrust / (0.5 * RHO * V0**2 * C *B)
        self.Cpow = power /  (0.5 * RHO * V0**3 * C *B)