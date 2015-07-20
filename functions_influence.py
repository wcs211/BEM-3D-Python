# -*- coding: utf-8 -*-
"""Module for flow field functions that include all swimmers."""
import numpy as np
from functions_general import panel_vectors, transformation

def inf_sourcepanel(xp1, xp2, zp):
    """Returns a matrix of source panel influence coefficients."""
    return((xp1 * np.log(xp1**2 + zp**2) - xp2 * np.log(xp2**2 + zp**2) \
          + 2*zp*(np.arctan2(zp,xp2) - np.arctan2(zp,xp1)))/(4*np.pi))

def inf_doubletpanel(xp1, xp2, zp):
    """Returns a matrix of doublet panel influence coefficients."""
    return(-(np.arctan2(zp,xp2) - np.arctan2(zp,xp1))/(2*np.pi))

def quilt(Swimmers, influence_type, NT, NI, i):
    """Constructs a full transformation matrix that includes all Swimmers.

    The target is always the Bodies' collocation points, but the influence
    could be the Bodies, Edges, or Wakes.

    Args:
        Swimmers: List of Swimmer objects being simulated.
        influence_type: Type of influencing panel (Body, Edge, or Wake).
        NT: Total number of target panels (across all Swimmers).
        NI: Total number of influence panels (across all Swimmers).
        i: Time step number.

    Returns:
        xp1: Transformed coordinate matrix from panels' left endpoints.
        xp2: Transformed coordinate matrix from panels' right endpoints.
        zp: Transformed coordinate matrix from panels' z-levels.
    """
    xp1 = np.empty((NT,NI))
    xp2 = np.empty((NT,NI))
    zp = np.empty((NT,NI))

    for SwimT in Swimmers: # Target Swimmer (rows)
        (r0, rn) = (SwimT.i_b, SwimT.i_b+SwimT.Body.N) # Insertion row range

        for SwimI in Swimmers: # Influencing Swimmer (columns)
            if influence_type == 'Body':
                (c0, cn) = (SwimI.i_b, SwimI.i_b+SwimI.Body.N) # Insertion column range
                (xi, zi) = (SwimI.Body.AF.x, SwimI.Body.AF.z) # Coordinates of influences
            elif influence_type == 'Edge':
                (c0, cn) = (SwimI.i_e, SwimI.i_e+SwimI.Edge.N)
                (xi, zi) = (SwimI.Edge.x, SwimI.Edge.z)
            elif influence_type == 'Wake':
                (c0, cn) = (SwimI.i_w, SwimI.i_w+i)
                (xi, zi) = (SwimI.Wake.x[:i+1], SwimI.Wake.z[:i+1])
            else:
                print 'ERROR! Invalid influence type.'

            (xp1[r0:rn, c0:cn], xp2[r0:rn, c0:cn], zp[r0:rn, c0:cn]) = \
                transformation(SwimT.Body.AF.x_col, SwimT.Body.AF.z_col, xi, zi)

    return(xp1, xp2, zp)


def influence_matrices(Swimmers, i):
    """Constructs the influence coefficient matrices.

    Args:
        Swimmers: List of Swimmer objects being simulated.
        i: Time step number.

    Returns:
        sigma_all: Array containing all Swimmers' body source strengths.
        a_bodydoublet: Body doublet panels' influence matrix.
        b_bodysource: Body source panels' influence matrix.
        b_edgedoublet: Edge panels' influence matrix.
        b_wakedoublet: Wake panels' influence matrix.
        a_explicit: Augments the a_bodydoublet matrix when doing explicit Kutta.
    """
    n_b = 0
    n_e = 0
    n_w = 0
    for Swim in Swimmers:
        Swim.i_b = n_b
        Swim.i_e = n_e
        Swim.i_w = n_w
        n_b += Swim.Body.N
        n_e += 1
        n_w += i

    sigma_all = np.empty(n_b)
    mu_w_all = np.empty(n_w)
    for Swim in Swimmers:
        (r0, rn) = (Swim.i_b, Swim.i_b+Swim.Body.N)
        sigma_all[r0:rn] = Swim.Body.sigma[:]
        (r0, rn) = (Swim.i_w, Swim.i_w+i)
        mu_w_all[r0:rn] = Swim.Wake.mu[:i]

    (xp1, xp2, zp) = quilt(Swimmers, 'Body', n_b, n_b, i)
    # Body source singularities influencing the bodies (part of RHS)
    b_bodysource = inf_sourcepanel(xp1, xp2, zp)
    # Body doublet singularities influencing bodies themselves (the A matrix)
    a_bodydoublet = inf_doubletpanel(xp1, xp2, zp)

    (xp1, xp2, zp) = quilt(Swimmers, 'Edge', n_b, n_e, i)
    # Edge doublet singularities influencing the bodies (part of RHS)
    b_edgedoublet = inf_doubletpanel(xp1, xp2, zp)

    if i==0: # There are no wake panels until i==1
        b_wakedoublet = 0
    else:
        (xp1, xp2, zp) = quilt(Swimmers, 'Wake', n_b, n_w, i)
        # Wake doublet singularities influencing the bodies (part of RHS)
        b_wakedoublet = inf_doubletpanel(xp1, xp2, zp)

    a_explicit = np.zeros((n_b,n_b))

    return(sigma_all, mu_w_all, a_bodydoublet, a_explicit,
           b_bodysource, b_edgedoublet, b_wakedoublet)

def solve_phi(Swimmers, RHO, DEL_T, i, outerCorr):
    """Solves the boundary integral equation using a Kutta condition.

    Args:
        Swimmers: List of Swimmer objects being simulated.
        RHO: Fluid density.
        DEL_T: Time step length.
        i: Time step number.
    """
    for Swim in Swimmers:  
        if (outerCorr <= 1):
            # mu_past used in differencing for pressure
            Swim.Body.mu_past[1,:] = Swim.Body.mu_past[0,:]
            Swim.Body.mu_past[0,:] = Swim.Body.mu
    
    (sigma_all, mu_w_all, a_b, a_e, b_b, b_e, b_w) = influence_matrices(Swimmers, i)

    n_iter = 0
    while True:
        n_iter += 1

        if n_iter == 1:
            # Begin with explicit Kutta condition as first guess
            # Construct the augmented body matrix by combining body and trailing edge panel arrays
            for SwimI in Swimmers:
                a_e[:,SwimI.i_b] = -b_e[:, SwimI.i_e]
                a_e[:,SwimI.i_b+SwimI.Body.N-1] = b_e[:, SwimI.i_e]
            a_inv = np.linalg.inv(a_b + a_e)
            # Get right-hand side
            if i == 0:
                b = -np.dot(b_b, sigma_all)
            else:
                b = -np.dot(b_b, sigma_all) - np.dot(b_w, mu_w_all)
            # Solve for bodies' doublet strengths using explicit Kutta
            mu_b_all = np.dot(a_inv, b)
            # First mu_guess (from explicit Kutta)
            for Swim in Swimmers:
                Swim.mu_guess = np.empty(2) # [0] is current guess, [1] is previous
                Swim.delta_cp = np.empty(2) # [0] is current delta_cp, [1] is previous
                Swim.Body.mu[:] = mu_b_all[Swim.i_b:Swim.i_b+Swim.Body.N]
                Swim.mu_guess[0] = Swim.Body.mu[-1]-Swim.Body.mu[0]

        else:
            if n_iter == 2: # Make a second initial guess
                # Update phi_dinv so it no longer includes explicit Kutta condition
                a_inv = np.linalg.inv(a_b)

                Swimmers[0].mu_guess[1] = Swimmers[0].mu_guess[0]
                Swimmers[0].delta_cp[1] = Swimmers[0].delta_cp[0]
                Swimmers[0].mu_guess[0] = 0.8*Swimmers[0].mu_guess[1] # Multiply first (explicit) guess by arbitrary constant to get second guess

            else: # Newton method to get delta_cp == 0
                # Get slope, which is change in delta_cp divided by change in mu_guess
                slope = (Swimmers[0].delta_cp[0]-Swimmers[0].delta_cp[1])/(Swimmers[0].mu_guess[0]-Swimmers[0].mu_guess[1])
                Swimmers[0].mu_guess[1] = Swimmers[0].mu_guess[0]
                Swimmers[0].delta_cp[1] = Swimmers[0].delta_cp[0]
                Swimmers[0].mu_guess[0] = Swimmers[0].mu_guess[1] - Swimmers[0].delta_cp[0]/slope

            # Form right-hand side including mu_guess as an influence
            if i == 0:
                rhs = -np.dot(b_b, sigma_all) - np.squeeze(np.dot(b_e, Swimmers[0].mu_guess[0]))
            else:
                rhs = -np.dot(b_b, sigma_all) - np.dot(np.insert(b_w, 0, b_e[:,0], axis=1), np.insert(Swimmers[0].Wake.mu[:i], 0, Swimmers[0].mu_guess[0]))

            Swimmers[0].Body.mu = np.dot(a_inv, rhs)

        for Swim in Swimmers:
            Swim.Body.pressure(RHO, DEL_T, i)
        if len(Swimmers) > 1 or Swimmers[0].SW_KUTTA == 0:
            break
        Swimmers[0].delta_cp[0] = np.absolute(Swimmers[0].Body.cp[-1]-Swimmers[0].Body.cp[0])

        # wcs211: Added a max iteration break for the implicit Kutta loop
        if Swimmers[0].delta_cp[0] < 0.0001 or n_iter >= 1000:
            break

    for Swim in Swimmers:  
        if (outerCorr <= 1):
#            # mu_past used in differencing for pressure
#            Swim.Body.mu_past[1,:] = Swim.Body.mu_past[0,:]
#            Swim.Body.mu_past[0,:] = Swim.Body.mu

            Swim.Edge.mu = Swim.mu_guess[0]
            Swim.Edge.gamma[0] = -Swim.Edge.mu
            Swim.Edge.gamma[1] = Swim.Edge.mu

            # Get gamma of body panels for use in wake rollup
            Swim.Body.gamma[0] = -Swim.Body.mu[0]
            Swim.Body.gamma[1:-1] = Swim.Body.mu[:-1]-Swim.Body.mu[1:]
            Swim.Body.gamma[-1] = Swim.Body.mu[-1]

def wake_rollup(Swimmers, DEL_T, i):
    """Performs wake rollup on the swimmers' wake panels.

    Args:
        Swimmers: List of Swimmer objects being simulated.
        DEL_T: Time step length.
        i: Time step number.
    """
    # Wake panels initialize when i==1
    if i == 0:
        pass

    else:
        NT = i # Number of targets (wake panel points that are rolling up)
        for SwimT in Swimmers:
            SwimT.Wake.vx = np.zeros(NT)
            SwimT.Wake.vz = np.zeros(NT)
            DELTA_CORE = SwimT.DELTA_CORE
            for SwimI in Swimmers:
                # Coordinate transformation for body panels influencing wake
                (xp1, xp2, zp) = transformation(SwimT.Wake.x[1:i+1], SwimT.Wake.z[1:i+1], SwimI.Body.AF.x, SwimI.Body.AF.z)

                # Angle of normal vector with respect to global z-axis
                (nx, nz) = panel_vectors(SwimI.Body.AF.x, SwimI.Body.AF.z)[2:4]
                beta = np.arctan2(-nx, nz)

                # Katz-Plotkin eqns 10.20 and 10.21 for body source influence
                dummy1 = np.log((xp1**2+zp**2)/(xp2**2+zp**2))/(4*np.pi)
                dummy2 = (np.arctan2(zp,xp2)-np.arctan2(zp,xp1))/(2*np.pi)

                # Rotate back to global coordinates
                dummy3 = dummy1*np.cos(beta) - dummy2*np.sin(beta)
                dummy4 = dummy1*np.sin(beta) + dummy2*np.cos(beta)

                # Finish eqns 10.20 and 10.21 for induced velocity by multiplying with sigma
                SwimT.Wake.vx += np.dot(dummy3, SwimI.Body.sigma)
                SwimT.Wake.vz += np.dot(dummy4, SwimI.Body.sigma)

                # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
                NI = SwimI.Body.N+1
                xp = np.repeat(SwimT.Wake.x[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Body.AF.x[:,np.newaxis].T, NT, 0)
                zp = np.repeat(SwimT.Wake.z[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Body.AF.z[:,np.newaxis].T, NT, 0)

                # Find distance r_b between each influence/target
                r_b = np.sqrt(xp**2+zp**2)

                # Katz-Plotkin eqns 10.9 and 10.10 for body doublet (represented as point vortices) influence
                dummy1 = zp/(2*np.pi*(r_b**2+DELTA_CORE**2))
                dummy2 = -xp/(2*np.pi*(r_b**2+DELTA_CORE**2))

                # Finish eqns 10.9 and 10.10 by multiplying with Body.gamma, add to induced velocity
                SwimT.Wake.vx += np.dot(dummy1, SwimI.Body.gamma)
                SwimT.Wake.vz += np.dot(dummy2, SwimI.Body.gamma)

                # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
                NI = SwimI.Edge.N+1
                xp = np.repeat(SwimT.Wake.x[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Edge.x[:,np.newaxis].T, NT, 0)
                zp = np.repeat(SwimT.Wake.z[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Edge.z[:,np.newaxis].T, NT, 0)

                # Find distance r_e between each influence/target
                r_e = np.sqrt(xp**2+zp**2)

                # Katz-Plotkin eqns 10.9 and 10.10 for edge (as point vortices) influence
                dummy1 = zp/(2*np.pi*(r_e**2+DELTA_CORE**2))
                dummy2 = -xp/(2*np.pi*(r_e**2+DELTA_CORE**2))

                # Finish eqns 10.9 and 10.10 by multiplying with Edge.gamma, add to induced velocity
                SwimT.Wake.vx += np.dot(dummy1, SwimI.Edge.gamma)
                SwimT.Wake.vz += np.dot(dummy2, SwimI.Edge.gamma)

                # Formation of (x-x0) and (z-z0) matrices, similar to xp1/xp2/zp but coordinate transformation is not necessary
                NI = i+1
                xp = np.repeat(SwimT.Wake.x[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Wake.x[:i+1,np.newaxis].T, NT, 0)
                zp = np.repeat(SwimT.Wake.z[1:i+1,np.newaxis], NI, 1) - np.repeat(SwimI.Wake.z[:i+1,np.newaxis].T, NT, 0)

                # Find distance r_w between each influence/target
                r_w = np.sqrt(xp**2+zp**2)

                # Katz-Plotkin eqns 10.9 and 10.10 for wake (as point vortices) influence
                dummy1 = zp/(2*np.pi*(r_w**2+DELTA_CORE**2))
                dummy2 = -xp/(2*np.pi*(r_w**2+DELTA_CORE**2))

                # Finish eqns 10.9 and 10.10 by multiplying with Wake.gamma, add to induced velocity
                SwimT.Wake.vx += np.dot(dummy1, SwimI.Wake.gamma[:i+1])
                SwimT.Wake.vz += np.dot(dummy2, SwimI.Wake.gamma[:i+1])

        for Swim in Swimmers:
            # Modify wake with the total induced velocity
            Swim.Wake.x[1:i+1] += Swim.Wake.vx*DEL_T
            Swim.Wake.z[1:i+1] += Swim.Wake.vz*DEL_T