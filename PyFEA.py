#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import numpy as np
import scipy.io as sio

class PyFEA(object):
    def __init__(self, Solid, FRAC_DELT, endTime, E, RHO_S):
        """
        Iniitalizes object related variables needed for other class methods.
        
        Args:
            Solid (object): A solid object created from the solid class
            FRAC_DELT (float): Fraction of the fluid solver time-step to break 
                the structural solver time-step up into.
            endTime (float): Total elapsed time for the structural solver.
            E (float): Young's Modulus of thesolid object.
            RHO_S (float): Solid object's density.
        """
        
        self.Nelements = Solid.Nelements
        self.M = np.zeros((3 * (Solid.Nnodes), 3 * (Solid.Nnodes)))
        self.K = np.zeros((3 * (Solid.Nnodes), 3 * (Solid.Nnodes)))
        self.deltaT = FRAC_DELT * endTime
        self.endTime = endTime
        self.E = E
        self.I = np.zeros((Solid.Nelements,1))
        self.A = np.zeros((Solid.Nelements,1))
        self.l = np.zeros((Solid.Nelements,1))
        self.RHO_S = RHO_S
        self.Fload = np.zeros((3*Solid.Nnodes,1))
        self.Fext_n = np.zeros((3*Solid.Nnodes,1))
        self.Fext_nPlus = np.zeros((3*Solid.Nnodes,1))
        
        # Initial Displacements
        temp = 3 * Solid.fixedCounter
        self.U_n = np.zeros((3*Solid.Nnodes,1))
        self.Udot_n = np.zeros((3*Solid.Nnodes,1))
        self.UdotDot_n = np.zeros((3*Solid.Nnodes-temp,1))
        
        # Final Displacements
        self.U_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
        self.Udot_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
        self.UdotDot_nPlus = np.zeros((3*Solid.Nnodes-temp,1))
        
        self.initU = np.zeros((3*Solid.Nnodes,1))
        self.initUdot = np.zeros((3*Solid.Nnodes,1))
        
    def elementStiffnessMatrix(self, E, I, A, l):
        """
        Calculates the element siffness matrix for bending and axial loads.

        Args:
            E (float): Young's Modulus of thesolid object.
            I (float): Element's area moment of inertia.
            A (float): Element's cross-sectional area.
            l (float): Length of the element.

        Return:
            k_e (float): NumPy 2D array of the element stiffness matrix.
        """
        C1 = (E * A / l)
        C2 = (E * I / l**3)
        k_e = np.array(
                       [[1.*C1,          0.,           0.,     -1.*C1,           0.,           0.],
                        [0.,         12.*C2,      6.*l*C2,         0.,      -12.*C2,      6.*l*C2],
                        [0.,        6.*l*C2,    4*l**2*C2,         0.,     -6.*l*C2,   2.*l**2*C2],
                        [-1.*C1,         0.,           0.,      1.*C1,           0.,           0.],
                        [0.,        -12.*C2,     -6.*l*C2,         0.,       12.*C2,     -6.*l*C2],
                        [0.,        6.*l*C2,   2.*l**2*C2,         0.,     -6.*l*C2,   4.*l**2*C2]]        
                      )
        return k_e
        
    def elementMassMatrix(self, RHO_S, A, l, mType):
        """
        Calculates the element mass matrix for bending and axial loads. This can
        return either a 'consistent' or 'lumped' mass matrix

        Args:
            RHO_S (float): Solid object's density.
            A (float): Element's cross-sectional area.
            l (float): Length of the element.
            mType (str): Type of Mass Matrix. must be 'consistent' or 'lumped'.

        Returns:
            m_e (float): NumPy 2D array of the element mass matrix.
            
        Raises:
            ValueError: If 'mType' is not defined as 'consistent' or 'lumped'.
        """
        
        if (mType == 'consistent'):
            C1 = RHO_S * A * l / 420
            C2 = RHO_S * A * l / 6
            m_e = np.array(
                           [[ 2*C2,         0,          0,       1.*C2,         0,          0],
                            [    0,    156*C1,    22*l*C1,           0,     54*C1,   -13*l*C1],
                            [    0,   22*l*C1,  4*l**2*C1,           0,   13*l*C1, -3*l**2*C1],
                            [1.*C2,         0,          0,        2*C2,         0,          0],
                            [    0,     54*C1,    13*l*C1,           0,    156*C1,   -22*l*C1],
                            [    0,  -13*l*C1, -3*l**2*C1,           0,  -22*l*C1,  4*l**2*C1]]
                          )
        elif (mType == 'lumped'):
            C1 = RHO_S * A * l / 420
            C2 = RHO_S * A * l / 6
            m_e = np.array(
                           [[2*C2,         0,          0,          C2,         0,          0],
                            [   0,        C1,          0,           0,         0,          0],
                            [   0,         0,          0,           0,         0,          0],
                            [  C2,         0,          0,        2*C2,         0,          0],
                            [   0,         0,          0,           0,        C1,          0],
                            [   0,         0,          0,           0,         0,          0]]
                          )           
        else:
            #TODO: Figure out how to throw an exception and hault exec.
            # An exception should be thrown and execuition haulted
            print 'ERROR: Invalid mass matrix type "%s"' % mType
            print 'Valid types are:'
            print '    "consistent"'
            print '    "lumped"'
            
        return m_e
        
    def elementConnectivityMatrix(self, element, theta):
        """
        Calculates the element connectivity matrix. This is used to formulate 
        the global mass and stiffness matricies.
        
        Args:
            elememnt (int): The current global element number.
            theta (float): The initial theta displacement.
            
        Returns:
            l_e (float): Element's local to global connectivity matrix.
        """
        element += 1
        l_e = np.zeros((6,3*(self.Nelements+1)))
        C = np.cos(theta)
        S = np.sin(theta)
        temp = np.array(
                        [   [ C,  S,  0,  0,  0,  0],
                            [-S,  C,  0,  0,  0,  0],
                            [ 0,  0,  1,  0,  0,  0],
                            [ 0,  0,  0,  C,  S,  0],
                            [ 0,  0,  0, -S,  C,  0],
                            [ 0,  0,  0,  0,  0,  1]    ]                       
                       )             
        l_e[:,3*element-2-1:5+3*element-2] = np.copy(temp)
        
        return l_e
        
    def HHT(self, alpha, beta, gamma, Fext_n, Fext_nPlus, fixedNodes, U_n, Udot_n, UdotDot_n):
        """
        Solves a dynamic system of equations using the Hilber-Hughes-Taylor 
        (HHT) Method. This is a transient, implicit method with numerical
        dissipation in the high frequency domain. This method has second-order
        accuracy.
        
        Args:
            alpha (float): Integration constant.
            beta (float): Integration constant.
            gamma (float): Integration constant.
            Fext_n (float): Force exterted at the begining of the time-step.
            Fext_nPlus (float): Force exterted at the end of the time-step.
            fixedNodes (int): Number of nodes with a zero dispacement condition.
            U_n (float): NumPy array of initial displacements.
            Udot_n (float): NumPy array of initial velocities.
            UdotDot_n (float): NumPy array of initial accelerations.
        
        Returns:
            U_nPlus (float): NumPy array of final displacements.
            Udot_nPlus (float): NumPy array of final velocities.
            UdotDot_nPlus (float): NumPy array of final accelerations.
        """
        
        temp = 3 * fixedNodes
        
        # Form the 'A' matrix
        A = (self.M[temp:,temp:] + beta * self.deltaT**2 * (1 - alpha) * self.K[temp:,temp:]) / (beta * self.deltaT**2)
               
        # Form the 'B' matrix
        B = (1 - alpha) * Fext_nPlus[temp:,:] + 1 / (beta * self.deltaT**2) * \
            np.dot(self.M[temp:,temp:], U_n[temp:,:] + self.deltaT * Udot_n[temp:,:] + self.deltaT**2 * \
            (0.5 - beta) * UdotDot_n) + \
            alpha * Fext_n[temp:,:] - alpha * np.dot(self.K[temp:,temp:], U_n[temp:,:])
            
        # Solve the system to get the displacements
        U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        UdotDot_nPlus = (U_nPlus - (U_n[temp:,:] + self.deltaT * Udot_n[temp:,:] + self.deltaT**2 * (0.5 - beta) * UdotDot_n)) / (beta * self.deltaT**2)
        
        # Solve for the velocities
        Udot_nPlus = (Udot_n[temp:,:] + self.deltaT * (1 - gamma) * UdotDot_n) + gamma * self.deltaT * UdotDot_nPlus
        
        return (U_nPlus, Udot_nPlus, UdotDot_nPlus)
        
    def NEWMARK(self, beta, gamma, Fext_n, Fext_nPlus, fixedNodes, U_n, Udot_n, UdotDot_n):
        """
        Solves a dynamic system of equations using the NEWMARK Method.
        This is a transient, implicit method with numerical dissipation in the 
        high frequency domain. This method has first-order accuracy.
        
        Args:
            beta (float): Integration constant.
            gamma (float): Integration constant.
            Fext_n (float): Force exterted at the begining of the time-step.
            Fext_nPlus (float): Force exterted at the end of the time-step.
            fixedNodes (int): Number of nodes with a zero dispacement condition.
            U_n (float): NumPy array of initial displacements.
            Udot_n (float): NumPy array of initial velocities.
            UdotDot_n (float): NumPy array of initial accelerations.
        
        Returns:
            U_nPlus (float): NumPy array of final displacements.
            Udot_nPlus (float): NumPy array of final velocities.
            UdotDot_nPlus (float): NumPy array of final accelerations.
        """
        
        # Form the 'A' matrix
        A = (self.M + beta * self.deltaT**2 * self.K) / (beta * self.deltaT**2)
        
        # Form the 'B' matrix
        B = Fext_nPlus + 1 / (beta * self.deltaT**2) * self.M * (U_n + \
            self.deltaT * Udot_n + \
            self.deltaT**2 * (0.5 - beta) * UdotDot_n)
            
        # Solve the system to get the displacements
        U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the accelerations
        UdotDot_nPlus = (U_nPlus - (U_n + self.deltaT * Udot_n + self.deltaT**2 * (0.5 - beta) * UdotDot_n)) / (beta * self.deltaT**2)
        
        # Solve for the velocities
        Udot_nPlus = (Udot_n + self.deltaT * (1 - gamma) * UdotDot_n) + gamma * self.deltaT * UdotDot_nPlus

        return (U_nPlus, Udot_nPlus, UdotDot_nPlus)
        
    def TRAPEZOIDAL(self, Fext_nPlus, fixedNodes, U_n, Udot_n, UdotDot_n):
        """
        Solves for the system dynamics using the trapezoidal rule.
        
        Args:
            Fext_nPlus (float): Force exterted at the end of the time-step.
            fixedNodes (int): Number of nodes with a zero dispacement condition.
            U_n (float): NumPy array of initial displacements.
            Udot_n (float): NumPy array of initial velocities.
            UdotDot_n (float): NumPy array of initial accelerations.
        
        Returns:
            U_nPlus (float): NumPy array of final displacements.
            Udot_nPlus (float): NumPy array of final velocities.
            UdotDot_nPlus (float): NumPy array of final accelerations.
        """
        temp = 3 * fixedNodes
        
        # Form the 'A' matrix
        A = (self.K[temp:,temp:] + (2 / self.deltaT)**2 * self.M[temp:,temp:])
        
        # Form the 'B' matrix

        B = (Fext_nPlus[temp:,:] + np.dot(self.M[temp:,temp:],((2 / self.deltaT)**2 * U_n[temp:,:] + \
            (4 / self.deltaT) * Udot_n[temp:,:] + UdotDot_n)))
            
        # Solve the system to get the displacements
        U_nPlus = np.linalg.solve(A, B)
        
        # Solve for the velocities
        Udot_nPlus = 2 * (U_nPlus - U_n[temp:,:]) / self.deltaT - Udot_n[temp:,:]
        
        # Solve for the accelerations
        UdotDot_nPlus = 2 * (Udot_nPlus - Udot_n[temp:,:]) / self.deltaT - UdotDot_n
        
        return (U_nPlus, Udot_nPlus, UdotDot_nPlus)
        
    def solve(self, Body, Solid, outerCorr, mType, method, alpha, beta, gamma):
        """
        Solves an unsteady finite element system of equations.
        
        Args: 
            Body (object): A body object created from the swimmer class.
            Solid (object): A solid object created from the solid class.
            outerCorr (int): Current FSI subiteration number.
            t (float): Current simulation time.
            mType (str): Type of Mass Matrix. must be 'consistent' or 'lumped'.
            method (str): Time integration method to use.
            alpha (float): Integration constant.
            beta (float): Integration constant.
            gamma (float): Integration constant.

        Raises:
            ValueError: If 'method' is not defined as 'HHT', 'NEWMARK', or 'TRAPEZOIDAL'.           
        """
        # Set the zero displacement constraints
        temp = 3 * Solid.fixedCounter
        
        U_n = np.copy(self.U_n)
        Udot_n = np.copy(self.Udot_n)
        UdotDot_n = np.copy(self.UdotDot_n)   
        
        # Reset mass and stiffness matrix to include all nodes
        self.M.fill(0.)
        self.K.fill(0.)
        
        # Assemble global mass and stiffness matricies
        for i in xrange(self.Nelements):
            # Clear/initialize old element matrix
            k_e = np.zeros((6, 6))
            m_e = np.zeros((6, 6))
            l_e = np.zeros((6, 2 * (self.Nelements + 1)))
            
            # Determine element stiffness, mass, and connectivity matricies
            k_e = self.elementStiffnessMatrix(self.E, self.I[i], self.A[i], self.l[i])
            m_e = self.elementMassMatrix(self.RHO_S, self.A[i], self.l[i], mType)
            l_e = self.elementConnectivityMatrix(i, -U_n[3*i-1,0])
            
            # Add element matricies to the global matricies
            self.M = self.M + np.dot(np.dot(np.transpose(l_e), m_e), l_e)
            self.K = self.K + np.dot(np.dot(np.transpose(l_e), k_e), l_e)
        
        # Solve for the initial acceleration matrix
        if (outerCorr == 1):
            self.Fext_n = np.copy(self.Fext_nPlus)
        Fext_n = np.copy(self.Fext_n)
#        RHS = Fext_n[temp:,:] - np.dot(self.K[temp:, temp:], U_n[temp:])
#        theright = Fext_n[temp:,:] - np.dot(self.K[temp:, temp:], U_n[temp:])
#        RHS = np.copy(Fext_n[temp:,:])
#        UdotDot_n = np.linalg.solve(self.M[temp:,temp:], RHS)
        
        # March through time until the total simulated time has elapsed
        j = np.size(np.arange(self.deltaT,self.endTime+self.deltaT,self.deltaT))
        for i in xrange(j):
            # Assume the acting force is constant through the time-step and set
            # the initial input force as the final input force.
            Fext_nPlus = np.copy(self.Fload)
            
            # Determine which integration method to use
            if (method == 'HHT'):
                (U_nPlus, Udot_nPlus, UdotDot_nPlus) = self.HHT(alpha, beta, gamma, Fext_n, Fext_nPlus, Solid.fixedCounter, U_n, Udot_n, UdotDot_n)
            elif (method == 'NEWMARK'):
                (U_nPlus, Udot_nPlus, UdotDot_nPlus) = self.NEWMARK(beta, gamma, Fext_n, Fext_nPlus, Solid.fixedCounter, U_n, Udot_n, UdotDot_n)
            elif (method == 'TRAPEZOIDAL'):
                (U_nPlus, Udot_nPlus, UdotDot_nPlus) = self.TRAPEZOIDAL(Fext_nPlus, Solid.fixedCounter, U_n, Udot_n, UdotDot_n)
            else:
                #TODO: Figure out how to throw an exception and hault exec.
                # Throw exception and hault execuition
                print 'ERROR! Invalid integration scheme "%s".' % method
                print 'Valid schemes are:'
                print '    HHT'
                print '    NEWMARK'
                print '    TRAPEZOIDAL'
                
            # Store the final displacmeents, velocities, and accelerations as 
            # the new intial values if the structural simulation has not reached its end time.
            if (i != j):
                U_n[temp:,:] = np.copy(U_nPlus)
                Udot_n[temp:,:] = np.copy(Udot_nPlus)
                UdotDot_n = np.copy(UdotDot_nPlus)
                
        # Store the final displacements, velocities, and accelerations
        
        self.U_nPlus = np.copy(U_nPlus)
        self.Udot_nPlus = np.copy(Udot_nPlus)
        self.UdotDot_nPlus = np.copy(UdotDot_nPlus)
        self.Fext_nPlus = np.copy(Fext_nPlus)