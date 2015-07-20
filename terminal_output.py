# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import os
import time
import socket

class print_output(object):
    'Toolkit for formated terminal output'
    def prog_title(self, VERSION):
        """
        Prints the program title to STDOUT.
        
        Args:
            VERSION (str): Program version number.
        """
        print "/*---------------------------------------------------------------------------*\\"
        print "| .   )\\      Py thon       |                                                 |"
        print "| \\`.-' `-oo                | Written by:  Lehigh Biofluids Group             |"
        print "|  ) _  __,0) B  oundary    | Version:     %s                              |" % VERSION
        print "| /.' )/      E  lement     | Web:         http://www.lehigh.edu/             |"
        print "| '           M  ethod      |                                                 |"
        print "\\*---------------------------------------------------------------------------*/"
        print 'Date     : %s' % time.strftime("%B %d, %Y")
        print 'Time     : %s' % time.strftime("%I:%M:%S %p")
        print 'Host     : %s' % socket.gethostname()
        print 'PID      : %i' % os.getpid()
        print 'Case     : %s' % os.getcwd()
        print ''
        print '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
        
    def calc_input(self,ALPHA_MAX,RE,AOA_MAX,DEL_T):
        """
        Prints calculated simulation inputs to STDOUT.
        
        Args:
            ALPHA_MAX (float): Maximum pitching angle.
            RE (float): Reynolds number based on fluid speed and chord length.
            AOA_MAX (float): Maximum angle of attack.
            DEL_T (float): Time-step size.
        """
        print 'Calculated inputs:'
        print '     ALPHA_MAX = %f' % ALPHA_MAX
        print '     RE = %e'        % RE
        print '     AOA_MAX = %f'   % AOA_MAX
        print '     DEL_T = %f'     % DEL_T
        
    def initialize_output(self,t):
        """
        Informs the user that a flow solution is being generated for t = 0
        by printing to STDOUT.
        
        Args:
            t (float): Simulation time.
        """
        print '\nInitializing the flow solution for FLOW TIME = %f\n' % t
        
    def timestep_header(self,i_t,t):
        """
        Prints the time-step number and flow time to STDOUT.
        
        Args:
            i_t (int): Time-step number.
            t (float): Flow time.
        """
        print '==============================================================================='
        print ' TIME-STEP NUMBER = %i, FLOW TIME = %f' % (i_t,t)
        print '-------------------------------------------------------------------------------'
        
    def fsi_header(self):
        """
        Prints the fluid-structure interaction header to STDOUT.
        """
        print '|  Iter.  |  Relax.    |  Max Displ |  Max Res  | RMS Res Norm | Max Res Norm |'
        print '+---------+------------+------------+-----------+--------------+--------------+'
        
    def fsi_iter_out(self,outer_corr,fsi_relaxation_factor,max_displ,max_mag_fsi_residual,fsi_residual_norm,max_fsi_residual_norm):
        """
        Prints the fluid structure ineraction information to STDOUT.
        
        Args:
            outer_corr (int): subiteration number for the FSI solver.
            fsi_relaxation_factor (float): Under-relaxation factor for the FSI coupling.
            max_displ (float): Maximum nodal displacment.
            max_mag_fsi_residual (float): maximum difference betweeen the solid and fluid domain.
            fsi_residual_norm (float): Scaled L2 norm of the FSI residual.
            max_fsi_residual_norm (float): Scaled Infinity norm of the FSI residual.
        """        
        print '| %7i |   %.2E |   %.2E |  %.2E |     %.2E |     %.2E |' % (
            outer_corr,fsi_relaxation_factor,max_displ,max_mag_fsi_residual,
            fsi_residual_norm,max_fsi_residual_norm )
        
    def fsi_converged(self):
        """Prints to STDOUT that the FSI solver has converged."""
        print '+-----------------------------------------------------------------------------+'
        print '| SOLUTION CONVERGED!                                                         |'
        print '+-----------------------------------------------------------------------------+'
        
    def fsi_not_converged(self):
        """
        Prints to STDOUT that the FSI solver has reached the maximum allotted 
        subiterations without meeting the convergence tolerance.
        """
        print '+-----------------------------------------------------------------------------+'
        print '| WARNING! MAX INNER-LOOP ITERATIONS REACHED                                  |'
        print '+-----------------------------------------------------------------------------+'
        
    def solution_output(self,cf,cl,ct,cpow):
        """
        Prints time-step solution information to STDOUT.
        
        Args:
            d_visc (float): Viscous body drag.
            cf (float): Force coefficient.
            cl (float): Lift coefficient.
            ct (float): Thrust coefficient.
            cpow (float): Power coefficient.
            gamma (float): Body circulation.
        """
        print '| Solution Information:                                                       |'
        print '|     cf       = %13e                                                |' % cf    
        print '|     cl       = %13e                                                |' % cl    
        print '|     ct       = %13e                                                |' % ct    
        print '|     cpow     = %13e                                                |' % cpow  
        
    def solution_avg_output(self,cl_avg,ct_avg,tnet_avg,d_avg,pow_avg,cpow_avg):
        """
        Prints current time-averaged solution information to STDOUT.
        
        Args:           
            cl_avg (float): Time averaged lift coefficient.
            ct_avg (float): Time averaged thrust coefficient.
            tnet_avg (float):Time averaged thrust.
            d_avg (float): Time averaged drag.
            pow_avg (float): Time averaged power.
            cpow (float): Time averaged power coefficient.
        """
        print '|                                                                             |'
        print '| Solution Average Cycle Information:                                         |'
        print '|     cl_avg   = %13e                                                |' % cl_avg  
        print '|     ct_avg   = %13e                                                |' % ct_avg  
        print '|     tnet_avg = %13e                                                |' % tnet_avg
        print '|     d_avg    = %13e                                                |' % d_avg   
        print '|     pow_avg  = %13e                                                |' % pow_avg 
        print '|     cpow_avg = %13e                                                |' % cpow_avg
        
    def solution_complete_output(self, perDone):
        """
        Prints to STDOUT the completion percentage of the simulation.
        
        Args:
            perDone (float): Percent of the simulation complete
        """
        print '|     %3i%% Complete                                                           |' % perDone
        print '+-----------------------------------------------------------------------------+\n'