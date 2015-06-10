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
        print 'Calculated inputs:'
        print '     ALPHA_MAX = %f' % ALPHA_MAX
        print '     RE = %e'        % RE
        print '     AOA_MAX = %f'   % AOA_MAX
        print '     DEL_T = %f'     % DEL_T
        
    def initialize_output(self,t):
        print '\nInitializing the flow solution for FLOW TIME = %f\n' % t
        
    def timestep_header(self,i_t,t):
        print '==============================================================================='
        print ' TIME-STEP NUMBER = %i, FLOW TIME = %f' % (i_t,t)
        print '-------------------------------------------------------------------------------'
        
    def fsi_header(self):
        print '|  Iter.  |  Relax.    |  Max Displ |  Max Res  | RMS Res Norm | Max Res Norm |'
        print '+---------+------------+------------+-----------+--------------+--------------+'
        
    def fsi_iter_out(self,outer_corr,fsi_relaxation_factor,max_displ,max_mag_fsi_residual,fsi_residual_norm,max_fsi_residual_norm):
        print '| %7i |   %.2E |   %.2E |  %.2E |     %.2E |     %.2E |' % (
            outer_corr,fsi_relaxation_factor,max_displ,max_mag_fsi_residual,
            fsi_residual_norm,max_fsi_residual_norm )
        
    def fsi_converged(self):
        print '+-----------------------------------------------------------------------------+'
        print '| SOLUTION CONVERGED!                                                         |'
        print '+-----------------------------------------------------------------------------+'
        
    def fsi_not_converged(self):
        print '+-----------------------------------------------------------------------------+'
        print '| WARNING! MAX INNER-LOOP ITERATIONS REACHED                                  |'
        print '+-----------------------------------------------------------------------------+'
        
    def solution_output(self,d_visc,cf,cl,ct,cpow,gamma):
        print '| Solution Information:                                                       |'
        print '|     d_visc   = %13e                                                |' % d_visc
        print '|     cf       = %13e                                                |' % cf    
        print '|     cl       = %13e                                                |' % cl    
        print '|     ct       = %13e                                                |' % ct    
        print '|     cpow     = %13e                                                |' % cpow  
        print '|     gamma    = %13e                                                |' % gamma 
        
    def solution_avg_output(self,cl_avg,ct_avg,tnet_avg,d_avg,pow_avg,cpow_avg):
        print '|                                                                             |'
        print '| Solution Average Cycle Information:                                         |'
        print '|     cl_avg   = %13e                                                |' % cl_avg  
        print '|     ct_avg   = %13e                                                |' % ct_avg  
        print '|     tnet_avg = %13e                                                |' % tnet_avg
        print '|     d_avg    = %13e                                                |' % d_avg   
        print '|     pow_avg  = %13e                                                |' % pow_avg 
        print '|     cpow_avg = %13e                                                |' % cpow_avg
        
    def solution_complete_output(self, perDone):
       print '|     %3i%% Complete                                                           |' % perDone
       print '+-----------------------------------------------------------------------------+\n'