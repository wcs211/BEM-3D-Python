#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import time
import numpy as np
from data_IO_class import DataIO
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from functions_influence import solve_phi, wake_rollup
from terminal_output import print_output as po
import functions_graphics as graph
from functions_general import archive, simulation_startup

po().prog_title('1.0.0')
DIO = DataIO(P)
start_time = time.time()

DEL_T = P['DEL_T']
DSTEP = P['DSTEP']
TSTEP = P['TSTEP']
T = P['T']
RHO = P['RHO']
RE = P['RE']

(START_COUNTER, COUNTER, SwiP, GeoP, MotP, Swimmers) = simulation_startup(P, DIO, PC, Swimmer)[0:6]

po().calc_input(MotP[0].THETA_MAX/np.pi*180.,RE,MotP[0].THETA_MAX/np.pi*180.,DEL_T)
po().initialize_output((START_COUNTER-1)*DEL_T)

# Data points per cycle == 1/(F*DEL_T)
for i in xrange(START_COUNTER, COUNTER):
    if i == 0:
        for Swim in Swimmers:
                Swim.Body.panel_positions(DSTEP, T[i], P['THETA'][i], P['HEAVE'][i])
                Swim.Body.surface_kinematics(DSTEP, TSTEP, P['THETA_MINUS'][i], P['THETA_PLUS'][i], P['HEAVE_MINUS'][i], P['HEAVE_PLUS'][i], DEL_T, T[i], i)
                Swim.edge_shed(DEL_T, i)
                Swim.wake_shed(DEL_T, i)
        solve_phi(Swimmers, RHO, DEL_T, i)
        for Swim in Swimmers:
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.body_plot(Swimmers[0].Edge, Swimmers[0].Body)
        DIO.write_data(P, i, DEL_T, SwiP, GeoP, MotP, Swimmers)

    else:
        if np.fmod(i,P['VERBOSITY']) == 0:
            po().timestep_header(i,T[i])

        for Swim in Swimmers:
            Swim.Body.panel_positions(DSTEP, T[i], P['THETA'][i], P['HEAVE'][i])
            Swim.Body.surface_kinematics(DSTEP, TSTEP, P['THETA_MINUS'][i], P['THETA_PLUS'][i], P['HEAVE_MINUS'][i], P['HEAVE_PLUS'][i], DEL_T, T[i], i)
            Swim.edge_shed(DEL_T, i)
            Swim.wake_shed(DEL_T, i)
        solve_phi(Swimmers, RHO, DEL_T, i)

        if np.fmod(i,P['VERBOSITY']) == 0:
            po().solution_output(0,0,0,0,0,0)
            po().solution_complete_output(i/float(COUNTER-1)*100.)
        wake_rollup(Swimmers, DEL_T, i)
        for Swim in Swimmers:
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.body_plot(Swimmers[0].Edge, Swimmers[0].Body)
        DIO.write_data(P, i, DEL_T, SwiP, GeoP, MotP, Swimmers)

total_time = time.time()-start_time
print "Simulation time:", np.round(total_time, 3), "seconds"

graph.body_wake_plot(Swimmers)