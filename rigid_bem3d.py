#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-3D
A 3D boundary element method code

"""
import time
import numpy as np
from data_IO_class import DataIO
from input_parameters import PARAMETERS as P
from swimmer_class import Swimmer
import parameter_classes as PC
from terminal_output import print_output as po
from functions_general import simulation_startup, archive
import functions_graphics as graph

po().prog_title('1.0.082515a')
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
po().initialize_output(P['T'][START_COUNTER])

outerCorr = 1
for i in xrange(START_COUNTER, COUNTER):
    if i == 0:
        for Swim in Swimmers:
                Swim.Body.panel_positions(DSTEP, T[i], P['THETA'][i], P['HEAVE'][i])
                Swim.Body.surface_kinematics(DSTEP, TSTEP, P['THETA_MINUS'][i], P['THETA_PLUS'][i], P['HEAVE_MINUS'][i], P['HEAVE_PLUS'][i], DEL_T, T[i], i)
                Swim.edge_shed(DEL_T, i)
                Swim.wake_shed(DEL_T, i)
                Swim.Body.force(P['THETA'][i], RHO, P['V0'], P['C'], 1.0, i, P['SW_SV_FORCES'])
#        solve_phi(Swimmers, RHO, DEL_T, i, outerCorr)
        for Swim in Swimmers:
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.plot_n_go(Swimmers, P['V0'], P['T'][i], P['HEAVE'][i], i, P['SW_PLOT_FIG'])
        DIO.write_data(P, i, DEL_T, SwiP, GeoP, MotP, Swimmers)

    else:
        if np.fmod(i,P['VERBOSITY']) == 0:
            po().timestep_header(i,T[i])

        for Swim in Swimmers:
            Swim.Body.panel_positions(DSTEP, T[i], P['THETA'][i], P['HEAVE'][i])
            Swim.Body.surface_kinematics(DSTEP, TSTEP, P['THETA_MINUS'][i], P['THETA_PLUS'][i], P['HEAVE_MINUS'][i], P['HEAVE_PLUS'][i], DEL_T, T[i], i)
            Swim.edge_shed(DEL_T, i)
            Swim.wake_shed(DEL_T, i)
            Swim.Body.force(P['THETA'][i], RHO, P['V0'], P['C'], 1.0, i, P['SW_SV_FORCES'])
#        solve_phi(Swimmers, RHO, DEL_T, i, outerCorr)
#        wake_rollup(Swimmers, DEL_T, i, P)
        for Swim in Swimmers:
            if np.fmod(i,P['VERBOSITY']) == 0:
                po().solution_output(Swim.Body.Cf, Swim.Body. Cl,Swim.Body.Ct,Swim.Body.Cpow)
                po().solution_complete_output(i/float(COUNTER-1)*100.)
            archive(Swim.Body.AF.x_mid)
            archive(Swim.Body.AF.z_mid)
        graph.plot_n_go(Swimmers, P['V0'], P['T'][i], P['HEAVE'][i], i, P['SW_PLOT_FIG'])
        DIO.write_data(P, i, DEL_T, SwiP, GeoP, MotP, Swimmers)

graph.basic_xyz(Swimmers[0].Body.BF.x,Swimmers[0].Body.BF.y,Swimmers[0].Body.BF.z)