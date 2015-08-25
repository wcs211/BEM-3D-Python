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
from functions_general import simulation_startup
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



graph.basic_xyz(Swimmers[0].Body.BF.x,Swimmers[0].Body.BF.y,Swimmers[0].Body.BF.z)