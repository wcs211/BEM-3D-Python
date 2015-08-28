#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-3D
A 3D boundary element method code

"""
from input_parameters import PARAMETERS as P
from terminal_output import print_output as po

po().prog_title('1.0.082815a')

if (P['SW_FSI'] == True):
    # Run the Fluid Structure Interaction Solver with the BEM Solver
    execfile("FSI_bem3d.py")
else:
    # Run the Boundary Element Solver
    execfile("rigid_bem3d.py")
