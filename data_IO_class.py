#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
BEM-2D
A 2D boundary element method code

"""
import os
import pickle
import numpy as np

class DataIO(object):
    def __init__(self, P):
        self.OUTPUT_DIR = P['OUTPUT_DIR']
        
        # Determine if the output directory exists. If not, create the directory.
        if not os.path.exists(self.OUTPUT_DIR):
            os.makedirs(self.OUTPUT_DIR)
            
    def write_data(self, P, i, DEL_T, SwiP, GeoP, MotP, Swimmers, solid=None, FSI=None, PyFEA=None):
        """
        Writes the simulation state to a binary file named with the current flow time.
        
        Args:
            P (dict):
            i (int):
            DEL_T (float):
            SwiP (list):
            GeoP (list):
            MotP (list):
            Swimmers (list):
            solid (Optional [list]):
            FSI (Optional [list]):
            PyFEA (Optional [list]):     
        """
        if (np.fmod(i,P['SAVE_EVERY']) == 0 and P['SW_SAVE_DATA'] == True):
            outfile = ''.join((self.OUTPUT_DIR, "/%.8f" % (i * DEL_T)))
            with open(outfile, 'wb') as f:
                pickle.dump([P, i, i*DEL_T, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA], f, 1)
            
    def read_data(self, INPUT_FILE):
        with open(INPUT_FILE, 'rb') as f:
            P, i, FLOWTIME, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA = pickle.load(f)
            
        return (P, i, FLOWTIME, SwiP, GeoP, MotP, Swimmers, solid, FSI, PyFEA)