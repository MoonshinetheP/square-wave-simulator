'''
===================================================================================================
Copyright (C) 2023 Steven Linfield

This file is part of the electrochemistry-simulations package. This package is free software: you 
can redistribute it and/or modify it under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, or (at your option) any later 
version. This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details. You should have received a copy of the GNU General 
Public License along with electrochemistry-simulations. If not, see https://www.gnu.org/licenses/
===================================================================================================

Package title:      electrochemistry-simulations
Repository:         https://github.com/MoonshinetheP/electrochemistry-simulations
Date of creation:   09/03/2023
Main author:        Steven Linfield (MoonshinetheP)
Collaborators:      None
Acknowledgements:   Oliver Rodriguez (oliverrdz)
    
Filename:           training.py

===================================================================================================
How to use this file:
    

===================================================================================================
'''
import sys
import os
import time

import numpy as np
import waveforms as wf
import simulations as sim
import capacitance as cap
import noise as noise

from errno import EEXIST


if __name__ == '__main__':

    '''1. MAKE A /DATA FOLDER'''
    cwd = os.getcwd()    
    try:
        os.makedirs(cwd + '/data')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/data'):
            pass
        else: 
            raise
    
    '''2. DEFINE THE START TIME'''
    start = time.time()

    '''3. DESCRIBE THE WAVEFORM'''
    shape = wf.DPV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.02, pt = 0.5, rt = 1.5, st = 0.01, detailed = True)
    
    '''4. RUN THE SIMULATION'''
    instance = sim.E(input = shape, E0 = 0.25, k0 = 0.1, a = 0.5, cR = 0.005, cO = 0.000, DR = 5E-6, DO = 5E-6, r = 0.15, expansion = 1.05, Nernstian = False, BV = True, MH = False)
    
    '''5. DEFINE THE END TIME'''
    end = time.time()
    print(end-start)

    '''6. SAVE THE DATA'''
    filepath = cwd + '/data/' + f'{shape.subtype}' + '.txt'
    with open(filepath, 'w') as file:
        for ix, iy, iz in instance.results():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')