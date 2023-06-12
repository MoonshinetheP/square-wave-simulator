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
    
Filename:           noise.py

===================================================================================================
How to use this file:
    

===================================================================================================
'''



import sys
import os
import time

import numpy as np
import waveforms as wf

from errno import EEXIST
from scipy.stats import boltzmann




class Noise:
    def __init__(self, input, electrical = False, shot = False, thermal = False):
        self.input = input
        self.electrical = electrical
        self.shot = shot
        self.thermal = thermal

        self.t = self.input.tPLOT

        if self.electricalnoise == True:
            self.electrical()

        if self.shotnoise == True:
            self.shot()

        if self.thermalnoise == True:
            self.thermal()
        
        #self.report(self, self.ie, self.ir, self.it)

    def electricalnoise(self):
        self.i = (1E-6)*np.sin(2*np.pi*50*self.t)

    def shotnoise(self):
        lambda_, N = 1.4, 19
        mean, var, skew, kurt = boltzmann.stats(lambda_, N, moments='mvsk')

        x = np.arange(boltzmann.ppf(0.01, lambda_, N),
        boltzmann.ppf(0.99, lambda_, N))
        self.ir = np.zeros((self.t))

    def thermalnoise(self):
        self.it = np.zeros((self.t))

    def report(self, ie, ir, it):
        i = ie + ir + it
        return i


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
    shape = wf.CV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.5, ns = 1)    
    
    '''4. RUN THE SIMULATION'''
    instance = Noise(input = shape, electrical = True)

    '''5. DEFINE THE END TIME'''
    end = time.time()
    print(end-start)

    '''6. SAVE THE DATA'''
    filepath = cwd + '/data/' + 'noise' + '.txt'
    with open(filepath, 'w') as file:
        for ix, iy in instance.electricalnoise():
            file.write(str(ix) + ',' + str(iy) + '\n')