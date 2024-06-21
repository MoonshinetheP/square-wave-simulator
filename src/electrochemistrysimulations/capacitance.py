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
    
Filename:           capacitance.py

===================================================================================================
How to use this file:
    

===================================================================================================
'''



import sys
import os
import time
import numpy as np

from errno import EEXIST

import waveforms as wf


class Capacitance:
    '''Simulates the charging of a double layer when using CV or CSV \n 
    Also takes input parameters and uses them to create potential waveforms for CV and CSV'''
    
    def __init__(self, input, Cd = 0.000250, Ru = 500):
        '''Defines the parameters of the simulation, checks for errors, and makes a potential waveform to be used by other functions'''

        self.input = input
        self.Cd = Cd            # Double layer capacitance in F
        self.Ru = Ru            # Uncompensated resistance in Ohms
        
        if self.input.detailed == False:
            self.simple()
        if self.input.detailed == True:
            self.detailed() 


    def simple(self):
        '''Returns E vs. i for a CV performed on a capacitor with parameters derived from the Capacitance() class\n
        Uses equation 1.6.23 from the 3rd edition of Electrochemical Methods:\n
        i = sr*Cd*(1-np.exp(-t/(Ru*Cd)))'''
        self.Eini = self.input.Eini        # Start potential in V
        self.Eupp = self.input.Eupp        # Upper vertex potential in V
        self.Elow = self.input.Elow        # Lower vertex potential in V
        self.dE = self.input.dE            # Step size in V
        self.sr = self.input.sr            # Scan rate in V/s
        self.ns = self.input.ns            # Number of scans (no unit)

        if self.Eini == self.Elow:
            self.i = np.array([])
            for iy in range(0, self.ns):
                self.i = np.append(self.i, self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                self.i = np.append(self.i, -self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
        
        if self.Eini == self.Eupp:
            self.i = np.array([])
            for iy in range(0, self.ns):
                self.i = np.append(self.i, -self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                self.i = np.append(self.i, self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
        
        if self.Elow < self.Eini < self.Eupp:  
            self.i = np.array([])
            if self.dE > 0:
                for iy in range(0, self.ns):
                    if iy == 0:
                        self.i = np.append(self.i, self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.uppdp]) / (self.Ru * self.Cd))))
                        self.i = np.append(self.i, -self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                    elif iy < self.ns:
                        self.i = np.append(self.i, self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                        self.i = np.append(self.i, -self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                
                self.i = np.append(self.i, self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.lowdp]) / (self.Ru * self.Cd))))

            if self.dE < 0:
                for iy in range(0, self.ns):
                    if iy == 0:
                        self.i = np.append(self.i, -self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.lowdp]) / (self.Ru * self.Cd))))
                        self.i = np.append(self.i, self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                    elif iy < self.ns:
                        self.i = np.append(self.i, -self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                        self.i = np.append(self.i, self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.dp]) / (self.Ru * self.Cd))))
                
                self.i = np.append(self.i, -self.sr * self.Cd * (1 - np.exp((-self.input.t[:self.input.uppdp]) / (self.Ru * self.Cd))))


    def detailed(self):
        '''Returns E vs. i for a CSV performed on a capacitor with parameters derived from the Capacitance() class\n
        Uses equation 1.6.17 from the 3rd edition of Electrochemical Methods:\n
        i = (dE/Ru)*np.exp(-t/(Ru*Cd))'''
        self.Eini = self.input.Eini        # Start potential in V
        self.Eupp = self.input.Eupp        # Upper vertex potential in V
        self.Elow = self.input.Elow        # Lower vertex potential in V
        self.dE = np.abs(self.input.dE)            # Step size in V
        self.sr = self.input.sr            # Scan rate in V/s
        self.ns = self.input.ns            # Number of scans (no unit)
        self.st = self.input.st
        self.detailed = self.input.detailed
        

        if self.Eini == self.Elow:
            self.i = np.array([])
            self.iplus = np.zeros(self.input.sp * self.input.dp)
            self.iminus = np.zeros(self.input.sp * self.input.dp)
            for ix in range(0, self.ns):
                for iy in range(0, self.input.dp):
                    space = int(iy * self.input.sp)
                    self.iplus[space:] = np.add(self.iplus[space:], (self.dE/self.Ru) * np.exp((-self.input.tWF[:self.input.sp * self.input.dp - space]) / (self.Ru * self.Cd)))
                    if self.input.sampled == True:
                        self.i = np.append(self.i, np.average(self.iplus[space + self.input.sp - round(self.input.sp * (1 - self.input.alpha)): space + self.input.sp]))
                for iy in range(0, self.input.dp):
                    space = int(iy * self.input.sp)
                    self.iminus[space:] = np.add(self.iminus[space:], (-self.dE/self.Ru) * np.exp((-self.input.tWF[:self.input.sp * self.input.dp - space]) / (self.Ru * self.Cd)))
                    if self.input.sampled == True:
                        self.i = np.append(self.i, -np.average(self.iplus[space + self.input.sp - round(self.input.sp * (1 - self.input.alpha)): space + self.input.sp]))
                if self.detailed == True and self.input.sampled == False:
                    self.i = np.append(self.i, self.iplus)
                    self.i = np.append(self.i, self.iminus)
                    self.iplus = np.zeros(self.input.sp * self.input.dp)
                    self.iminus = np.zeros(self.input.sp * self.input.dp)

        if self.Eini == self.Eupp:
            self.i = np.array([])
            self.iplus = np.zeros(self.input.sp * self.input.dp)
            self.iminus = np.zeros(self.input.sp * self.input.dp)
            for ix in range(0, self.ns):
                for iy in range(0, self.input.dp):
                    space = int(iy * self.input.sp)
                    self.iminus[space:] = np.add(self.iminus[space:], (-self.dE/self.Ru) * np.exp((-self.input.tWF[:self.input.sp * self.input.dp - space]) / (self.Ru * self.Cd)))
                for iy in range(0, self.input.dp):
                    space = int(iy * self.input.sp)
                    self.iplus[space:] = np.add(self.iplus[space:], (self.dE/self.Ru) * np.exp((-self.input.tWF[:self.input.sp * self.input.dp - space]) / (self.Ru * self.Cd)))
                self.i = np.append(self.i, self.iminus)
                self.i = np.append(self.i, self.iplus)
                self.iplus = np.zeros(self.input.sp * self.input.dp)
                self.iminus = np.zeros(self.input.sp * self.input.dp)

        if self.Elow < self.Eini < self.Eupp:  
            self.i = np.array([])
            self.iupp = np.zeros(self.input.sp * self.input.uppdp)
            self.ilow = np.zeros(self.input.sp * self.input.lowdp)
            self.iplus = np.zeros(self.input.sp * self.input.dp)
            self.iminus = np.zeros(self.input.sp * self.input.dp)
            if self.dE > 0:
                for ix in range(0, self.ns):
                    for iy in range(0, self.input.dp):
                        space = int(iy * self.input.sp)
                        self.iplus[space:] = np.add(self.iplus[space:], (self.dE/self.Ru) * np.exp((-self.input.tWF[:self.input.sp * self.input.dp - space]) / (self.Ru * self.Cd)))
                    for iy in range(0, self.input.dp):
                        space = int(iy * self.input.sp)
                        self.iminus[space:] = np.add(self.iminus[space:], (-self.dE/self.Ru) * np.exp((-self.input.tWF[:self.input.sp * self.input.dp - space]) / (self.Ru * self.Cd)))
                    self.i = np.append(self.i, self.iplus)
                    self.i = np.append(self.i, self.iminus)
                    self.iplus = np.zeros(self.input.sp * self.input.dp)
                    self.iminus = np.zeros(self.input.sp * self.input.dp)
            if self.dE < 0:
                pass

    
    def results(self):
        '''Returns the output'''
        self.output = zip(self.input.tPLOT, self.input.EPLOT, self.i)
        return self.output