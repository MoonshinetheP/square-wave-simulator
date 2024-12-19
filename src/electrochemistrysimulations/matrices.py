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
Acknowledgements:   Oliver Rodriguez (oliverrdz), Guy Denuault

Filename:           mechanism.py

===================================================================================================

Description:

This is the 

===================================================================================================

How to use this file:
    
    1.  

===================================================================================================

Note:


===================================================================================================
'''

import sys
import string
import collections as col
from scipy.sparse import diags as diagonals

import numpy as np


class Matrices:
    def __init__(self, geometry, spatial, coordinates, mechanism):
        
        self.geometry = geometry 
        self.spatial = spatial
        self.coordinates = coordinates
        self.n = self.coordinates[0]
        self.m = self.coordinates[1]
        self.mechanism = mechanism

        if self.geometry == '1D':
            self.x = self.spatial[0]
            for iy in self.mechanism.markers:
                iy.update({'Concentration array': np.ones((self.n, self.m)) * iy['Concentration']})
                
                self.alpha = np.ones(self.n -1)
                self.beta = np.ones(self.n)
                self.gamma = np.ones(self.n - 1)

                for iz in range(1, self.n):

                    try: 
                        self.xplus = self.x[iz + 1] - self.x[iz]
                    except: pass

                    self.xminus = self.x[iz] - self.x[iz - 1]
                    self.denominator = 1 / (self.xminus * (self.xplus ** 2) + self.xplus * (self.xminus **2))
                    
                    self.alpha[iz - 1] *= 2 * iy['Diffusion coefficient'] * self.xplus * self.denominator
                    self.beta[iz] *= 1 - (2 * iy['Diffusion coefficient'] * (self.xminus + self.xplus) * self.denominator)
                                        
                    try:
                        self.gamma[iz] *= 2 * iy['Diffusion coefficient'] * self.xminus * self.denominator
                    except: pass


                matrix = diagonals([self.alpha, self.beta, self.gamma], [-1,0,1]).toarray()
                matrix[0,:] = np.zeros(self.n)        
                matrix[-1,:] = np.zeros(self.n)
                matrix[0,0] = 1
                matrix[-1,-1] = 1
                iy.update({'Matrix' : matrix})
        




if __name__ == "__main__":

    pass