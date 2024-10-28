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

Filename:           simulations.py

===================================================================================================

Description:

This is the main simulation module of the electrochemistry-simulations package and can be used to 
generate simulations of diffusion-based electrochemical reactions. The module uses the output of 
the waveforms module to prepare a .txt file containing the time, potential, and current. 

===================================================================================================

How to use this file:
    
    1.  Scroll down to the section titled 'RUNNING THE SIMULATION FROM MAIN' (near the bottom)
    2.  Describe the waveform using the appropriate class from the waveforms (wf) module
    3.  Define the parameters of the simulation in the 'Diffusive' class
    4.  Run the file in Python
    
The programme will attempt to make a new folder in the current working directory to keep all the 
data. The name of this folder is /data by default, but this can be edited. If the /data folder 
already exists, then this programme will save data into this folder.

===================================================================================================

Note:

The code within this file is split into several sections. Each section starts with a title and a 
brief description of what the code within it does. Comments are placed throughout the code, with 
one hashtag (#) indicating a description of what the code below does and two hashtags (##) 
indicating the formal definition of a variable.

===================================================================================================
'''

import sys

class Reaction:
    def __init__(self, left, right, c, D, z):
        
        self.left = left
        self.right = right
        self.c = c
        self.D = D
        self.z = z
        
        '''DATA TYPE ERRORS'''
        if len(self.left) == 0:
            print('\n' + 'No species were entered on the left hand side' + '\n')
            sys.exit()
        
        if len(self.left) == 1:
            if isinstance(self.ledt,(str)) is False:
                print('\n' + ' numerical value' + '\n')
                sys.exit()
        

        '''DATA VALUE ERRORS'''
        if self.Eupp == self.Elow:
            print('\n' + 'Upper and lower vertex potentials must be different values' + '\n')
            sys.exit()




class E(Reaction):
    def __init__(self, left, right, c, D, z, E0, k0, a):
        super().__init__(left, right, c, D, z,)

        self.E0 = E0
        self.k0 = k0
        self.a = a

        

class C:
    def __init__(self, left, right, c, D, z, k1):
        super().__init__(left, right, c, D, z,)

        self.k1 = k1
    


class Reactions:
    """Packs the specified """
    def __init__(self, **kwargs):
        
        if kwargs == None:
            print('\n' + 'No reactions were entered' + '\n')
            sys.exit()
        for reaction in kwargs:
            if isinstance(reaction,(E, C)) is False:
                print('\n' + 'One or more of the entered reactions is not an E or a C type' + '\n')
                sys.exit() 
            pass
        
        # Make potential specific to electroactive species (maybe in E)
        # Make diagonals
        # Prepare boundary condition?
        # Define oxidation or reduction
        # probably I need two classes, one to sum up the reactions and put them in context, one to make matrices once the distance array and time array are known


class Matrices:
    def __init__(self, dimensions, **kwargs):
        
        self.dimensions = dimensions 
        pass

if __name__ == "__main__":

    instance = Reactions()