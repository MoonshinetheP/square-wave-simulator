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



class Reaction:
    """Parent class for all reaction types."""
    def __init__(self, species, molarities, charges, concentrations, diffusion):
        
        self.species = species
        self.molarities = molarities
        self.charges = charges
        self.concentrations = concentrations
        self.diffusion = diffusion
        

        '''DATA TYPE ERRORS'''
        if isinstance(self.species, (tuple)) is False:
            print('\n' + 'Information about species needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit() 
        
        for ix in self.species:
            if isinstance(ix, (list)) is False:
                print('\n' + 'Information about species needs to be entered as lists inside a (2,1) tuple.' + '\n')
                sys.exit()

            if len(ix) ==0:
                print('\n' + 'There is a reaction with no species entered on one or more sides.' + '\n')
                sys.exit() 

        for iy,iz in zip(self.species[0], self.species[1]):
            if isinstance(iy, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()
            
            if isinstance(iz, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()


        if isinstance(self.molarities, (tuple)) is False:
            print('\n' + 'Information about molarities needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit() 

        for ix in self.molarities:
            if isinstance(ix, (list)) is False:
                print('\n' + 'Information about molarities needs to be entered as lists inside a (2,1) tuple.' + '\n')
                sys.exit()
            if len(ix) ==0:
                print('\n' + 'There is a reaction with no molarities entered on the left side.' + '\n')
                sys.exit()

        for iy,iz in zip(self.molarities[0], self.molarities[1]):
            if isinstance(iy, (int)) is False:
                print('\n' + 'One or more of the entered molarities is not an integer' + '\n')
                sys.exit()
            if isinstance(iz, (int)) is False:
                print('\n' + 'One or more of the entered molarities is not an integer' + '\n')
                sys.exit()

        if isinstance(self.charges, (tuple)) is False:
            print('\n' + 'Information about charges needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit()

        for ix in self.charges:
            if isinstance(ix, (list)) is False:
                print('\n' + 'Information about charges needs to be entered as lists inside a (2,1) tuple.' + '\n')
                sys.exit()
            if len(ix) ==0:
                print('\n' + 'There is a reaction with no species entered on the left side.' + '\n')
                sys.exit() 

        for iy,iz in zip(self.species[0], self.species[1]):
            if isinstance(iy, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()
            if isinstance(iz, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()

        if isinstance(self.concentrations, (tuple)) is False:
            print('\n' + 'Information about concentrations needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit() 

        for ix in self.concentrations:
            if isinstance(ix, (list)) is False:
                print('\n' + 'Information about concentrations needs to be entered as lists inside a (2,1) tuple.' + '\n')
                sys.exit()
            if len(ix) ==0:
                print('\n' + 'There is a reaction with no species entered on the left side.' + '\n')
                sys.exit() 

        for iy,iz in zip(self.species[0], self.species[1]):
            if isinstance(iy, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()
            if isinstance(iz, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()

        if isinstance(self.diffusion, (tuple)) is False:
            print('\n' + 'Information about diffusion coefficients needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit() 

        for ix in self.diffusion:
            if isinstance(ix, (list)) is False:
                print('\n' + 'Information about diffusion coefficients needs to be entered as lists inside a (2,1) tuple.' + '\n')
                sys.exit()
            if len(iy) ==0:
                print('\n' + 'There is a reaction with no species entered on the left side.' + '\n')
                sys.exit()

        for iy,iz in zip(self.species[0], self.species[1]):
            if isinstance(iy, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()
            if isinstance(iz, (str)) is False:
                print('\n' + 'One or more of the entered species is not a string' + '\n')
                sys.exit()

       

        '''DATA VALUE ERRORS'''

        #molarities 0
        #all linkers 0
        #diffusion 0
        #molarities, diffusion or concentration negative
        #charges don't balance



class E(Reaction):
    """Child class for a generic electrochemical reaction. \n
    \n
    species is a (2, 1) tuple representing the two sides of a reaction
    """
    def __init__(self, species, molarities, charges, concentrations, diffusion, E0, k0, a):
        super().__init__(species, molarities, charges, concentrations, diffusion)

        self.E0 = E0
        self.k0 = k0
        self.a = a


        '''DATA TYPE ERRORS'''
        if isinstance(self.E0, (float)) is False:
            print('\n' + 'Information about species needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit() 

        if isinstance(self.k0, (float)) is False:
            print('\n' + 'Information about species needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit()

        if isinstance(self.a, (float)) is False:
            print('\n' + 'Information about species needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit() 


        '''DATA VALUE ERRORS'''
        if self.k0 < 0:
            print('\n' + 'Upper and lower vertex potentials must be different values' + '\n')
            sys.exit()
        
        if self.a < 0 or self.a > 1:
            print('\n' + 'Upper and lower vertex potentials must be different values' + '\n')
            sys.exit()
        


class C(Reaction):
    """Child class for a generic chemical reaction """
    def __init__(self, species, molarities, charges, concentrations, diffusion, k1):
        super().__init__(species, molarities, charges, concentrations, diffusion)

        self.k1 = k1
        
        '''DATA TYPE ERRORS'''
        if isinstance(self.k1, (float)) is False:
            print('\n' + 'Information about species needs to be entered as lists inside a (2,1) tuple.' + '\n')
            sys.exit() 


        '''DATA VALUE ERRORS'''
        if self.k1 < 0:
            print('\n' + 'Upper and lower vertex potentials must be different values' + '\n')
            sys.exit()
    


class Reactions:
    """Reads chemical information from the inputted reaction instances and packages into a relevant format for each species present."""
    def __init__(self, *args):
        
        if args == None:
            print('\n' + 'No reactions were entered' + '\n')
            sys.exit()
        else:
            species = []
            concs = []
            diffs = []
            for reaction in args:
                if isinstance(reaction,(E, C)) is False:
                    print('\n' + 'One or more of the entered reactions is not an E or a C type' + '\n')
                    sys.exit()
                for m,n,o in zip(reaction.species, reaction.concentrations, reaction.diffusion):
                    if species.count(m[0]) == 0:
                        species.append(m[0])
                        concs.append(n[0])
                        diffs.append(o[0])   
                    else: pass
                

        markerstrings = list(string.ascii_uppercase)[0:len(species)]
        self.markers = col.namedtuple('Markers', markerstrings)._make({} for _ in markerstrings)

        self.cmax = max(concs)
        self.Dmax = max(diffs)

        for n in range(0, len(species)):
            self.markers[n].update({'Marker': markerstrings[n]})            
            self.markers[n].update({'Species':species[n]})
            self.markers[n].update({'Concentration': concs[n]/self.cmax})
            self.markers[n].update({'Diffusion coefficient': diffs[n]/self.Dmax})
            self.markers[n].update({'Charge': None})
            self.markers[n].update({'Molarity': None})
            self.markers[n].update({'Oxidised': []})
            self.markers[n].update({'Reduced': []})
            self.markers[n].update({'Consumed': []})
            self.markers[n].update({'Produced': []})

        for reaction in args:
            sides = [index for index, char in enumerate(reaction.species)]
            if isinstance(reaction, (E)):
                for ix in sides:
                    
                    if ix == 0:
                        iy = 1
                    else: iy = 0

                    participant = [char for index, char in enumerate(reaction.species[ix])]
                    otherparticipant = [char for index, char in enumerate(reaction.species[iy])]
                    for iz in reaction.species[ix]:
                        for n in range(0, len(species)):
                            if self.markers[n]['Species'] == iz:
                                if reaction.charges[ix][0] < reaction.charges[iy][0]:
                                    entry = [(participant.remove(iz) ,otherparticipant), reaction.E0, reaction.k0, reaction.a]
                                    if entry in self.markers[n]['Oxidised'] is True:
                                        pass#indexing 0 each time not a good idea, need to distinguish between left and right
                                    else:
                                        self.markers[n].update({'Oxidised': entry})
                                if reaction.charges[ix][0] > reaction.charges[iy][0]:
                                    entry = [(participant.remove(iz), otherparticipant), reaction.E0, reaction.k0, reaction.a]
                                    if entry in self.markers[n]['Reduced'] is True:
                                        pass
                                    else:
                                        self.markers[n].update({'Reduced': entry})
                        else: pass

            
            if isinstance(reaction, (C)):
                for ix in sides:
                    
                    if ix == 0:
                        iy = 1
                    else: iy = 0

                    participant = [char for index, char in enumerate(reaction.species[ix])]
                    otherparticipant = [char for index, char in enumerate(reaction.species[iy])]
                    for iz in reaction.species[ix]:
                        for n in range(0, len(species)):
                            if self.markers[n]['Species'] == iz:
                                if ix == 0:
                                    entry = [(participant.remove(iz) ,otherparticipant), reaction.k1]
                                    if entry in self.markers[n]['Consumed'] is True:
                                        pass
                                    else:
                                        self.markers[n].update({'Consumed': entry})
                                if ix == 1:
                                    entry = [(participant.remove(iz), otherparticipant), reaction.k1]
                                    if entry in self.markers[n]['Produced'] is True:
                                        pass
                                    else:
                                        self.markers[n].update({'Produced': entry})




if __name__ == "__main__":

    E1 = E((['G'], ['H']), ([1],[1]), ([2],[3]), ([0.005],[0]), ([5E-6],[5E-6]), E0 = 0.2, k0 = 0.1, a = 0.5)
    E2 = E((['J'], ['I']), ([1],[1]), ([2],[1]), ([0],[0]), ([6E-6],[7E-6]), E0 = 0.4, k0 = 0.1, a = 0.5)
    C1 = C((['H'], ['J']), ([1],[1]), ([3],[2]), ([0],[0]), ([5E-6],[5E-6]), k1 = 0.1)

    instance = Reactions(E1, C1, E2)