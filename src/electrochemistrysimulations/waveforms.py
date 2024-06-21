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
    
Filename:           waveforms.py

===================================================================================================
How to use this file:
    
    1. Scroll down to the section titled 'MAKING WAVEFORMS FROM MAIN'
===================================================================================================
'''




import sys
import os
import time
from errno import EEXIST
import numpy as np
from scipy.signal import square



"""PARENT CLASSES"""
class sweep:
    '''Parent class for all sweep type waveforms'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):
        
        self.type = 'sweep'
        self.detailed = False

        self.Eini = Eini        # Start potential
        self.Eupp = Eupp        # Upper vertex potential
        self.Elow = Elow        # Lower vertex potential
        self.dE = dE            # Step size (i.e. the number of data points)
        self.sr = sr            # Scan rate
        self.ns = ns            # Number of scans for cyclic voltammetry

        '''DATATYPE ERRORS'''
        if isinstance(self.Eini, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the start potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.Eupp, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the upper vertex potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()        
        if isinstance(self.Elow, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the lower vertex potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()      
        if isinstance(self.dE, (float)) is False:
            print('\n' + 'An invalid datatype was used for the step potential. Enter a float value corresponding to a potential in V.' + '\n')
            sys.exit()    
        if isinstance(self.sr, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the scan rate. Enter a float or an integer value corresponding to the scan rate in V/s.' + '\n')
            sys.exit() 
        if isinstance(self.ns, (int)) is False:
            print('\n' + 'An invalid datatype was used for the number of scans. Enter an integer value corresponding to the scan rate in V/s.' + '\n')
            sys.exit() 

        '''DATA VALUE ERRORS'''
        if self.Eupp == self.Elow:
            print('\n' + 'Upper and lower vertex potentials must be different values' + '\n')
            sys.exit()
        if self.Eupp < self.Elow:
            print('\n' + 'Upper vertex potential must be greater than lower vertex potential' + '\n')
            sys.exit()  
        if self.Eini < self.Elow:
            print('\n' + 'Start potential must be higher than or equal to the lower vertex potential' + '\n')
            sys.exit()
        if self.Eini > self.Eupp:
            print('\n' + 'Start potential must be lower than or equal to the upper vertex potential' + '\n')
            sys.exit()
        if self.dE == 0:
            print('\n' + 'Step potential must be a non-zero value' + '\n')
            sys.exit()
        if np.abs(self.dE) > (self.Eupp - self.Elow):
            print('\n' + 'Step potential cannot be greater than potential window' + '\n')
            sys.exit()
        if self.Eini == self.Elow and self.dE < 0:
            print('\n' + 'Step potential must be a positive value for a positive scan direction' + '\n')
            sys.exit()
        if self.Eini == self.Eupp and self.dE > 0:
            print('\n' + 'Step potential must be a negative value for a negative scan direction' + '\n')
            sys.exit()
        if self.sr <= 0:
            print('\n' + 'Scan rate must be a positive non-zero value' + '\n')
            sys.exit()
        if self.ns <=0:
            print('\n' + 'Number of scans must be a positive non-zero value' + '\n')
            sys.exit()

    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped
            
class step:
    '''Parent class for all step type waveforms'''
    def __init__(self, dE, dt, st):
        
        self.type = 'step'
        
        self.dE = dE            # Potential for single step chronoamperommetry
        self.dt = dt            # Step period
        self.st = st
        
        '''DATATYPE ERRORS'''
        if isinstance(self.dE, (list)) is True and isinstance(self.dt, (list)) is True:
            self.multiple = True
            if len(dE) != len(dt):
                print('\n' + 'The list of potentials, the list of times, and the list of data points were not equal lengths.' + '\n')
                sys.exit() 
            for ix in self.dE:
                if isinstance(ix, (float, int)) is False:
                    print('\n' + 'An invalid datatype was used for one or more potentials. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
                    sys.exit()
            for iy in self.dt:
                if isinstance(iy, (float, int)) is False:
                    print('\n' + 'An invalid datatype was used for one or more step times. Enter either a float or an integer value corresponding to a time in s.' + '\n')
                    sys.exit()   

        elif isinstance(self.dE, (list)) is False and isinstance(self.dt, (list)) is False:
            self.multiple = False
            if isinstance(self.dE, (tuple, dict)) is True and isinstance(self.dt, (tuple, dict)) is True:
                print('\n' + 'When entering multiple steps, make sure to enter potentials and times as lists and not tuples or dictionaries' + '\n')
                sys.exit()
            if isinstance(self.dE, (float, int)) is False:
                print('\n' + 'An invalid datatype was used for the potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
                sys.exit()
            if isinstance(self.dt, (float, int)) is False:
                print('\n' + 'An invalid datatype was used for the step time. Enter either a float or an integer value corresponding to a time in s.' + '\n')
                sys.exit()        

        else:
            print('\n' + 'When entering multiple steps, make sure to enter a separate lists of potentials and times' + '\n')
            sys.exit()
        
        if isinstance(self.st, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the sampling time. Enter either a float or an integer value corresponding to a time in s.' + '\n')
            sys.exit()

        '''DATA VALUE ERRORS'''
        if self.multiple == True:
            for ix in range(0, len(self.dE) -1):
                if self.dE[ix] == self.dE[ix + 1]:
                    print('\n' + 'Adjacent potentials must not be equal' + '\n')
                    sys.exit()        
            for iy in self.dt:
                if iy <=0:
                    print('\n' + 'All step times must have a positive non-zero value' + '\n')
                    sys.exit()

        else:      
            if self.dt <= 0:
                print('\n' + 'Step time must be a positive non-zero value' + '\n')
                sys.exit()
        
        if self.st <= 0:
            print('\n' + 'Sampling time must be a positive non-zero value' + '\n')
            sys.exit()  
    
    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped

class pulse:
    '''Parent class for all pulse type waveforms'''
    def __init__(self, Eini, Efin, dEs, dEp, pt, rt, st, detailed, sampled, alpha):

        self.type = 'pulse'

        self.Eini = Eini        # Start potential for pulsed techniques
        self.Efin = Efin        # End potential for pulsed techniques
        self.dEs = dEs          # Step size for pulsed techniques
        self.dEp = dEp          # Pulse size for pulsed techniques
        self.pt = pt            # Pulse period
        self.rt = rt            # Rest period
        self.st = st            # Sampling period for more detailed waveforms

        self.detailed = detailed
        self.sampled = sampled
        self.alpha = alpha
               
        '''DATATYPE ERRORS'''
        if isinstance(self.Eini, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the start potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.Efin, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the end potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()        
        if isinstance(self.dEs, (float)) is False:
            print('\n' + 'An invalid datatype was used for the step size. Enter a float value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.dEp, (float)) is False:
            print('\n' + 'An invalid datatype was used for the pulse size. Enter a float value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.pt, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the pulse period. Enter either a float or an integer value corresponding to a time in s.' + '\n')
            sys.exit()
        if isinstance(self.rt, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the rest period. Enter either a float or an integer value corresponding to a time in s.' + '\n')
            sys.exit()
        if isinstance(self.st, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the sampling time. Enter either a float or an integer value corresponding to a time in s.' + '\n')
            sys.exit()
        if isinstance(self.detailed, (bool)) is False:
            print('\n' + 'An invalid datatype was used for the detailed argument. Enter either True or False.' + '\n')
            sys.exit()
        if isinstance(self.sampled, (bool)) is False:
            print('\n' + 'An invalid datatype was used for the sampled argument. Enter either True or False.' + '\n')
            sys.exit()
        if isinstance(self.detailed, (bool)) is False:
            print('\n' + 'An invalid datatype was used for the detailed argument. Enter either True or False.' + '\n')
            sys.exit()
        if isinstance(self.alpha,(int,float)) is False:
            print('\n' + 'The alpha parameter must be a numerical value' + '\n')
            sys.exit()

        
                  

        '''DATA VALUE ERRORS'''
        if self.Eini == self.Efin:
            print('\n' + 'Start and end potentials must be different values' + '\n')
            sys.exit()  
        if self.dEs == 0:
            print('\n' + 'Step potential must be a non-zero value' + '\n')
            sys.exit()
        if np.abs(self.dEs) > (self.Efin - self.Eini):
            print('\n' + 'Step potential cannot be greater than potential window' + '\n')
            sys.exit()
        if self.dEp == 0:
            print('\n' + 'Pulse potential must be a non-zero value' + '\n')
            sys.exit()
        if self.Eini < self.Efin and self.dEs < 0:
            print('\n' + 'Step potential must be a positive value for a positive scan direction' + '\n')
            sys.exit()
        if self.Eini > self.Efin and self.dEs > 0:
            print('\n' + 'Step potential must be a negative value for a negative scan direction' + '\n')
            sys.exit()          
        if self.pt <= 0:
            print('\n' + 'Pulse period must be a positive non-zero value' + '\n')
            sys.exit()
        if self.rt <= 0:
            print('\n' + 'Rest period must be a positive non-zero value' + '\n')
            sys.exit() 
        if self.st <= 0:
            print('\n' + 'Sampling time must be a positive non-zero value' + '\n')
            sys.exit()
        if self.st >= self.pt:
            print('\n' + 'Sampling time must be less than the pulse period' + '\n')
            sys.exit()
        if self.st >= self.rt:
            print('\n' + 'Sampling time must be less than the rest period' + '\n')
            sys.exit()
        if (0 <= self.alpha < 1) == False:
            print('\n' + 'The alpha parameter must be between 0 and 1' + '\n')
            sys.exit()


    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped
    
class hybrid:
    '''Parent class for all waveforms composed of both steps and sweeps'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns, st, detailed, sampled, alpha):

        self.type = 'hybrid'

        self.Eini = Eini        # Start potential for sweeping step techniques
        self.Eupp = Eupp        # Upper vertex potential for sweeping step techniques
        self.Elow = Elow        # Lower vertex potential for sweeping step techniques
        self.dE = dE            # Step size for sweeping step techniques
        self.sr = sr            # Scan rate for sweeping step techniques
        self.ns = ns            # Number of scans for cyclic staircase voltammetry
        self.st = st            # Sampling time
        

        self.detailed = detailed
        self.sampled = sampled
        self.alpha = alpha

        '''DATATYPE ERRORS'''
        if isinstance(self.Eini, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the start potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.Eupp, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the upper vertex potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()        
        if isinstance(self.Elow, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the lower vertex potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()      
        if isinstance(self.dE, (float)) is False:
            print('\n' + 'An invalid datatype was used for the step potential. Enter a float value corresponding to a potential in V.' + '\n')
            sys.exit()    
        if isinstance(self.sr, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the scan rate. Enter a float or an integer value corresponding to the scan rate in V/s.' + '\n')
            sys.exit() 
        if isinstance(self.ns, (int)) is False:
            print('\n' + 'An invalid datatype was used for the number of scans. Enter an integer value corresponding to the scan rate in V/s.' + '\n')
            sys.exit() 
        if isinstance(self.st, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the sampling time. Enter either a float or an integer value corresponding to a time in s.' + '\n')
            sys.exit()
        if isinstance(self.detailed, (bool)) is False:
            print('\n' + 'An invalid datatype was used for the detailed argument. Enter either True or False.' + '\n')
            sys.exit()
        if isinstance(self.detailed, (bool)) is False:
            print('\n' + 'An invalid datatype was used for the detailed argument. Enter either True or False.' + '\n')
            sys.exit()
        if isinstance(self.alpha,(int,float)) is False:
            print('\n' + 'The alpha parameter must be a numerical value' + '\n')
            sys.exit()

        '''DATA VALUE ERRORS'''
        if self.Eupp == self.Elow:
            print('\n' + 'Upper and lower vertex potentials must be different values' + '\n')
            sys.exit()
        if self.Eupp < self.Elow:
            print('\n' + 'Upper vertex potential must be greater than lower vertex potential' + '\n')
            sys.exit()  
        if self.Eini < self.Elow:
            print('\n' + 'Start potential must be higher than or equal to the lower vertex potential' + '\n')
            sys.exit()
        if self.Eini > self.Eupp:
            print('\n' + 'Start potential must be lower than or equal to the upper vertex potential' + '\n')
            sys.exit()
        if self.dE == 0:
            print('\n' + 'Step potential must be a non-zero value' + '\n')
            sys.exit()
        if np.abs(self.dE) > (self.Eupp - self.Elow):
            print('\n' + 'Step potential cannot be greater than potential window' + '\n')
            sys.exit()   
        if self.Eini == self.Elow and self.dE < 0:
            print('\n' + 'Step potential must be a positive value for a positive scan direction' + '\n')
            sys.exit()
        if self.Eini == self.Eupp and self.dE > 0:
            print('\n' + 'Step potential must be a negative value for a negative scan direction' + '\n')
            sys.exit()
        if self.sr <= 0:
            print('\n' + 'Scan rate must be a positive non-zero value' + '\n')
            sys.exit()    
        if self.ns <=0:
            print('\n' + 'Number of scans must be a positive non-zero value' + '\n')
            sys.exit()
        if self.st <= 0:
            print('\n' + 'Sampling time must be a positive non-zero value' + '\n')
            sys.exit()
        if (0 <= self.alpha < 1) == False:
            print('\n' + 'The alpha parameter must be between 0 and 1' + '\n')
            sys.exit()

    
    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped

class impedance:
    '''Parent class for impedance waveforms'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):

        self.type = 'impedance'

        self.Eini = Eini        # Start potential for sweeping step techniques
        self.Eupp = Eupp        # Upper vertex potential for sweeping step techniques
        self.Elow = Elow        # Lower vertex potential for sweeping step techniques
        self.dE = dE            # Step size for sweeping step techniques
        self.sr = sr            # Scan rate for sweeping step techniques
        self.ns = ns            # Number of scans for cyclic staircase voltammetry

        '''DATATYPE ERRORS'''
        if isinstance(self.Eini, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the start potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.Eupp, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the upper vertex potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()        
        if isinstance(self.Elow, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the lower vertex potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()      
        if isinstance(self.dE, (float)) is False:
            print('\n' + 'An invalid datatype was used for the step potential. Enter a float value corresponding to a potential in V.' + '\n')
            sys.exit()    
        if isinstance(self.sp, (int)) is False:
            print('\n' + 'An invalid datatype was used for the number of data points in a step. Enter an integer value.' + '\n')
            sys.exit()
        if isinstance(self.sr, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the scan rate. Enter a float or an integer value corresponding to the scan rate in V/s.' + '\n')
            sys.exit() 
        if isinstance(self.ns, (int)) is False:
            print('\n' + 'An invalid datatype was used for the number of scans. Enter an integer value corresponding to the scan rate in V/s.' + '\n')
            sys.exit() 

        '''DATA VALUE ERRORS'''
        if self.Eupp == self.Elow:
            print('\n' + 'Upper and lower vertex potentials must be different values' + '\n')
            sys.exit()
        if self.Eupp < self.Elow:
            print('\n' + 'Upper vertex potential must be greater than lower vertex potential' + '\n')
            sys.exit()  
        if self.Eini < self.Elow:
            print('\n' + 'Start potential must be higher than or equal to the lower vertex potential' + '\n')
            sys.exit()
        if self.Eini > self.Eupp:
            print('\n' + 'Start potential must be lower than or equal to the upper vertex potential' + '\n')
            sys.exit()
        if self.dE == 0:
            print('\n' + 'Step potential must be a non-zero value' + '\n')
            sys.exit()
        if np.abs(self.dE) > (self.Eupp - self.Elow):
            print('\n' + 'Step potential cannot be greater than potential window' + '\n')
            sys.exit()   
        if self.Eini == self.Elow and self.dE < 0:
            print('\n' + 'Step potential must be a positive value for a positive scan direction' + '\n')
            sys.exit()
        if self.Eini == self.Eupp and self.dE > 0:
            print('\n' + 'Step potential must be a negative value for a negative scan direction' + '\n')
            sys.exit()
        if self.sr <= 0:
            print('\n' + 'Scan rate must be a positive non-zero value' + '\n')
            sys.exit()    
        if self.ns <=0:
            print('\n' + 'Number of scans must be a positive non-zero value' + '\n')
            sys.exit()
    
    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped


"""SWEEP CLASSES"""
class LSV(sweep):
    '''Linear Sweep Voltammetry (LSV) waveform\n
    \n
    Eini - start potential \n
    Eupp - upper vertex potential \n
    Elow - lower vertex potential \n
    dE   - step size (in this case, number of data points) \n
    sr   - scan rate \n
    ns   - number of scans (ignored)'''
    def __init__(self,Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        self.subtype = 'LSV'
        self.detailed = False

        if self.Elow < self.Eini < self.Eupp:
            print('\n' + 'Initial potential should be equal to either upper vertex or lower vertex potential in LSV' + '\n')
            sys.exit()
        
        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                 
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
        
            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.round(np.linspace(self.Eini, self.Eupp, self.dp + 1, endpoint = True, dtype = np.float32), 3)
        
        '''STARTING FROM UPPER VERTEX POTENTIAL''' 
        if self.Eini == self.Eupp:                 
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
        
            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.round(np.linspace(self.Eini, self.Elow, self.dp + 1, endpoint = True, dtype = np.float32), 3)
        
        '''PLOTTING WAVEFORM'''
        self.tPLOT = self.t
        self.EPLOT = self.E

        '''EXPORTED WAVEFORM'''
        self.indexWF = self.index
        self.tWF = self.t        
        self.EWF = self.E

class CV(sweep):
    '''Waveform for cyclic voltammetry \n
    \n
    Eini - start potential \n
    Eupp - upper vertex potential \n
    Elow - lower vertex potential \n
    dE   - step size (in this case, number of data points) \n
    sr   - scan rate \n
    ns   - number of scans'''

    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        self.subtype = 'CV'
        self.detailed = False

        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(2 * self.ns * self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.round(np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32), 3))
                self.E = np.append(self.E, np.round(np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32), 3))

    
        '''STARTING FROM UPPER VERTEX POTENTIAL'''
        if self.Eini == self.Eupp:     
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(2 * self.ns * self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.round(np.linspace(self.Eini + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32), 3))
                self.E = np.append(self.E, np.round(np.linspace(self.Elow - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32), 3))


        '''STARTING IN BETWEEN VERTEX POTENTIALS'''
        if self.Elow < self.Eini < self.Eupp:        
            self.uppwindow = round(self.Eupp - self.Eini, 3)
            self.window = round(self.Eupp - self.Elow, 3)
            self.lowwindow = round(self.Eini - self.Elow, 3)
            self.uppdp = round(np.abs(self.uppwindow / self.dE))
            self.dp = round(np.abs(self.window / self.dE))
            self.lowdp = round(np.abs(self.lowwindow / self.dE))
            self.tmax = round(self.ns * (self.uppwindow + self.window + self.lowwindow) / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL WITH POSITIVE SCAN DIRECTION'''
            if self.dE > 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.round(np.linspace(self.Eini + self.dE, self.Eupp, self.uppdp, endpoint = True, dtype = np.float32), 3))
                    self.E = np.append(self.E, np.round(np.linspace(self.Eupp - self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32), 3))
                    self.E = np.append(self.E, np.round(np.linspace(self.Elow + self.dE, self.Eini, self.lowdp, endpoint = True, dtype = np.float32), 3))

            '''POTENTIAL WITH NEGATIVE SCAN DIRECTION'''
            if self.dE < 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.round(np.linspace(self.Eini - self.dE, self.Elow, self.lowdp, endpoint = True, dtype = np.float32), 3))
                    self.E = np.append(self.E, np.round(np.linspace(self.Elow + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32), 3))
                    self.E = np.append(self.E, np.round(np.linspace(self.Eupp + self.dE, self.Eini, self.uppdp, endpoint = True, dtype = np.float32), 3))
            
        '''PLOTTING WAVEFORM'''
        self.tPLOT = self.t
        self.EPLOT = self.E

        '''EXPORTED WAVEFORM'''
        self.indexWF = self.index
        self.tWF = self.t        
        self.EWF = self.E


"""STEP CLASSES"""
class CA(step):
    '''Chronoamperommetry (CA) waveform\n   
    \n
    dE  =   step potential \n
    dt  =   step duration \n
    st  =   sampling time'''

    def __init__(self, dE, dt, st):
        super().__init__(dE, dt, st)

        self.subtype = 'CA'
        self.detailed = True

        '''SINGLE STEP CHRONOAMPEROMMETRY'''
        if self.multiple is False:
            '''INDEX'''
            self.index = 0
        
            '''TIME'''
            self.t = np.array([0])
            self.t = np.append(self.t, self.dt)
            
            '''POTENTIAL'''
            self.E = self.dE
        
        '''MULTIPLE STEP CHRONOAMPEROMMETRY'''
        if self.multiple is True:
            '''INDEX'''
            self.index = np.arange(0, len(self.dE), 1, dtype = np.int32)
     
            '''TIME'''
            self.t = np.array([0])
            for ix in np.asarray(self.dt):
                self.t = np.append(self.t, self.t[-1] + ix)
            
            '''POTENTIAL'''
            self.E = np.asarray(self.dE)
        
        
        '''PLOTTING WAVEFORM'''
        self.tPLOT = np.array([])
        for ix in range(0, self.t.size):
            try:
                self.tPLOT = np.append(self.tPLOT, np.arange(self.t[ix], self.t[ix + 1], self.st))
            except:
                pass

        '''EXPORTED WAVEFORM'''
        self.indexWF = np.arange(0, self.tPLOT.size, 1)
        self.tWF = self.tPLOT        
        
        if self.multiple == False:
            self.EWF = self.dE * np.ones(round(self.dt/self.st))
        if self.multiple == True:
            self.EWF = np.array([])
            for ix in range(0, self.E.size):
                self.EWF = np.append(self.EWF, self.E[ix] * np.ones(round(self.dt[ix]/self.st)))


"""PULSE CLASSES"""
class DPV(pulse):
    '''Waveform for differential pulse voltammetry \n
    \n
    Eini - start potential \n
    Efin - end potential \n
    dEs  - step potential \n
    dEp  - pulse potential \n
    pt   - pulse time \n
    rt   - rest time \n
    st   - sampling time'''
    def __init__(self, Eini, Efin, dEs, dEp, pt, rt, st, detailed, sampled, alpha):
        super().__init__(Eini, Efin, dEs, dEp, pt, rt, st, detailed, sampled, alpha)

        self.subtype = 'DPV'

        if self.pt == self.rt:
            print('\n' + 'Wouldn\'t you rather be using square wave voltammetry?' + '\n')
            sys.exit()       

        '''TIME VARIABLES'''
        self.window = round(self.Efin - self.Eini, 3)
        self.dp = round(np.abs(self.window/self.dEs))
        self.dt = round(self.pt + self.rt, 6)
        self.tmax = round(self.dp * self.dt, 6)

        '''SAMPLING VARIABLES'''
        self.sp = round(self.dt / self.st)
        self.pp = round(self.pt / self.st)
        self.rp = round(self.rt / self.st)
        
        '''INDEX'''
        self.index = np.arange(0, round((2 * (self.tmax + self.dt)  / self.dt) + 1, 9), 1, dtype = np.int32)
        
        '''TIME'''
        self.t = np.array([0])
        for ix in range(0, self.dp + 1):
            self.t = np.append(self.t, round(ix * self.dt + self.pt, 9))
            self.t = np.append(self.t, round(ix * self.dt + self.dt, 9))
        
        '''POTENTIAL'''
        self.E = np.array([0])
        for ix in range(0, self.dp + 1):
            self.E = np.append(self.E, round(ix * self.dEs + self.dEp, 9))
            self.E = np.append(self.E, round(ix * self.dEs - self.dEp, 9))


        '''PLOTTING WAVEFORM'''
        if self.detailed == False:
            self.tPLOT = np.array([])
            self.EPLOT = np.array([])
            for ix in self.index:
                if ix % 2 == 0:
                    self.tPLOT = np.append(self.tPLOT, self.t[ix])
                    self.EPLOT = np.append(self.EPLOT, self.E[ix])
                else:
                    pass

        if self.detailed == True:
            self.tPLOT = np.array([])
            for ix in range (1, self.dp + 1):
                self.tPLOT = np.append(self.tPLOT, np.linspace((ix - 1) * self.dt, ix * self.dt, self.sp, endpoint = False))
        
            self.EPLOT = np.array([])
            for ix in range(1, self.dp + 1):
                self.EPLOT = np.append(self.EPLOT, np.linspace(self.Eini + ((ix - 1) * self.dEs), self.Eini + (ix * self.dEs), self.sp, endpoint = False))


        '''EXPORTED WAVEFORM'''
        if self.detailed == False:
            self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp        
            self.EWF = np.array([0])
            for ix in range(1, self.E.size - 1, 2):
                self.EWF = np.append(self.EWF, np.ones((self.pp)) * self.E[ix])
                self.EWF = np.append(self.EWF, np.ones((self.rp)) * self.E[ix + 1])
        
        if self.detailed == True:
            self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp        
            self.EWF = np.array([0])
            for ix in range(1, self.E.size - 1, 2):
                self.EWF = np.append(self.EWF, np.ones((self.pp)) * self.E[ix])
                self.EWF = np.append(self.EWF, np.ones((self.rp)) * self.E[ix + 1])

class SWV(pulse):
    '''Waveform for square wave voltammetry \n   
    \n
    Eini - start potential \n
    Eupp - upper vertex potential \n
    Elow - lower vertex potential \n
    dE   - step size \n
    sr   - scan rate \n
    ns   - number of scans \n
    st   - sampling time'''
    def __init__(self, Eini, Efin, dEs, dEp, pt, rt, st, detailed):
        super().__init__(Eini, Efin, dEs, dEp, pt, rt, st, detailed)

        self.subtype = 'SWV'

        if self.pt != self.rt:
            print('\n' + 'Wouldn\'t you rather be using differential pulse voltammetry?' + '\n')
            sys.exit()       

        '''TIME VARIABLES'''
        self.window = round(self.Efin - self.Eini, 3)
        self.dp = round(np.abs(self.window/self.dEs))
        self.dt = round(self.pt + self.rt, 6)
        self.tmax = round(self.dp * self.dt, 6)

        '''SAMPLING VARIABLES'''
        self.sp = round(self.dt / self.st)
        self.pp = round(self.pt / self.st)
        self.rp = round(self.rt / self.st)
        
        '''INDEX'''
        self.index = np.arange(0, round((2 * (self.tmax + self.dt)  / self.dt) + 1, 9), 1, dtype = np.int32)
        
        '''TIME'''
        self.t = np.array([0])
        for ix in range(0, self.dp + 1):
            self.t = np.append(self.t, round(ix * self.dt + self.pt, 9))
            self.t = np.append(self.t, round(ix * self.dt + self.dt, 9))
        
        '''POTENTIAL'''
        self.E = np.array([0])
        for ix in range(0, self.dp + 1):
            self.E = np.append(self.E, round(ix * self.dEs + self.dEp, 9))
            self.E = np.append(self.E, round(ix * self.dEs - self.dEp, 9))


        '''PLOTTING WAVEFORM'''
        if self.detailed == False:
            self.tPLOT = np.array([])
            self.EPLOT = np.array([])
            for ix in self.index:
                if ix % 2 == 0:
                    self.tPLOT = np.append(self.tPLOT, self.t[ix])
                    self.EPLOT = np.append(self.EPLOT, self.E[ix])
                else:
                    pass

        if self.detailed == True:
            self.tPLOT = np.array([])
            for ix in range (1, self.dp + 1):
                self.tPLOT = np.append(self.tPLOT, np.linspace((ix - 1) * self.dt, ix * self.dt, self.sp, endpoint = False))
        
            self.EPLOT = np.array([])
            for ix in range(1, self.dp + 1):
                self.EPLOT = np.append(self.EPLOT, np.linspace(self.Eini + ((ix - 1) * self.dEs), self.Eini + (ix * self.dEs), self.sp, endpoint = False))


        '''EXPORTED WAVEFORM'''
        if self.detailed == False:
            self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp        
            self.EWF = np.array([0])
            for ix in range(1, self.E.size - 1, 2):
                self.EWF = np.append(self.EWF, np.ones((self.pp)) * self.E[ix])
                self.EWF = np.append(self.EWF, np.ones((self.rp)) * self.E[ix + 1])
        
        if self.detailed == True:
            self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp        
            self.EWF = np.array([0])
            for ix in range(1, self.E.size - 1, 2):
                self.EWF = np.append(self.EWF, np.ones((self.pp)) * self.E[ix])
                self.EWF = np.append(self.EWF, np.ones((self.rp)) * self.E[ix + 1])

class NPV(pulse):
    '''Waveform for normal pulse voltammetry\n   
    \n
    Eini - start potential \n
    Efin - upper vertex potential \n
    dEs   - step size \n   
    dEp   - step size \n
    sr   - scan rate \n
    ns   - number of scans \n
    st   - sampling time'''
    def __init__(self, Eini, Efin, dEs, dEp, pt, rt, st, detailed, sampled, alpha):
        super().__init__(Eini, Efin, dEs, dEp, pt, rt, st, detailed, sampled, alpha)

        self.subtype = 'NPV'

        '''TIME VARIABLES'''
        self.window = round(self.Efin - self.Eini, 3)
        self.dp = round(np.abs(self.window/self.dEs))
        self.dt = round(self.pt + self.rt, 6)
        self.tmax = round(self.dp * self.dt, 6)

        '''SAMPLING VARIABLES'''
        self.sp = round(self.dt / self.st)
        self.pp = round(self.pt / self.st)
        self.rp = round(self.rt / self.st)
        
        '''INDEX'''
        self.index = np.arange(0, round((2 * (self.tmax + self.dt)  / self.dt) + 1, 9), 1, dtype = np.int32)
        
        '''TIME'''
        self.t = np.array([0])
        for ix in range(0, self.dp + 1):
            self.t = np.append(self.t, round(ix * self.dt + self.pt, 9))
            self.t = np.append(self.t, round(ix * self.dt + self.dt, 9))
        
        '''POTENTIAL'''
        self.E = np.array([self.Eini])
        for ix in range(1, self.dp + 1):
            if ix != self.dp:
                self.E = np.append(self.E, square(2 * np.pi * (1/self.dt) * self.t[ix:ix + 1], duty = self.pt/self.dt) *(0.5 * ix * self.dEs) + (0.5 * ix * self.dEs))
                self.E = np.append(self.E, self.Eini)
            else:
                self.E = np.append(self.E, square(2 * np.pi * (1/self.dt) * self.t[ix - 1:self.dp], duty = self.pt/self.dt) *(0.5 * self.dp * self.dEs) + (0.5 * self.dp * self.dEs))
                self.E = np.append(self.E, self.Eini)

        '''PLOTTING WAVEFORM'''
        if self.detailed == False:
            self.tPLOT = np.array([])
            self.EPLOT = np.array([])
            for ix in self.index:
                if ix % 2 == 0:
                    self.tPLOT = np.append(self.tPLOT, self.t[ix])
                    self.EPLOT = np.append(self.EPLOT, self.E[ix])
                else:
                    pass

        if self.detailed == True:
            self.tPLOT = np.array([])
            for ix in range (1, self.dp + 1):
                self.tPLOT = np.append(self.tPLOT, np.linspace((ix - 1) * self.dt, ix * self.dt, self.sp, endpoint = False))
        
            self.EPLOT = np.array([])
            for ix in range(1, self.dp + 1):
                self.EPLOT = np.append(self.EPLOT, np.linspace(self.Eini + ((ix - 1) * self.dEs), self.Eini + (ix * self.dEs), self.sp, endpoint = False))


        '''EXPORTED WAVEFORM'''
        if self.detailed == False:
            self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp        
            self.EWF = np.array([0])
            for ix in range(1, self.E.size - 1, 2):
                self.EWF = np.append(self.EWF, np.ones((self.pp)) * self.E[ix])
                self.EWF = np.append(self.EWF, np.ones((self.rp)) * self.E[ix + 1])
        
        if self.detailed == True:
            self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp        
            self.EWF = np.array([0])
            for ix in range(1, self.E.size - 1, 2):
                self.EWF = np.append(self.EWF, np.ones((self.pp)) * self.E[ix])
                self.EWF = np.append(self.EWF, np.ones((self.rp)) * self.E[ix + 1])
                

"""HYBRID CLASSES"""
class CSV(hybrid):
    '''Waveform for cyclic staircase voltammetry \n   
    \n
    Eini - start potential \n
    Eupp - upper vertex potential \n
    Elow - lower vertex potential \n
    dE   - step size \n
    sr   - scan rate \n
    ns   - number of scans \n
    st   - sampling time'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns, st, detailed, sampled, alpha):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns, st, detailed, sampled, alpha)
        
        self.subtype = 'CSV'

        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(2 * self.ns * self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
            self.sp = round(self.dt / self.st)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))
            self.E = np.round(self.E, 6)

    
        '''STARTING FROM UPPER VERTEX POTENTIAL'''
        if self.Eini == self.Eupp:     
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(2 * self.ns * self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
            self.sp = round(self.dt / self.st)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Elow - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))
            self.E = np.round(self.E, 6)


        '''STARTING IN BETWEEN VERTEX POTENTIALS'''
        if self.Elow < self.Eini < self.Eupp:        
            self.uppwindow = round(self.Eupp - self.Eini, 3)
            self.window = round(self.Eupp - self.Elow, 3)
            self.lowwindow = round(self.Eini - self.Elow, 3)
            self.uppdp = round(np.abs(self.uppwindow / self.dE))
            self.dp = round(np.abs(self.window / self.dE))
            self.lowdp = round(np.abs(self.lowwindow / self.dE))
            self.tmax = round(self.ns * (self.uppwindow + self.window + self.lowwindow) / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
            self.sp = round(self.dt / self.st)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL WITH POSITIVE SCAN DIRECTION'''
            if self.dE > 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.uppdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.lowdp, endpoint = True, dtype = np.float32))

            '''POTENTIAL WITH NEGATIVE SCAN DIRECTION'''
            if self.dE < 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini - self.dE, self.Elow, self.lowdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Eupp + self.dE, self.Eini, self.uppdp, endpoint = True, dtype = np.float32))
            
            self.E = np.round(self.E, 6)
        

        '''PLOTTING WAVEFORM'''
        if self.detailed == False:
            self.tPLOT = self.t
            self.EPLOT = self.E

        if self.detailed == True:
            if self.sampled == True:
                self.tPLOT = self.t
            else:
                self.tPLOT = np.array([])
                for ix in range(0, self.t.size):
                    try:
                        self.tPLOT = np.append(self.tPLOT, np.linspace(self.t[ix], self.t[ix + 1], self.sp, endpoint = False))
                    except:
                        self.tPLOT = np.append(self.tPLOT, np.linspace(self.t[ix], self.t[ix] + self.dt, self.sp, endpoint = False))
                self.tPLOT = np.round(self.tPLOT, 9)

            if self.sampled == True:
                self.EPLOT = self.E
            else:
                self.EPLOT = np.array([])
                for iy in range(0, self.E.size):
                    try:
                        self.EPLOT = np.append(self.EPLOT, np.linspace(self.E[iy], self.E[iy + 1], self.sp, endpoint = False))
                    except:
                        self.EPLOT = np.append(self.EPLOT, np.linspace(self.E[iy], self.E[iy] + self.dE, self.sp, endpoint = False))
                self.EPLOT = np.round(self.EPLOT, 9)

    
        '''EXPORTED WAVEFORM'''
        if self.detailed == False:
            self.indexWF = np.arange(0, 2 * round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
            self.tWF = np.array(0)
            for ix in range(1, round(self.indexWF.size / 2)):
                self.tWF = np.append(self.tWF, np.ones(2) * self.indexWF[ix] * (self.dt))
            self.tWF = np.append(self.tWF, self.tWF[-1] + self.dt)
            self.EWF = np.array([])
            for ix in self.E:
                self.EWF = np.append(self.EWF, np.ones((2)) * ix)

        if self.detailed == True:
            self.indexWF = np.arange(0, self.sp * round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp
            self.EWF = np.array([])
            for ix in self.E:
                self.EWF = np.append(self.EWF, np.ones((self.sp)) * ix)

class AC(hybrid):
    '''Waveform for AC voltammetry \n   
    \n
    Eini - start potential \n
    Eupp - upper vertex potential \n
    Elow - lower vertex potential \n
    dE   - step size \n
    sr   - scan rate \n
    ns   - number of scans \n
    st   - sampling time'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns, st, detailed):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns, st, detailed)
        
        self.subtype = 'AC'
        
        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(2 * self.ns * self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
            self.sp = round(self.dt / self.st)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))
            self.E = np.round(self.E, 6)

    
        '''STARTING FROM UPPER VERTEX POTENTIAL'''
        if self.Eini == self.Eupp:     
            self.window = round(self.Eupp - self.Elow, 3)
            self.dp = round(np.abs(self.window / self.dE))
            self.tmax = round(2 * self.ns * self.window / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
            self.sp = round(self.dt / self.st)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Elow - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))
            self.E = np.round(self.E, 6)


        '''STARTING IN BETWEEN VERTEX POTENTIALS'''
        if self.Elow < self.Eini < self.Eupp:        
            self.uppwindow = round(self.Eupp - self.Eini, 3)
            self.window = round(self.Eupp - self.Elow, 3)
            self.lowwindow = round(self.Eini - self.Elow, 3)
            self.uppdp = round(np.abs(self.uppwindow / self.dE))
            self.dp = round(np.abs(self.window / self.dE))
            self.lowdp = round(np.abs(self.lowwindow / self.dE))
            self.tmax = round(self.ns * (self.uppwindow + self.window + self.lowwindow) / self.sr, 6)
            self.dt = round(np.abs(self.dE / self.sr), 9)
            self.sp = round(self.dt / self.st)

            '''INDEX'''
            self.index = np.arange(0, round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL WITH POSITIVE SCAN DIRECTION'''
            if self.dE > 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.uppdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.lowdp, endpoint = True, dtype = np.float32))

            '''POTENTIAL WITH NEGATIVE SCAN DIRECTION'''
            if self.dE < 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini - self.dE, self.Elow, self.lowdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Eupp + self.dE, self.Eini, self.uppdp, endpoint = True, dtype = np.float32))
            
            self.E = np.round(self.E, 6)
        

        '''PLOTTING WAVEFORM'''
        if self.detailed == False:
            self.tPLOT = self.t
            self.EPLOT = self.E

        if self.detailed == True:
            self.tPLOT = np.array([])
            for ix in range(0, self.t.size):
                try:
                    self.tPLOT = np.append(self.tPLOT, np.linspace(self.t[ix], self.t[ix + 1], self.sp, endpoint = False))
                except:
                    self.tPLOT = np.append(self.tPLOT, np.linspace(self.t[ix], self.t[ix] + self.dt, self.sp, endpoint = False))
            self.tPLOT = np.round(self.tPLOT, 9)

            self.EPLOT = np.array([])
            for iy in range(0, self.E.size):
                try:
                    self.EPLOT = np.append(self.EPLOT, np.linspace(self.E[iy], self.E[iy + 1], self.sp, endpoint = False))
                except:
                    self.EPLOT = np.append(self.EPLOT, np.linspace(self.E[iy], self.E[iy] + self.dE, self.sp, endpoint = False))
            self.EPLOT = np.round(self.EPLOT, 9)

    
        '''EXPORTED WAVEFORM'''
        if self.detailed == False:
            self.indexWF = np.arange(0, 2 * round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
            self.tWF = np.array(0)
            for ix in range(1, round(self.indexWF.size / 2)):
                self.tWF = np.append(self.tWF, np.ones(2) * self.indexWF[ix] * (self.dt))
            pass
            self.tWF = np.append(self.tWF, self.tWF[-1] + self.dt)
            self.EWF = np.array([])
            for ix in self.E:
                self.EWF = np.append(self.EWF, np.ones((2)) * ix)

        if self.detailed == True:
            self.indexWF = np.arange(0, self.sp * round((self.tmax + self.dt) / self.dt, 9), 1, dtype = np.int32)
            self.tWF = (self.indexWF * self.dt) / self.sp
            self.EWF = np.array([])
            for ix in self.E:
                self.EWF = np.append(self.EWF, np.ones((self.sp)) * ix)


"""IMPEDANCE CLASSES"""
class EIS(impedance):
    pass







"""
===================================================================================================
GENERATING WAVEFORMS FROM MAIN
===================================================================================================
"""


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
    '''Sweeps'''
    #wf = LSV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1)
    #wf = CV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1)
    
    '''STEPS'''
    #wf = CA(dE = [0.5], dt = [1], st = 0.001)
    
    '''PULSES'''
    wf = DPV(Eini = 0, Efin = 0.25, dEs = 0.005, dEp = 0.05, pt = 0.2, rt = 0.4, st = 0.001, detailed = False, sampled = True, alpha = 0.25)
    #wf = SWV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.02, pt = 0.1, rt = 0.1, st = 0.001, detailed = True, sampled = True, alpha = 0.25)
    #wf = NPV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.02, pt = 0.05, rt = 0.15, st = 0.001, detailed = True, sampled = True, alpha = 0.25)
    
    '''HYBRID'''
    #wf = CSV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.005, sr = 0.1, ns = 1, st = 0.0001, detailed = False, sampled = True, alpha = 0.25)
    #wf = AC(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1, st = 0.001, detailed = True, sampled = True, alpha = 0.25)
    
    '''5. DEFINE THE END TIME'''
    end = time.time()
    print(f'The waveform took {end-start} seconds to generate')

    '''6. SAVE THE DATA'''
    filepath = f'{cwd}/data/{time.strftime("%Y-%m-%d %H-%M-%S")} {wf.subtype} waveform.txt'
    with open(filepath, 'w') as file:
        for ix, iy, iz in wf.output():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')