import sys
import os
from errno import EEXIST
import numpy as np
from scipy.signal import square



"""PARENT CLASSES"""
class sweep:
    '''Parent class for all sweep type waveforms'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):
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
        zipped = zip(self.index, self.t, self.E)
        return zipped
            
class step:
    '''Parent class for all step type waveforms'''
    def __init__(self, dE, dt, sp):
        self.dE = dE              # Potential for single step chronoamperommetry
        self.dt = dt            # Step period
        self.sp = sp            # Sample points in a step
        
        '''DATATYPE ERRORS'''
        if isinstance(self.dE, (list)) is True and isinstance(self.dt, (list)) is True and isinstance(self.sp, (list)) is True:
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
            for iz in self.sp:
                if isinstance(iz, (int)) is False:
                    print('\n' + 'An invalid datatype was used for the number of data points in a step. Enter an integer value.' + '\n')
                    sys.exit()       

        elif isinstance(self.dE, (list)) is False and isinstance(self.dt, (list)) is False:
            self.multiple = False
            if isinstance(self.dE, (float, int)) is False:
                print('\n' + 'An invalid datatype was used for the potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
                sys.exit()
            if isinstance(self.dt, (float, int)) is False:
                print('\n' + 'An invalid datatype was used for the step time. Enter either a float or an integer value corresponding to a time in s.' + '\n')
                sys.exit()        
            if isinstance(self.sp, (int)) is False:
                print('\n' + 'An invalid datatype was used for the number of data points in a step. Enter an integer value.' + '\n')
                sys.exit()
        else:
            print('\n' + 'When entering multiple steps, make sure to enter a list of potentials, a list of times, and a list of data points.' + '\n')
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
            for iz in self.sp:
                if iz <= 0:
                    print('\n' + 'All steps must have a number of data points with a positive non-zero value' + '\n')
                    sys.exit() 
        else:      
            if self.dt <= 0:
                print('\n' + 'Step time must be a positive non-zero value' + '\n')
                sys.exit() 
            if self.sp <= 0:
                print('\n' + 'Number of data points in a step must be a positive non-zero value' + '\n')
                sys.exit()        
    
    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.index, self.t, self.E)
        return zipped
            
class hybrid:
    '''Parent class for all waveforms composed of both steps and sweeps'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, sp, ns):
        self.Eini = Eini        # Start potential for sweeping step techniques
        self.Eupp = Eupp        # Upper vertex potential for sweeping step techniques
        self.Elow = Elow        # Lower vertex potential for sweeping step techniques
        self.dE = dE            # Step size for sweeping step techniques
        self.sr = sr            # Scan rate for sweeping step techniques
        self.sp = sp            # Sample points in a step
        self.ns = ns            # Number of scans for cyclic staircase voltammetry

        '''DATATYPE ERRORS'''
        if isinstance(self.E, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.E1, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the 1st potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()        
        if isinstance(self.E2, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the 2nd potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
            sys.exit()
        if isinstance(self.dt, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the step time. Enter either a float or an integer value corresponding to a time in s.' + '\n')
            sys.exit()        
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
        if isinstance(self.ns, (int)) is False:
            print('\n' + 'An invalid datatype was used for the number of scans. Enter an integer value corresponding to the scan rate in V/s.' + '\n')
            sys.exit() 

        '''DATA VALUE ERRORS'''
        if self.E1 == self.E2:
            print('\n' + 'First and second potentials must be different values' + '\n')
            sys.exit()  
        if self.dt <= 0:
            print('\n' + 'Step time must be a positive non-zero value' + '\n')
            sys.exit()    
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
        if self.Eini == self.Elow and self.dE < 0:
            print('\n' + 'Step potential must be a positive value for a positive scan direction' + '\n')
            sys.exit()
        if self.Eini == self.Eupp and self.dE > 0:
            print('\n' + 'Step potential must be a negative value for a negative scan direction' + '\n')
            sys.exit()
        if self.sp <= 0:
            print('\n' + 'Number of data points in a step must be a positive non-zero value' + '\n')
            sys.exit()        
        if self.ns <=0:
            print('\n' + 'Number of scans must be a positive non-zero value' + '\n')
            sys.exit()
    
    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.index, self.t, self.E)
        return zipped
    
"""SWEEP CLASSES"""
class LSV(sweep):
    '''Waveform for linear sweep voltammetry'''
    def __init__(self,Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        if self.Elow < self.Eini < self.Eupp:
            print('\n' + 'Initial potential should be equal to either upper vertex or lower vertex potential in LSV' + '\n')
            sys.exit()
        
        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                 
            self.window = self.Eupp - self.Elow
            self.dp = int(self.window / self.dE)
            self.tmax = self.window / self.sr
            self.dt = self.dE / self.sr
        
            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.linspace(self.Eini, self.Eupp, self.dp + 1, endpoint = True, dtype = np.float32)
        
        '''STARTING FROM UPPER VERTEX POTENTIAL''' 
        if self.Eini == self.Eupp:                 
            self.window = self.Eupp - self.Elow
            self.dp = int(self.window / np.abs(self.dE))
            self.tmax = self.window / self.sr
            self.dt = np.abs(self.dE) / self.sr
        
            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.linspace(self.Eini, self.Elow, self.dp + 1, endpoint = True, dtype = np.float32)

class CV(sweep):
    '''Waveform for cyclic voltammetry'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                
            self.window = self.Eupp - self.Elow
            self.dp = int(self.window / self.dE)
            self.tmax = (2 * self.ns * self.window) / self.sr
            self.dt = self.dE / self.sr

            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))

        '''STARTING FROM UPPER VERTEX POTENTIAL'''
        if self.Eini == self.Eupp:     
            self.window = self.Eupp - self.Elow
            self.dp = int(self.window / np.abs(self.dE))
            self.tmax = (2 * self.ns * self.window) / self.sr
            self.dt = np.abs(self.dE) / self.sr

            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Elow - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))


        '''STARTING IN BETWEEN VERTEX POTENTIALS'''
        if self.Elow < self.Eini < self.Eupp:        
            self.uppwindow = self.Eupp - self.Eini
            self.window = self.Eupp - self.Elow
            self.lowwindow = self.Eini - self.Elow
            self.uppdp = int(self.uppwindow / self.dE)
            self.dp = int(self.window / self.dE)
            self.lowdp = int(self.lowwindow / self.dE)
            self.tmax = self.ns * (self.uppwindow + self.window + self.lowwindow) / self.sr
            self.dt = self.dE / self.sr

            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL WITH POSITIVE SCAN DIRECTION'''
            if dE > 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.uppdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.lowdp, endpoint = True, dtype = np.float32))

            '''POTENTIAL WITH NEGATIVE SCAN DIRECTION'''
            if dE < 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini - self.dE, self.Elow, self.lowdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.uppdp, endpoint = True, dtype = np.float32))



"""STEP CLASSES"""
class CA(step):
    '''Waveform for chronoamperommetry'''
    def __init__(self, dE, dt, sp):
        super().__init__(dE, dt, sp)

        '''SINGLE STEP CHRONOAMPEROMMETRY'''
        if self.multiple is False:
            '''INDEX'''
            self.index = np.arange(0, (self.sp + 1), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = (self.index/self.sp) * self.dt
            
            '''POTENTIAL'''
            self.E = np.ones((self.sp + 1)) * self.dE
        
        '''MULTIPLE STEP CHRONOAMPEROMMETRY'''
        if self.multiple is True:
            '''INDEX'''
            self.index = np.arange(0, (sum(self.sp) + 1), 1, dtype = np.int32)

            #1000 - 5, 400 - 2, 500- 5
            '''TIME''' # This doesn't work
            self.t = np.array([])
            for ix in range(0, len(self.dt)):
                if ix == 0:
                    self.t = np.append(self.t, (self.index[0:self.sp[ix] + 1])/(self.sp[ix]) * self.dt[ix])
                else:
                    self.t = np.append(self.t, (self.index[sum(self.sp[:ix]):sum(self.sp[:ix + 1]) + 1] - sum(self.sp[:ix]))/(sum(self.sp[:ix + 1])) * self.dt[ix])
            
            '''POTENTIAL'''
            self.E = np.array(self.dE[0])
            for iy in range(0, len(self.dE)):
                self.E = np.append(self.E, np.ones(self.sp[iy]) * self.dE[ix])

class CSV(hybrid):
    '''Waveform for cyclic staircase voltammetry'''
    def __init__(self, E, E1, E2, dt, Eini, Eupp, Elow, dE, sp, ns):
        super().__init__(E, E1, E2, dt, Eini, Eupp, Elow, dE, sp, ns)
        
        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                
            self.window = self.Eupp - self.Elow
            self.dp = int(self.window / self.dE)
            self.tmax = (2 * self.ns * self.window) / self.sr

            '''INDEX'''
            self.index = np.arange(0, self.sp * ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))

        '''STARTING FROM UPPER VERTEX POTENTIAL'''
        if self.Eini == self.Eupp:     
            self.window = self.Eupp - self.Elow
            self.dp = int(self.window / np.abs(self.dE))
            self.tmax = (2 * self.ns * self.window) / self.sr
            self.dt = np.abs(self.dE) / self.sr

            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Elow - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))


        '''STARTING IN BETWEEN VERTEX POTENTIALS'''
        if self.Elow < self.Eini < self.Eupp:        
            self.uppwindow = self.Eupp - self.Eini
            self.window = self.Eupp - self.Elow
            self.lowwindow = self.Eini - self.Elow
            self.uppdp = int(self.uppwindow / self.dE)
            self.dp = int(self.window / self.dE)
            self.lowdp = int(self.lowwindow / self.dE)
            self.tmax = self.ns * (self.uppwindow + self.window + self.lowwindow) / self.sr
            self.dt = self.dE / self.sr

            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL WITH POSITIVE SCAN DIRECTION'''
            if dE > 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.uppdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.lowdp, endpoint = True, dtype = np.float32))

            '''POTENTIAL WITH NEGATIVE SCAN DIRECTION'''
            if dE < 0:
                self.E = np.array([self.Eini])
                for ix in range(0, self.ns):
                    self.E = np.append(self.E, np.linspace(self.Eini - self.dE, self.Elow, self.lowdp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.uppdp, endpoint = True, dtype = np.float32))

class NPV(hybrid):
    '''Waveform for normal pulse voltammetry'''
    def __init__(self):
        pass

class DPV(hybrid):
    '''Waveform for differential pulse voltammetry'''
    def __init__(self):
        pass

class RPV(hybrid):
    '''Waveform for reverse pulse voltammetry'''
    def __init__(self):
        pass

class SWV(hybrid):
    """Waveform for square wave voltammetry"""
    
    def __init__(self, Eini = 0, Efin = 1, dEs = 0.002, dEp = 0.05, f = 25, sp = 1000):
        """Makes an instance of the SWV class and produces sweep, step, and square wave voltammetry waveforms"""
        
        '''USER-DEFINED VARIABLES'''
        self.Eini = Eini                # Initial potential
        self.Efin = Efin                # Final potential
        self.dEs = dEs                  # Step potential
        self.dEp = dEp                  # Pulse amplitude
        self.f = f                      # Pulse frequency
        self.sp = sp                    # Data points in a step

        '''DERIVED VARIABLES'''
        self.window = self.Efin - self.Eini             # Potential window             
        self.dp = int(self.window / self.dEs)           # Number of steps 
        self.step_period = 1 / self.f                   # Step interval

        '''TIME INDEX'''
        self.index = np.arange(0, int(1000 * self.step_period * (self.dp + 1)), int(self.step_period * 1000))
        
        '''SWEEP WAVEFORM'''
        self.sweep_time = self.index / 1000
        self.sweep_waveform = np.linspace(self.Eini, self.Efin, self.dp + 1, endpoint = True)

        '''STEP WAVEFORM'''
        self.step_time = np.array([])
        for ix in range(0, self.index.size):
            try:
                self.step_time = np.append(self.step_time, np.arange(self.index[ix], self.index[ix + 1], (1000 * self.step_period / self.sp))/1000)
            except:
                self.step_time = np.append(self.step_time, np.arange(self.index[ix], self.index[ix] + self.index[1], (1000 * self.step_period / self.sp))/1000)

        self.step_waveform = np.array([])
        for iy in range(0, self.sweep_waveform.size):
            self.step_waveform = np.append(self.step_waveform, np.ones((self.sp)) * self.sweep_waveform[iy])

        '''SQUARE WAVE VOLTAMMETRY WAVEFORM'''
        self.square_waveform = square(2 * np.pi * self.f * self.step_time, duty = 0.5) * self.dEp
        self.combined = self.step_waveform + self.square_waveform

        '''INPUT POTENTIALS FOR SIMULATION'''
        self.E = np.array([])
        for ix in range(self.sweep_waveform.size):
            self.E = np.append(self.E, self.sweep_waveform[ix] + self.dEp)
            self.E = np.append(self.E, self.sweep_waveform[ix] - self.dEp)


if __name__ == '__main__':
    
        
    cwd = os.getcwd()

    try:
        os.makedirs(cwd + '/data')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/data'):
            pass
        else: 
            raise
    filepath = cwd + '/data/' + 'waveform.txt'

    wf = CA(dE = [0.5, 0.7, -0.2], dt = [5, 2, 5], sp = [1000, 400, 500])

    with open(filepath, 'w') as file:
        for ix, iy, iz in wf.output():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')