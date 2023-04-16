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
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped
            
class step:
    '''Parent class for all step type waveforms'''
    def __init__(self, dE, dt):
        
        self.type = 'step'
        
        self.dE = dE              # Potential for single step chronoamperommetry
        self.dt = dt            # Step period
        
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

        elif isinstance(self.dE, (list)) is False and isinstance(self.dt, (list)) is False:
            self.multiple = False
            if isinstance(self.dE, (float, int)) is False:
                print('\n' + 'An invalid datatype was used for the potential. Enter either a float or an integer value corresponding to a potential in V.' + '\n')
                sys.exit()
            if isinstance(self.dt, (float, int)) is False:
                print('\n' + 'An invalid datatype was used for the step time. Enter either a float or an integer value corresponding to a time in s.' + '\n')
                sys.exit()        

        else:
            print('\n' + 'When entering multiple steps, make sure to enter a separate lists of potentials and times' + '\n')
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
    
    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped

class pulse:
    '''Parent class for all pulse type waveforms'''
    def __init__(self, Eini, Efin, dEs, dEp, pt, rt, st, detailed):

        self.type = 'pulse'

        self.Eini = Eini        # Start potential for pulsed techniques
        self.Efin = Efin        # End potential for pulsed techniques
        self.dEs = dEs          # Step size for pulsed techniques
        self.dEp = dEp          # Pulse size for pulsed techniques
        self.pt = pt            # Pulse period
        self.rt = rt            # Rest period
        self.st = st            # Sampling period for more detailed waveforms

        self.detailed = detailed
               
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
        
                  

        '''DATA VALUE ERRORS'''
        if self.Eini == self.Efin:
            print('\n' + 'Start and end potentials must be different values' + '\n')
            sys.exit()  
        if self.dEs == 0:
            print('\n' + 'Step potential must be a non-zero value' + '\n')
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


    def output(self):
        '''Function that returns the waveform for checking or data processing purposes'''
        zipped = zip(self.indexWF, self.tWF, self.EWF)
        return zipped
    
class hybrid:
    '''Parent class for all waveforms composed of both steps and sweeps'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns, st, detailed):

        self.type = 'hybrid'

        self.Eini = Eini        # Start potential for sweeping step techniques
        self.Eupp = Eupp        # Upper vertex potential for sweeping step techniques
        self.Elow = Elow        # Lower vertex potential for sweeping step techniques
        self.dE = dE            # Step size for sweeping step techniques
        self.sr = sr            # Scan rate for sweeping step techniques
        self.ns = ns            # Number of scans for cyclic staircase voltammetry
        self.st = st            # Sampling time

        self.detailed = detailed

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
        if self.st <= 0:
            print('\n' + 'Sampling time must be a positive non-zero value' + '\n')
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
    '''Waveform for linear sweep voltammetry'''
    def __init__(self,Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        self.subtype = 'LSV'

        if self.Elow < self.Eini < self.Eupp:
            print('\n' + 'Initial potential should be equal to either upper vertex or lower vertex potential in LSV' + '\n')
            sys.exit()
        
        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                 
            self.window = np.abs(self.Eupp - self.Elow)
            self.dp = int(self.window / np.abs(self.dE))
            self.tmax = self.window / self.sr
            self.dt = np.abs(self.dE) / self.sr
        
            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.linspace(self.Eini, self.Eupp, self.dp + 1, endpoint = True, dtype = np.float32)
        
        '''STARTING FROM UPPER VERTEX POTENTIAL''' 
        if self.Eini == self.Eupp:                 
            self.window = np.abs(self.Eupp - self.Elow)
            self.dp = int(self.window / np.abs(self.dE))
            self.tmax = self.window / self.sr
            self.dt = np.abs(self.dE) / self.sr
        
            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
            '''TIME'''
            self.t = self.index * self.dt
            
            '''POTENTIAL'''
            self.E = np.linspace(self.Eini, self.Elow, self.dp + 1, endpoint = True, dtype = np.float32)
        
        '''PLOTTING WAVEFORM'''
        self.tPLOT = self.t
        self.EPLOT = self.E

        '''EXPORTED WAVEFORM'''
        self.indexWF = self.index
        self.tWF = self.t        
        self.EWF = self.E

class CV(sweep):
    '''Waveform for cyclic voltammetry'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        self.subtype = 'CV'

        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                
            self.window = np.abs(self.Eupp - self.Elow)
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
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))

    
        '''STARTING FROM UPPER VERTEX POTENTIAL'''
        if self.Eini == self.Eupp:     
            self.window = np.abs(self.Eupp - self.Elow)
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
            self.uppwindow = np.abs(self.Eupp - self.Eini)
            self.window = np.abs(self.Eupp - self.Elow)
            self.lowwindow = np.abs(self.Eini - self.Elow)
            self.uppdp = int(self.uppwindow / np.abs(self.dE))
            self.dp = int(self.window / np.abs(self.dE))
            self.lowdp = int(self.lowwindow / np.abs(self.dE))
            self.tmax = self.ns * (self.uppwindow + self.window + self.lowwindow) / self.sr
            self.dt = np.abs(self.dE) / self.sr

            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
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
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.uppdp, endpoint = True, dtype = np.float32))
            
        '''PLOTTING WAVEFORM'''
        self.tPLOT = self.t
        self.EPLOT = self.E

        '''EXPORTED WAVEFORM'''
        self.indexWF = self.index
        self.tWF = self.t        
        self.EWF = self.E


"""STEP CLASSES"""
class CA(step):
    '''Waveform for chronoamperommetry'''
    def __init__(self, dE, dt):
        super().__init__(dE, dt)

        self.subtype = 'CA'

        '''SINGLE STEP CHRONOAMPEROMMETRY'''
        if self.multiple is False:
            '''INDEX'''
            self.index = 0
        
            '''TIME'''
            self.t = self.dt
            
            '''POTENTIAL'''
            self.E = self.dE
        
        '''MULTIPLE STEP CHRONOAMPEROMMETRY'''
        if self.multiple is True:
            '''INDEX'''
            self.index = np.arange(0, len(self.dE), 1, dtype = np.int32)
     
            '''TIME'''
            self.t = self.dt.toarray()
            
            '''POTENTIAL'''
            self.E = self.dE.toarray()
        
        
        '''PLOTTING WAVEFORM'''
        self.tPLOT = self.t
        self.EPLOT = self.E

        '''EXPORTED WAVEFORM'''
        self.indexWF = self.index
        self.tWF = self.t        
        self.EWF = self.E


"""PULSE CLASSES"""
class DPV(pulse):
    '''Waveform for differential pulse voltammetry'''
    def __init__(self, Eini, Efin, dEs, dEp, pt, rt, st, detailed):
        super().__init__(Eini, Efin, dEs, dEp, pt, rt, st, detailed)

        self.subtype = 'DPV'

        if 2 * self.pt == self.rt:
            print('\n' + 'Wouldn\'t you rather be using square wave voltammetry?' + '\n')
            sys.exit()       

        '''TIME VARIABLES'''
        self.window = np.abs(self.Efin - self.Eini)
        self.dp = int(self.window / np.abs(self.dEs))
        self.dt = self.pt + self.rt
        self.tmax = self.dp * self.dt

        '''SAMPLING VARIABLES'''
        self.sp = int(self.dt / self.st)
        self.pp = int(self.pt / self.st)
        self.rp = int(self.rt / self.st)
        
        '''INDEX'''
        self.index = np.arange(0, (2 * (self.tmax + self.dt)  / self.dt) + 1, 1, dtype = np.int32)
        
        '''TIME'''
        self.t = np.array([0])
        for ix in range(0, self.dp + 1):
            self.t = np.append(self.t, ix * self.dt + self.pt)
            self.t = np.append(self.t, ix * self.dt + self.dt)
        
        '''POTENTIAL'''
        self.E = np.array([0])
        for ix in range(0, self.dp + 1):
            self.E = np.append(self.E, ix * self.dEs + self.dEp )
            self.E = np.append(self.E, ix * self.dEs - self.dEp)

        '''PLOTTING WAVEFORM'''
        self.tPLOT = np.array([0])
        for ix in range (1, self.dp + 1):
            self.tPLOT = np.append(self.tPLOT, ix * self.dt)
        
        self.EPLOT = np.array([0])
        for ix in range(1, self.dp + 1):
            self.EPLOT = np.append(self.EPLOT, ix * self.dEs)

        '''EXPORTED WAVEFORM'''
        self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
        self.tWF = (self.indexWF * self.dt) / self.sp        
        self.EWF = np.array([0])
        for ix in range(1, self.E.size - 1, 2):
            self.EWF = np.append(self.EWF, np.ones((self.pp)) * self.E[ix])
            self.EWF = np.append(self.EWF, np.ones((self.rp)) * self.E[ix + 1])

class SWV(pulse):
    """Waveform for square wave voltammetry"""
    def __init__(self, Eini, Efin, dEs, dEp, pt, rt, st, detailed):
        super().__init__(Eini, Efin, dEs, dEp, pt, rt, st, detailed)

        self.subtype = 'SWV'

        if 2 * self.pt != self.dt:
            print('\n' + 'Wouldn\'t you rather be using differential pulse voltammetry?' + '\n')
            sys.exit()       

        self.window = np.abs(self.Efin - self.Eini)
        self.dp = int(self.window / np.abs(self.dEs))
        self.tmax = self.dp * self.dt

        '''INDEX'''
        self.index = np.arange(0, ((self.sp * self.tmax + self.dt) / self.dt), 1, dtype = np.int32)

        '''TIME'''
        self.t = (self.index / self.sp) * self.dt

        '''POTENTIAL'''
        self.step = np.array([])
        for ix in range(0, self.dp):
            self.step = np.append(self.step, np.ones((self.sp)) * (ix * self.dEs))
        self.square = square(2 * np.pi * (1/self.dt) * self.t[:-1], duty = self.pt/self.dt) * self.dEp
        
        self.E = np.array([self.Eini])
        self.E = np.append(self.E, (self.step + self.square))

class NPV(pulse):
    '''Waveform for normal pulse voltammetry'''
    def __init__(self, Eini, Efin, dEs, dEp, dt, pt):
        super().__init__(Eini, Efin, dEs, dEp, dt, pt)

        self.subtype = 'NPV'

        self.window = np.abs(self.Efin - self.Eini)
        self.dp = int(self.window / np.abs(self.dEs))
        self.tmax = self.dp * self.dt

        '''INDEX'''
        self.index = np.arange(0, (self.tmax + self.dt) / self.dt, 1, dtype = np.int32)

        '''TIME'''
        self.t = self.index * self.dt

        '''POTENTIAL'''
        self.E = np.array([self.Eini])
        for ix in range(1, self.dp + 1):
            if ix != self.dp:
                self.E = np.append(self.E, square(2 * np.pi * (1/self.dt) * self.t[ix:ix + 1], duty = self.pt/self.dt) *(0.5 * ix * self.dEs) + (0.5 * ix * self.dEs))
                self.E = np.append(self.E, self.Eini)
            else:
                self.E = np.append(self.E, square(2 * np.pi * (1/self.dt) * self.t[ix - 1:self.dp], duty = self.pt/self.dt) *(0.5 * self.dp * self.dEs) + (0.5 * self.dp * self.dEs))
                self.E = np.append(self.E, self.Eini)
                

"""HYBRID CLASSES"""
class CSV(hybrid):
    '''Waveform for cyclic staircase voltammetry'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns, st, detailed):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns, st, detailed)
        
        self.subtype = 'CSV'

        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                
            self.window = np.abs(self.Eupp - self.Elow)
            self.dp = int(self.window / np.abs(self.dE))
            self.tmax = (2 * self.ns * self.window) / self.sr
            self.dt = np.abs(self.dE) / self.sr
            self.sp = int(self.dt / self.st)

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
            self.window = np.abs(self.Eupp - self.Elow)
            self.dp = int(self.window / np.abs(self.dE))
            self.tmax = (2 * self.ns * self.window) / self.sr
            self.dt = np.abs(self.dE) / self.sr
            self.sp = int(self.dt / self.st)

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
            self.uppwindow = np.abs(self.Eupp - self.Eini)
            self.window = np.abs(self.Eupp - self.Elow)
            self.lowwindow = np.abs(self.Eini - self.Elow)
            self.uppdp = int(self.uppwindow / np.abs(self.dE))
            self.dp = int(self.window / np.abs(self.dE))
            self.lowdp = int(self.lowwindow / np.abs(self.dE))
            self.tmax = self.ns * (self.uppwindow + self.window + self.lowwindow) / self.sr
            self.dt = np.abs(self.dE) / self.sr
            self.sp = int(self.dt / self.st)

            '''INDEX'''
            self.index = np.arange(0, ((self.tmax + self.dt) / self.dt), 1, dtype = np.int32)
        
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
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                    self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.uppdp, endpoint = True, dtype = np.float32))
        
        '''PLOTTING WAVEFORM'''
        self.tPLOT = self.t
        self.EPLOT = self.E

        '''EXPORTED WAVEFORM'''
        self.indexWF = np.arange(0, ((self.sp * (self.tmax + self.dt))  / self.dt) + 1, 1, dtype = np.int32)
        self.tWF = (self.indexWF * self.dt) / self.sp        
        self.EWF = np.array([0])
        for ix in range(0, self.E.size, 1):
            self.EWF = np.append(self.EWF, np.ones((self.sp)) * self.E[ix])

class AC(hybrid):
    pass

class EIS(impedance):
    pass


if __name__ == '__main__': 
        
    start = time.time()    
    cwd = os.getcwd()

    try:
        os.makedirs(cwd + '/data')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/data'):
            pass
        else: 
            raise
    filepath = cwd + '/data/' + 'waveform.txt'

    wf = CSV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.005, sr = 0.1, ns = 1, st = 0.0001, detailed = False)

    with open(filepath, 'w') as file:
        for ix, iy, iz in wf.output():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')

    end = time.time()
    print(f'The simulation took {end-start} seconds to complete')