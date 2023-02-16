import sys
import numpy as np
from scipy.signal import square



"""PARENT CLASSES"""
class sweep:
    '''Parent class for all sweep type waveforms'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):
        self.Eini = Eini
        self.Eupp = Eupp
        self.Elow = Elow
        self.dE = dE
        self.sr = sr
        self.ns = ns

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

            
class step:
    '''Parent class for all step type waveforms'''
    def __init__(self, Eini, Eupp, Elow, dE, sp, sr, ns):
        self.Eini = Eini
        self.Eupp = Eupp
        self.Elow = Elow
        self.dE = dE
        self.sp = sp
        self.sr = sr
        self.ns = ns

        '''DATATYPE ERRORS'''
        if isinstance(self.Eini, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the start potential. Enter either a float or an integer.' + '\n')
            sys.exit()
        if isinstance(self.Eupp, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the upper vertex potential. Enter either a float or an integer.' + '\n')
            sys.exit()        
        if isinstance(self.Elow, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the lower vertex potential. Enter either a float or an integer.' + '\n')
            sys.exit()        
        if isinstance(self.dE, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the step potential. Enter either a float or an integer.' + '\n')
            sys.exit()
        if isinstance(self.sr, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the scan rate. Enter either a float or an integer.' + '\n')
            sys.exit()
        if isinstance(self.ns, (int)) is False:
            print('\n' + 'An invalid datatype was used for the number of scans. Enter an integer.' + '\n')
            sys.exit()
        if isinstance(self.Cd, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the double layer capacitance. Enter either a float or an integer.' + '\n')
            sys.exit()
        if isinstance(self.Ru, (float, int)) is False:
            print('\n' + 'An invalid datatype was used for the uncompensated resistance. Enter either a float or an integer.' + '\n')
            sys.exit()
        if isinstance(self.sp, (int)) is False:
            print('\n' + 'An invalid datatype was used for the number of data points in a step. Enter an integer.' + '\n')
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
        if self.Cd <= 0:
            print('\n' + 'Double layer capacitance must be a positive non-zero value' + '\n')
            sys.exit()
        if self.Ru <= 0:
            print('\n' + 'Uncompensated resistance must be a positive non-zero value' + '\n')
            sys.exit()
        if self.sp <= 0:
            print('\n' + 'Number of data points in a step must be a positive non-zero value' + '\n')
            sys.exit()


        '''WAVEFORM GENERATION'''        
        if self.Eini == self.Elow:          
            self.segments = 2 * self.ns          # Number of segments expected          
            self.window = self.Eupp - self.Elow         # Potential window of each segment
            self.dp = int(self.window/self.dE)       # Number of data points in each potential window
            
            self.index = np.arange(0, (1000 * self.segments * self.window / self.sr), (1000))

            self.sweeptime = np.array([0])
            for iw in range(1, self.segments + 1):
                self.sweeptime = np.append(self.sweeptime, np.linspace((self.sweeptime[-1] + self.dE/self.sr), (iw*self.window/self.sr), self.dp, endpoint = True, dtype = np.float32))
            
            self.sweep = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.sweep = np.append(self.sweep, np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                self.sweep = np.append(self.sweep, np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))

            self.steptime = np.array([])
            for iy in range(0, self.sweeptime.size - 1):
                self.steptime = np.append(self.steptime, np.linspace(self.sweeptime[iy], self.sweeptime[iy + 1], self.sp, dtype = np.float64))

            self.step = np.array([])
            for iz in range(0, self.segments * self.dp):
                self.step = np.append(self.step, np.linspace(self.sweep[iz], self.sweep[iz + 1], self.sp, dtype = np.float64))


        if self.Eini == self.Eupp:
            self.segments = 2 * self.ns          # Number of segments expected          
            self.window = self.Eupp - self.Elow         # Potential window of each segment
            self.dp = int(self.window/self.dE)       # Number of data points in each potential window
            
            self.sweeptime = np.array([0])
            for iw in range(1, self.segments + 1):
                self.sweeptime = np.append(self.sweeptime, np.round(np.linspace((self.sweeptime[-1] + self.dE/self.sr), (iw*self.window/self.sr), self.dp, endpoint = True), decimals = 3))
            
            self.sweep = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.sweep = np.append(self.sweep, np.round(np.linspace(self.Eini + self.dE, self.Elow, self.dp, endpoint = True), decimals = 4))
                self.sweep = np.append(self.sweep, np.round(np.linspace(self.Elow - self.dE, self.Eini, self.dp, endpoint = True), decimals = 4))

        
            self.step = np.array([])
            for iz in range(0, self.segments * self.dp):
                try:
                    self.step = np.append(self.step, np.round(np.linspace(self.sweep[iz], self.sweep[iz + 1] - self.dE/self.sp, self.sp), decimals = 7))
                except:
                    self.step = np.append(self.step, np.round(np.linspace(self.sweep[iz], self.sweep[-1], self.sp + 1), decimals = 7))


        if self.Elow < self.Eini < self.Eupp:
            self.segments = 3 * self.ns          # Number of segments expected          
            self.uppwindow = self.Eupp - self.Eini         # Potential window of each segment
            self.window = self.Eupp - self.Elow
            self.lowwindow = self.Eini - self.Elow
            self.uppdp = int(self.uppwindow/self.dE)
            self.dp = int(self.window/self.dE)       # Number of data points in each potential window
            self.lowdp = int(self.lowwindow/self.dE)

            self.time = np.array([])
            for ix in range(0, self.ns):
                self.time = np.append(self.time, np.round(np.linspace(0, (self.window - self.dE)/self.sr, self.dp), decimals = 3))
                self.time = np.append(self.time, np.round(np.linspace(0, (self.window - self.dE)/self.sr, self.dp), decimals = 3))
                self.time = np.append(self.time, np.round(np.linspace(0, (self.window - self.dE)/self.sr, self.dp), decimals = 3))

            self.sweep = np.array([])
            for iy in range(0, self.ns):
                self.sweep = np.append(self.sweep, np.round(np.linspace(self.Eini, self.Eupp - self.dE, self.uppdp), decimals = 4))
                self.sweep = np.append(self.sweep, np.round(np.linspace(self.Eupp, self.Elow + self.dE, self.dp), decimals = 4))
                self.sweep = np.append(self.sweep, np.round(np.linspace(self.Elow, self.Eini - self.dE, self.lowdp), decimals = 4))
            

"""SWEEP CLASSES"""
class LSV(sweep):
    '''Waveform for linear sweep voltammetry'''
    def __init__(self,Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        if self.ns != 1:
            print('\n' + 'Number of scans should be 1 for LSV' + '\n')
            sys.exit()

        if self.Elow < self.Eini < self.Eupp:
            print('\n' + 'Initial potential should be equal to either upper vertex or lower vertex potential in LSV' + '\n')
            sys.exit()
        
        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:                 
            self.window = self.Eupp - self.Elow
            self.dp = self.window / self.dE
            self.tmax = self.window / self.sr
        
            '''INDEX'''
            self.index = np.arange(0, (1000 * self.tmax), (self.tmax/self.dp))
        
            '''TIME'''
            self.t = self.index / 1000
            
            '''POTENTIAL'''
            self.E = np.linspace(self.Eini, self.Eupp, self.dp + 1, endpoint = True, dtype = np.float32)
        
        '''STARTING FROM UPPER VERTEX POTENTIAL''' 
        if self.Eini == self.Eupp:                 
            self.window = self.Eupp - self.Elow
            self.dp = self.window / self.dE
            self.tmax = self.window / self.sr
        
            '''INDEX'''
            self.index = np.arange(0, (1000 * self.tmax), (self.tmax/self.dp))
        
            '''TIME'''
            self.t = self.index / 1000
            
            '''POTENTIAL'''
            self.E = np.linspace(self.Eini, self.Elow, self.dp + 1, endpoint = True, dtype = np.float32)


class CV(sweep):
    '''Waveform for cyclic voltammetry'''
    def __init__(self, Eini, Eupp, Elow, dE, sr, ns):
        super().__init__(Eini, Eupp, Elow, dE, sr, ns)

        '''STARTING FROM LOWER VERTEX POTENTIAL''' 
        if self.Eini == self.Elow:          
            self.segments = 2 * self.ns        
            self.window = self.Eupp - self.Elow
            self.dp = (2 * self.window) / self.dE
            self.tmax = (2 * self.window) / self.sr

            '''INDEX'''
            self.index = np.arange(0, (1000 * self.segments * self.T), (self.T/self.dp))
        
            '''TIME'''
            self.t = self.index / 1000
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini + self.dE, self.Eupp, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Eupp - self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))

        '''STARTING FROM UPPER VERTEX POTENTIAL'''
        if self.Eini == self.Eupp:
            self.segments = 2 * self.ns        
            self.window = self.Eupp - self.Elow
            self.dp = (2 * self.window) / self.dE
            self.tmax = (2 * self.window) / self.sr

            '''INDEX'''
            self.index = np.arange(0, (1000 * self.segments * self.T), (self.T/self.dp))
        
            '''TIME'''
            self.t = self.index / 1000
            
            '''POTENTIAL'''
            self.E = np.array([self.Eini])
            for ix in range(0, self.ns):
                self.E = np.append(self.E, np.linspace(self.Eini - self.dE, self.Elow, self.dp, endpoint = True, dtype = np.float32))
                self.E = np.append(self.E, np.linspace(self.Elow + self.dE, self.Eini, self.dp, endpoint = True, dtype = np.float32))


        '''STARTING IN BETWEEN VERTEX POTENTIALS'''
        if self.Elow < self.Eini < self.Eupp:
            self.segments = 3 * self.ns         
            self.uppwindow = self.Eupp - self.Eini
            self.window = self.Eupp - self.Elow
            self.lowwindow = self.Eini - self.Elow
            self.uppdp = self.uppwindow / self.dE
            self.dp = self.window / self.dE
            self.lowdp = self.lowwindow / self.dE
            self.tmax = (self.uppwindown + self.window + self.lowwindow) / self.sr

            '''INDEX'''
            self.index = np.arange(0, (1000 * self.segments * self.tmax), (self.tmax/self.dp))
        
            '''TIME'''
            self.t = self.index / 1000
            
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
    def __init__(self):
        pass

class CSV(step):
    '''Waveform for cyclic staircase voltammetry'''
    def __init__(self):
        pass

class NPV(step):
    '''Waveform for normal pulse voltammetry'''
    def __init__(self):
        pass

class DPV(step):
    '''Waveform for differential pulse voltammetry'''
    def __init__(self):
        pass

class RPV(step):
    '''Waveform for reverse pulse voltammetry'''
    def __init__(self):
        pass

class SWV(step):
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
    
    instance = CV(Eini = 0, Eupp = 'a', Elow = -0.5, sr = 0.1, ns = 1)
    data = 'C:/Users/SLinf/Documents/data.txt'
    '''with open(data, 'w') as file:
        for ix in instance.combined:
            file.write(str(ix) + '\n')'''