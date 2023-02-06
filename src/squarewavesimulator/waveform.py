import numpy as np
from scipy.signal import square

class SWV:
    """Simulation of square wave voltammetry"""
    
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
    
    instance = SWV()
    data = 'C:/Users/SLinf/Documents/data.txt'
    with open(data, 'w') as file:
        for ix in instance.combined:
            file.write(str(ix) + '\n')