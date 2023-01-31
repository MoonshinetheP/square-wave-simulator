import numpy as np
from scipy.signal import square

class SWV:
    '''Makes a square wave voltammetry waveform'''
    
    def __init__(self, Eini = 0, Efin = 1, dEs = 0.002, dEp = 0.05, f = 25, sp = 1000):

        self.Eini = Eini
        self.Efin = Efin
        self.dEs = dEs
        self.dEp = dEp
        self.f = f
        self.sp = sp

        self.window = self.Efin - self.Eini
        self.dp = int(self.window / self.dEs)
        self.step_period = 1 / self.f
        self.pulse_period = self.step_period / 2


        self.sweep_time = np.linspace(0 , self.step_period * self.dp, self.dp + 1, endpoint = True)
        self.sweep_waveform = np.linspace(self.Eini, self.Efin, self.dp + 1, endpoint = True)

        self.step_time = np.array([])
        for ix in range(0, self.sweep_time.size):
            try:
                self.step_time = np.append(self.step_time, np.linspace(self.sweep_time[ix], self.sweep_time[ix + 1], self.sp))
            except:
                pass

        self.step_waveform = np.array([])
        for iy in range(1, self.sweep_waveform.size):
            self.step_waveform = np.append(self.step_waveform, np.ones((self.sp)) * self.sweep_waveform[iy])

        self.square_waveform = square(2 * np.pi * self.f * self.step_time, duty = 0.5) * self.dEp
        
        self.combined = self.step_waveform + self.square_waveform

if __name__ == '__main__':
    
    instance = SWV()
    data = 'C:/Users/SLinf/Documents/data.txt'
    with open(data, 'w') as file:
        for ix in instance.combined:
            file.write(str(ix) + '\n')