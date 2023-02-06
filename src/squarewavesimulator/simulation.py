import numpy as np
import os
from errno import EEXIST

import waveform

class Conversion:
    """Converts the inputs of an electrochemical simulation into dimensionless equivalents"""

    def __init__(self):
        pass

class Deconversion:
    """Converts the outputs of an electrochemical simulation into quantities with dimensionality"""
    
    def __init__(self):
        pass


class EEC:
    """Simulation"""

    def __init__(self, theta, tmax = 1, omegax = 1.05, h = 1E-4, omegat = 1.015):
        
        self.theta = theta
        self.tmax = tmax
        self.omegax = omegax
        self.h = h
        self.omegat = omegat
        
        self.deltat = self.tmax / 4000
        self.ts = self.tmax / 2

        self.maxx = 6 * np.sqrt(self.tmax)
        
        '''Expanding grid'''
        self.x = np.array([0])
        while self.x[-1] < self.maxx:
            self.x = np.append(self.x, self.x[-1] + self.h)
            self.h *= self.omegax

        self.n = int(self.x.size)
        
        '''Expanding time'''
        self.t = np.array([0])
        while self.t[-1] < self.tmax:
            self.t = np.append(self.t, self.t[-1] + self.deltat)
            if self.t[-1] > self.ts:
                self.deltat *= self.omegat

        self.m = int(self.t.size)
        self.t[self.m - 1] = self.tmax
        
        '''Containing arrays'''
        self.alpha = np.zeros((self.n - 1,))
        self.beta = np.zeros((self.n - 1,))
        self.gamma = np.zeros((self.n - 1,))
        self.gmod = np.zeros((self.n - 1,))
        self.delta = np.ones((self.n - 1,))
        self.dmod = np.zeros((self.n - 1,))
        self.C = np.ones((self.n,))

        
        '''Thomas coefficients for T < Ts'''
        self.delt = self.t[1] - self.t[0]
        self.delx = np.zeros((self.n,))
        self.delx[0] = self.x[1] - self.x[0]

        for i in range(1, self.n - 1):
            self.delx[i] = self.x[i + 1] - self.x[i]

            self.alpha[i] = -(2 * self.delt) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
            self.gamma[i] = -(2 * self.delt) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
            self.beta[i] = 1 - self.alpha[i] - self.gamma[i]
        

        '''Modified gamma coefficients for T < Ts'''   
        self.gmod[0] = 0
        for i in range(1, self.n -1, 1):
            self.gmod[i] = self.gamma[i] / (self.beta[i] - self.gmod[i - 1] * self.alpha[i]) 


        self.time = np.array([])
        self.flux = np.array([])
        self.qj = np.array([])

        for k in range(1, self.m, 1):
            if self.t[k] > self.ts:
                for i in range(1, self.n - 1):
                    self.delt = self.t[k] - self.t[k-1]
                    self.alpha[i] = -(2 * self.delt) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                    self.gamma[i] = -(2 * self.delt) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                    self.beta[i] = 1 - self.alpha[i] - self.gamma[i]
            
                self.gmod[0] = 0
                for i in range(1, self.n - 1):
                    self.gmod[i] = self.gamma[i] / (self.beta[i] - self.gmod[i - 1] * self.alpha[i])
                
            self.dmod[0] = 1 / (1 + np.exp(-self.theta))
            for i in range(1, self.n - 1):
                self.dmod[i] = (self.delta[i] - self.dmod[i - 1] * self.alpha[i]) / (self.beta[i] - self.gmod[i - 1] * self.alpha[i])

            self.C[self.n - 1] = 1
            for i in range(self.n - 2, -1, -1):
                self.C[i] = self.dmod[i] - self.gmod[i] * self.C[i + 1]

                self.delta[i] = self.C[i]

            self.time = np.append(self.time, self.t[k])
            self.flux = np.append(self.flux, -(self.C[1] - self.C[0]) / (self.x[1] - self.x[0]))
            self.ref = -1 / np.sqrt(np.pi * self.t[k])
            self.qj = np.append(self.qj, (self.flux[-1] - self.ref) / (self.ref) * 100)

        self.dist = np.array([])
        self.conc = np.array([])    
        for i in range(0, self.n):
            self.dist = np.append(self.dist, self.x[i])
            self.conc = np.append(self.conc, self.C[i])
        
        self.output = zip(self.time, self.flux, self.qj)
        self.output2 = zip(self.dist, self.conc)
        
    def results(self):
        return(self.output)


if __name__ == '__main__':
        
    cwd = os.getcwd()

    try:
        os.makedirs(cwd + '/data')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/data'):
            pass
        else: 
            raise
    filepath = cwd + '/data/' + 'test.txt'

    instance = EEC(theta = -15)

    with open(filepath, 'w') as file:
        for ix, iy, iz in instance.results():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')