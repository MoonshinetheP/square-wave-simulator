import numpy as np
import sys
import os
from errno import EEXIST

import waveforms as wf

class E:
    '''Simulation of single electron transfer reaction: A + e- = B'''
    def __init__(self, input, E0, k0, a, cA, DA, r, h, expansion, Nernstian = False, BV = False, MH = False, Explicit = False, Implicit = False):
        '''Waveform variables'''
        self.input = input
        self.index = self.input.index
        self.t = self.input.t
        self.E = self.input.E
        
        self.Eini = self.input.Eini
        self.Eupp = self.input.Eupp
        self.Elow = self.input.Elow
        self.dE = self.input.dE
        self.sr = self.input.sr
        self.ns = self.input.ns
        
        '''Mechanism variables'''
        self.E0 = E0
        self.k0 = k0
        self.a = a
        self.cA = cA
        self.DA = DA

        self.F = 96485
        self.R = 8.314
        self.Temp = 298

        '''Spatial grid variables'''
        self.r = r
        self.h = h
        self.expansion = expansion

        '''Kinetics'''
        self.Nernstian = Nernstian
        self.BV = BV
        self.MH = MH

        self.methods = 0
        for ix in [self.Nernstian, self.BV, self.MH]:
            if ix == True:
                self.methods += 1

        if self.methods == 0:
            print('\n' + 'No kinetic model was chosen' + '\n')
            sys.exit()

        if self.methods >= 2:
            print('\n' + 'More than one kinetic model was chosen' + '\n')
            sys.exit()

        '''Simulation style'''
        self.Explicit = Explicit
        self.Implicit = Implicit

        self.style = 0
        for ix in [self.Explicit, self.Implicit]:
            if ix == True:
                self.style += 1
        
        if self.style == 0:
            print('\n' + 'No simulation style was chosen' + '\n')
            sys.exit()

        if self.style == 2:
            print('\n' + 'More than one simulation style was chosen' + '\n')
            sys.exit()

        '''Dimensionless variables''' 
        self.CA = self.cA / self.cA
        
        self.dA = self.DA / self.DA

        self.dX = self.h / self.r

        self.sigma = ((self.r ** 2) / self.DA) * (self.F / (self.R * self.Temp)) * self.sr

        self.T = (self.DA * self.t) / (self.r ** 2)
        self.dT = self.dX / self.sigma        
        self.Tmax = self.T[-1] 

        self.Xmax = 6 * np.sqrt(self.dA * self.Tmax)

        self.theta = (self.F / (self.R * self.Temp)) * (self.E - self.E0)

        self.K0 = (self.k0 * self.r) / self.DA
        
        if self.Implicit == True:
            '''Expanding spatial grid'''         
            self.n = 1 + int(self.Xmax / self.dX)
            self.m = int(self.theta.size)          
            
            '''Containing arrays'''
            self.gmod = np.zeros((self.n - 1,))
            self.delta = np.ones((self.n - 1,))
            self.dmod = np.zeros((self.n - 1,))
            self.C = np.ones((self.n))
        
            if self.Nernstian == True:

                self.l = self.dT / (self.dX * self.dX)
                self.alpha = -self.l
                self.beta = 2 * self.l + 1
                self.gamma = -self.l

                self.gmod[0] = 0
                for i in range(1, self.n -1, 1): # same
                    self.gmod[i] = self.gamma / (self.beta - self.gmod[i-1] * self.alpha) 

                '''Output arrays'''
                self.potential = np.array([])
                self.flux = np.array([])

                '''Solving'''
                for k in range(1, self.theta.size + 1, 1):

                    '''Forward sweep'''
                    self.dmod[0] = 1 / (1 + np.exp(-self.theta[k - 1]))
                    for x in range(1, self.n - 1):
                        self.dmod[x] = (self.delta[x] - self.dmod[x - 1] * self.alpha) / (self.beta - self.gmod[x - 1] * self.alpha)

                    self.C[self.n - 1] = 1
                    for y in range(self.n - 2, -1, -1):
                        self.C[y] = self.dmod[y] - self.gmod[y] * self.C[y + 1]

                        self.delta[y] = self.C[y]

                    '''Appending results'''
                    
                    self.potential = np.append(self.potential, self.theta[k-1])
                    self.flux = np.append(self.flux, -(self.C[1] - self.C[0]) / (self.dX))

                '''Finalise results'''
                self.output = zip(self.potential, self.flux)
    
    def results(self):
        return self.output


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

    shape = wf.CV(Eini = 0.5, Eupp = 0.5, Elow = 0, dE = -0.001, sr = 0.1, ns = 1)
    instance = E(input = shape, E0 = 0.25, k0 = 1, a = 0.5, cA = 1, DA = 5E-6, r = 0.15, h = 5E-6, expansion = 1.05, Implicit = True, Nernstian = True)

    with open(filepath, 'w') as file:
        for ix, iy in instance.results():
            file.write(str(ix) + ',' + str(iy) + '\n')