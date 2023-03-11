import sys
import os
import time

import numpy as np
import waveforms as wf

from errno import EEXIST
from scipy.sparse import diags as diagonals
from scipy.integrate import solve_ivp as solver

class E:
    """Simulation of an E mechanism using solving from scipy \n
    E: A -> B + e"""

    def __init__(self, input, E0, k0, a, cA, cB, DA, DB, r, expansion):
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
        self.cB = cB
        self.DA = DA
        self.DB = DB

        self.F = 96485
        self.R = 8.314
        self.Temp = 298

        '''Spatial grid variables'''
        self.r = r
        self.expansion = expansion

        '''Dimensionless variables''' 
        if self.cA >= self.cB:
            self.cmax = self.cA
        elif self.cA < self.cB:
            self.cmax = self.cB
        
        self.CA = self.cA / self.cmax
        self.CB = self.cB / self.cmax


        if self.DA >= self.DB:
            self.Dmax = self.DA
        elif self.DA < self.DB:
            self.Dmax = self.DB

        self.dA = self.DA / self.Dmax
        self.dB = self.DB / self.Dmax     
        self.d = self.Dmax / self.Dmax

        self.T = (self.Dmax * self.t) / (self.r ** 2)
        self.dT = self.T[1] - self.T[0]      
        self.Tmax = self.T[-1] 

        self.Xmax = 6 * np.sqrt(self.d * self.Tmax)

        self.theta = (self.F / (self.R * self.Temp)) * (self.E - self.E0)

        self.K0 = (self.k0 * self.r) / self.Dmax
        
        '''Expanding spatial grid'''                 
        self.dX = np.sqrt(2 * self.dT)
        self.x = np.array([0])
        while self.x[-1] < self.Xmax:
            self.x = np.append(self.x, self.x[-1] + self.dX)
            self.dX *= self.expansion
        
        self.n = int(self.x.size) 
        self.m = int(self.theta.size)

        self.C = np.ones((self.n, self.m)) * self.CA

        self.alpha = np.ones(self.n - 1)
        self.beta = np.ones(self.n)
        self.gamma = np.ones(self.n - 1)
            
          
        for ix in range(1, self.n - 1):
            self.xplus = self.x[ix + 1] - self.x[ix]
            self.xminus = self.x[ix] - self.x[ix - 1]
            self.denominator = self.xminus * (self.xplus ** 2) + self.xplus * (self.xminus **2)
            self.alpha[ix - 1] *= 2 * self.xplus / self.denominator
            self.beta[ix] *= -2 * (self.xminus + self.xplus) / self.denominator
            self.gamma[ix] *= 2 * self.xminus / self.denominator
            

        A = diagonals([self.alpha, self.beta, self.gamma], [-1,0,1]).toarray()
        A[0,:] = np.zeros(self.n)
        A[0,0] = 1     

        def function(t,y):
            return np.dot(A,y)

        self.flux = np.array([])

        for k in range(1,self.m):
            self.C[0, k] = (self.C[1, k - 1] + (self.x[1] - self.x[0]) * self.K0 * np.exp(-self.a * self.theta[k - 1]))/(1 + (self.x[1] - self.x[0]) * self.K0 * (np.exp((1 - self.a) * self.theta[k - 1]) + np.exp(-self.a * self.theta[k - 1])))
            
            integrator = solver(function, [0, self.dT], self.C[:,k - 1], t_eval=[self.dT], method='RK45')
            
            self.C[1:-1, k] = integrator.y[1:-1, 0]

                            
            self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cA * self.DA) * ((self.C[1, k] - self.C[0, k]) / (self.x[1] - self.x[0]))) # did I save flux properly, or is something still wrong with BV?

        self.output = zip(self.E, self.flux)
        self.simend = time.time()
    
    def results(self):
        return self.output


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
    
    shape = wf.CV(Eini = 0, Eupp = 0.5, Elow = 0.0, dE = 0.002,sr = 0.1, ns = 1)
    instance = E(input = shape, E0 = 0.25, k0 = 0.01, a = 0.5, cA = 0.005, cB = 0.000, DA = 5E-6, DB = 5E-6, r = 0.15, expansion = 1.05)
    
    filepath = cwd + '/data/' + 'test K1 ' + '.txt'
    with open(filepath, 'w') as file:
        for ix, iy in instance.results():
            file.write(str(ix) + ',' + str(iy) + '\n')

    end = time.time()
    print(end-start)

