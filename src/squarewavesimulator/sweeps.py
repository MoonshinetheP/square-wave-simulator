import numpy as np
import sys
import os
from errno import EEXIST

import waveforms as wf

class E:
    '''Simulation of single electron transfer reaction: A + e- = B'''
    def __init__(self, input, E0, k0, a, cA, cB, DA, DB, r, h, expansion, Nernstian = False, BV = False, MH = False, Explicit = False, Implicit = False):
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

        self.dX = self.h / self.r

        self.sigma = ((self.r ** 2) / self.Dmax) * (self.F / (self.R * self.Temp)) * self.sr

        self.T = (self.Dmax * self.t) / (self.r ** 2)
        self.dT = self.T[1] - self.T[0]      
        self.Tmax = self.T[-1] 

        self.Xmax = 6 * np.sqrt(self.d * self.Tmax)

        self.theta = (self.F / (self.R * self.Temp)) * (self.E - self.E0)

        self.K0 = (self.k0 * self.r) / self.Dmax
        
        if self.Explicit == True:
            self.lamb = 0.45
            self.dX = np.sqrt(self.dT / self.lamb)
            
            self.x = np.array([0])
            while self.x[-1] < self.Xmax:
                self.x = np.append(self.x, self.x[-1] + self.dX)
                
            self.n = int(self.x.size)
            self.m = int(self.theta.size)

            self.dmod_A = np.zeros((self.n - 1,))
            self.C_A = np.ones((self.n)) * self.CA

            self.dmod_B = np.zeros((self.n - 1,))
            self.C_B = np.ones((self.n)) * self.CB
            
            self.potential = np.array([])
            self.flux = np.array([])

            for k in range(1, self.m + 1, 1):
                '''Forward sweep'''
                if self.Nernstian == True:
                    self.dmod_A[0] = 1 / (1 + np.exp(-self.theta[k - 1])) #self.C_A[0] + self.C_B[0] * np.exp(self.theta[k - 1]) 
                    self.dmod_B[0] =  1 / (1 + np.exp(self.theta[k - 1])) #self.C_B[0] + (self.C_A[0] * np.exp(-self.theta[k - 1]))
                
                if self.BV == True:
                    self.dmod_A[0] = (self.C_A[1] + self.dX  * self.K0 * np.exp((1-self.a) * self.theta[k-1])) / (1 + self.dX * self.K0 * (np.exp((-self.a) * self.theta[k-1]) + np.exp((1-self.a) * self.theta[k - 1])))
                    #CR[k,0] = (CR1kb + dX*K0*np.exp(-alpha*eps[k])*(CO1kb + CR1kb/DOR))/(1 + dX*K0*(np.exp((1-alpha)*eps[k]) + np.exp(-alpha*eps[k])/DOR))
                    # CO[k,0] = CO1kb + (CR1kb - CR[k,0])/DOR

                for x in range(1, self.n - 1):
                    
                    self.dmod_A[x] = self.lamb*self.C_A[x - 1] + (1 - 2*self.lamb)*self.C_A[x] + self.lamb*self.C_A[x + 1]
                
                    #self.dmod_B[x] = self.lamb*self.C_B[x - 1] + (1 - 2*self.lamb)*self.C_B[x] + self.lamb*self.C_B[x + 1]

                self.C_A[self.n - 1] = self.CA
                #self.C_B[self.n - 1] = self.CB

                for y in range(self.n - 2, -1, -1):
                    
                    self.C_A[y] = self.dmod_A[y]
            
                    #self.C_B[y] = self.dmod_B[y]
                self.potential = np.append(self.potential, (self.E0 + ((self.R * self.Temp) / self.F) * self.theta[k-1]))
                self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cA * self.DA) * (-(-self.C_A[2] + 4*self.C_A[1] - 3*self.C_A[0]) / (2 * self.dX)))# - (self.F * np.pi * self.r * self.cB * self.DB) * (-(-self.C_B[2] + 4*self.C_B[1]- 3*self.C_B[0]) / (2*self.dX)))
            self.output = zip(self.potential, self.flux)

        if self.Implicit == True:
            '''Expanding spatial grid'''         
            self.x = np.array([0])
            while self.x[-1] < self.Xmax:
                self.x = np.append(self.x, self.x[-1] + self.dX)
                self.h *= self.expansion

            self.n = int(self.x.size)
            self.m = int(self.theta.size)          
            
            '''Containing arrays'''
            self.alpha_A = np.zeros((self.n - 1,))
            self.beta_A = np.zeros((self.n - 1,))
            self.gamma_A = np.zeros((self.n - 1,))
            self.gmod_A = np.zeros((self.n - 1,))
            self.delta_A = np.ones((self.n - 1,))
            self.dmod_A = np.zeros((self.n - 1,))
            self.C_A = np.ones((self.n)) * self.CA

            self.alpha_B = np.zeros((self.n - 1,))
            self.beta_B = np.zeros((self.n - 1,))
            self.gamma_B = np.zeros((self.n - 1,))
            self.gmod_B = np.zeros((self.n - 1,))
            self.delta_B = np.ones((self.n - 1,))
            self.dmod_B = np.zeros((self.n - 1,))
            self.C_B = np.ones((self.n)) * self.CB
        
            self.delx = np.zeros((self.n,))
            self.delx[0] = self.x[1] - self.x[0]

            for i in range(1, self.n - 1):
                self.delx[i] = self.x[i + 1] - self.x[i]

                self.alpha_A[i] = -(2 * self.dA * self.dT) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                self.gamma_A[i] = -(2 * self.dA * self.dT) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                self.beta_A[i] = 1 - self.alpha_A[i] - self.gamma_A[i]

                self.alpha_B[i] = -(2 * self.dB * self.dT) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                self.gamma_B[i] = -(2 * self.dB * self.dT) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                self.beta_B[i] = 1 - self.alpha_B[i] - self.gamma_B[i]
                    
            self.gmod_A[0] = 0
            for i in range(1, self.n -1, 1):
                self.gmod_A[i] = self.gamma_A[i] / (self.beta_A[i] - self.gmod_A[i-1] * self.alpha_A[i]) 

            self.gmod_B[0] = 0
            for i in range(1, self.n -1, 1):
                self.gmod_B[i] = self.gamma_B[i] / (self.beta_B[i] - self.gmod_B[i-1] * self.alpha_B[i])

            '''Output arrays'''
            self.potential = np.array([])
            self.flux = np.array([])

            

                

            '''Solving'''
            for k in range(1, self.theta.size + 1, 1):
                '''Forward sweep'''
                if self.Nernstian == True:
                    self.dmod_A[0] = 1 / (1 + np.exp(-self.theta[k - 1]))
                    self.dmod_B[0] = 1 / (1 + np.exp(self.theta[k - 1]))

                if self.BV == True:
                    self.dmod_A[0] = self.h *np.exp(-1 * self.a * self.theta[k-1]) * self.K0 * (self.CA + self.CB) * np.exp(self.theta[k-1])
                    self.beta_A[0] = 1 + self.h * np.exp(-1 * self.a * self.theta[k-1]) * self.K0 * (1 + np.exp(self.theta[k-1]))
                    self.gamma_A[0] = -1

                    self.dmod_B[0] = self.h *np.exp(1 - self.a * self.theta[k-1]) * self.K0 * (self.CA + self.CB) * np.exp(self.theta[k-1])
                    self.beta_B[0] = 1 + self.h * np.exp(1 - self.a * self.theta[k-1]) * self.K0 * (1 + np.exp(self.theta[k-1]))
                    self.gamma_A[0] = -1

                for x in range(1, self.n - 1):
                    self.dmod_A[x] = (self.delta_A[x] - self.dmod_A[x - 1] * self.alpha_A[x]) / (self.beta_A[x] - self.gmod_A[x - 1] * self.alpha_A[x])
                
                    self.dmod_B[x] = (self.delta_B[x] - self.dmod_B[x - 1] * self.alpha_B[x]) / (self.beta_B[x] - self.gmod_B[x - 1] * self.alpha_B[x])

                self.C_A[self.n - 1] = self.CA
                self.C_B[self.n - 1] = self.CB
                for y in range(self.n - 2, -1, -1):
                    
                    self.C_A[y] = self.dmod_A[y] - self.gmod_A[y] * self.C_A[y + 1]
                    self.delta_A[y] = self.C_A[y]
            
                    self.C_B[y] = self.dmod_B[y] - self.gmod_B[y] * self.C_B[y + 1]
                    self.delta_B[y] = self.C_B[y]
                    
                '''Appending results'''
                self.potential = np.append(self.potential, (self.E0 + ((self.R * self.Temp) / self.F) * self.theta[k-1]))
                self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cA * self.DA) * (-(self.C_A[1] - self.C_A[0]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cB * self.DB) * (-(self.C_B[1] - self.C_B[0]) / (self.x[1] - self.x[0])))

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
    
    iterables = [0.01]
    for ix in iterables:
        shape = wf.CV(Eini = 0.5, Eupp = 0.5, Elow = 0, dE = -0.001, sr = 0.1, ns = 1)
        instance = E(input = shape, E0 = 0.25, k0 = ix, a = 0.5, cA = 0.005, cB = 0.000, DA = 5E-6, DB = 5E-6, r = 0.15, h = 1E-4, expansion = 1.05, Explicit = True, BV = True)
        filepath = cwd + '/data/' + 'test BV' + str(ix) + '.txt'
        with open(filepath, 'w') as file:
            for ix, iy in instance.results():
                file.write(str(ix) + ',' + str(iy) + '\n')