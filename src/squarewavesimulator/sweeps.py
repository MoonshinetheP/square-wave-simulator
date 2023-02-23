import numpy as np
import sys
import os
from errno import EEXIST

import waveforms as wf

class E:
    '''Simulation of single electron transfer reaction: A + e- = B'''
    def __init__(self, input, E0, k0, a, cA, cB, DA, DB, r, h, omegax, Nernstian = False, BV = False, MH = False):
        '''Imported variables'''
        self.input = input
        self.E = self.input.E
        self.t = self.input.t
        self.dE = self.input.dE
        self.sr = self.input.sr

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

        '''Spatial grid variables'''
        self.r = r
        self.h = h
        self.omegax = omegax

        '''Dimensionless variables''' # I have some problem here, dimensionless variables seem to mess up the simulation
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

        self.dX = self.h / self.r

        self.T = (self.Dmax * self.t) / (self.r ** 2)
        self.Tmax = self.T[-1]
        
        self.sigma = ((self.r ** 2) / self.Dmax) * (self.F / (self.R * self.Temp)) * self.sr

        self.theta = (self.F / (self.R * self.Temp)) * (self.E - self.E0)

        self.K0 = (self.k0 * self.r) / self.Dmax

        '''Expanding spatial grid'''
        self.Xmax = 6 * np.sqrt(self.Dmax * self.Tmax)
        
        self.X = np.array([0])
        
        while self.X[-1] < self.Xmax:
            self.X = np.append(self.X, self.X[-1] + self.dX)
            self.dX *= self.omegax

        self.n = int(self.X.size)

        '''Static temporal grid'''
        self.dT = self.dE / self.sigma
        
        self.m = int(self.input.index.size)
        #int(self.t[-1] / self.dT)

        if self.Nernstian == True:
            '''Containing arrays for A'''
            self.alpha_A = np.zeros((self.n - 1,))
            self.beta_A = np.zeros((self.n - 1,))
            self.gamma_A = np.zeros((self.n - 1,))
            self.gmod_A = np.zeros((self.n - 1,))
            self.delta_A = np.ones((self.n - 1,)) * self.CA
            self.dmod_A = np.zeros((self.n - 1,))
            self.C_A = np.ones((self.n,)) * self.CA # i think shoud be multiplied by dimensionless conc of a
        
            '''Containing arrays for B'''
            self.alpha_B = np.zeros((self.n - 1,))
            self.beta_B = np.zeros((self.n - 1,))
            self.gamma_B = np.zeros((self.n - 1,))
            self.gmod_B = np.zeros((self.n - 1,))
            self.delta_B = np.ones((self.n - 1,)) * self.CB
            self.dmod_B = np.zeros((self.n - 1,))
            self.C_B = np.zeros((self.n,)) * self.CB 
                
            self.Xspace = np.zeros((self.n - 1,)) #should this be n?
            self.Xspace[0] = self.X[1] - self.X[0]

            for i in range(1, self.n - 1):
                self.Xspace[i] = self.X[i + 1] - self.X[i]

                self.alpha_A[i] = (-2 * self.dA * self.dT) / (self.Xspace[i - 1] * (self.Xspace[i - 1] + self.Xspace[i]))
                self.gamma_A[i] = (-2 * self.dA * self.dT) / (self.Xspace[i] * (self.Xspace[i - 1] + self.Xspace[i]))
                self.beta_A[i] = 1 - self.dA * (self.alpha_A[i] + self.gamma_A[i])
                
                self.alpha_B[i] = (-2 * self.dB * self.dT) / (self.Xspace[i - 1] * (self.Xspace[i - 1] + self.Xspace[i]))
                self.gamma_B[i] = (-2 * self.dB * self.dT) / (self.Xspace[i] * (self.Xspace[i - 1] + self.Xspace[i]))
                self.beta_B[i] = 1 - self.dB * (self.alpha_B[i] + self.gamma_B[i])

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
            for k in range(0, self.m, 1):
                '''Forward sweep'''
                self.dmod_A[0] = 1 / (1 + np.exp(-self.theta[k]))
                for x in range(1, self.n - 1):
                    self.dmod_A[x] = (self.delta_A[x] - self.dmod_A[x - 1] * self.alpha_A[x]) / (self.beta_A[x] - self.gmod_A[x - 1] * self.alpha_A[x])
                
                self.dmod_B[0] = 1 / (1 + np.exp(self.theta[k]))
                for x in range(1, self.n - 1):
                    self.dmod_B[x] = (self.delta_B[x] - self.dmod_B[x - 1] * self.alpha_B[x]) / (self.beta_B[x] - self.gmod_B[x - 1] * self.alpha_B[x])


                '''Back substitution'''
                self.C_A[self.n - 1] = 1
                for y in range(self.n - 2, -1, -1):
                    self.C_A[y] = self.dmod_A[y] - self.gmod_A[y] * self.C_A[y + 1]

                    self.delta_A[y] = self.C_A[y]

                self.C_B[self.n - 1] = 1
                for y in range(self.n - 2, -1, -1):
                    self.C_B[y] = self.dmod_B[y] - self.gmod_B[y] * self.C_B[y + 1]

                    self.delta_B[y] = self.C_B[y]

                '''Appending results'''
                self.potential = np.append(self.potential, self.theta[k])
                self.flux = np.append(self.flux, (-self.C_A[1] + self.C_A[0]) / (self.X[1] - self.X[0]))

            '''Finalise results'''
            self.output = zip(self.potential, self.flux)
    
    def results(self):
        return self.output

class EC:

    def __init__(self):
        pass

class CE:

    def __init__(self):
        pass

class ECP:

    def __init__(self):
        pass

class ECE:

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
                    self.alpha_A[i] = -(2 * self.delt) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                    self.gamma_A[i] = -(2 * self.delt) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                    self.beta_A[i] = 1 - self.alpha[i] - self.gamma[i]

                    self.alpha_AP[i] = -(2 * self.delt) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                    self.gamma_AP[i] = -(2 * self.delt) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                    self.beta_AP[i] = 1 - self.alpha[i] - self.gamma[i]

                    self.alpha_B[i] = -(2 * self.delt) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                    self.gamma_B[i] = -(2 * self.delt) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                    self.beta_B[i] = 1 - self.alpha[i] - self.gamma[i]

                    self.alpha_C[i] = -(2 * self.delt) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                    self.gamma_C[i] = -(2 * self.delt) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                    self.beta_C[i] = 1 - self.alpha[i] - self.gamma[i]

                    self.alpha_D[i] = -(2 * self.delt) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
                    self.gamma_D[i] = -(2 * self.delt) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
                    self.beta_D[i] = 1 - self.alpha[i] - self.gamma[i]
                    
            #forward sweep
                self.gmod[0] = 0 # probably will change with new boundary condition
                for i in range(1, self.n - 1):
                    self.gmod[i] = self.gamma[i] / (self.beta[i] - self.gmod[i - 1] * self.alpha[i])
                
            self.dmod[0] = 1 / (1 + np.exp(-self.theta)) #This is the boundary condition
            for i in range(1, self.n - 1):
                self.dmod[i] = (self.delta[i] - self.dmod[i - 1] * self.alpha[i]) / (self.beta[i] - self.gmod[i - 1] * self.alpha[i])

            #back sweep
            self.C[self.n - 1] = 1 #This is the other boundary condition
            for i in range(self.n - 2, -1, -1):
                self.C[i] = self.dmod[i] - self.gmod[i] * self.C[i + 1]

                self.delta[i] = self.C[i] # maybe this has a part about chemical reaction

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

    shape = wf.CV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1)
    instance = E(input = shape, E0 = 0.25, k0 = 1, a = 0.5, cA = 1, cB = 0, DA = 5E-6, DB = 5E-6, r = 0.15, h = 2E-4, omegax = 1.05, Nernstian = True)

    with open(filepath, 'w') as file:
        for ix, iy in instance.results():
            file.write(str(ix) + ',' + str(iy) + '\n')