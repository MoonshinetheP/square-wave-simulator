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
    E: R -> O + e"""

    def __init__(self, input, E0, k0, a, cR, cO, DR, DO, r, expansion, Nernstian, BV, MH):
        '''Waveform variables'''
        self.input = input
        self.index = self.input.index
        self.t = self.input.t
        self.E = self.input.E
        self.tPLOT = self.input.tPLOT
        self.EPLOT = self.input.EPLOT

        if self.input.type == 'sweep':
            self.Eini = self.input.Eini
            self.Eupp = self.input.Eupp
            self.Elow = self.input.Elow
            self.dE = self.input.dE
            self.sr = self.input.sr
            self.ns = self.input.ns

            self.detailed = False

        if self.input.type == 'pulse':
            self.Eini = self.input.Eini
            self.Efin = self.input.Efin
            self.dEs = self.input.dEs
            self.dEp = self.input.dEp
            self.pt = self.input.pt
            self.rt = self.input.rt
            self.st = self.input.st
            self.detailed = self.input.detailed


        
        '''Mechanism variables'''
        self.E0 = E0
        self.k0 = k0
        self.a = a
        self.cR = cR
        self.cO = cO
        self.DR = DR
        self.DO = DO

        self.F = 96485
        self.R = 8.314
        self.Temp = 298
        
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
        self.expansion = expansion

        '''Dimensionless variables''' 
        concentrations = np.array([self.cR, self.cO])
        self.cmax = np.amax(concentrations)       
        self.CR = self.cR / self.cmax
        self.CO = self.cO / self.cmax

        diffusions = np.array([self.DR, self.DO])
        self.Dmax = np.amax(diffusions)

        self.dR = self.DR / self.Dmax
        self.dO = self.DO / self.Dmax  
        self.d = self.Dmax / self.Dmax

        self.T = (self.Dmax * self.t) / (self.r ** 2)
        self.dT = np.diff(self.T)
        
        if self.detailed == False:
            self.sT = self.dT
        else:
            self.sT = np.ones((self.dT.size)) * (self.Dmax * self.st) / (self.r ** 2)
        self.Tmax = self.T[-1] 
        

        self.Xmax = 6 * np.sqrt(self.d * self.Tmax)

        self.theta = (self.F / (self.R * self.Temp)) * (self.E - self.E0)

        self.K0 = (self.k0 * self.r) / self.Dmax
        
        '''Expanding spatial grid'''                 
        self.dX = np.sqrt(2.05 * self.dT[0])
        self.x = np.array([0])
        while self.x[-1] < self.Xmax:
            self.x = np.append(self.x, self.x[-1] + self.dX)
            self.dX *= self.expansion
        
        self.n = int(self.x.size) 
        self.m = int(self.theta.size)

        self.C_R = np.ones((self.n, self.m)) * self.CR
        self.C_O = np.ones((self.n, self.m)) * self.CO

        self.alpha_R = np.ones(self.n - 1)
        self.beta_R = np.ones(self.n)
        self.gamma_R = np.ones(self.n - 1)
        
        self.alpha_O = np.ones(self.n - 1)
        self.beta_O = np.ones(self.n)
        self.gamma_O = np.ones(self.n - 1)            
          
        for ix in range(1, self.n - 1):
            self.xplus = self.x[ix + 1] - self.x[ix]
            self.xminus = self.x[ix] - self.x[ix - 1]
            self.denominator = self.xminus * (self.xplus ** 2) + self.xplus * (self.xminus **2)
            
            self.alpha_R[ix - 1] *= 2 * self.dR * self.xplus / self.denominator
            self.beta_R[ix] *= -2 * self.dR * (self.xminus + self.xplus) / self.denominator
            self.gamma_R[ix] *= 2 * self.dR * self.xminus / self.denominator
            
            self.alpha_O[ix - 1] *= 2 * self.dO * self.xplus / self.denominator
            self.beta_O[ix] *= -2 * self.dO * (self.xminus + self.xplus) / self.denominator
            self.gamma_O[ix] *= 2 * self.dO * self.xminus / self.denominator
        
            

        R = diagonals([self.alpha_R, self.beta_R, self.gamma_R], [-1,0,1]).toarray()
        R[0,:] = np.zeros(self.n)
        R[0,0] = 1     
        
        O = diagonals([self.alpha_O, self.beta_O, self.gamma_O], [-1,0,1]).toarray()
        O[0,:] = np.zeros(self.n)
        O[0,0] = 1   

        def reduced(t,y):
            return np.dot(R,y)
        
        def oxidised(t,y):
            return np.dot(O,y)

        self.flux = np.array([0])
        self.fluxplus = np.array([])
        self.fluxminus = np.array([])
        for k in range(1,self.m):
            '''Boundary conditions'''
            if self.Nernstian == True:
                '''Nernstian'''
                self.C_R[0, k] = (self.C_R[1, k - 1] + (self.dR/self.dO) * self.C_O[1, k - 1])/(1 + (self.dR/self.dO) * np.exp(self.theta[k-1]))

                self.C_O[0, k] = (self.C_R[1, k - 1] + (self.dR/self.dO) * self.C_O[1, k - 1])/((self.dR/self.dO) + np.exp(-self.theta[k-1]))
            
            if self.BV == True:
                '''Butler-Volmer'''
                self.C_R[0, k] = (-self.C_R[1, k - 1] + (self.x[1] - self.x[0]) * self.K0 * np.exp(-self.a * self.theta[k - 1]) * (self.C_O[1, k - 1] + (self.dR/self.dO) * self.C_R[1, k - 1]))/((self.x[1] - self.x[0]) * self.K0 * (np.exp((1 - self.a) * self.theta[k - 1]) + (self.dR/self.dO) * np.exp((-self.a) * self.theta[k - 1])) - 1)

                self.C_O[0, k] = (-self.C_O[1, k - 1] + (self.x[1] - self.x[0]) * self.K0 * np.exp((1 - self.a) * self.theta[k - 1]) * (self.C_R[1, k - 1] + (self.dO/self.dR) * self.C_O[1, k - 1]))/((self.x[1] - self.x[0]) * self.K0 * (np.exp(-self.a * self.theta[k - 1]) + (self.dO/self.dR) * np.exp((1 - self.a) * self.theta[k - 1])) - 1)
           
            oxidation = solver(reduced, [0, self.dT[k - 1]], self.C_R[:,k - 1], t_eval=[self.sT[k - 1]], method='RK45')
            self.C_R[1:-1, k] = oxidation.y[1:-1, 0]
            
            reduction = solver(oxidised, [0, self.dT[k - 1]], self.C_O[:,k - 1], t_eval=[self.sT[k - 1]], method='RK45')
            self.C_O[1:-1, k] = reduction.y[1:-1, 0]
            
            if self.input.type == 'sweep':
                self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_R[1, k] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_O[1, k] - self.C_O[0, k]) / (self.x[1] - self.x[0])))
            
            if self.input.type == 'pulse':
                if k % 2 != 0:
                    self.fluxplus = np.append(self.fluxplus, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_R[1, k] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_O[1, k] - self.C_O[0, k]) / (self.x[1] - self.x[0])))
                if k % 2 == 0:
                    self.fluxminus = np.append(self.fluxminus, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_R[1, k] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_O[1, k] - self.C_O[0, k]) / (self.x[1] - self.x[0])))

        if self.input.type == 'pulse':
            self.flux = self.fluxplus - self.fluxminus

        self.output = zip(self.tPLOT, self.EPLOT, self.flux)
    
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
    
    shape = wf.DPV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.01, pt = 0.01, rt = 0.03, st = 0.000001, detailed = False)
    instance = E(input = shape, E0 = 0.25, k0 = 0.1, a = 0.5, cR = 0.005, cO = 0.000, DR = 5E-6, DO = 5E-6, r = 0.15, expansion = 1.05, Nernstian = False, BV = True, MH = False)
    
    filepath = cwd + '/data/' + 'pulse' + '.txt'
    with open(filepath, 'w') as file:
        for ix, iy, iz in instance.results():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')

    end = time.time()
    print(end-start)

