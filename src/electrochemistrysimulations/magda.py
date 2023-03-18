import numpy as np
import sys
import os
from errno import EEXIST
import time
import waveforms as wf

class EEC:
    """Simulation of E(E)C mechanism using explicit solving of differential equations \n
    E1: A -> B + e \n 
    E2: B + e = C \n 
    C1: B -> D """

    def __init__(self, input, E0A, E0B, k1, k2, k3, a1, a2, cA, cB, cC, cD, DA, DB, DC, DD, r, h, expansion):
        '''Waveform variables'''
        self.input = input
        self.index = self.input.index
        self.t = self.input.t
        self.E = self.input.simulation()
        
        self.Eini = self.input.Eini
        self.Efin = self.input.Efin
        self.dEs = self.input.dEs
        self.dEp = self.input.dEp
        self.dt = self.input.dt
        self.pt = self.input.pt
        self.sp = self.input.sp
        
        '''Mechanism variables'''
        self.E0A = E0A
        self.E0B = E0B
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.a1 = a1
        self.a2 = a2
        self.cA = cA
        self.cB = cB
        self.cC = cC
        self.cD = cD
        self.DA = DA
        self.DB = DB
        self.DC = DC
        self.DD = DD

        self.F = 96485
        self.R = 8.314
        self.Temp = 298

        '''Spatial grid variables'''
        self.r = r
        self.h = h
        self.expansion = expansion

        '''Dimensionless variables''' 
        if self.cA >= self.cB:
            self.cmax = self.cA
        elif self.cA < self.cB:
            self.cmax = self.cB
        
        self.CA = self.cA / self.cmax
        self.CB = self.cB / self.cmax
        self.CC = self.cC / self.cmax
        self.CD = self.cD / self.cmax

        if self.DA >= self.DB:
            self.Dmax = self.DA
        elif self.DA < self.DB:
            self.Dmax = self.DB

        self.dA = self.DA / self.Dmax
        self.dB = self.DB / self.Dmax
        self.dC = self.DC / self.Dmax
        self.dD = self.DD / self.Dmax        
        self.d = self.Dmax / self.Dmax

        self.dX = self.h / self.r

        self.T = (self.Dmax * self.t) / (self.r ** 2)
        self.dT = self.T[1] - self.T[0]      
        self.Tmax = self.T[-1] 

        self.Xmax = 6 * np.sqrt(self.d * self.Tmax)

        self.theta = (self.F / (self.R * self.Temp)) * (self.E - self.E0A)
        self.gap = (self.F / (self.R * self.Temp)) * (self.E0B - self.E0A)

        self.K1 = (self.k1 * self.r) / self.Dmax
        self.K2 = (self.k2 * self.r) / self.Dmax
        self.K3 = (self.k3 * self.r) / self.Dmax
        
        '''Expanding spatial grid'''         
        self.lamb = 0.45
        self.dX = np.sqrt(self.dT / self.lamb)
        
        self.x = np.array([0])
        while self.x[-1] < self.Xmax:
            self.x = np.append(self.x, self.x[-1] + self.dX)
            self.dX *= self.expansion

                
        self.n = int(self.x.size)
        self.m = int(self.theta.size)
            
        self.sT = (self.Dmax * self.dt) / (self.r ** 2)
        self.pT = (self.Dmax * self.pt) / (self.r ** 2) # and maybe here I need to define some expanding time

        self.sn = int(abs(self.sT - self.pT) / self.dT)
        self.pn = int(self.pT / self.dT)
        
            
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

        self.alpha_C = np.zeros((self.n - 1,))
        self.beta_C = np.zeros((self.n - 1,))
        self.gamma_C = np.zeros((self.n - 1,))
        self.gmod_C = np.zeros((self.n - 1,))
        self.delta_C = np.ones((self.n - 1,))
        self.dmod_C = np.zeros((self.n - 1,))
        self.C_C = np.ones((self.n)) * self.CC

        self.alpha_D = np.zeros((self.n - 1,))
        self.beta_D = np.zeros((self.n - 1,))
        self.gamma_D = np.zeros((self.n - 1,))
        self.gmod_D = np.zeros((self.n - 1,))
        self.delta_D = np.ones((self.n - 1,))
        self.dmod_D = np.zeros((self.n - 1,))
        self.C_D = np.ones((self.n)) * self.CD
        
        self.delx = np.zeros((self.n,))
        self.delx[0] = self.x[1] - self.x[0]

        for i in range(1, self.n - 1):
            self.delx[i] = self.x[i + 1] - self.x[i]

            self.alpha_A[i] = -(2 * self.dA * self.dT) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
            self.gamma_A[i] = -(2 * self.dA * self.dT) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
            self.beta_A[i] = 1 - self.alpha_A[i] - self.gamma_A[i]

            self.alpha_B[i] = -(2 * self.dB * self.dT) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
            self.gamma_B[i] = -(2 * self.dB * self.dT) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
            self.beta_B[i] = 1 - self.alpha_B[i] - self.gamma_B[i] + self.K3*self.C_B[i]
            
            self.alpha_C[i] = -(2 * self.dC * self.dT) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
            self.gamma_C[i] = -(2 * self.dC * self.dT) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
            self.beta_C[i] = 1 - self.alpha_C[i] - self.gamma_C[i]
            
            self.alpha_D[i] = -(2 * self.dD * self.dT) / (self.delx[i - 1] * (self.delx[i - 1] + self.delx[i]))
            self.gamma_D[i] = -(2 * self.dD * self.dT) / (self.delx[i] * (self.delx[i - 1] + self.delx[i]))
            self.beta_D[i] = 1 - self.alpha_D[i] - self.gamma_D[i] - self.K3*self.C_B[i]
                    
        
        self.potential = np.linspace(self.Eini, self.Efin, self.T.size, endpoint = True)
        self.flux = np.array([])
        self.possteppotential = np.array([])
        self.posstepflux = np.array([])
        self.negsteppotential = np.array([])
        self.negstepflux = np.array([])


        self.waveforms = time.time()
        for k in range(1, self.m + 1, 1):
            
            #if k ==2: print(f'Expected completion time: {(self.m - 1) * (time.time() - self.waveforms)} seconds')

            if k % 2 != 0:
                for j in range(0, self.sn):
            
                    self.dmod_A[0] = (self.C_A[1] + self.dX  * self.K1 * np.exp((1-self.a1) * self.theta[k-1])) / (1 + self.dX * self.K1 * (np.exp((-self.a1) * self.theta[k-1]) + np.exp((1-self.a1) * self.theta[k - 1])))
                    
                    self.dmod_B[0] = (self.C_B[1] + self.dX  * self.K1 * np.exp((-self.a1) * self.theta[k-1])) / (1 + self.dX * self.K1 * (np.exp((1-self.a1) * self.theta[k-1]) + np.exp((-self.a1) * self.theta[k - 1]))) + (self.C_B[1] + self.dX  * self.K2 * np.exp((1-self.a2) * (self.theta[k-1] - self.gap))) / (1 + self.dX * self.K2 * (np.exp((-self.a2) * (self.theta[k-1] - self.gap)) + np.exp((1-self.a2) * (self.theta[k-1] - self.gap))))
                    
                    self.dmod_C[0] = (self.C_C[1] + self.dX  * self.K2 * np.exp((1-self.a2) * (self.theta[k-1] - self.gap))) / (1 + self.dX * self.K2 * (np.exp((-self.a2) * (self.theta[k-1] - self.gap)) + np.exp((1-self.a2) * (self.theta[k - 1]) - self.gap)))
                    
                    self.dmod_D[0] = 0

                    for x in range(1, self.n - 1):  #I probably need to edit these to deal with expanding grid
                        self.dmod_A[x] = self.lamb*self.C_A[x - 1] + (1 - 2*self.lamb)*self.C_A[x] + self.lamb*self.C_A[x + 1]
                        self.dmod_B[x] = self.lamb*self.C_B[x - 1] + (1 - 2*self.lamb)*self.C_B[x] + self.lamb*self.C_B[x + 1]
                        self.dmod_C[x] = self.lamb*self.C_C[x - 1] + (1 - 2*self.lamb)*self.C_C[x] + self.lamb*self.C_C[x + 1]
                        self.dmod_D[x] = self.lamb*self.C_D[x - 1] + (1 - 2*self.lamb)*self.C_D[x] + self.lamb*self.C_D[x + 1]

                    self.C_A[self.n - 1] = self.CA
                    self.C_B[self.n - 1] = self.CB
                    self.C_C[self.n - 1] = self.CC
                    self.C_D[self.n - 1] = self.CD

                    for y in range(self.n - 2, -1, -1): 
                        self.C_A[y] = self.dmod_A[y]
                        self.C_B[y] = self.dmod_B[y]
                        self.C_C[y] = self.dmod_C[y]
                        self.C_D[y] = self.dmod_D[y]
                                      
                self.negsteppotential = np.append(self.negsteppotential, (self.E0A + ((self.R * self.Temp) / self.F) * self.theta[k-1]))
                self.negstepflux = np.append(self.negstepflux, (self.F * np.pi * self.r * self.cA * self.DA) * (-(self.C_A[1] - self.C_A[0]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cB * self.DB) * (-(self.C_B[1] - self.C_B[0]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cC * self.DC) * (-(self.C_C[1] - self.C_C[0]) / (self.x[1] - self.x[0])))
                 
            self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cA * self.DA) * (-(-self.C_A[2] + 4*self.C_A[1] - 3*self.C_A[0]) / (2 * self.dX)) - (self.F * np.pi * self.r * self.cB * self.DB) * (-(-self.C_B[2] + 4*self.C_B[1]- 3*self.C_B[0]) / (2*self.dX)) - (self.F * np.pi * self.r * self.cC * self.DC) * (-(-self.C_C[2] + 4*self.C_C[1]- 3*self.C_C[0]) / (2*self.dX)))
            

            if k % 2 == 0:
                for l in range(0, self.pn):
            
                    self.dmod_A[0] = (self.C_A[1] + self.dX  * self.K1 * np.exp((1-self.a1) * self.theta[k-1])) / (1 + self.dX * self.K1 * (np.exp((-self.a1) * self.theta[k-1]) + np.exp((1-self.a1) * self.theta[k - 1])))
                    
                    self.dmod_B[0] = (self.C_B[1] + self.dX  * self.K1 * np.exp((-self.a1) * self.theta[k-1])) / (1 + self.dX * self.K1 * (np.exp((1-self.a1) * self.theta[k-1]) + np.exp((-self.a1) * self.theta[k - 1])))
                    
                    self.dmod_C[0] = (self.C_C[1] + self.dX  * self.K2 * np.exp((1-self.a2) * (self.theta[k-1] - self.gap))) / (1 + self.dX * self.K2 * (np.exp((-self.a2) * (self.theta[k-1] - self.gap)) + np.exp((1-self.a2) * (self.theta[k - 1]) - self.gap)))
                    
                    self.dmod_D[0] = 0

                    for x in range(1, self.n - 1):  
                        self.dmod_A[x] = self.lamb*self.C_A[x - 1] + (1 - 2*self.lamb)*self.C_A[x] + self.lamb*self.C_A[x + 1]
                        self.dmod_B[x] = self.lamb*self.C_B[x - 1] + (1 - 2*self.lamb)*self.C_B[x] + self.lamb*self.C_B[x + 1]
                        self.dmod_C[x] = self.lamb*self.C_C[x - 1] + (1 - 2*self.lamb)*self.C_C[x] + self.lamb*self.C_C[x + 1]
                        self.dmod_D[x] = self.lamb*self.C_D[x - 1] + (1 - 2*self.lamb)*self.C_D[x] + self.lamb*self.C_D[x + 1]

                    self.C_A[self.n - 1] = self.CA
                    self.C_B[self.n - 1] = self.CB
                    self.C_C[self.n - 1] = self.CC
                    self.C_D[self.n - 1] = self.CD

                    for y in range(self.n - 2, -1, -1): 
                        self.C_A[y] = self.dmod_A[y]
                        self.C_B[y] = self.dmod_B[y]
                        self.C_C[y] = self.dmod_C[y]
                        self.C_D[y] = self.dmod_D[y]
                
                self.possteppotential = np.append(self.possteppotential, (self.E0A + ((self.R * self.Temp) / self.F) * self.theta[k-1]))
                self.posstepflux = np.append(self.posstepflux, (self.F * np.pi * self.r * self.cA * self.DA) * (-(self.C_A[1] - self.C_A[0]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cB * self.DB) * (-(self.C_B[1] - self.C_B[0]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cC * self.DC) * (-(self.C_C[1] - self.C_C[0]) / (self.x[1] - self.x[0])))
            
            self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cA * self.DA) * (-(-self.C_A[2] + 4*self.C_A[1] - 3*self.C_A[0]) / (2 * self.dX)) - (self.F * np.pi * self.r * self.cB * self.DB) * (-(-self.C_B[2] + 4*self.C_B[1]- 3*self.C_B[0]) / (2*self.dX)) - (self.F * np.pi * self.r * self.cC * self.DC) * (-(-self.C_C[2] + 4*self.C_C[1]- 3*self.C_C[0]) / (2*self.dX)))
            
            #print(f'{k} out of {self.theta.size}')

        self.negsimplified = zip(self.negsteppotential, self.negstepflux)
        self.possimplified = zip(self.possteppotential, self.posstepflux)
        self.detailed = zip(self.potential, self.flux)
        self.all = zip(self.potential, self.flux, self.possteppotential, self.posstepflux, self.negsteppotential, self.negstepflux)
        self.simend = time.time()
    def results(self):
        return self.all
    def time(self):
        return self.simend - self.waveforms
    

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
    
    k1_values = [0.1,1,10]
    k2_values = [0.1,1,10]
    k3_values = [0.1, 1, 10]
    number = 0
    
    shape = wf.SWV(Eini = 0, Efin = 0.5, dEs = 0.002, dEp = 0.05, dt = 0.04, pt = 0.02, sp = 350)

    for ix in k1_values:
        for iy in k2_values:
            for iz in k3_values:
                
                instance = EEC(input = shape, E0A = 0.25, E0B = 0.25, k1 = ix, k2 = iy, k3 = iz, a1 = 0.5, a2 = 0.5, cA = 0.005, cB = 0.000, cC = 0.000, cD = 0.000, DA = 5E-6, DB = 5E-6, DC = 5E-6, DD = 5E-6, r = 0.15, h = 1E-4, expansion = 1.05)
                filepath = cwd + '/data/' + 'test K1 ' + str(ix) + 'test K2 ' + str(iy)+ 'test K3 ' + str(iz) + '.txt'
                with open(filepath, 'w') as file:
                    for ix, iy, iz, jx, jy, jz in instance.results():
                        file.write(str(ix) + ',' + str(iy) + ',' + str(iz) +',' + str(jx) +',' + str(jy) +',' + str(jz) + '\n')
                number += 1
                print(f'Completed: {number} out of {len(k1_values) * len(k2_values) * len(k3_values)}. Remaining time estimated: {(len(k1_values) * len(k2_values) * len(k3_values) - number) * (instance.time())} seconds')
    end = time.time()
    print(end-start)