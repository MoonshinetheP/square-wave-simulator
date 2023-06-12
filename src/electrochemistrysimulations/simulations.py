'''
===================================================================================================
Copyright (C) 2023 Steven Linfield

This file is part of the electrochemistry-simulations package. This package is free software: you 
can redistribute it and/or modify it under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, or (at your option) any later 
version. This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details. You should have received a copy of the GNU General 
Public License along with electrochemistry-simulations. If not, see https://www.gnu.org/licenses/
===================================================================================================

Package title:      electrochemistry-simulations
Repository:         https://github.com/MoonshinetheP/electrochemistry-simulations
Date of creation:   09/03/2023
Main author:        Steven Linfield (MoonshinetheP)
Collaborators:      None
Acknowledgements:   Oliver Rodriguez (oliverrdz)
    
Filename:           simulations.py

===================================================================================================
How to use this file:
    
    1.  Scroll down to the section titled 'RUNNING THE SIMULATION FROM MAIN' (near the bottom)
    2.  Describe the waveform using the appropriate class from the waveforms (wf) module
    3.  Define the parameters of the simulation in the E class
    4.  Run the file in Python
    
The programme will attempt to make a new folder in the current working directory to keep all the 
data. The name of this folder is /data by default, but this can be edited. If the /data folder 
already exists, then this programme will save data into this folder.
===================================================================================================
'''




'''MODULES'''
import sys
import os
import time

import numpy as np
import waveforms as wf
import capacitance as cap
import noise as noise

from errno import EEXIST
from scipy.sparse import diags as diagonals
from scipy.integrate import solve_ivp as solver


'''SIMULATION CLASS'''
class Diffusive:
    """Simulation of an E mechanism using solving from scipy \n
    E: R -> O + e"""

    def __init__(self, input, E0, k0, a, cR, cO, DR, DO, eqE, eqt, Cd, Ru, Nernstian, BV, MH, electrical, shot, thermal, r, expansion):
        '''
        ------------------
        WAVEFORM VARIABLES
        ------------------

        In this section, several variables are imported from the waveforms module through the input argument.
        Specific waveform types are defined in order to differentiate between variables that are available in
        some waveforms, but not others. For sweeps, the detailed variable is set to False by default. For steps,
        the detailed variable is set to True by default. For all other waveforms, this variable is set according
        to the value entered in the waveform generation.
        '''

        self.input = input                                      # The actual waveform object passed into the simulation class
        self.index = self.input.index                           # Index (i.e. the list of data points) derived from the waveform object
        self.t = self.input.t                                   # Time for simulations derived from the waveform object
        self.E = self.input.E                                   # Potential for simulations derived from the waveform object
        self.tPLOT = self.input.tPLOT                           # Time for plotting simulation data
        self.EPLOT = self.input.EPLOT                           # Potential for plotting simulation data


        if self.input.type == 'sweep':
            self.Eini = self.input.Eini                         # Initial potential for a sweep waveform
            self.Eupp = self.input.Eupp                         # Upper vertex potential for a sweep waveform (ignored in LSV if dE is negative)
            self.Elow = self.input.Elow                         # Lower vertex potential for a sweep waveform (ignored in LSV if dE is positive)
            self.dE = self.input.dE                             # Step potential for a sweep waveform (equates to number of data points)
            self.sr = self.input.sr                             # Scan rate for a sweep waveform
            self.ns = self.input.ns                             # Number of scans for a sweep waveform

            self.detailed = False                               # Whether the waveform should be detailed or not (set to False for a sweep waveform)


        if self.input.type == 'step':
            self.dE = self.input.dE                             # Applied potential for a step waveform
            self.dt = self.input.dt                             # Duration for which the potential is applied in a step waveform

            self.detailed = True                                # Whether the waveform should be detailed or not (set to True for a step waveform)


        if self.input.type == 'pulse':
            self.Eini = self.input.Eini                         # Initial potential for a pulse waveform
            self.Efin = self.input.Efin                         # Final potential for a pulse waveform
            self.dEs = self.input.dEs                           # Step potential for a pulse waveform
            self.dEp = self.input.dEp                           # Pulse amplitude (not peak to peak) for a pulse waveform
            self.pt = self.input.pt                             # Duration of the pulse in a pulse waveform
            self.rt = self.input.rt                             # Duration of the rest period in a pulse waveform
            self.st = self.input.st                             # Sampling time for a detailed pulse waveform (ignored if detailed == False)

            self.detailed = self.input.detailed                 # Whether the pulse waveform should be detailed or not


        if self.input.type == 'hybrid':
            self.Eini = self.input.Eini                         # Initial potential for a hybrid waveform
            self.Eupp = self.input.Eupp                         # Upper vertex potential for a hybrid waveform
            self.Elow = self.input.Elow                         # Lower vertex potential for a hybrid waveform
            self.dE = self.input.dE                             # Step potential for a hybrid waveform
            self.sr = self.input.sr                             # Scan rate for a hybrid waveform
            self.ns = self.input.ns                             # Number of scans for a hybrid waveform
            self.st = self.input.st                             # Sampling time for a detailed hybrid waveform (ignored if detailed == False)

            self.detailed = self.input.detailed                 # Whether the hybrid waveform should be detailed or not

        
        '''
        -------------------
        MECHANISM VARIABLES
        -------------------
            
        '''
        self.E0 = E0                                            # Standard redox potential (in V)
        self.k0 = k0                                            # 
        self.a = a                                              # Transfer coefficient (no units)
        self.cR = cR                                            # Concentration of the reduced species (in mol/cm^3)
        self.cO = cO                                            # Concentration of the oxidised species (in mol/cm^3)
        self.DR = DR                                            #
        self.DO = DO                                            #

        self.F = 96485                                          # Faraday constant (in C/mol)
        self.R = 8.314                                          # Ideal gas constant (in J/Kmol)
        self.Temp = 298                                         # Standar temperature (in K)
        
        


        '''
        -------------------
        KINETIC MODELS
        -------------------
            
        '''

        self.Nernstian = Nernstian                              #
        self.BV = BV                                            #
        self.MH = MH                                            #

        self.kinetics = 0
        for ix in [self.Nernstian, self.BV, self.MH]:
            if ix == True:
                self.kinetics += 1

        if self.kinetics == 0:
            print('\n' + 'No kinetic model was chosen' + '\n')
            sys.exit()

        if self.kinetics >= 2:
            print('\n' + 'More than one kinetic model was chosen' + '\n')
            sys.exit()


        '''
        -------------------
        NOISE MODELS
        -------------------
        
        '''
        self.electrical = electrical
        self.shot = shot
        self.thermal = thermal

        self.noise = 0
        for ix in [self.electrical, self.shot, self.thermal]:
            if ix == True:
                self.noise += 1


        '''
        -------------------
        SPATIAL GRID VARIABLES
        -------------------
            
        '''
        self.r = r                                              #
        self.expansion = expansion                              #

        
        '''
        -------------------
        DIMENSIONLESS VARIABLES
        -------------------
            
        '''
        concentrations = np.array([self.cR, self.cO])           #
        self.cmax = np.amax(concentrations)                     #
        self.CR = self.cR / self.cmax                           #
        self.CO = self.cO / self.cmax                           #

        diffusions = np.array([self.DR, self.DO])               #
        self.Dmax = np.amax(diffusions)                         #

        self.dR = self.DR / self.Dmax                           #
        self.dO = self.DO / self.Dmax                           #
        self.d = self.Dmax / self.Dmax                          #

        self.T = (self.Dmax * self.t) / (self.r ** 2)           # Dimensionless time
        self.dT = np.diff(self.T)                               #

        if self.input.type == 'pulse':
            self.pT = np.array([])
            self.rT = np.array([])
            for ix in range (0, self.dT.size):
                if ix % 2 == 0:
                    self.pT = np.append(self.pT, self.dT[ix])
                else:
                    self.rT = np.append(self.rT, self.dT[ix])
        
        if self.detailed == False:
            if self.input.type == 'sweep' or self.input.type == 'hybrid':
                self.sT = np.array([])
                self.sT = np.append(self.sT, np.amin(self.dT))
            
            if self.input.type == 'pulse':
                self.psT = np.empty((self.pT.size, 1))
                for ix in range(0, self.pT.size):
                    self.psT[ix,:] = np.amin(self.pT)
                
                self.rsT = np.empty((self.rT.size, 1))
                for iy in range(0, self.rT.size):
                    self.rsT[iy,:] = np.amin(self.rT)
                
        else:
            if self.input.type == 'pulse':
                self.psT = np.empty((self.pT.size, self.input.pp))
                for ix in range(0, self.pT.size):
                    self.psT[ix,:] = np.arange(0, np.amin(self.pT), (self.Dmax * self.st) / (self.r ** 2))
                
                self.rsT = np.empty((self.rT.size, self.input.rp))
                for iy in range(0, self.rT.size):
                    self.rsT[iy,:] = np.arange(0, np.amin(self.rT), (self.Dmax * self.st) / (self.r ** 2))

            if self.input.type == 'hybrid':
                self.sT = np.array([])
                self.sT = np.append(self.sT, np.arange(0, np.amin(self.dT), (self.Dmax * self.st) / (self.r ** 2)))
                if self.sT[-1] > self.dT[1]:
                    self.sT = self.sT[:-1]
        
        self.Tmax = self.T[-1]                                  #
        

        self.Xmax = 6 * np.sqrt(self.d * self.Tmax)             #

        self.theta = (self.F / (self.R * self.Temp)) * (self.E - self.E0)

        self.K0 = (self.k0 * self.r) / self.Dmax                #
        
        '''
        -------------------
        CELL CHARACTERISTICS
        -------------------
        
        '''
        self.eqE = eqE
        self.eqtheta = (self.F / (self.R * self.Temp)) * self.eqE
        
        self.eqt = eqt
        self.eqT = (self.Dmax * self.eqt) / (self.r ** 2)         
        
        self.Cd = Cd
        self.Ru = Ru
        '''
        -------------------
        EXPANDING SPATIAL GRID
        -------------------
            
        '''                
        self.dX = np.sqrt(2.05 * self.dT[0])                    #
        self.x = np.array([0])                                  #
        while self.x[-1] < self.Xmax:
            self.x = np.append(self.x, self.x[-1] + self.dX)
            self.dX *= self.expansion
        

        '''
        -------------------
        DIFFUSION MATRICES
        -------------------
            
        '''
        self.n = int(self.x.size)                               # 
        self.m = int(self.theta.size)                           #

        self.C_R = np.ones((self.n, self.m)) * self.CR          #
        self.C_O = np.ones((self.n, self.m)) * self.CO          #

        self.alpha_R = np.ones(self.n - 1)                      #
        self.beta_R = np.ones(self.n)                           #
        self.gamma_R = np.ones(self.n - 1)                      #
        
        self.alpha_O = np.ones(self.n - 1)                      #
        self.beta_O = np.ones(self.n)                           #
        self.gamma_O = np.ones(self.n - 1)                      #    
          
        for ix in range(1, self.n - 1):
            self.xplus = self.x[ix + 1] - self.x[ix]            #
            self.xminus = self.x[ix] - self.x[ix - 1]           #
            self.denominator = self.xminus * (self.xplus ** 2) + self.xplus * (self.xminus **2)
            
            self.alpha_R[ix - 1] *= 2 * self.dR * self.xplus / self.denominator
            self.beta_R[ix] *= -2 * self.dR * (self.xminus + self.xplus) / self.denominator
            self.gamma_R[ix] *= 2 * self.dR * self.xminus / self.denominator
            
            self.alpha_O[ix - 1] *= 2 * self.dO * self.xplus / self.denominator
            self.beta_O[ix] *= -2 * self.dO * (self.xminus + self.xplus) / self.denominator
            self.gamma_O[ix] *= 2 * self.dO * self.xminus / self.denominator
        
            

        R = diagonals([self.alpha_R, self.beta_R, self.gamma_R], [-1,0,1]).toarray()
        R[0,:] = np.zeros(self.n)
        R[0,0] = self.CR    
        
        O = diagonals([self.alpha_O, self.beta_O, self.gamma_O], [-1,0,1]).toarray()
        O[0,:] = np.zeros(self.n)
        O[0,0] = self.CO   

        def reduced(t,y):
            return np.dot(R,y)
        
        def oxidised(t,y):
            return np.dot(O,y)

        '''
        -------------------
        CONTAINERS
        -------------------
            
        '''
        self.flux = np.array([])
        self.pulseflux = np.array([])
        self.restflux = np.array([])

        '''
        -------------------
        MAIN LOOP
        -------------------
            
        '''
        if self.eqT > 0:
            if self.Nernstian == True:
                '''Nernstian'''
                self.C_R[0, 1] = (self.C_R[1, 0] + (self.dR/self.dO) * self.C_O[1, 0])/(1 + (self.dR/self.dO) * np.exp(self.eqtheta))

                self.C_O[0, 1] = (self.C_R[1, 0] + (self.dR/self.dO) * self.C_O[1, 0])/((self.dR/self.dO) + np.exp(-self.eqtheta))
                
            if self.BV == True:
                '''Butler-Volmer'''
                self.C_R[0, 1] = (-self.C_R[1, 0] + (self.x[1] - self.x[0]) * self.K0 * np.exp(-self.a * self.eqtheta) * (self.C_O[1, 0] + (self.dR/self.dO) * self.C_R[1, 0]))/((self.x[1] - self.x[0]) * self.K0 * (np.exp((1 - self.a) * self.eqtheta) + (self.dR/self.dO) * np.exp((-self.a) * self.eqtheta)) - 1)

                self.C_O[0, 1] = (-self.C_O[1, 0] + (self.x[1] - self.x[0]) * self.K0 * np.exp((1 - self.a) * self.eqtheta) * (self.C_R[1, 0] + (self.dO/self.dR) * self.C_O[1, 0]))/((self.x[1] - self.x[0]) * self.K0 * (np.exp(-self.a * self.eqtheta) + (self.dO/self.dR) * np.exp((1 - self.a) * self.eqtheta)) - 1)
            
            oxidation = solver(reduced, [0, self.eqT], self.C_R[:,0], t_eval=[self.eqT], method='RK45')
            self.C_R[1:-1, 1] = oxidation.y[1:-1, -1]
        
            reduction = solver(oxidised, [0, self.eqT], self.C_O[:,0], t_eval=[self.eqT], method='RK45')
            self.C_O[1:-1, 1] = reduction.y[1:-1, -1]
        else:
            pass


        for k in range(1,self.m):
            try:
                self.offset = self.flux[-1] * self.Ru
                #self.theta[k - 1] = self.theta[k - 1] - (self.F / (self.R * self.Temp)) * (self.offset)
            except:
                pass
            
            '''Boundary conditions'''
            if self.Nernstian == True:
                '''Nernstian'''
                self.C_R[0, k] = (self.C_R[1, k - 1] + (self.dR/self.dO) * self.C_O[1, k - 1])/(1 + (self.dR/self.dO) * np.exp(self.theta[k-1]))

                self.C_O[0, k] = (self.C_R[1, k - 1] + (self.dR/self.dO) * self.C_O[1, k - 1])/((self.dR/self.dO) + np.exp(-self.theta[k-1]))
            
            if self.BV == True:
                '''Butler-Volmer'''
                self.C_R[0, k] = (-self.C_R[1, k - 1] + (self.x[1] - self.x[0]) * self.K0 * np.exp(-self.a * self.theta[k - 1]) * (self.C_O[1, k - 1] + (self.dR/self.dO) * self.C_R[1, k - 1]))/((self.x[1] - self.x[0]) * self.K0 * (np.exp((1 - self.a) * self.theta[k - 1]) + (self.dR/self.dO) * np.exp((-self.a) * self.theta[k - 1])) - 1)

                self.C_O[0, k] = (-self.C_O[1, k - 1] + (self.x[1] - self.x[0]) * self.K0 * np.exp((1 - self.a) * self.theta[k - 1]) * (self.C_R[1, k - 1] + (self.dO/self.dR) * self.C_O[1, k - 1]))/((self.x[1] - self.x[0]) * self.K0 * (np.exp(-self.a * self.theta[k - 1]) + (self.dO/self.dR) * np.exp((1 - self.a) * self.theta[k - 1])) - 1)
           
            if self.input.type == 'sweep' or self.input.type == 'hybrid':
                oxidation = solver(reduced, [0, self.dT[k - 1]], self.C_R[:,k - 1], t_eval=self.sT, method='RK45')
                self.C_R[1:-1, k] = oxidation.y[1:-1, -1]
                if self.detailed == True:
                    self.C_Rdet = oxidation.y[1:-1, :]
                
                reduction = solver(oxidised, [0, self.dT[k - 1]], self.C_O[:,k - 1], t_eval=self.sT, method='RK45')
                self.C_O[1:-1, k] = reduction.y[1:-1, -1]
                if self.detailed == True:
                    self.C_Odet = reduction.y[1:-1, :]

            if self.input.type == 'pulse':
                if k % 2 != 0:
                    oxidation = solver(reduced, [0, self.dT[k - 1]], self.C_R[:,k - 1], t_eval=self.psT[(k-1)//2], method='RK45')
                if k % 2 == 0:
                    oxidation = solver(reduced, [0, self.dT[k - 1]], self.C_R[:,k - 1], t_eval=self.rsT[(k-2)//2], method='RK45')
            
                self.C_R[1:-1, k] = oxidation.y[1:-1, -1]

                if self.detailed == True:
                    self.C_Rdet = oxidation.y[1:-1, :]
           
                if k % 2 != 0:
                    reduction = solver(oxidised, [0, self.dT[k - 1]], self.C_O[:,k - 1], t_eval=self.psT[(k-1)//2], method='RK45')
                if k % 2 == 0:
                    reduction = solver(oxidised, [0, self.dT[k - 1]], self.C_O[:,k - 1], t_eval=self.rsT[(k-2)//2], method='RK45')
                
                self.C_O[1:-1, k] = reduction.y[1:-1, -1]

                if self.detailed == True:
                    self.C_Odet = reduction.y[1:-1, :]
            
            if self.input.type == 'sweep':
                self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_R[1, k] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_O[1, k] - self.C_O[0, k]) / (self.x[1] - self.x[0])))

            if self.input.type == 'pulse':
                if self.detailed == False:
                    if k % 2 != 0:
                        self.pulseflux = np.append(self.pulseflux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_R[1, k] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_O[1, k] - self.C_O[0, k]) / (self.x[1] - self.x[0])))
                    if k % 2 == 0:
                        self.restflux = np.append(self.restflux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_R[1, k] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_O[1, k] - self.C_O[0, k]) / (self.x[1] - self.x[0])))
                if self.detailed == True:
                    if k % 2 != 0:
                        for ix in range(0, self.input.pp):
                            self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_Rdet[0, ix] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_Odet[0, ix] - self.C_O[0, k]) / (self.x[1] - self.x[0])))
                    if k % 2 == 0:
                        for ix in range(0, self.input.rp):
                            self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_Rdet[0, ix] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_Odet[0, ix] - self.C_O[0, k]) / (self.x[1] - self.x[0])))

            if self.input.type == 'hybrid':
                if self.detailed == False:
                    self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_R[1, k] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_O[1, k] - self.C_O[0, k]) / (self.x[1] - self.x[0])))
                if self.detailed == True:
                    if self.input.sampled == False:
                        for ix in range(0, self.input.sp):
                            self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cR * self.DR) * ((self.C_Rdet[0, ix] - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((self.C_Odet[0, ix] - self.C_O[0, k]) / (self.x[1] - self.x[0])))
                    if self.input.sampled == True:
                        self.flux = np.append(self.flux, (self.F * np.pi * self.r * self.cR * self.DR) * ((np.average(self.C_Rdet[0, self.input.sp - round(self.input.sp * (1 - self.input.alpha)):]) - self.C_R[0, k]) / (self.x[1] - self.x[0])) - (self.F * np.pi * self.r * self.cO * self.DO) * ((np.average(self.C_Odet[0, self.input.sp - round(self.input.sp * (1 - self.input.alpha)):]) - self.C_O[0, k]) / (self.x[1] - self.x[0])))
       
        
        '''
        -------------------
        POST SIMULATION CALCULATIONS
        -------------------
            
        '''
        self.fluxcapacitance = cap.Capacitance(input = self.input, Cd = self.Cd, Ru = self.Ru)
        self.fluxnoise = noise.Noise(input = self.input, electrical = self.electrical, shot = self.shot, thermal = self.thermal)
        if self.input.type == 'pulse':
            if self.detailed == False:
                self.flux = self.pulseflux - self.restflux
            
            else:
                pass
        
        '''
        -------------------
        OUTPUTS
        -------------------
            
        '''
        self.flux += self.fluxcapacitance.i
        self.flux += self.fluxnoise.i
        self.output = zip(self.tPLOT, self.EPLOT, self.flux)
        self.concentrations = zip(self.C_O, self.C_R)
    

    '''
    -------------------
    CLASS FUNCTIONS
    -------------------
            
    '''
    def results(self):
        return self.output
    
    def profiles(self):
        return self.concentrations
    




'''RUNNING THE SIMULATION FROM MAIN'''

if __name__ == '__main__':
    
    '''1. MAKE A /DATA FOLDER'''
    cwd = os.getcwd()    
    try:
        os.makedirs(cwd + '/data')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/data'):
            pass
        else: 
            raise
    

    '''2. DEFINE THE START TIME'''
    start = time.time()


    '''3. DESCRIBE THE WAVEFORM'''
    '''Sweeps'''
    #shape = wf.LSV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1)
    #shape = wf.CV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1)
    
    '''STEPS'''
    #shape = wf.CA(dE = [0.5], dt = [1], st = 0.001)
    
    '''PULSES'''
    #shape = wf.DPV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.02, pt = 0.05, rt = 0.15, st = 0.001, detailed = True, sampled = True, alpha = 0.25)
    #shape = wf.SWV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.02, pt = 0.1, rt = 0.1, st = 0.001, detailed = True, sampled = True, alpha = 0.25)
    #shape = wf.NPV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.02, pt = 0.05, rt = 0.15, st = 0.001, detailed = True, sampled = True, alpha = 0.25)
    
    '''HYBRID'''
    shape = wf.CSV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1, st = 0.001, detailed = True, sampled = True, alpha = 0.05)
    #shape = wf.AC(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1, st = 0.001, detailed = True, sampled = True, alpha = 0.25)
    

    '''4. RUN THE SIMULATION'''
    instance = Diffusive(input = shape, E0 = 0.25, k0 = 10, a = 0.5, cR = 0.000005, cO = 0.000000, DR = 5E-6, DO = 5E-6, eqE = 0, eqt = 5, Cd = 0.000050, Ru = 200, Nernstian = False, BV = True, MH = False, electrical = False, shot = False, thermal = False, r = 0.1, expansion = 1.05)
    
    '''5. DEFINE THE END TIME'''
    end = time.time()
    print(end-start)

    '''6. SAVE THE DATA'''
    filepath = cwd + '/data/' + f'{shape.subtype}' + '.txt'
    with open(filepath, 'w') as file:
        for ix, iy, iz in instance.results():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')