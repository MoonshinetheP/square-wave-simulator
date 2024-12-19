import numpy as np
from scipy.integrate import solve_ivp as solver

class BoundaryCondition:
    """Placing it here, although it could conceivably be done in Solver. In general, calling from Solver itself, but could also be called from simulations.py in order to cut down on arguments to Solver (maybe)"""
    def __init__(self, k, mechanism, species, geometry, spatial, model):
        
        self.k = k
        self.mechanism = mechanism
        self.geometry = geometry
        self.spatial = spatial
        self.x = self.spatial[0]
        self.model = model
        self.species = species
        
        self.BC = 0

        '''Butler-Volmer'''
        if len(self.mechanism['Oxidised']) != 0:
            for ix in self.species:
                if ix['Species'] == self.mechanism['Oxidised'][0][1][0]:
                    self.comp = ix['Concentration array'][1, k - 1]
                    self.dOther = ix['Diffusion coefficient']
            self.theta = self.mechanism['Oxidised'][1]
            self.K0 = self.mechanism['Oxidised'][2]
            self.a = self.mechanism['Oxidised'][3]

            self.BC += (self.mechanism['Concentration array'][1, k - 1] + self.x[1] * self.K0 * np.exp(-self.a * self.theta[k - 1]) * (self.comp + (self.mechanism['Diffusion coefficient']/self.dOther) * self.mechanism['Concentration array'][1, k - 1]))/(self.x[1] *  self.K0 * (np.exp((1 - self.a) * self.theta[k - 1]) + (self.mechanism['Diffusion coefficient']/self.dOther) * np.exp((-self.a) * self.theta[k - 1])) + 1)

               
        if len(self.mechanism['Reduced']) != 0:
            for ix in self.species:
                if ix['Species'] == self.mechanism['Reduced'][0][1][0]:
                    self.comp = ix['Concentration array'][1, k - 1] 
                    self.dOther = ix['Diffusion coefficient']
            self.theta = self.mechanism['Reduced'][1]
            self.K0 = self.mechanism['Reduced'][2]
            self.a = self.mechanism['Reduced'][3]

            self.BC += (self.mechanism['Concentration array'][1, k - 1] + self.x[1] * self.K0 * np.exp((1-self.a) * self.theta[k - 1]) * (self.mechanism['Concentration array'][1, k - 1] + (self.dOther/self.mechanism['Diffusion coefficient']) * self.comp))/(self.x[1] *  self.K0 * (np.exp((1 - self.a) * self.theta[k - 1]) + (self.mechanism['Diffusion coefficient']/self.dOther) * np.exp((-self.a) * self.theta[k - 1])) + 1)

        if len(self.mechanism['Consumed']) != 0 or len(self.mechanism['Produced']) != 0:
            if self.BC == 0:
                self.BC += self.mechanism['Concentration array'][0, k -1] #need to check concentration profiles to see if this worked
            else: pass
        
        if self.BC < 0:
            self.BC = 0
        #if self.BC > 1:
            #self.BC = 1
        #still need to import K0 theta for each reaction

class Solver:
    """Needs: current timestep, dT information, concnetration matrix for species of interest and complementary species, diagonal arrays, theta and K0 for reaction. Also need to import geometry and kinetic model at some point. Needs to output concentration and current for that timestep. ALso need to know if I have oxidation or reduction. Also need to enact different schemes depending on geometry, two dimensions requires solving twice, at least in solve_ivp."""
    def __init__(self, k, solvetime, mechanism, species, geometry, spatial, model):
        
        #timestep
        #dT and sT
        #species (oxidation)
        #all concentration arrays (how can I predict which ones will be needed? No clear way to update Reaction class with arrays, which are made after Reactions is called)
        #geometry
        #kinetics
        #
        self.k = k
        self.dT = solvetime[0]
        self.sT = solvetime[1]
        self.mechanism = mechanism
        self.species = species
        self.geometry = geometry
        self.spatial = spatial
        self.model = model

        array = self.mechanism['Matrix']

              
        def update(t,y):
            return np.dot(array,y)
            
        self.mechanism['Concentration array'][0,k] = BoundaryCondition(self.k, self.mechanism, self.species, self.geometry, self.spatial, self.model).BC
        #somehow need to think about how to import information regarding related species, because need complimentary diffusion coefficient and concentration
        #also need to import time. maybe should try boundary condition in separate class?
        #the idea from martin was grouping each array in one, then dotting it to produce 4,n array (although I still get stuck on how to input this into the solve_ivp, unless I reshape). he also suggested matlab specific solvers pdepd or something like that. also when it comes to solver, it supposedly uses t-1 data and transposes it to new time then solves until correct
        solved = solver(update, [0, self.dT[k - 1]], self.mechanism['Concentration array'][:,k - 1], t_eval=self.sT, method='RK45')
        self.mechanism['Concentration array'][1:-1, k] = solved.y[1:-1, -1]
        
        #for ix in range(0, len(self.spatial[0])):
        #    if self.mechanism['Concentration array'][ix,k] > 1:
        #        self.mechanism['Concentration array'][ix,k] = 1
                
        if len(self.mechanism['Consumed']) != 0:
            for iw in self.species:
                if iw['Species'] == self.mechanism['Consumed'][0][1][0]:
                    for pos in range(0, len(self.spatial[0])):
                        self.mechanism['Concentration array'][pos, k] -= self.mechanism['Concentration array'][pos, self.k-1] * self.mechanism['Consumed'][1] # iw['Concentration array'][pos, self.k-1] *
                        if self.mechanism['Concentration array'][pos, k] < 0:
                            self.mechanism['Concentration array'][pos, k] = 0
                            
        if len(self.mechanism['Produced']) != 0:
            for iw in self.species:
                if iw['Species'] == self.mechanism['Produced'][0][1][0]:
                    for pos in range(0, len(self.spatial[0])):
                        self.mechanism['Concentration array'][pos, k] += iw['Concentration array'][pos, self.k-1] * self.mechanism['Produced'][1]
