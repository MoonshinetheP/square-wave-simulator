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
    
Filename:           animate.py

===================================================================================================
How to use this file:
    

===================================================================================================
'''




import os
from errno import EEXIST
import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import waveforms as wf
import simulations as sim


if __name__ == '__main__':
    
    cwd = os.getcwd()

    try:
        os.makedirs(cwd + '/data')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/data'):
            pass
        else: 
            raise
    
    try:
        os.makedirs(cwd + '/animations')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/animations'):
            pass
        else: 
            raise
    
    '''SIMULATION'''
    start = time.time()
    
    shape = wf.CSV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.001, sr = 0.1, ns = 1, st = 0.001, detailed = True)
    instance = sim.E(input = shape, E0 = 0.25, k0 = 10, a = 0.5, cR = 0.005, cO = 0.000, DR = 5E-6, DO = 5E-6, r = 0.15, expansion = 1.05, Nernstian = False, BV = True, MH = False)
    
    end = time.time()
    print(f'The simulation took {end-start} seconds to complete')


    '''SAVE DATA'''
    filepath = f'{cwd}/data/{shape.type} {shape.subtype} 0.005.txt'
    with open(filepath, 'w') as file:
        for ix, iy, iz in instance.results():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')
    

    '''PLOT GENERATION'''
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 5))
    left, = ax1.plot([], [], linewidth = 1, linestyle = '-', color = 'blue', marker = None, label = None, visible = True)
    right, = ax2.plot([], [], linewidth = 1, linestyle = '-', color = 'red', marker = None, label = None, visible = True)


    '''PLOT SETTINGS'''
    ax1.set_xlim(np.amin(shape.tWF) - (0.1 * (np.amax(shape.tWF) - np.amin(shape.tWF))), np.amax(shape.tWF) + (0.1 * (np.amax(shape.tWF) - np.amin(shape.tWF))))
    ax1.set_ylim(np.amin(shape.EWF) - (0.1 * (np.amax(shape.EWF) - np.amin(shape.EWF))), np.amax(shape.EWF) + (0.1 * (np.amax(shape.EWF) - np.amin(shape.EWF))))
    ax1.set_title('E vs. t', pad = 15, fontsize = 20)
    ax1.set_xlabel('t / s', labelpad = 5, fontsize = 15)
    ax1.set_ylabel('E / V', labelpad = 5, fontsize = 15)

    ax2.set_xlim(np.amin(instance.EPLOT) - (0.1 * (np.amax(instance.EPLOT) - np.amin(instance.EPLOT))), np.amax(instance.EPLOT) + (0.1 * (np.amax(instance.EPLOT) - np.amin(instance.EPLOT))))
    ax2.set_ylim(np.amin(instance.flux) - (0.1 * (np.amax(instance.flux) - np.amin(instance.flux))), np.amax(instance.flux) + (0.1 * (np.amax(instance.flux) - np.amin(instance.flux)))) 
    ax2.set_title('i vs. E', pad = 15, fontsize = 20)
    ax2.set_xlabel('E / V', labelpad = 5, fontsize = 15)
    ax2.set_ylabel('i / A', labelpad = 5, fontsize = 15)


    '''ANIMATION FUNCTIONS'''
    def potential(i):
        left.set_data(shape.tWF[:i], shape.EWF[:i])
        return left, 

    def current(i):
        right.set_data(instance.EPLOT[:i], instance.flux[:i])
        return right,


    '''ANIMATION'''
    Evt = animation.FuncAnimation(fig, potential, frames = shape.indexWF.size, interval = 4, repeat = False, blit = True) 
    if shape.type == 'hybrid' and instance.detailed == False:
        ivE = animation.FuncAnimation(fig, current, frames = instance.EPLOT.size, interval = 10, repeat = False, blit = True)
    if shape.type == 'hybrid' and instance.detailed == True:
        ivE = animation.FuncAnimation(fig, current, frames = instance.EPLOT.size, interval = 4, repeat = False, blit = True)
    
    plt.show()
    plt.close()
    
    
    
    
   
    
 

