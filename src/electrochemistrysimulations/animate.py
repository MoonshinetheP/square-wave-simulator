import os
from errno import EEXIST
import time

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import waveforms as wf
import E


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
    
    shape = wf.CV(Eini = 0, Eupp = 0.5, Elow = 0, dE = 0.005, sr = 0.1, ns = 1)
    #shape = wf.DPV(Eini = 0, Efin = 0.5, dEs = 0.005, dEp = 0.01, pt = 0.01, rt = 0.03, st = 0.005, detailed = False)
    instance = E.E(input = shape, E0 = 0.25, k0 = 0.1, a = 0.5, cR = 0.005, cO = 0.000, DR = 5E-6, DO = 5E-6, r = 0.15, expansion = 1.05, Nernstian = False, BV = True, MH = False)

    filepath = f'{cwd}/data/{shape.type} {shape.subtype}.txt'
    with open(filepath, 'w') as file:
        for ix, iy, iz in instance.results():
            file.write(str(ix) + ',' + str(iy) + ',' + str(iz) + '\n')
    
    end = time.time()
    print(f'The simulation took {end-start} seconds to complete')


    '''ANIMATION'''
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12, 5))

    def potential(i):
        ax1.plot(shape.tWF[:i], shape.EWF[:i], linewidth = 1, linestyle = '-', color = 'red', marker = None, label = None, visible = True)
        ax1.set_xlim(np.amin(shape.tWF) - (0.1 * (np.amax(shape.tWF) - np.amin(shape.tWF))), np.amax(shape.tWF) + (0.1 * (np.amax(shape.tWF) - np.amin(shape.tWF))))
        ax1.set_ylim(np.amin(shape.EWF) - (0.1 * (np.amax(shape.EWF) - np.amin(shape.EWF))), np.amax(shape.EWF) + (0.1 * (np.amax(shape.EWF) - np.amin(shape.EWF))))
        ax1.set_title('E vs. t', pad = 15, fontsize = 20)
        ax1.set_xlabel('t / s', labelpad = 5, fontsize = 15)
        ax1.set_ylabel('E / V', labelpad = 5, fontsize = 15)

    def current(i):
        ax2.plot(instance.EPLOT[:i], instance.flux[:i], linewidth = 1, linestyle = '-', color = 'blue', marker = None, label = None, visible = True)
        ax2.set_xlim(np.amin(instance.EPLOT) - (0.1 * (np.amax(instance.EPLOT) - np.amin(instance.EPLOT))), np.amax(instance.EPLOT) + (0.1 * (np.amax(instance.EPLOT) - np.amin(instance.EPLOT))))
        ax2.set_ylim(np.amin(instance.flux) - (0.1 * (np.amax(instance.flux) - np.amin(instance.flux))), np.amax(instance.flux) + (0.1 * (np.amax(instance.flux) - np.amin(instance.flux)))) 
        ax2.set_title('i vs. E', pad = 15, fontsize = 20)
        ax2.set_xlabel('E / V', labelpad = 5, fontsize = 15)
        ax2.set_ylabel('i / A', labelpad = 5, fontsize = 15)

    Evt = animation.FuncAnimation(fig, potential, frames = shape.indexWF.size, interval = 1, repeat = False, blit = True) 
    ivE = animation.FuncAnimation(fig, current, frames = instance.tPLOT.size, interval = 1, repeat = False, blit = True)
    
    try:
        os.makedirs(cwd + '/animations')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/animations'):
            pass
        else: 
            raise
    
    plt.show()
    filepath = f'{cwd}/animations/{shape.type} {shape.subtype} i vs. E.gif'
    #ivE.save(filepath, writer = 'imagemagick', fps = 30)
    
    
    plt.close()
    
    
    
    
   
    
 

