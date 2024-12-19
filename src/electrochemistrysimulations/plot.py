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
    
Filename:           plot.py

===================================================================================================
How to use this file:
    

===================================================================================================
'''




import os
import time
import collections as col
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
from matplotlib.widgets import Button, RadioButtons, Slider, TextBox

class Plotter:
    '''Plots potential waveforms and analysed oscilloscope data generated from other files in the package\n
    
    Requires:\n
'''

    def __init__(self, shape, data, display = True, save = True, waveform = True, plot = True, profile = True, map = True):

        '''PARAMETER INITIALISATION'''
        self.shape = shape      #
        self.data = data        #
        self.display = display      #
        self.save = save        #
        self.waveform = waveform
        self.plot = plot
        self.profile = profile
        self.map = map
        
        self.mosaic = [
            ['A','A','A','A','B','B','B','B','C','C','C','C','C','C'],
            ['A','A','A','A','B','B','B','B','C','C','C','C','C','C'],
            ['A','A','A','A','B','B','B','B','C','C','C','C','C','C'],
            ['A','A','A','A','B','B','B','B','C','C','C','C','C','C'],
            ['D','D','D','D','D','D','D','D','C','C','C','C','C','C'],
            ['D','D','D','D','D','D','D','D','C','C','C','C','C','C'],
            ['D','D','D','D','D','D','D','D','E','F','G','H','H','I'],
            ['D','D','D','D','D','D','D','D','J','J','J','K','K','I']
                  ]
        
        '''PLOT DEFINITION'''
        fig, axes = plt.subplot_mosaic((self.mosaic), figsize=(16, 10))        # defines a matplotlib figure with two horizontally arranged subplots
        
        if waveform == True:       
            left, = axes['A'].plot(self.shape.tWF, self.shape.EWF, linewidth = 1, figsize=(4, 4), linestyle = '-', color = 'blue', marker = None, label = None, visible = True)       # plots the potential waveform from waveforms.py on the left-hand subplot
            axes['A'].set_xlim(np.amin(self.shape.tWF) - (0.1 * (np.amax(self.shape.tWF) - np.amin(self.shape.tWF))), np.amax(self.shape.tWF) + (0.1 * (np.amax(self.shape.tWF) - np.amin(self.shape.tWF))))      # sets the x-axis limits of the left-hand subplot to +/- 10% of the waveform's time range
            axes['A'].set_ylim(np.amin(self.shape.EWF) - (0.1 * (np.amax(self.shape.EWF) - np.amin(self.shape.EWF))), np.amax(self.shape.EWF) + (0.1 * (np.amax(self.shape.EWF) - np.amin(self.shape.EWF))))      # sets the y-axis limits of the left-hand subplot to +/- 10% of the waveform's potential range
            axes['A'].set_title('E vs. t', pad = 15, fontsize = 20)       # defines the title and settings of the left-hand subplot
            axes['A'].set_xlabel('t / s', labelpad = 5, fontsize = 15)        # defines the x-axis label and settings of the left-hand subplot
            axes['A'].set_ylabel('E / V', labelpad = 5, fontsize = 15)        # defines the y-axis labe and settings of the left-hand subplot
        
        if plot == True:
            right, = axes['B'].plot(self.data.E[1:], self.data.i, linewidth = 1, linestyle = '-', color = 'red', marker = None, label = None, visible = True)     # plots the oscilloscope data from operations.py on the right-hand subplot
            axes['B'].set_xlim(np.amin(self.data.E) - (0.1 * (np.amax(self.data.E) - np.amin(self.data.E))), np.amax(self.data.E) + (0.1 * (np.amax(self.data.E) - np.amin(self.data.E))))        # sets the x-axis of the right-hand subplot to +/- 10% of the oscilloscope data's potential range
            axes['B'].set_ylim(np.nanmin(self.data.i) - (0.1 * (np.nanmax(self.data.i) - np.nanmin(self.data.i))), np.nanmax(self.data.i) + (0.1 * (np.nanmax(self.data.i) - np.nanmin(self.data.i))))        # sets the y-axis of the right-hand subplot to +/- 10% of the oscilloscope data's current range
            axes['B'].set_title('i vs. E', pad = 15, fontsize = 20)       # defines the title and settings of the right-hand subplot
            axes['B'].set_xlabel('E / V', labelpad = 5, fontsize = 15)        # defines the x-axis label and settings of the right-hand subplot 
            axes['B'].set_ylabel('i / A', labelpad = 5, fontsize = 15)        # defines the y-axis label and settings of the right-hand subplot
        
        cprofilelist = []
        markerlist = []
        colours = mcolors.TABLEAU_COLORS
        colourkeys = list(colours.keys())

        for n in range(1, len(self.data.mechanism.markers) + 1):
            cprofilelist.append(f'bottom{n}')
            markerlist.append(self.data.mechanism.markers[n - 1]['Species'])
        self.cprofile = col.namedtuple('Lines', cprofilelist)._make({} for _ in cprofilelist)
        for n in range(0, len(self.cprofile)):
            self.cprofile[n].update({'Axes': n})
            self.cprofile[n].update({'Line': 0})

        for n in self.cprofile:
            n['Line'], = axes['D'].plot(self.data.x, self.data.mechanism.markers[n['Axes']]['Concentration array'][:,1], linewidth = 1, linestyle = '-', color = colours[colourkeys[n['Axes']]], marker = None, label = markerlist[n['Axes']], visible = True)
        
        axes['D'].set_xlim(0, self.data.x[-1])
        axes['D'].set_ylim(-0.1 , 1.1)        # sets the y-axis of the right-hand subplot to +/- 10% of the oscilloscope data's current range
        axes['D'].legend()
        
        bw = Button(axes['E'], '<<')

        pp = Button(axes['F'], '>')

        fw = Button(axes['G'], '>>')

        speed = Slider(axes['H'], '', 0, 2, valinit = 1, valstep = 0.25)

        selector = RadioButtons(
            axes['I'], activecolor = 'red', labels = markerlist
        )


        samp = Slider(
        axes['J'], '', 0, self.data.m,
        valinit = 0, valstep=np.arange(0, self.data.m, 1),
        color="green"
        )

        sp = TextBox(axes['K'], 'Stop points')
        
        def play(val):
            pass
        def pause(val):
            pass
        def rewind(val):
            pass
        def forward(val):
            pass
        def speed(val):
            pass

        def update(val):
            amp = samp.val
            for n in self.cprofile:
                n['Line'].set_ydata(self.data.mechanism.markers[n['Axes']]['Concentration array'][:,int(amp)])

            #axes['C'].set_ylim(np.amin(np.array([self.data.C_R[:,int(amp)], self.data.C_O[:,int(amp)]])) - (0.1 * (np.amax(np.array([self.data.C_R[:,int(amp)], self.data.C_O[:,int(amp)]])) - np.amin(np.array([self.data.C_R[:,int(amp)], self.data.C_O[:,int(amp)]]))))), np.amax(np.array([self.data.C_R[:,int(amp)], self.data.C_O[:,int(amp)]])) + (0.1 * (np.amax(np.array([self.data.C_R[:,int(amp)], self.data.C_O[:,int(amp)]])) - np.amin(np.array([self.data.C_R[:,int(amp)], self.data.C_O[:,int(amp)]]))))

            fig.canvas.draw_idle()


        samp.on_changed(update)

        '''PLOT GENERATION & MANAGEMENT '''
        if self.save == True:
            plt.savefig(f'{os.getcwd()}/plots/{time.strftime("%Y-%m-%d %H-%M-%S")} {shape.subtype} data.png')      # saves the figure as a .png image
        
        if self.display == True:
            plt.show()      # displays the figure

        self.animate == False
        if self.animate == True:
        
            '''ANIMATION FUNCTIONS'''
            def potential(i):
                left.set_data(shape.tWF[:i], shape.EWF[:i])
                return left, 

            def current(i):
                right.set_data(shape.EPLOT[:i], data.i[:i])
                return right,


            '''ANIMATION'''
            Evt = animation.FuncAnimation(fig, potential, frames = shape.indexWF.size, interval = 4, repeat = False, blit = True) 
            ivE = animation.FuncAnimation(fig, current, frames = data.EPLOT.size, interval = 10, repeat = False, blit = True)
            plt.show()

        '''PLOT CLOSURE'''
        plt.close()     # closes the figure