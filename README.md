# electrochemistry-simulations
Scripts for simulating various electrochemical techniques and mechanisms

### Description
The aim of this project is to develop an electrochemical simulation platform 
which can accomodate different waveforms and various kinetic models whilst 
also allowing for both simplified (sampled) and detailed (analogue) responses
to be calculated.

### Explanation
This project has been split into several modules which are used in sequence to 
prepare the simulation. The workflow looks a little like this:
    1.  waveforms.py is used to generate a range of waveforms including the main
        potential waveform used to perform the simulation, a plotting waveform used
        to plot the output of the simulation (often differs to the main wavefrom,
        especially for more detailed simulations), and an exported waveform for 
        visualising the applied potential
    2.  simulations.py takes the main waveform and uses it to calculate the current
        given a set of reaction parameters
    3.  plot.py can be used to visualise the waveform and the resulting current profile