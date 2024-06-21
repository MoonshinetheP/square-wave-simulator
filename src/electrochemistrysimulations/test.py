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
    
Filename:           test_waveforms.py

===================================================================================================
How to use this file:
    

===================================================================================================
'''

import unittest
import waveforms

class TestCV(unittest.TestCase):
    def setUp(self):
        CV1 = waveforms.CV(Eini = 0.0, Eupp = 0.5, Elow = 0.0, dE = 0.001, sr = 0.1, ns  = 1)
        CV2 = waveforms.CV(Eini = 0.2, Eupp = 0.8, Elow = 0.2, dE = 0.002, sr = 0.05, ns  = 2)
        CV3 = waveforms.CV(Eini = -0.3, Eupp = 0.3, Elow = -0.3, dE = 0.001, sr = 0.1, ns  = 1)
        CV4 = waveforms.CV(Eini = 0.5, Eupp = 0.5, Elow = -0.1, dE = -0.001, sr = 0.2, ns  = 1)
        
        self.CV = [CV1, CV2, CV3, CV4]

        self.CV_window = [0.5, 0.6, 0.6, 0.6]
        self.CV_dp = [500, 300, 600, 600]    
        self.CV_tmax = [10, 48, 12, 6]
        self.CV_dt = [0.01, 0.04, 0.01, 0.005]

        self.CV_index = [1001, 1201, 1201, 1201]
        self.CV_t = [1001, 1201, 1201, 1201]
        self.CV_E = [1001, 1201, 1201, 1201]


    def test_CV_window(self):
        for ix in range(0, len(self.CV)):
            self.assertEqual(self.CV[ix].window, self.CV_window[ix], msg = f'CV {ix + 1}')

    def test_CV_dp(self):
        for ix in range(0, len(self.CV)):
            self.assertEqual(self.CV[ix].dp, self.CV_dp[ix], msg = f'CV {ix + 1}')

    def test_CV_tmax(self):
        for ix in range(0, len(self.CV)):
            self.assertEqual(self.CV[ix].tmax, self.CV_tmax[ix], msg = f'CV {ix + 1}')

    def test_CV_dt(self):
        for ix in range(0, len(self.CV)):
            self.assertEqual(self.CV[ix].dt, self.CV_dt[ix], msg = f'CV {ix + 1}')
    
    def test_CV_index(self):
        for ix in range(0, len(self.CV)):
            self.assertEqual(self.CV[ix].index.size, self.CV_index[ix], msg = f'CV {ix + 1}')
    
    def test_CV_t(self):
        for ix in range(0, len(self.CV)):
            self.assertEqual(self.CV[ix].t.size, self.CV_t[ix], msg = f'CV {ix + 1}')

    def test_CV_E(self):
        for ix in range(0, len(self.CV)):
            self.assertEqual(self.CV[ix].index.size, self.CV_index[ix], msg = f'CV {ix + 1}')


if __name__ == '__main__':
    
    unittest.main()