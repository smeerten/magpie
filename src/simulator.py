#!/usr/bin/env python3

# Copyright 2021 Bas van Meerten and Wouter Franssen

# This file is part of magpie.
#
# magpie is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# magpie is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ssNake. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from scipy.integrate import odeint


class Simulator():
    
    def __init__(self, pulseSeq=None, settings=None, sample=None):
        self.pulseSeq = pulseSeq
        self.settings = settings
        self.sample = sample
        
    def simulate(self, startingInten=None):
        """
        Simulate using the defined pulseSeq, settings, and sample.

        Parameters
        ----------
        startingInten : array_like, optional
            Starting intensities of the spins. Should have the same length as those returned by sample.
            If None, all spins will have relative intensity 1.

        Returns
        -------
        Array
            The datapoints in the FID blocks of the pulseSeq.
        """
        # spinInfo = self.sample.getSpins() # obtain a 2D array with for each spin [intensity, frequency, T1, T2]
        
        

