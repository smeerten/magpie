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
# along with magpie. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt

import sample

def diffMat(B, T1, T2):
    Bx, By, Bz = B
    return np.array([[-1/T2, Bz,  -By],
                     [-Bz,   -1/T2, Bx],
                     [By,   -Bx,  -1/T1]])

def diffEq(t, M, B, T1, T2, M0):
    mat = diffMat(B, T1, T2)
    dMdt = np.matmul(mat, M)
    dMdt[2] += M0/T1
    return dMdt
    

class Simulator():
    
    def __init__(self, pulseSeq=None, settings=None, sample=None):
        self.pulseSeq = pulseSeq
        self.settings = settings
        self.sample = sample
        
    def simulate(self, scans=1):
        """
        Simulate using the defined pulseSeq, settings, and sample.

        Parameters
        ----------
        scans : int, optional
            Number of scans for which the pulse sequence will be simulated.
            Default is 1.

        Returns
        -------
        Array
            The datapoints in the FID blocks of the pulseSeq.
        """
        allSpins = self.sample.expandSystems(self.settings['B0'], self.settings['observe'])
        FID = None
        for spinInfo in allSpins:
            Frequency, Intensity, T1, T2, T2prime = spinInfo
            spinFID = self.simulateSpin(Intensity, Frequency, T1, T2)
            if FID is None:
                FID = spinFID
            else:
                FID += spinFID
        return FID
                

    def simulateSpin(self, amp, freq, T1, T2, scans=1):
        M = np.array([0, 0, amp])
        fullResults = None
        for i in range(scans):
            M[0] = 0.0
            M[1] = 0.0
            scanResults = []
            for pulseStep in self.pulseSeq:
                if pulseStep[0] == 'pulse':
                    rf = pulseStep[2]
                    phase = pulseStep[3]
                    B = 2 * np.pi * np.array([rf*np.cos(phase), rf*np.sin(phase), freq])
                else:
                    B = 2 * np.pi * np.array([0, 0, freq])
                if pulseStep[0] == 'FID':
                    npoints = pulseStep[2]
                    t_eval = np.linspace(0, pulseStep[1], pulseStep[2])
                else:
                    t_eval = None
                sol = solve_ivp(diffEq, (0, pulseStep[1]), M, t_eval=t_eval, args=(B, T1, T2, amp))
                M = sol.y[:,-1]
                if pulseStep[0] == 'FID':
                    scanResults.append(sol.y[0] + 1j*sol.y[1])
            scanResults = np.concatenate(scanResults)
            if fullResults is None:
                fullResults = scanResults
            else:
                fullResults += scanResults
        return fullResults
            

if __name__ == '__main__':
    tube = sample.sample()
    tube.addMolecule([('1H',0,1),('1H',1,1),('13C',1,1)],np.array([[0,10,5],[10,0,0],[5,0,0]]),1,3,0.3,1)

    settings = {'B0':1.0, 'observe':'1H'}
    pulseSeq = [('delay', 1),
                ('pulse', 1e-6, 100e3, -np.pi/2.0),
                ('FID', 10, 2048)]
    
    sim = Simulator(pulseSeq, settings, tube)
    fid = sim.simulate()
    spectrum = np.fft.fftshift(np.fft.fft(fid))
    plt.plot(np.real(spectrum))
    plt.plot(np.imag(spectrum))
    plt.show()
