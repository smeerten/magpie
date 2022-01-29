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
import pandas as pd
from scipy.integrate import solve_ivp
import ast
import matplotlib.pyplot as plt
import scipy.io
import time

import sample
import helperFunctions as helpFie


BITS = 8

NVALS = 2**BITS
HRANGE = NVALS // 2 
DIGITPOINTS = np.arange(-HRANGE+1, HRANGE+1, dtype=float)
DIGITBINS = DIGITPOINTS[:-1] + 0.5
DIGITPOINTS /= HRANGE
DIGITBINS /= HRANGE


def ADC(fid):
    fidReal = DIGITPOINTS[np.digitize(np.real(fid), DIGITBINS)]
    fidImag = DIGITPOINTS[np.digitize(np.imag(fid), DIGITBINS)]
    return fidReal + 1j*fidImag

def diffEq(t, M, B, T1, T2, M0):
    Bx, By, Bz = B
    
    mat = np.array([[-1/T2, Bz,  -By],
                    [-Bz,   -1/T2, Bx],
                    [By,   -Bx,  -1/T1]])
    dMdt = np.matmul(mat, M)
    dMdt[2] += M0/T1
    return dMdt
    

class Simulator():
    
    def __init__(self, pulseSeq=None, settings=None, sample=None):
        if settings is None:
            self.settings = {'B0':7.0, 'observe':'1H', 'offset':0.0, 'gain':1.0}
        else:
            self.settings = settings
        self.sample = sample
        self.allSpins = None
        self.pulseSeq = None
        if pulseSeq is not None:
            self.loadPulseSeq(pulseSeq)
        self.reset()
        
    def reset(self):
        if self.sample is not None:
            self.allSpins = np.array(list(self.sample.expandBroadening(self.sample.expandSystems(self.settings['B0'], self.settings['observe']))))
            if len(self.allSpins) == 0: # When there are no spins, include a 'zero' spin
                self.allSpins = np.array([[0.0, 0.0, 1.0, 1.0, 1.0]])
            self.allSpinsCurrentAmp = np.copy(self.allSpins[:,1])
        else:
            self.allSpins = None
            self.allSpinsCurrentAmp = None
        self.phaseIter = 0
        self.arrayIter = 0
        self.FID = []
        self.scaledFID = []
        self.FIDtime = 0.0

    def step(self):
        self.phaseIter = 0
        self.arrayIter += 1

    def setSettings(self, field=None, nuclei=None, offset=None, gain=None):
        if field is not None:
            self.settings['B0'] = field
        if nuclei is not None:
            self.settings['observe'] = nuclei
        if offset is not None:
            self.settings['offset'] = offset
        if gain is not None:
            self.settings['gain'] = gain
        
    def setSample(self, sample):
        self.sample = sample
        
    def loadPulseSeq(self, name):
        """
        Loads the pulse sequence from a given csv file
        """
        pulseSeq = pd.read_csv(name, index_col=['name','type'])
        pulseSeq = pulseSeq.applymap(lambda x: ast.literal_eval(x) if isinstance(x,str) else x)
        pulseSeq = pulseSeq.reset_index()
        self.pulseSeq = pulseSeq

    def setPulseSeq(self, pulseSeq):
        self.pulseSeq = pulseSeq

    def getData(self):
        if len(self.scaledFID) == 0:
            return [], []
        t = np.linspace(0, self.FIDtime, len(self.scaledFID[0]))
        return t, np.array(self.scaledFID)

    def scan(self, realTime=False):
        """
        Simulate using the defined pulseSeq, settings, and sample.
        """
        fid = 0
        if realTime:
            listFilter = self.pulseSeq['time'].apply(lambda x: isinstance(x, list))
            timeLists = self.pulseSeq['time'][listFilter]
            timeLists = timeLists.apply(lambda x: x[self.arrayIter % len(x)])
            time.sleep(timeLists.sum() + self.pulseSeq['time'][~listFilter].sum())
        for i, spinInfo in enumerate(self.allSpins):
            Frequency, Intensity, T1, T2 = spinInfo
            spinFID, self.allSpinsCurrentAmp[i], self.FIDtime = self.simulateSpin(self.allSpinsCurrentAmp[i], Frequency, T1, T2, len(self.allSpins), Intensity)
            fid += spinFID
        fid = ADC(self.settings['gain']*fid)
        if len(self.FID) == self.arrayIter:
            self.FID.append(fid)
            self.scaledFID.append(self.FID[-1])
        else:
            self.FID[self.arrayIter] += fid
            self.scaledFID[self.arrayIter] = self.FID[self.arrayIter] / (self.phaseIter + 1)
        self.phaseIter += 1

    def simulateSpin(self, amp, freq, T1, T2, numSpins, M0):
        M = np.array([0, 0, amp])
        scanResults = []
        FIDtime = 0
        freq -= self.settings['offset']
        for ind, pulseStep in self.pulseSeq.iterrows():
            if pulseStep['type'] == 'pulse':
                rf = pulseStep['amp']
                if hasattr(rf, '__iter__'):
                    rf = rf[self.arrayIter % len(rf)]
                phase = pulseStep['phase']
                if hasattr(phase, '__iter__'):
                    phase = phase[self.phaseIter % len(phase)]
                phase = np.deg2rad(phase)
                B = 2 * np.pi * np.array([rf*np.cos(phase), rf*np.sin(phase), freq])
            else:
                #B = 2 * np.pi * np.array([0, 0, freq])
                B = 2 * np.pi * np.array([0, 0, 0]) # rotating frame per spin
            timeStep = pulseStep['time']
            if hasattr(timeStep, '__iter__'):
                timeStep = timeStep[self.arrayIter % len(timeStep)]
            if pulseStep['type'] == 'FID':
                npoints = int(pulseStep['amp'])
                t_eval = np.linspace(0, timeStep, npoints)
                FIDtime += timeStep
            else:
                t_eval = None
            sol = solve_ivp(diffEq, (0, timeStep), M, t_eval=t_eval, args=(B, T1, T2, M0), vectorized=True)
            if pulseStep['type'] == 'pulse':
                data = np.copy(sol.y)
            else:
                # Convert spin rot-frame to global rot-frame
                data = np.copy(sol.y)
                phi = 2 * np.pi * freq * sol.t
                data[0] = np.cos(phi) * sol.y[0] + np.sin(phi) * sol.y[1]
                data[1] = -1 * np.sin(phi) * sol.y[0] + np.cos(phi) * sol.y[1]
            M = data[:,-1]
            if pulseStep['type'] == 'sat':
                M *= 0
            if pulseStep['type'] == 'FID':
                # TODO: proper scaling factor for noise factor
                SNR = helpFie.getGamma(self.settings['observe']) * np.sqrt(helpFie.getGamma(self.settings['observe'])**3 * self.settings['B0']**3)
                SNR *= (sol.t[1]-sol.t[0]) / 10.0
                noise = np.random.normal(0, 1, len(data[0])) + 1j*np.random.normal(0, 1, len(data[0]))
                fid = SNR * (data[0] - 1j*data[1]) + 0*noise/np.sqrt(float(numSpins))
                fid *= 0.01 # Arbitrary scaling of the signal
                scanResults.append(fid)
        if len(scanResults) > 0:
            scanResults = np.concatenate(scanResults)
        else:
            scanResults = np.array(scanResults)
        return scanResults, M[2], FIDtime

    def saveMatlabFile(self, filePath, name='FID'):
        if len(self.scaledFID) == 0:
            print('No simulated FID available')
            return
        freq = 1e6 * self.settings['B0'] * helpFie.getGamma(self.settings['observe'])
        sw = len(self.scaledFID[0]) / self.FIDtime
        struct = {}
        if len(self.scaledFID) > 1:
            struct['dim'] = 2
            struct['data'] = np.array([self.scaledFID])
            struct['spec'] = [False, False]
            struct['freq'] = [freq-self.settings['offset'], freq-self.settings['offset']]
            struct['sw'] = [sw, sw]
            struct['ref'] = [freq, freq]
            struct['xaxArray'] = [np.arange(len(self.scaledFID)), np.linspace(0, self.FIDtime, len(self.scaledFID[0]))]
            struct['wholeEcho'] = [False, False]
            struct['hyper'] = [0]
        else:
            struct['dim'] = 1
            struct['spec'] = [False]
            struct['data'] = np.array(self.scaledFID)
            struct['freq'] = [freq-self.settings['offset']]
            struct['sw'] = [sw]
            struct['ref'] = [freq]
            struct['xaxArray'] = [np.linspace(0, self.FIDtime, len(self.scaledFID[0]))]
            struct['wholeEcho'] = [False]
            struct['hyper'] = [0]
        matlabStruct = {name: struct}
        scipy.io.savemat(filePath, matlabStruct)
    

if __name__ == '__main__':
    tube = sample.sample()
    tube.addMolecule([('1H',0,1),('1H',1,1),('13C',1,1)],np.array([[0,10,5],[10,0,0],[5,0,0]]), 1, 0.3, 0.3, 1)

    settings = None # {'B0':1.0, 'observe':'1H'}
    pulseSeq = '../pulseSeqs/onePulse.csv'
    
    sim = Simulator(pulseSeq, settings, tube)
    sim.reset()
    for i in range(10):
        sim.scan()
        fid = sim.FID[0] / sim.phaseIter
        fid[0] *= 0.5
        spectrum = np.fft.fftshift(np.fft.fft(fid))
        plt.plot(np.real(spectrum))
        #plt.plot(np.imag(spectrum))
        #sim.step()
    plt.show()
