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
from scipy.interpolate import interp1d
import ast
import matplotlib.pyplot as plt
import scipy.io
import pathlib
import time

import sample
import helperFunctions as helpFie


BITS = 10

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

def diffEq(t, M, totalTime, B, T1, T2, M0, shape=None):
    Bx, By, Bz = B
    if shape is not None and totalTime != 0:
        if t > totalTime:
            amp = 0
        else:
            amp = shape(t/totalTime)
    else:
        amp = 1.0
    mat = np.array([[-1/T2, Bz,  -amp*By],
                    [-Bz,   -1/T2, Bx],
                    [amp*By,   -amp*Bx,  -1/T1]])
    dMdt = np.matmul(mat, M)
    dMdt[2] += M0/T1
    return dMdt

def diffEqExc(t, M, totalTime, B, T1, T2, M0, k, shape=None):
    Bx0, By0, Bz0, Bx1, By1, Bz1 = B
    relk0 = k * M0[0] / np.mean(M0)
    relk1 = k * M0[1] / np.mean(M0)
    if shape is not None and totalTime != 0:
        amp = shape(t/totalTime)
    else:
        amp = 1.0
    mat = np.array([[-1/T2[0]-relk0, Bz0, -amp*By0, relk0, 0, 0],
                    [-Bz0, -1/T2[0]-relk0, amp*Bx0, 0, relk0, 0],
                    [amp*By0, -amp*Bx0, -1/T1[0]-relk0, 0, 0, relk0],
                    [relk1, 0, 0, -1/T2[1]-relk1, Bz1, -amp*By1],
                    [0, relk1, 0, -Bz1, -1/T2[1]-relk1, amp*Bx1],
                    [0, 0, relk1, amp*By1, -amp*Bx1, -1/T1[1]-relk1]])
    dMdt = np.matmul(mat, M)
    dMdt[2] += M0[0]/T1[0]
    dMdt[5] += M0[1]/T1[1]
    return dMdt


class Simulator():
    
    def __init__(self, pulseSeq=None, settings=None, sample=None):
        if settings is None:
            self.settings = {'B0':7.0, 'observe':'1H', 'decouple':None, 'offset':0.0, 'gain':1.0}
        else:
            self.settings = settings
        self.sample = sample
        self.allSpins = None # array with spin information, column order: [Frequency, Intensity, T1, T2]
        self.allPairs = None # array with exchanging pairs information, column order: [Frequency_0, Frequency_1, Intensity_0, Intensity_1, T1_0, T1_1, T2_0, T2_1, k]
        self.pulseSeq = None
        if pulseSeq is not None:
            self.loadPulseSeq(pulseSeq)
        self.reset()
        
    def reset(self):
        if self.sample is not None:
            self.allSpins = np.array(list(self.sample.expandBroadening(self.sample.expandSystems(self.settings['B0'], self.settings['observe'], self.settings['decouple']))))
            if len(self.allSpins) == 0: # When there are no spins, include a 'zero' spin
                self.allSpins = np.array([[0.0, 0.0, 1.0, 1.0]])
            self.allSpinsCurrentAmp = np.copy(self.allSpins[:,1])
            self.allPairs = np.array(list(self.sample.expandPairs(self.settings['B0'], self.settings['observe'], self.settings['decouple'])))
            if self.allPairs is not None and len(self.allPairs) > 0:
                self.allPairsCurrentAmp = np.copy(self.allPairs[:,2:4])
            else:
                self.allPairsCurrentAmp = None
        else:
            self.allSpins = None
            self.allPairs = None
            self.allSpinsCurrentAmp = None
            self.allPairsCurrentAmp = None
        self.phaseIter = 0
        self.arrayIter = 0
        self.FID = []
        self.scaledFID = []
        self.FIDtime = 0.0

    def step(self):
        self.phaseIter = 0
        self.arrayIter += 1

    def setSettings(self, field=None, nuclei=None, decouple=None, offset=None, gain=None):
        if field is not None:
            self.settings['B0'] = field
        if nuclei is not None:
            self.settings['observe'] = nuclei
        if offset is not None:
            self.settings['offset'] = offset
        if gain is not None:
            self.settings['gain'] = gain
        self.settings['decouple'] = decouple                
        
    def setSample(self, sample):
        self.sample = sample
        
    def loadPulseSeq(self, name):
        """
        Loads the pulse sequence from a given csv file
        """
        pulseSeq = pd.read_csv(name, index_col=['name','type'])
        pulseSeq[['time', 'amp', 'phase']] = pulseSeq[['time', 'amp', 'phase']].applymap(lambda x: ast.literal_eval(x) if isinstance(x,str) else x)
        pulseSeq = pulseSeq.reset_index()
        if 'shapefile' in pulseSeq.keys():
            self.pulseshapes = self.loadShapeFiles(pulseSeq['shapefile'].dropna().unique())
        else:
            self.pulseshapes = {}
        self.pulseSeq = pulseSeq

    def setPulseSeq(self, pulseSeq):
        self.pulseSeq = pulseSeq

    def loadShapeFiles(self, filenames):
        pulseshapes = {}
        for f in filenames:
            shapeFile = pathlib.Path(__file__).parents[1] / 'shapefiles' / f # shapefile in library
            if not shapeFile.is_file():
                # If not in library, the filename is the filepath
                shapeFile = pathlib.Path(f)
            shape = np.loadtxt(shapeFile)
            shapefunc = interp1d(np.linspace(0,1,len(shape)), shape)
            pulseshapes[f] = shapefunc
        return pulseshapes

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
            Freq, Intensity, T1, T2 = spinInfo
            spinFID, self.allSpinsCurrentAmp[i], self.FIDtime = self.simulateSpin(self.allSpinsCurrentAmp[i], Freq, T1, T2, len(self.allSpins), Intensity)
            fid += spinFID
        for i, spinInfo in enumerate(self.allPairs):
            Freq0, Freq1, Intensity0, Intensity1, T1_0, T1_1, T2_0, T2_1, k = spinInfo
            spinFID, self.allPairsCurrentAmp[i], self.FIDtime = self.simulateSpin(self.allPairsCurrentAmp[i], np.array([Freq0, Freq1]), np.array([T1_0, T1_1]), np.array([T2_0, T2_1]), None, np.array([Intensity0, Intensity1]), k)
            fid += spinFID
        fid = ADC(self.settings['gain']*fid)
        if len(self.FID) == self.arrayIter:
            self.FID.append(fid)
            self.scaledFID.append(self.FID[-1])
        else:
            self.FID[self.arrayIter] += fid
            self.scaledFID[self.arrayIter] = self.FID[self.arrayIter] / (self.phaseIter + 1)
        self.phaseIter += 1

    def simulateSpin(self, amp, freq, T1, T2, numSpins, M0, k=None):
        if k is None:
            M = np.array([0, 0, amp])
        else:
            M = np.array([0, 0, amp[0], 0, 0, amp[1]])
        scanResults = []
        FIDtime = 0
        freq -= self.settings['offset']
        for ind, pulseStep in self.pulseSeq.iterrows():
            if pulseStep['type'] in ['pulse', 'shapedPulse']:
                rf = pulseStep['amp']
                if hasattr(rf, '__iter__'):
                    rf = rf[self.arrayIter % len(rf)]
                phase = pulseStep['phase']
                if hasattr(phase, '__iter__'):
                    phase = phase[self.phaseIter % len(phase)]
                phase = np.deg2rad(phase)
                if k is None:
                    B = 2 * np.pi * np.array([rf*np.cos(phase), rf*np.sin(phase), freq])
                else:
                    B = 2 * np.pi * np.array([rf*np.cos(phase), rf*np.sin(phase), freq[0], rf*np.cos(phase), rf*np.sin(phase), freq[1]])
            else:
                if k is None:
                    rotFreq = freq
                    B = 2 * np.pi * np.array([0, 0, 0]) # rotating frame per spin
                else:
                    rotFreq = np.mean(freq)
                    B = 2 * np.pi * np.array([0, 0, freq[0]-rotFreq, 0, 0, freq[1]-rotFreq]) # rotating frame for spin 1
            timeStep = pulseStep['time']
            if hasattr(timeStep, '__iter__'):
                timeStep = timeStep[self.arrayIter % len(timeStep)]
            if pulseStep['type'] == 'FID':
                npoints = int(pulseStep['amp'])
                t_eval = np.linspace(0, timeStep, npoints)
                FIDtime += timeStep
            else:
                t_eval = None
            if pulseStep['type'] == 'shapedPulse':
                shapeFunc = self.pulseshapes[pulseStep['shapefile']]
            else:
                shapeFunc = None
            if k is None:
                sol = solve_ivp(diffEq, (0, timeStep), M, t_eval=t_eval, args=(timeStep, B, T1, T2, M0, shapeFunc), vectorized=True)
            else:
                sol = solve_ivp(diffEqExc, (0, timeStep), M, t_eval=t_eval, args=(timeStep, B, T1, T2, M0, k, shapeFunc), vectorized=True)
            if pulseStep['type'] in ['pulse', 'shapedPulse']:
                data = np.copy(sol.y)
            else:
                # Convert spin rot-frame to global rot-frame
                data = np.copy(sol.y)
                phi = 2 * np.pi * rotFreq * sol.t
                if k is None:
                    data[0] = np.cos(phi) * sol.y[0] + np.sin(phi) * sol.y[1]
                    data[1] = -1 * np.sin(phi) * sol.y[0] + np.cos(phi) * sol.y[1]
                else:
                    data[0] = np.cos(phi) * sol.y[0] + np.sin(phi) * sol.y[1]
                    data[1] = -1 * np.sin(phi) * sol.y[0] + np.cos(phi) * sol.y[1]
                    data[3] = np.cos(phi) * sol.y[3] + np.sin(phi) * sol.y[4]
                    data[4] = -1 * np.sin(phi) * sol.y[3] + np.cos(phi) * sol.y[4]
            M = data[:,-1]
            if pulseStep['type'] == 'sat':
                M *= 0
            if pulseStep['type'] == 'FID':
                # TODO: proper scaling factor for noise factor
                SNR = 0.5 * helpFie.getGamma(self.settings['observe']) * np.sqrt(helpFie.getGamma(self.settings['observe'])**3 * self.settings['B0']**3)
                SNR *= (sol.t[1]-sol.t[0])
                noise = np.random.normal(0, 1, len(data[0])) + 1j*np.random.normal(0, 1, len(data[0]))
                if k is None:
                    fid = SNR * (data[0] - 1j*data[1]) + noise/np.sqrt(float(numSpins))
                else:
                    # No noise in the pairs!
                    fid = SNR * (data[0] - 1j*data[1] + data[3] - 1j*data[4])
                fid *= 0.001 # Arbitrary scaling of the signal
                scanResults.append(fid)
        if len(scanResults) > 0:
            scanResults = np.concatenate(scanResults)
        else:
            scanResults = np.array(scanResults)
        if k is None:
            finalMz = M[2]
        else:
            finalMz = [M[2], M[5]]
        return scanResults, finalMz, FIDtime

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
