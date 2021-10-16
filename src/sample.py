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
import loadIsotopes

ISOTOPES = loadIsotopes.getIsotopes('IsotopeProperties')
GAMMASCALE = 42.576e6/100

def getSplittingPattern(I,Multi):
    Kernel = np.ones((int(2*I+1)))
    IntenFull = [1]
    for _ in range(Multi):
        IntenFull = np.convolve(IntenFull, Kernel)
    IntenFull /= np.sum(IntenFull) 
    Split = np.arange(len(IntenFull)) - len(IntenFull)/2 + 0.5
    return Split, IntenFull

def kronAdd(A,B):
    """
    Performs a Kronecker addition.

    Parameters
    -------
    A: array_like
    B: array_like

    Returns
    -------
    array_like: result of the Kronecker addition
    """
    Final = np.array([])
    for val in A:
        Final = np.append(Final,B + val)
    return Final

def scaleScalarCouplings(J,spins,decouple,B0):
    """
    Applies J coupling scaling based on decoupling settings.
    Only spins with the same type as the decoupling are affected.

    Jscaling is cos(theta)
    with theta = atan(B1/offset)
    According to Ernst, Bodenhausen, Wokaun. p 235-236.

    Parameters
    -------
    J: 2D array, J coupling matrix
    spins: list of list, for each spin [Nucleus,shift,multi,T1,T2,T2prime]
    decouple: list, decoupling settings [Nucleus,offset [Hz],strength [Hz]]
    B0: float, magnetic field strength [T]

    Returns
    -------
    2D array: scaled J matrix
    """
    for pos, spin in enumerate(spins):
        if spin[0] == decouple[0]: # If nucleus is decoupled
            shift = spin[1] * B0 * ISOTOPES[spin[0]][1] * GAMMASCALE * 1e-6
            offset = decouple[1] - shift
            scale = np.cos(np.arctan(decouple[2]/offset))
            # Scale row and column of Jmatrix with this value
            J[pos,:] = J[pos,:] * scale
            J[:,pos] = J[pos,:] * scale
    return J

def lorentz(T2,axis):
    """
    Calculates a Lorentz distribution for a T2 and an sampling axis.
    No shape offset is used. The values are normalized to the sum.

    Parameters
    -------
    T2: float, T2 lifetime (1/width)
    axis: 1D array, 'frequency' axis values.

    Returns
    -------
    1D array: intensity values.
    """
    w = 1/(T2*np.pi)
    x = axis/(w/2)
    L = 1/(1+x**2)
    L /= np.sum(L)
    return L

class sample():
    def __init__(self):
        self.moleculeList = []

    def addMolecule(self,nucl,Jmatrix,amp,T1=None,T2=None,T2prime=None):
        """
        nuclei: list of tuples with (type,shift,multi,relax=None)
            Type: string value of the nucleus ('1H')
            shift: float, chemical shift in ppm
            multi: int, multiplicity of the nucleus (e.g. 3 for a CH3 group)
            relax: list of [T1,T2,T2prime] in s. Is optional. If not given, molecule master
                   values are used. Within the tuple, None can be used to default to master value.
        The relaxation times are optional, as the can also be input for the whole molecule
        as master values. 

        Parameters
        ----------
        Jmatrix: 2D array holding the J value between nuclei in Hz. (Symmatric)
        amp: amplitude (concentration) of the molecule
        T1: float, master T1 value is seconds
        T2: float, master T2 value is seconds
        T2prime: float, master T2prime value is seconds

        Note that for every spin, all three relaxation values need to be known one way or the other.
        A message in printed in case there are missing values.
        """
        if Jmatrix.ndim != 2 or Jmatrix.shape[0] != Jmatrix.shape[1] or Jmatrix.shape[0] != len(nucl):
            print('Jmatrix dimensions not correct')
            return

        Jmatrix = np.array(Jmatrix,dtype=float)

        molecule = dict()
        # Create a tuple for each nucleus
        spins = []
        for pos, s in enumerate(nucl):
            if len(s) == 4:
                Type, shift, multi, relax = s
                if relax[0] is None:
                    relax[0] = T1
                if relax[1] is None:
                    relax[1] = T2
                if relax[2] is None:
                    relax[2] = T2prime
            elif len(s) == 3:
                Type, shift, multi = s
                relax = [T1,T2,T2prime]
            if relax[0] is None or relax[1] is None or relax[2] is None:
                print(f'Spin #{pos} is missing a relaxation value')

            spins.append([Type,shift,multi,relax[0],relax[1],relax[2]])
        molecule['spins'] = spins
        molecule['J'] = Jmatrix
        molecule['amp'] = amp
        self.moleculeList.append(molecule)

    def removeMolecule(self,index):
        self.moleculeList.pop(index)

    def expandSystems(self,B0,observe,decouple = None):
        """
        Expands the spin systems in individual 'lines', with frequency, intensity,
        relaxation values.
        Only spins relating to the Observed nucleus are returned.

        Parameters
        ----------
        B0: float, magnetic field strength in Tesla
        observe: string, observe nucleus
        decouple [optional]: nucleus, offset [Hz] and CW decouple power [Hz]

        Returns
        -------
        List of list with the individual 'lines', having:
            [Frequency, Intensity, T1, T2, T2prime]
        """

        # Scale the J matrix before, with the decoupling

        Expanded = []
        for system in self.moleculeList:
            spins = system['spins']
            concentration = system['amp']
            J = system['J']

            # Splitting patterns
            FreqSplits = []
            IntSplits = []
            for spin in spins:
                Multi = spin[2]
                I = ISOTOPES[spin[0]][0]
                Split, IntenFull = getSplittingPattern(I,Multi)
                FreqSplits.append(Split)
                IntSplits.append(IntenFull)

            # Jscaling
            if decouple is not None:
                J = scaleScalarCouplings(J,spins,decouple,B0)
               
            for pos, spin in enumerate(spins):
                if spin[0] == observe:
                    Jlist = J[pos,:]
                    Shift = spin[1]
                    Multi = spin[2]
                    freqList = np.array([Shift * B0 * ISOTOPES[observe][1] * GAMMASCALE * 1e-6]) # Freq of spin in Hz
                    intList = np.array([Multi])
                    for pos2, Jval in enumerate(Jlist):
                        if pos2 != pos and Jval != 0:
                            Jmatr = FreqSplits[pos2] * Jval
                            IntMatr = IntSplits[pos2]
                            freqList = kronAdd(freqList,Jmatr)
                            intList = np.kron(intList,IntMatr)
                    
                    for place, _ in enumerate(freqList):
                        yield [freqList[place],intList[place] * concentration,spin[3],spin[4],spin[5]]


    def expandBroadening(self,spinList,widthMax = 3, factor = 50):
        """
        Expand each spin, to take T2' broadening (e.g. shimming)
        into account. Do this via subsampling with lorentz lines.

        Parameters
        ----------
        spinList: iterable of lists with [freq,intensity,T1,T2,T2prime]
        widthMax: float, indicates the maximum distance from the centre 
             of the T2prime lorentz that needs to be sampled in units of (1/T2prime)
             [default = 3]
        factor: int, indicates the amount of sampled points, which is multiplied
             with the T2:T2prime ratio.

        Returns
        -------
        Iterable of list with the individual 'lines', having:
            [Frequency, Intensity, T1, T2]
        """
        for spin in spinList:
            print(spin)
            freq, inten, T1, T2, T2prime = spin
            samples = int(T2/T2prime * factor)
            if samples%2 == 0:
                samples += 1
            print(samples)
            limit = widthMax / T2prime
            sel = np.linspace(-limit,limit,samples)
            distr = lorentz(T2prime,sel)
            for pos, shift in enumerate(sel):
                intenNew = inten * distr[pos]
                freqNew = freq + shift
                yield [freqNew,intenNew,T1,T2]


def loadSampleFile(loc):
    with open(loc,'r') as f:
        lines = [x.strip() for x in f.readlines()]

    if lines[0] != '###SAMPLE###':
        raise Exception('Incorrect first line') 
    
    # Initialize vars
    mainAmount = 1
    if lines[1].startswith('amount'):
        mainAmount = float(lines[1].split()[1])

    molPos = [x for x in range(len(lines)) if lines[x] == '###MOLECULE###']
    for val in molPos:
        lin = lines[val+1:]
        molecule = dict()
        molecule['spins'] = []
        for elem in lin:
            if elem.startswith('amount'):
                molecule['amount'] = float(elem.split()[1]) * mainAmount
            elif elem.startswith('T1 '):
                molecule['T1'] = float(elem.split()[1])
            elif elem.startswith('T2 '):
                molecule['T2'] = float(elem.split()[1])
            elif elem.startswith('T2prime '):
                molecule['T2prime'] = float(elem.split()[1])
            elif elem.startswith('spin '):
                spinelem = elem.split()[1:]
                if len(spinelem) == 3:
                    iso, shift, multi = spinelem
                    shift = float(shift)
                    multi = int(multi)
                    molecule['spins'].append([iso,shift,multi])



    print(molecule)

        



if __name__ == '__main__':

    # Some test code
    loadSampleFile(r'TestFiles/Ethanol.txt')

    #tube = sample()
    #tube.addMolecule([('1H',0,1),('1H',1,1),('13C',1,1)],np.array([[0,10,5],[10,0,0],[5,0,0]]),1,3,1,1)
    #elems = tube.expandSystems(1,'1H',['13C',0,10000])
    #elem2 = tube.expandBroadening(elems)


    

