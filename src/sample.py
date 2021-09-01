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
                """
                Jscaling is cos(theta)
                with theta = atan(B1/offset)
                According to Ernst, Bodenhausen, Wokaun. p 235-236
                """
                pass
               
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
                        Expanded.append([freqList[place],intList[place] * concentration,spin[3],spin[4],spin[5]])
        return Expanded

            

    def expandBroadening(self,nsteps=10):
        """
        Expand the list of spins further, making a series for the T2prime broadening
        nsteps: int, number of subdivisions per spin
        or: resolution of subdivision?
        """
        pass


if __name__ == '__main__':

    # Some test code
    tube = sample()
    tube.addMolecule([('1H',0,2),('1H',1,2),('1H',2,1)],np.array([[0,10,5],[10,0,0],[5,0,0]]),1,3,1,1)
    tube.expandSystems(7,'1H')

    

