#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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


import os.path
from PyQt5 import QtGui, QtCore, QtWidgets
import sys
import numpy as np

import matplotlib
# First import matplotlib and Qt
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import loadIsotopes
VERSION = '0.0'

ISOTOPES = loadIsotopes.getIsotopes('IsotopeProperties')
NUCLEI = [x for x in ISOTOPES.keys() if not x.startswith('-')]

class MainProgram(QtWidgets.QMainWindow):

    def __init__(self, root):
        super(MainProgram, self).__init__()
        #self.setAcceptDrops(True)
        self.main_widget = QtWidgets.QWidget(self)
        self.setCentralWidget(self.main_widget)
        
        self.mainFrame = QtWidgets.QGridLayout(self.main_widget)
        self.spectrometerFrame = SpectrometerFrame(self)
        self.mainFrame.addWidget(self.spectrometerFrame,0,0)
        self.plotFrame = PlotFrame(self)
        self.mainFrame.addWidget(self.plotFrame,0,1)
        self.paramFrame = ParameterFrame(self)
        self.mainFrame.addWidget(self.paramFrame,1,0,1,2)
        self.mainFrame.setColumnStretch(1, 1)
        self.mainFrame.setRowStretch(0, 1)
        self.exportFrame = ExportFrame(self)
        self.mainFrame.addWidget(self.exportFrame,2,0,1,2)


        self.resize(800, 700)
        self.show()


class ExportFrame(QtWidgets.QWidget):
    def __init__(self, parent):
        super(ExportFrame, self).__init__(parent)
        self.father = parent
        grid = QtWidgets.QGridLayout(self)
        self.setLayout(grid)
        self.exportFIDButton = QtWidgets.QPushButton('Export FID')
        self.exportFIDButton.setEnabled(False)
        grid.addWidget(self.exportFIDButton,0,10)
        self.exportSpectrumButton = QtWidgets.QPushButton('Export Spectrum')
        self.exportSpectrumButton.setEnabled(False)
        grid.addWidget(self.exportSpectrumButton,0,11)
        grid.setColumnStretch(0, 1)


class ParameterFrame(QtWidgets.QTabWidget):
    def __init__(self, parent):
        super(ParameterFrame, self).__init__(parent)
        self.father = parent
        self.addSaturation('Saturation')
        self.addDelay('D1')
        self.addPulse('P1')
        self.addShapedPulse('P2')
        self.addAcq('Acquisition')
        self.currentChanged.connect(self.tabChanged)
        #self.tabChanged(0)

    def tabChanged(self,index):
        self.father.plotFrame.sequenceFrame.highlightElement(index)
       
    def addDelay(self,title):
        self.addTab(DelayWidget(self),title)

    def addPulse(self,title):
        self.addTab(PulseWidget(self),title)

    def addShapedPulse(self,title):
        self.addTab(PulseWidget(self),title)

    def addAcq(self,title):
        self.addTab(AcqWidget(self),title)

    def addSaturation(self,title):
        self.addTab(ParameterWidget(self),title)

    def getSettings(self):
        """
        Returns all the settings in a dictionary or list format.
        To be used when acquiring
        """

def safeEval(inp, Type='All'):
    """
    Creates a more restricted eval environment.

    Parameters
    ----------
    inp : str
        String to evaluate.
    Type : {'All', 'FI', 'C'}, optional
        Type of expected output. 'All' will return all types, 'FI' will return a float or int, and 'C' will return a complex number.
        By default Type is set to 'All'

    Returns
    -------
    Object
        The result of the evaluated string.
    """
    env = vars(np).copy()
    try:
        val = eval(inp, env)
        if isinstance(val, str):
            return None
        if Type == 'All':
            return val
        if Type == 'FI':  #single float/int type
            if isinstance(val, (float, int)) and not np.isnan(val) and not np.isinf(val):
                return val
            return None
        if Type == 'C': #single complex number
            if isinstance(val, (float, int, complex)) and not np.isnan(val) and not np.isinf(val):
                return val
            return None
    except Exception:
        return None

class ParameterWidget(QtWidgets.QWidget):
    def __init__(self, parent):
        super(ParameterWidget, self).__init__(parent)
        self.grid = QtWidgets.QGridLayout(self)
        self.grid.setColumnStretch(10, 1)
        self.params = []

    def returnValues(self):
        values = []
        for elem in self.params:
            values.append(safeEval(elem.text())) #Needs safe eval, and array detection.
        print(values)
        return values

class DelayWidget(ParameterWidget):
    def __init__(self, parent):
        super(DelayWidget, self).__init__(parent)
        self.grid.addWidget(QtWidgets.QLabel('Duration [s]'),0,0)
        self.delay = QtWidgets.QLineEdit('1')
        self.params.append(self.delay)
        self.grid.addWidget(self.delay,0,1)

class PulseWidget(ParameterWidget):
    def __init__(self, parent):
        super(PulseWidget, self).__init__(parent)
        self.grid.addWidget(QtWidgets.QLabel('Duration [µs]'),0,0)
        self.delay = QtWidgets.QLineEdit('1')
        self.params.append(self.delay)
        self.grid.addWidget(self.delay,0,1)

        self.grid.addWidget(QtWidgets.QLabel('Amplitude [kHz]'),0,2)
        self.amplitude = QtWidgets.QLineEdit('1')
        self.params.append(self.amplitude)
        self.grid.addWidget(self.amplitude,0,3)

class AcqWidget(ParameterWidget):
    def __init__(self, parent):
        super(AcqWidget, self).__init__(parent)
        self.grid.addWidget(QtWidgets.QLabel('Spectral Width [kHz]'),0,0)
        self.sw = QtWidgets.QLineEdit('100')
        self.params.append(self.sw)
        self.grid.addWidget(self.sw,0,1)

        self.grid.addWidget(QtWidgets.QLabel('# of points'),0,2)
        self.np = QtWidgets.QLineEdit('1024')
        self.params.append(self.np)
        self.grid.addWidget(self.np,0,3)

        self.grid.addWidget(QtWidgets.QLabel('Acq. time [s]'),0,4)
        self.acqTime = QtWidgets.QLineEdit('1')
        self.params.append(self.acqTime)
        self.grid.addWidget(self.acqTime,0,5)


class PlotFrame(QtWidgets.QTabWidget):
    def __init__(self, parent):
        super(PlotFrame, self).__init__(parent)
        self.father = parent
        self.specFrame = SpecPlotFrame(self)
        self.specFrame.plot([-1,-2],[3,4])
        self.fidFrame = FidPlotFrame(self)
        self.sequenceFrame = SequenceDiagram(self)
        self.addTab(self.specFrame,'Spectrum')
        self.addTab(self.fidFrame,'FID')
        self.addTab(self.sequenceFrame,'Pulse Sequence')


class AbstractPlotFrame(QtWidgets.QWidget):
    SPEC = False
    def __init__(self, parent):
        super(AbstractPlotFrame, self).__init__(parent)
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(self.canvas, 0, 0)
        self.ax = self.fig.add_subplot(111)

    def plot(self,xdata,ydata):
        self.ax.plot(xdata,ydata)
        self.canvas.draw()
        if self.SPEC:
            self.ax.set_xlim(np.max(xdata), np.min(xdata))
        else:
            self.ax.set_xlim(np.min(xdata),np.max(xdata))

    def resetPlot(self):
        self.ax.cla()


class SequenceDiagram(AbstractPlotFrame):
    PULSEW = 0.75
    PULSEH = 1
    SHAPEW = 1.5
    DELAYW = 2
    FIDW = 3
    FIDH = 1
    SATH = 0.2
    SATW = 2
    LINEWIDTH = 2
    FONTSIZE = 12
    TEXTHEIGHT = 1.1

    def __init__(self, parent):
        super(SequenceDiagram, self).__init__(parent)
        self.father = parent
        self.ax.axis('off')
        self.xpos = 0
        self.elems = []
        self.backElems = [] #Holds transparent background images
        self.canvas.mpl_connect('pick_event', self.pickHandler)
        self.setIsotope('1H')
        self.drawSaturation('Sat.')
        self.drawDelay('RD')
        self.drawPulse('p1')
        self.drawShapedPulse('p2')
        self.drawFid('Acq')
        self.ax.set_xlim([-1,10])
        self.ax.set_ylim([-1,1.5])

    def pickHandler(self,pickEvent):
        artist = pickEvent.artist
        for pos, elem in enumerate(self.elems):
            if elem[0] is artist:
                self.highlightElement(pos)

        for pos, elem in enumerate(self.backElems):
            if elem is artist:
                self.highlightElement(pos)

    def setIsotope(self,iso):
        self.ax.text(-.2,0,iso,horizontalalignment='right',verticalalignment='center',fontsize=self.FONTSIZE)

    def drawPulse(self,text=None):
        elem = self.ax.add_patch(matplotlib.patches.Rectangle((self.xpos, 0), self.PULSEW, self.PULSEH,linewidth=self.LINEWIDTH,edgecolor='tab:blue',picker=True))
        self.elems.append([elem,'tab:blue'])
        self._addText(text,self.PULSEW)
        #Add transparent background
        self._addBackgroundRect(self.PULSEW)
        self.xpos += self.PULSEW

    def drawShapedPulse(self,text=None):
        samples = 100
        xdata = np.linspace(-2*np.pi,2*np.pi,samples)
        ydata = self.PULSEH * np.sin(xdata)/xdata 
        elem = self.ax.plot((xdata + 2*np.pi) * self.SHAPEW / (4*np.pi) + self.xpos,ydata,c='tab:blue',linewidth=self.LINEWIDTH,picker=True)
        self.elems.append([elem[0],'tab:blue'])
        self._addText(text,self.SHAPEW)
        self._addBackgroundRect(self.PULSEW)
        self.xpos += self.SHAPEW

    def drawSaturation(self,text=None):
        elem = self.ax.add_patch(matplotlib.patches.Rectangle((self.xpos, 0), self.SATW, self.SATH,linewidth=self.LINEWIDTH,edgecolor='tab:blue',picker=True))
        self.elems.append([elem,'tab:blue'])
        self._addText(text,self.SATW)
        self._addBackgroundRect(self.SATW)
        self.xpos += self.SATW

    def drawDelay(self,text=None):
        elem = self.ax.plot([self.xpos,self.xpos+self.DELAYW],[0,0],c='tab:blue',linewidth=self.LINEWIDTH,picker=True)
        self.elems.append([elem[0],'tab:blue'])
        self._addText(text,self.DELAYW)

        self._addBackgroundRect(self.DELAYW)
        self.xpos += self.DELAYW

    def drawFid(self,text=None):
        samples = 100
        T2 = 0.3
        xdata = np.linspace(0,1,samples)
        ydata = self.FIDH * np.cos(xdata * 30) * np.exp(-xdata / T2)
        self._addText(text,self.FIDW)
        elem = self.ax.plot(xdata * self.FIDW + self.xpos,ydata,c='tab:blue',linewidth=self.LINEWIDTH,picker=True)
        self.elems.append([elem[0],'tab:blue'])
        self._addBackgroundRect(self.FIDW)
        self.xpos += self.FIDW

    def _addText(self,text,xAdder):
        if text is not None:
            self.ax.text(self.xpos + 0.5*xAdder,self.TEXTHEIGHT,text,horizontalalignment='center',fontsize=self.FONTSIZE)

    def _addBackgroundRect(self,width):
        backelem = self.ax.add_patch(matplotlib.patches.Rectangle((self.xpos, -self.PULSEH), width, 2*self.PULSEH,linewidth=self.LINEWIDTH,picker=True,alpha=0))
        self.backElems.append(backelem)

    def resetColors(self):
        for elem in self.elems:
            elem[0].set_color(elem[1])

    def highlightElement(self,pos):
        for index, elem in enumerate(self.elems):
            if pos == index:
                elem[0].set_color('red')
            else:
                elem[0].set_color(elem[1])
        self.canvas.draw()
        self.father.father.paramFrame.setCurrentIndex(pos)

    def resetPlot(self):
        AbstractPlotFrame.resetPlot(self)
        self.ax.axis('off')
        self.xpos = 0
        self.elems = []
        self.backElems = []


class FidPlotFrame(AbstractPlotFrame):
    def __init__(self, parent):
        super(FidPlotFrame, self).__init__(parent)
        self.ax.set_xlabel('Time [s]')
        self.ax.set_ylabel('Intensity [arb. u.]')


class SpecPlotFrame(AbstractPlotFrame):
    SPEC = True
    def __init__(self, parent):
        super(SpecPlotFrame, self).__init__(parent)
        self.ax2 = self.ax.twiny()
        self.ax.set_xlabel('Shift [ppm]')
        self.ax.set_ylabel('Intensity [arb. u.]')
        self.ax2.set_xlabel('Frequency [kHz]')

    def plot(self,xdata,ydata):
        AbstractPlotFrame.plot(self,xdata,ydata)
        self.ax2.set_xlim([100,-100])


class SpectrometerFrame(QtWidgets.QFrame):
    def __init__(self, parent):
        super(SpectrometerFrame, self).__init__(parent)
        self.father = parent
        self.setFrameShape(1)
        grid = QtWidgets.QGridLayout(self)
        self.setLayout(grid)
        self.title = QtWidgets.QLabel('<b>Spectrometer settings:</b>')
        self.title.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(self.title,0,0,1,2)
        self.loadSample = QtWidgets.QPushButton('Load sample')
        grid.addWidget(self.loadSample,1,0,1,2)
        self.sampleLabel = QtWidgets.QLabel('<i>No sample</i>')
        self.sampleLabel.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(self.sampleLabel,2,0,1,2)
        self.loadSequence = QtWidgets.QPushButton('Load pulse sequence')
        grid.addWidget(self.loadSequence,3,0,1,2)
        self.sequenceLabel = QtWidgets.QLabel('<i>No sequence</i>')
        self.sequenceLabel.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(self.sequenceLabel,4,0,1,2)

        grid.addWidget(QtWidgets.QLabel('Detect:'),5,0)
        self.decoupleLabel = QtWidgets.QLabel('Decouple:')
        grid.addWidget(self.decoupleLabel,6,0)

        self.detectDrop = QtWidgets.QComboBox()
        self.detectDrop.addItems(NUCLEI)
        grid.addWidget(self.detectDrop,5,1)
        self.decoupleDrop = QtWidgets.QComboBox()
        self.decoupleDrop.addItems(NUCLEI)
        self.decoupleDrop.setEnabled(False)
        self.decoupleLabel.setEnabled(False)
        grid.addWidget(self.decoupleDrop,6,1)


        grid.addWidget(QtWidgets.QLabel('B<sub>0</sub>:'),7,0)
        self.b0Drop = QtWidgets.QComboBox()
        self.b0Drop.addItems(['7.0T','9.4T','11.7','14.1T','18.8T'])
        grid.addWidget(self.b0Drop,7,1)

        self.acquirePush = QtWidgets.QPushButton('Acquire')
        self.acquirePush.clicked.connect(self.simulate)
        grid.addWidget(self.acquirePush,9,0,1,2)

        #grid.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        grid.setRowStretch(20, 1)
        
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.grid = grid

    def simulate(self):
        parameters = self.father.paramFrame.getSettings()
        pass


if __name__ == '__main__':

    root = QtWidgets.QApplication(sys.argv)
    #root.setWindowIcon(QtGui.QIcon(os.path.dirname(os.path.realpath(__file__)) + '/Icons/Logo.png'))
    mainProgram = MainProgram(root)
    mainProgram.setWindowTitle(f"Magpie - {VERSION}")
    mainProgram.show()
    sys.exit(root.exec_())
