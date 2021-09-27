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


class ParameterFrame(QtWidgets.QFrame):
    def __init__(self, parent):
        super(ParameterFrame, self).__init__(parent)
        self.father = parent
        self.setFrameShape(1)
        grid = QtWidgets.QGridLayout(self)
        self.setLayout(grid)
        self.title = QtWidgets.QLabel('<b>Acquisition parameters:</b>')
        #self.title.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(self.title,0,0,1,10)

        self.swLabel = QtWidgets.QLabel('Spectral Width [kHz]:')
        self.swValue = QtWidgets.QLineEdit('100')
        grid.addWidget(self.swLabel,1,0)
        grid.addWidget(self.swValue,1,1)

        self.offsetLabel = QtWidgets.QLabel('Offset [kHz]:')
        self.offsetValue = QtWidgets.QLineEdit('0')
        grid.addWidget(self.offsetLabel,2,0)
        grid.addWidget(self.offsetValue,2,1)
        self.acqtimeLabel = QtWidgets.QLabel('Acquisition Time [s]:')
        self.acqtimeValue = QtWidgets.QLineEdit('1')
        grid.addWidget(self.acqtimeLabel,3,0)
        grid.addWidget(self.acqtimeValue,3,1)
        grid.setColumnStretch(20, 1)

    def getSettings(self):
        """
        Returns all the settings in a dictionary or list format.
        To be used when acquiring
        """

class PlotFrame(QtWidgets.QTabWidget):
    def __init__(self, parent):
        super(PlotFrame, self).__init__(parent)
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
    DELAYW = 2
    FIDW = 3
    FIDH = 1
    LINEWIDTH = 2
    FONTSIZE = 12
    TEXTHEIGHT = 1.1

    def __init__(self, parent):
        super(SequenceDiagram, self).__init__(parent)
        self.ax.axis('off')
        self.xpos = 0

        self.setIsotope('1H')
        self.drawDelay('RD')
        self.drawPulse('p1')
        self.drawFid('Acq')
        self.ax.set_xlim([-1,10])
        self.ax.set_ylim([-1,1.5])

    def setIsotope(self,iso):
        self.ax.text(-.2,0,iso,horizontalalignment='right',fontsize=self.FONTSIZE)

    def drawPulse(self,text=None):
        self.ax.add_patch(matplotlib.patches.Rectangle((self.xpos, 0), self.PULSEW, self.PULSEH,linewidth=self.LINEWIDTH,edgecolor='tab:blue'))
        self._addText(text,self.PULSEW)
        self.xpos += self.PULSEW

    def drawDelay(self,text=None):
        self.ax.plot([self.xpos,self.xpos+self.DELAYW],[0,0],c='tab:blue',linewidth=self.LINEWIDTH)
        self._addText(text,self.DELAYW)
        self.xpos += self.DELAYW

    def drawFid(self,text=None):
        samples = 100
        T2 = 0.3
        xdata = np.linspace(0,1,samples)
        ydata = self.FIDH * np.cos(xdata * 30) * np.exp(-xdata / T2)
        self._addText(text,self.FIDW)
        self.ax.plot(xdata * self.FIDW + self.xpos,ydata,c='tab:blue',linewidth=self.LINEWIDTH)
        self.xpos += self.FIDW

    def _addText(self,text,xAdder):
        if text is not None:
            self.ax.text(self.xpos + 0.5*xAdder,self.TEXTHEIGHT,text,horizontalalignment='center',fontsize=self.FONTSIZE)

    def resetPlot(self):
        AbstractPlotFrame.resetPlot(self)
        self.ax.axis('off')
        self.xpos = 0


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
