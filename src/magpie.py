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
import pandas as pd
from safeEval import safeEval
import helperFunctions as helpFie

import matplotlib
# First import matplotlib and Qt
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import loadIsotopes
import simulator
import sample
VERSION = '0.0'

TIMEDELAY = 500 # ms

ISOTOPES = loadIsotopes.getIsotopes('IsotopeProperties')
NUCLEI = [x for x in ISOTOPES.keys() if not x.startswith('-')]

class MainProgram(QtWidgets.QMainWindow):

    def __init__(self, root):
        super(MainProgram, self).__init__()
        #self.setAcceptDrops(True)
        self.main_widget = QtWidgets.QWidget(self)
        self.setCentralWidget(self.main_widget)

        self.simulator = simulator.Simulator()
        self.pulseSeqName = None
        self.sampleName = None
        self.numScans = 1
        self.timer = None
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

    def loadPulseSeq(self):
        fname, ftype = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '../pulseSeqs/', 'Pulse sequence (*.csv)')
        if fname:
            self.simulator.loadPulseSeq(fname)
            self.pulseSeqName = os.path.splitext(os.path.basename(fname))[0]
            self.drawPulseSeq()
            self.spectrometerFrame.upd()

    def loadSample(self):
        fname, ftype = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '../samples/', 'Sample (*.txt)')
        if fname:
            self.simulator.setSample(sample.loadSampleFile(fname))
            self.sampleName = os.path.splitext(os.path.basename(fname))[0]
            self.spectrometerFrame.upd()
            
    def drawPulseSeq(self):
        self.plotFrame.drawPulseSeq(self.simulator.settings['observe'], self.simulator.pulseSeq)
        self.paramFrame.setPulseSeq(self.simulator.pulseSeq)
        
    def getPulseSeq(self):
        return self.simulator.pulseSeq

    def getSettings(self):
        return self.simulator.settings
    
    def simulate(self):
        self.stop()
        self.spectrometerFrame.setRunning(True)
        parameters = self.paramFrame.getParameters()
        nuclei, field, self.numScans = self.spectrometerFrame.getSettings()
        self.simulator.setSettings(field, nuclei)
        self.simulator.setPulseSeq(parameters)
        self.simulator.reset()
        self.iScan = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.simulateScan)
        self.timer.start(TIMEDELAY)

    def simulateScan(self):
        self.iScan += 1
        self.simulator.scan()
        self.plotData()
        if self.iScan >= self.numScans:
            self.stop()

    def stop(self):
        self.spectrometerFrame.setRunning(False)
        if self.timer is not None:
            self.timer.stop()
            self.timer = None
            
    def plotData(self):
        time, FIDarray = self.simulator.getData()
        self.plotFrame.plotData(time, FIDarray)

        
class ExportFrame(QtWidgets.QWidget):
    def __init__(self, main):
        super(ExportFrame, self).__init__(main)
        self.main = main
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
    def __init__(self, main):
        super(ParameterFrame, self).__init__(main)
        self.main = main
        self.parWidgets = []
        self.currentChanged.connect(self.tabChanged)

    def setPulseSeq(self, pulseSeq, **kwargs):
        self.reset()
        if pulseSeq is None:
            return
        for ind, pulseStep in pulseSeq.iterrows():
            if pulseStep['type'] == 'sat':
                stepWidget = ParameterWidget(self, pulseStep)
            elif pulseStep['type'] == 'delay':
                stepWidget = DelayWidget(self, pulseStep)
            elif pulseStep['type'] in ['pulse', 'shapedPulse']:
                stepWidget = PulseWidget(self, pulseStep)
            elif pulseStep['type'] == 'FID':
                stepWidget = AcqWidget(self, pulseStep)
            self.parWidgets.append(stepWidget)
            self.addTab(self.parWidgets[-1], pulseStep['name'])
        
    def tabChanged(self, index):
        self.main.plotFrame.sequenceFrame.highlightElement(index)

    def getParameters(self):
        pars = [wid.returnValues() for wid in self.parWidgets]
        return pd.DataFrame(pars)

    def reset(self):
        self.parWidgets = []
        for _ in range(self.count()):
            self.removeTab(0)

class ParameterWidget(QtWidgets.QWidget):
    def __init__(self, parent, pulseStep):
        super(ParameterWidget, self).__init__(parent)
        self.pulseStep = pulseStep
        self.grid = QtWidgets.QGridLayout(self)
        self.grid.setColumnStretch(10, 1)

    def returnValues(self):
        return self.pulseStep

class DelayWidget(ParameterWidget):
    def __init__(self, parent, pulseStep):
        super(DelayWidget, self).__init__(parent, pulseStep)
        self.grid.addWidget(QtWidgets.QLabel('Duration [s]:'), 0, 0)
        self.time = QtWidgets.QLineEdit(str(self.pulseStep['time']))
        self.grid.addWidget(self.time, 0, 1)

    def returnValues(self):
        self.pulseStep['time'] = safeEval(self.time.text())
        return self.pulseStep

class PulseWidget(ParameterWidget):
    def __init__(self, parent, pulseStep):
        super(PulseWidget, self).__init__(parent, pulseStep)
        self.grid.addWidget(QtWidgets.QLabel('Duration [µs]:'), 0, 0)
        self.time = QtWidgets.QLineEdit(str(1e6*self.pulseStep['time']))
        self.grid.addWidget(self.time, 0, 1)

        self.grid.addWidget(QtWidgets.QLabel('Amplitude [kHz]:'), 0, 2)
        self.amplitude = QtWidgets.QLineEdit(str(1e-3*self.pulseStep['amp']))
        self.grid.addWidget(self.amplitude, 0, 3)

    def returnValues(self):
        timeVal = safeEval(self.time.text())
        if isinstance(timeVal, (int, float)):
            self.pulseStep['time'] = 1e-6 * timeVal
        else:
            self.pulseStep['time'] = [1e-6*i for i in timeVal]
        ampVal = safeEval(self.amplitude.text())
        if isinstance(ampVal, (int, float)):
            self.pulseStep['amp'] = 1e3 * ampVal
        else:
            self.pulseStep['amp'] = [1e3*i for i in ampVal]
        return self.pulseStep        

class AcqWidget(ParameterWidget):
    # TODO: recalculate the values after filling in a field
    def __init__(self, parent, pulseStep):
        super(AcqWidget, self).__init__(parent, pulseStep)
        self.grid.addWidget(QtWidgets.QLabel('Offset [kHz]:'), 0, 0)
        self.offset = QtWidgets.QLineEdit('0')
        self.grid.addWidget(self.offset, 0, 1)
                
        self.grid.addWidget(QtWidgets.QLabel('Spectral Width [kHz]:'), 0, 2)
        self.sw = QtWidgets.QLineEdit(str(1e-3*self.pulseStep['amp']/self.pulseStep['time']))
        self.sw.editingFinished.connect(self.swChanged)
        self.grid.addWidget(self.sw, 0, 3)

        self.grid.addWidget(QtWidgets.QLabel('# of points:'), 0, 4)
        self.np = QtWidgets.QLineEdit(str(int(self.pulseStep['amp'])))
        self.np.editingFinished.connect(self.npChanged)
        self.grid.addWidget(self.np, 0, 5)

        self.grid.addWidget(QtWidgets.QLabel('Acq. time [s]:'), 0, 6)
        self.time = QtWidgets.QLabel('{:.6g}'.format(self.pulseStep['time']))
        self.time.setAlignment(QtCore.Qt.AlignCenter)
        self.grid.addWidget(self.time, 0, 7)
        
    def swChanged(self):
        sw = safeEval(self.sw.text())
        self.sw.setText(str(sw))
        points = np.floor(safeEval(self.np.text()))
        self.time.setText('{:.6g}'.format(points / (sw*1000)))
        
    def npChanged(self):
        points = np.floor(safeEval(self.np.text()))
        self.np.setText(str(int(points)))
        sw = safeEval(self.sw.text())
        self.time.setText('{:.6g}'.format(points / (sw*1000)))
        
    def returnValues(self):
        self.pulseStep['time'] = float(self.time.text())
        self.pulseStep['amp'] = int(self.np.text())
        self.pulseStep['offset'] = float(self.offset.text())
        return self.pulseStep
    
class PlotFrame(QtWidgets.QTabWidget):
    
    def __init__(self, main):
        super(PlotFrame, self).__init__(main)
        self.main = main
        self.specFrame = SpecPlotFrame(self)
        self.fidFrame = FidPlotFrame(self)
        self.sequenceFrame = SequenceDiagram(self, self.main)
        self.addTab(self.sequenceFrame, 'Pulse Sequence')
        self.addTab(self.fidFrame, 'FID')
        self.addTab(self.specFrame, 'Spectrum')

    def drawPulseSeq(self, *args):
        self.sequenceFrame.drawPulseSeq(*args)

    def plotData(self, time, FIDarray):
        parameters = self.main.paramFrame.getParameters()
        Offset = float(parameters.loc[parameters['name'] == 'Acq']['offset']) * 1000
        ppmHz = self.main.simulator.settings['B0'] * helpFie.getGamma(self.main.simulator.settings['observe']) # Amount of Hz per ppm, e.g. mainFreq/1e6     
        
        self.fidFrame.clearPlot()
        self.specFrame.clearPlot()
        freq = np.fft.fftshift(np.fft.fftfreq(len(time), time[1]-time[0]))
        freqppm = freq/ppmHz + Offset/ppmHz
        spec = np.fft.fftshift(np.fft.fft(FIDarray, axis=1), axes=1)
        self.fidFrame.plot(time, np.real(FIDarray).T)
        self.specFrame.plot(freqppm,freq, np.real(spec).T)
        
class AbstractPlotFrame(QtWidgets.QWidget):
    MIRRORX = False
    
    def __init__(self, parent):
        super(AbstractPlotFrame, self).__init__(parent)
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(self.canvas, 0, 0)
        
        self.canvas.mpl_connect('button_press_event', self.buttonPress)
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()
        self.rect = [None, None, None, None]  # lines for zooming
        self.leftMouse = False  # is the left mouse button currently pressed
        self.rightMouse = False  # is the right mouse button currently pressed
        self.panX = None  # start position of dragging the spectrum
        self.panY = None  # start position of dragging the spectrum
        self.xmaxlim = None
        self.xminlim = None
        self.ymaxlim = None
        self.yminlim = None
        
        self.zoomX1 = None
        self.zoomX2 = None
        self.zoomY1 = None
        self.zoomY2 = None
        
        self.xdata = None
        self.ydata = None
        
        self.clearPlot()

    def scroll(self, event):
        if self.xdata is None or self.ydata is None:
            return
        modifiers = QtWidgets.QApplication.keyboardModifiers()
        if self.rightMouse:
            middle = (self.xmaxlim + self.xminlim) / 2.0
            width = self.xmaxlim - self.xminlim
            if modifiers == QtCore.Qt.ControlModifier:
                width = width * 0.6**event.step
            else:
                width = width * 0.9**event.step
            self.xmaxlim = middle + width / 2.0
            self.xminlim = middle - width / 2.0
            if self.MIRRORX:
                self.ax.set_xlim(self.xmaxlim, self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim, self.xmaxlim)
        else:
            if modifiers == QtCore.Qt.ControlModifier:
                self.ymaxlim *= 0.6**event.step
                self.yminlim *= 0.6**event.step
            else:
                self.ymaxlim *= 0.9**event.step
                self.yminlim *= 0.9**event.step
            self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.update()
        self.canvas.draw_idle()
            
    def buttonPress(self, event):
        if self.xdata is None or self.ydata is None:
            return
        inv = self.ax.transData.inverted()
        x, y = inv.transform((event.x, event.y))
        if event.button == 1:
            self.leftMouse = True
            self.zoomX1 = x
            self.zoomY1 = y
        elif (event.button == 3) and event.dblclick:
            self.plotReset()
        elif event.button == 3:
            self.rightMouse = True
            self.panX = x
            self.panY = y

    def buttonRelease(self, event):
        if self.xdata is None or self.ydata is None:
            return
        if event.button == 1:
            self.leftMouse = False
            try:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                if self.rect[1] is not None:
                    self.rect[1].remove()
                if self.rect[2] is not None:
                    self.rect[2].remove()
                if self.rect[3] is not None:
                    self.rect[3].remove()
            finally:
                self.rect = [None, None, None, None]
            if self.zoomX2 is not None and self.zoomY2 is not None:
                self.xminlim = min([self.zoomX1, self.zoomX2])
                self.xmaxlim = max([self.zoomX1, self.zoomX2])
                self.yminlim = min([self.zoomY1, self.zoomY2])
                self.ymaxlim = max([self.zoomY1, self.zoomY2])
                if self.MIRRORX:
                    self.ax.set_xlim(self.xmaxlim, self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim, self.xmaxlim)
                self.ax.set_ylim(self.yminlim, self.ymaxlim)
            self.zoomX1 = None
            self.zoomX2 = None
            self.zoomY1 = None
            self.zoomY2 = None
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw_idle()
     
    def pan(self, event):
        if self.xdata is None or self.ydata is None:
            return
        modifiers = QtWidgets.QApplication.keyboardModifiers()
        if self.rightMouse and self.panX is not None and self.panY is not None:
            inv = self.ax.transData.inverted()
            point = inv.transform((event.x, event.y))
            diffx = point[0] - self.panX
            diffy = point[1] - self.panY
            if modifiers == QtCore.Qt.ControlModifier:
                self.xmaxlim = self.xmaxlim - diffx
                self.xminlim = self.xminlim - diffx
            elif modifiers == QtCore.Qt.ShiftModifier:
                self.ymaxlim = self.ymaxlim - diffy
                self.yminlim = self.yminlim - diffy
            else:
                self.xmaxlim = self.xmaxlim - diffx
                self.xminlim = self.xminlim - diffx
                self.ymaxlim = self.ymaxlim - diffy
                self.yminlim = self.yminlim - diffy
            if self.MIRRORX:
                self.ax.set_xlim(self.xmaxlim, self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim, self.xmaxlim)
            self.ax.set_ylim(self.yminlim, self.ymaxlim)
            self.canvas.draw_idle()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            inv = self.ax.transData.inverted()
            point = inv.transform((event.x, event.y))
            self.zoomX2 = point[0]
            self.zoomY2 = point[1]
            if self.rect[0] is not None:
                try:
                    if self.rect[0] is not None:
                        self.rect[0].remove()
                    if self.rect[1] is not None:
                        self.rect[1].remove()
                    if self.rect[2] is not None:
                        self.rect[2].remove()
                    if self.rect[3] is not None:
                        self.rect[3].remove()
                finally:
                    self.rect = [None, None, None, None]
            self.rect[0], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY2, self.zoomY2], 'k', clip_on=False)
            self.rect[1], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY1, self.zoomY1], 'k', clip_on=False)
            self.rect[2], = self.ax.plot([self.zoomX1, self.zoomX1], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.rect[3], = self.ax.plot([self.zoomX2, self.zoomX2], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.canvas.draw_idle()
     
    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        if self.xdata is None or self.ydata is None:
            return
        miny = min(np.real(self.ydata))
        maxy = max(np.real(self.ydata))
        differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if xReset:
            self.xminlim = min(self.xdata)
            self.xmaxlim = max(self.xdata)
        if self.MIRRORX:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)

    def clearPlot(self):
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)

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

    def __init__(self, plotframe, main):
        super(SequenceDiagram, self).__init__(plotframe)
        self._setRelLinewidth()
        self.plotframe = plotframe
        self.main = main
        self.ax.axis('off')
        self.xpos = 0
        self.elems = []
        self.backElems = [] #Holds transparent background images
        self.canvas.mpl_connect('pick_event', self.pickHandler)

    def drawPulseSeq(self, isotope, pulseSeq):
        self.clearPlot()
        self.setIsotope(isotope)
        if pulseSeq is None:
            return
        for ind, pulseStep in pulseSeq.iterrows():
            if pulseStep['type'] == 'sat':
                self.drawSaturation(pulseStep['name'])
            elif pulseStep['type'] == 'delay':
                self.drawDelay(pulseStep['name'])
            elif pulseStep['type'] == 'pulse':
                self.drawPulse(pulseStep['name'])
            elif pulseStep['type'] == 'FID':
                self.drawAcq(pulseStep['name'])
        self.ax.set_xlim([-1,10])
        self.ax.set_ylim([-1,1.5])
        
    def _setRelLinewidth(self):
        linescaler = self.ax.transData.transform([(0, 0), (1, 1)])
        linescaler = linescaler[0][0] - linescaler[0][1]
        self.relLineWidth = self.LINEWIDTH/linescaler

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
        elem = self.ax.add_patch(matplotlib.patches.Rectangle((self.xpos, 0), self.PULSEW, self.PULSEH,linewidth=self.LINEWIDTH,color='black',picker=True))
        self.elems.append([elem,'black'])
        self._addText(text,self.PULSEW)
        #Add transparent background
        self._addBackgroundRect(self.PULSEW)
        self.xpos += self.PULSEW + self.relLineWidth
 
    def drawShapedPulse(self,text=None):
        samples = 100
        xdata = np.linspace(-2*np.pi,2*np.pi,samples)
        ydata = self.PULSEH * np.sin(xdata)/xdata 
        elem = self.ax.plot((xdata + 2*np.pi) * self.SHAPEW / (4*np.pi) + self.xpos,ydata,c='black',linewidth=self.LINEWIDTH,picker=True)
        self.elems.append([elem[0],'black'])
        self._addText(text,self.SHAPEW)
        self._addBackgroundRect(self.PULSEW)
        self.xpos += self.SHAPEW + self.relLineWidth

    def drawSaturation(self,text=None):
        elem = self.ax.add_patch(matplotlib.patches.Rectangle((self.xpos, 0), self.SATW, self.SATH,linewidth=self.LINEWIDTH,color='black',hatch='////',fill=False,picker=True))
        self.elems.append([elem,'black'])
        self._addText(text,self.SATW)
        self._addBackgroundRect(self.SATW)
        self.xpos += self.SATW + self.relLineWidth

    def drawDelay(self,text=None):
        elem = self.ax.plot([self.xpos,self.xpos+self.DELAYW],[0,0],c='black',linewidth=self.LINEWIDTH,picker=True)
        self.elems.append([elem[0],'black'])
        self._addText(text,self.DELAYW)

        self._addBackgroundRect(self.DELAYW)
        self.xpos += self.DELAYW + self.relLineWidth

    def drawAcq(self,text=None):
        samples = 100
        T2 = 0.3
        xdata = np.linspace(0,1,samples)
        ydata = self.FIDH * np.cos(xdata * 30) * np.exp(-xdata / T2)
        self._addText(text,self.FIDW)
        elem = self.ax.plot(xdata * self.FIDW + self.xpos,ydata,c='black',linewidth=self.LINEWIDTH,picker=True)
        self.elems.append([elem[0],'black'])
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
        self.main.paramFrame.setCurrentIndex(pos)

    def clearPlot(self):
        AbstractPlotFrame.clearPlot(self)
        self.ax.axis('off')
        self.xpos = 0
        self.elems = []
        self.backElems = []


class FidPlotFrame(AbstractPlotFrame):

    def plot(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata
        self.ax.plot(xdata, ydata)
        self.ax.set_xlim(np.min(xdata), np.max(xdata))
        self.ax.set_xlabel('Time [s]')
        self.ax.set_ylabel('Intensity [arb. u.]')
        
        self.xmaxlim = np.max(xdata)
        self.xminlim = np.min(xdata)
        self.ymaxlim = np.max(ydata)
        self.yminlim = np.min(ydata)
        self.plotReset()
        self.canvas.draw()


class SpecPlotFrame(AbstractPlotFrame):
    MIRRORX = True

    def plot(self, xdata, xdata2, ydata):
        self.xdata = xdata
        self.ydata = ydata
        self.ax.plot(xdata, ydata)
        self.ax2 = self.ax.twiny()
        self.ax.set_xlim(np.max(xdata), np.min(xdata))
        self.ax2.set_xlim([np.max(xdata2)/1000, np.min(xdata2)/1000])
        self.ax.set_xlabel('Shift [ppm]')
        self.ax.set_ylabel('Intensity [arb. u.]')
        self.ax2.set_xlabel('Frequency [kHz]')
        self.xmaxlim = np.max(xdata)
        self.xminlim = np.min(xdata)
        self.ymaxlim = np.max(ydata)
        self.yminlim = np.min(ydata)
        self.plotReset()
        self.canvas.draw()


class SpectrometerFrame(QtWidgets.QFrame):
    def __init__(self, main):
        super(SpectrometerFrame, self).__init__(main)
        self.main = main
        self.setFrameShape(1)
        grid = QtWidgets.QGridLayout(self)
        self.setLayout(grid)
        self.title = QtWidgets.QLabel('<b>Spectrometer settings:</b>')
        self.title.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(self.title,0,0,1,2)
        self.loadSample = QtWidgets.QPushButton('Load sample')
        self.loadSample.clicked.connect(self.main.loadSample)
        grid.addWidget(self.loadSample,1,0,1,2)
        self.sampleLabel = QtWidgets.QLabel()
        self.sampleLabel.setAlignment(QtCore.Qt.AlignCenter)
        grid.addWidget(self.sampleLabel,2,0,1,2)
        self.loadSequence = QtWidgets.QPushButton('Load pulse sequence')
        self.loadSequence.clicked.connect(self.main.loadPulseSeq)
        grid.addWidget(self.loadSequence,3,0,1,2)
        self.sequenceLabel = QtWidgets.QLabel()
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
        self.b0List = [7.0, 9.4, 11.7, 14.1, 18.8]
        self.b0Drop.addItems([f'{i:.1f}T' for i in self.b0List])
        grid.addWidget(self.b0Drop,7,1)

        grid.addWidget(QtWidgets.QLabel('# Scans:'),8,0)
        self.scanBox = QtWidgets.QSpinBox()
        self.scanBox.setMinimum(1)
        grid.addWidget(self.scanBox,8,1)
        
        self.acquirePush = QtWidgets.QPushButton('Acquire')
        self.acquirePush.clicked.connect(self.main.simulate)
        grid.addWidget(self.acquirePush,9,0,1,2)

        self.stopPush = QtWidgets.QPushButton('Stop')
        self.stopPush.clicked.connect(self.main.stop)
        grid.addWidget(self.stopPush,9,0,1,2)
        self.stopPush.setVisible(False)
        
        #grid.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        grid.setRowStretch(20, 1)
        
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.grid = grid
        self.upd()

    def setRunning(self, running):
        if running:
            self.acquirePush.setVisible(False)
            self.stopPush.setVisible(True)
        else:
            self.acquirePush.setVisible(True)
            self.stopPush.setVisible(False)
            
    def getSettings(self):
        nuclei = self.detectDrop.currentText()
        field = self.b0List[self.b0Drop.currentIndex()]
        nScans = self.scanBox.value()
        return nuclei, field, nScans
        
    def upd(self):
        if self.main.sampleName is None:
            self.sampleLabel.setText('<i>No sample</i>')
        else:
            self.sampleLabel.setText('<i>' + self.main.sampleName + '</i>')
        if self.main.pulseSeqName is None:
            self.sequenceLabel.setText('<i>No sequence</i>')
        else:
            self.sequenceLabel.setText('<i>' + self.main.pulseSeqName + '</i>')

        
if __name__ == '__main__':

    root = QtWidgets.QApplication(sys.argv)
    #root.setWindowIcon(QtGui.QIcon(os.path.dirname(os.path.realpath(__file__)) + '/Icons/Logo.png'))
    mainProgram = MainProgram(root)
    mainProgram.setWindowTitle(f"Magpie - {VERSION}")
    mainProgram.show()
    sys.exit(root.exec_())
