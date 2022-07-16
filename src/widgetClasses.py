#!/usr/bin/env python3

# Copyright 2021 - 2022 Bas van Meerten and Wouter Franssen

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

import os
import sys

from PyQt5 import QtGui, QtCore, QtWidgets


class ToolWindow(QtWidgets.QWidget):
    """
    The base class of the toolwindows.
    Implements the basic features shared by all toolwindows.
    Toolwindows inherit this class and alter its behaviour by the constants or by reimplementing functions.
    """

    NAME = ""              # The name displayed in the title of the window
    RESIZABLE = False      # should the window be resizable
    MENUDISABLE = True     # Should the window disable the menu of the main window
    APPLYANDCLOSE = True   # Should the window close after the ok button is pressed
    CANCELNAME = "&Cancel" # The name on the cancel button
    OKNAME = "&Ok"         # The name on the ok button

    def __init__(self, parent):
        """
        Initializes the ToolWindow.

        Parameters
        ----------
        parent : Main1DWindow or AbstractParamFrame
            Parent of the toolwindow.
        """
        super(ToolWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle(self.NAME)
        self.layout = QtWidgets.QGridLayout(self)
        self.grid = QtWidgets.QGridLayout()
        self.box = QtWidgets.QDialogButtonBox()
        self.layout.addLayout(self.grid, 0, 0, 1, 2)
        self.cancelButton = QtWidgets.QPushButton(self.CANCELNAME)
        self.cancelButton.clicked.connect(self.closeEvent)
        self.okButton = QtWidgets.QPushButton(self.OKNAME)
        self.okButton.clicked.connect(self.applyAndClose)
        self.okButton.setFocus()
        self.box.addButton(self.cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        self.box.addButton(self.okButton, QtWidgets.QDialogButtonBox.AcceptRole)
        self.show()
        self.layout.addWidget(self.box, 3, 0)
        if not self.RESIZABLE:
            self.setFixedSize(self.size())
        if self.MENUDISABLE:
            self.father.menuEnable(False)
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height(), 0, 0)


    def applyFunc(self):
        """
        Dummy function for the apply function.
        Should be reimplemented by the toolwindows.
        """
        pass

    def applyAndClose(self):
        """
        Runs the apply function and when APPLYANDCLOSE is set, closes the window.
        """
        self.applyFunc()
        if self.APPLYANDCLOSE:
            self.closeEvent()

    def closeEvent(self, *args):
        """
        Updates the view and closes the toolwindow.

        Parameters
        ----------
        *args
            Any arguments are ignored.
        """
        if self.MENUDISABLE:
            self.father.menuEnable(True)
        self.deleteLater()


