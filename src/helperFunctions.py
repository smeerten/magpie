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

import loadIsotopes

ISOTOPES = loadIsotopes.getIsotopes('IsotopeProperties')
GAMMASCALE = 42.576e6/100

def getGamma(isotope):
    return ISOTOPES[isotope][1] * GAMMASCALE * 1e-6 

def getSpinQuant(isotope):
    return ISOTOPES[isotope][0]
