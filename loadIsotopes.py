#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def fOrNone(inp):
    """Converts a string to a float and dashes to None"""
    if inp == '-':
        return None
    return float(inp)

def getIsotopes(isoPath):
    with open(isoPath) as isoFile:
        isoList = [line.strip().split('\t') for line in isoFile]
    isoList = isoList[1:] #Cut off header

    spins = dict()
    for i, _ in enumerate(isoList):
        quant = fOrNone(isoList[i][4])
        freqRatio = fOrNone(isoList[i][5])
        name = isoList[i][3]+isoList[i][1]
        spins[name] = [quant,freqRatio]
    return spins


