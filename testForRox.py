# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:55:12 2016

@author: Bonan
"""

from preset import *
from simClasses import Helix
#%% Investigate the effect of differnt divisions
pl.figure()
for div in np.arange(10,34,3):
    h = heli(CNC,150, 1000, d = div)
    s.setStructure([h])
    res = s.scanSpectrum(wlRange)
    pl.plot(res[0],res[1], label = 'Division = ' + str(div))
pl.legend()

#%% Here we don't speed up the calculation by repeating each pitch
pl.figure()
for div in np.arange(20,50,3):
    h = Helix()
    h.setPhyParas(CNC, 150, 1000, div, -1)
    s.setStructure([h])
    res = s.scanSpectrum(wlRange)
    pl.plot(res[0],res[1], label = 'Division = ' + str(div))
pl.legend()


#%%How the R-R was supressed?
pl.figure()
for div in np.arange(10,34,3):
    h = heli(CNC,150, 1000, d = div)
    s.setStructure([h])
    res = s.scanSpectrum(wlRange, coupling= 'RR')
    pl.plot(res[0],res[1], label = 'div = ' + str(div))
pl.legend()