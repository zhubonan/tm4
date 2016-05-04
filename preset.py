# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 20:09:06 2016

@author: Bonan
"""

import simClasses as sim
import matplotlib.pyplot as pl
import numpy as np
from colourTools import specToRGB
pl.rcParams['figure.figsize'] = (8,6)
pl.rcParams['savefig.dpi'] = 100
#%%
#define materials
CNC = sim.UniaxialMaterial(1.586,1.524) # this might not be right
air = sim.HomogeneousNondispersiveMaterial(1)
cellulose = sim.HomogeneousNondispersiveMaterial(1.55)
# Set the glass as a practical backend
glass = sim.HomogeneousNondispersiveMaterial(1.55)
front = sim.IsotropicHalfSpace(air)
airhalf= sim.IsotropicHalfSpace(air)
glasshalf = sim.IsotropicHalfSpace(glass)
s = sim.OptSystem()
s.setHalfSpaces(airhalf,glasshalf)
heli = sim.HeliCoidalStructure
h1 = heli(CNC,180,1000)
s.setStructure([h1])
wlRange = np.linspace(400,800,100)
"""
print('Followings are added to the scope', flush = True)
print('Materials: CNC, air, cellulosem, glass', flush = True)
print('HalfSpace: airhalf, glasshalf', flush = True)
print('OptSystem:s', flush = True)
print('heli as HeliCoidalStructure', flush = True)
print('wlRange as 400 to 800 nm', flush = True)
"""
def plotSpectrum(OptSys, wlRange):
    result = OptSys.scanSpectrum(wlRange)
    pl.plot(result[0],result[1], color = specToRGB(result))
    return