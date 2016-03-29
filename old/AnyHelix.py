# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 21:08:38 2016
Test AnyHelix Class operation 
@author: Bonan
"""

import simClasses as sim
import PitchProfile as pit
import matplotlib.pyplot as pl
import numpy as np
pl.rcParams['figure.figsize'] = (8,6)
pl.rcParams['savefig.dpi'] = 100
#%%
#define materials
cellulose = sim.UniaxialMaterial(1.586,1.524) # this might not be right
air = sim.HomogeneousNondispersiveMaterial(1)
# Set the glass as a practical backend
glass = sim.HomogeneousNondispersiveMaterial(1.55)
front = sim.IsotropicHalfSpace(air)
airback = sim.IsotropicHalfSpace(air)
glassback = sim.IsotropicHalfSpace(glass)
system = sim.OptSystem()
#%%
helix = sim.AnyHeliCoidalStructure(cellulose, 100, 1000)
helix.setPitchProfile(pit.PolyPitchProfile([-0.002,0.0001,0], 150, 1000))
system.setHalfSpaces(front, glassback)
system.setStructure([helix])
system.setIncidence(500,0,0)
#%% Plot spectrum
wlList = np.linspace(300,800,50)
result = system.scanSpectrum(wlList)
pl.plot(result[0],result[1])