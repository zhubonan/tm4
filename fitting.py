# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 20:07:53 2016
A script trys to fit the spectrum
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
helix.setPitchProfile(pit.PolyPitchProfile([0,1,0], 150, 1000))
system.setHalfSpaces(front, glassback)
system.setStructure([helix])
system.setIncidence(500,0,0)
#%% Plot spectrum
def merit(s, target):
    """
    A merit function to be minimised
    system: the initalised system to be used. Will try to call:
    s.setIncidence(i,self.Theta,self.Phi)
    s.updateStructurePartialTransfer()
    s.getTransferMatrix()
    
    Target: a 3xN array. First row is the wavelength, second is the measured value
    . The third row is for the tolerance
    """
    target = np.array(target)
    cV = np.array(s.scanSpectrum(target[0])[2]) # calculate the model values
    diff = cV - target[1]
    return np.sum((diff/target[2])**2)

x = merit(system, [[440,450,460],  [0.15, 0.14, 0.12]
, [0.01,0.01,0.01]])