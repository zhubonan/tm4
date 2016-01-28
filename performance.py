# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 21:08:38 2016
Performance test
@author: Bonan
"""

import simClasses as sim
import matplotlib.pyplot as pl
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
helix1= sim.customHeliCoidal(cellulose, 150, 100, 1000)
helix1.propagtor.setMethod('Pade')
system.setHalfSpaces(front, glassback)
system.setStructure([helix1])
system.setIncidence(500,0,0)
system.updateStructurePartialTransfer()
system.getTransferMatrix()