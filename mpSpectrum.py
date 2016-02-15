# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 09:48:41 2016
A test of the multiplrocessing 
@author: Bonan
"""

import numpy as np
import matplotlib.pyplot as pl
from multiprocessing import Pool
import simClasses as sim
#%%
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
h1 = heli(CNC,150,1000)
s.setStructure([h1])
#%%
def work(wlRange):
    return wlRange*2

#%%
p = Pool(processes = 4)
res = p.map(work, range(1,7))
print(res)