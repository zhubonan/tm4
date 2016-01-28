# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 14:38:34 2016
A script to reproduce the data in Digital Color in Cellulose Nanocrystal Films
Ahu Gümrah Dumanli,*,†,‡ Hanne M. van der Kooij,† Gen Kamita,† Erwin Reisner,‡ Jeremy J. Baumberg,† Ullrich Steiner,†,§ and Silvia Vignolini*,‡
†Cavendish
@author: Bonan
"""

import numpy as np
import simClasses as mul
import matplotlib.pyplot as pl
#%%new classes
material = mul.UniaxialMaterial(1.7,1.5)
helix1= mul.HeliCoidalStructure(material, 325, 25 ,325*5)
air = mul.HomogeneousNondispersiveMaterial(1.6)
front = mul.IsotropicHalfSpace(air)
back = mul.IsotropicHalfSpace(air)
system = mul.OptSystem()
system.setHalfSpaces(front,back)
system.setStructure([helix1])
#%%Test selfplotting
wlRange = np.linspace(600,1600,100)
result = system.scanSpectrum(wlRange)
pl.plot(result[1], result[2])