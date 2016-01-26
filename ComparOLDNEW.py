# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 14:38:34 2016
A script to reproduce the data in Digital Color in Cellulose Nanocrystal Films
Ahu Gümrah Dumanli,*,†,‡ Hanne M. van der Kooij,† Gen Kamita,† Erwin Reisner,‡ Jeremy J. Baumberg,† Ullrich Steiner,†,§ and Silvia Vignolini*,‡
†Cavendish
@author: Bonan
"""

import numpy as np
import multilayerNew as mul
import matplotlib.pyplot as pl
#%%new classes
cellulose = mul.UniaxialMaterial(1.586,1.524)
helix1= mul.HeliCoidalStructure(cellulose, 150,30,600)
helix2 = mul.HeliCoidalStructure(cellulose, 200,30,600)
air = mul.HomogeneousNondispersiveMaterial(1)
front = mul.IsotropicHalfSpace(air)
back = mul.IsotropicHalfSpace(air)
system = mul.OptSystem()
system.setHalfSpaces(front,back)
#%%Test selfplotting
wlRange = np.linspace(300,800,100)
for i in range(1,8):
    helix1= mul.HeliCoidalStructure(cellulose, 150,30,150*i)
    helix2 = mul.HeliCoidalStructure(cellulose, 200,30,200*i)
    system.setStructure([helix1,helix2])
    system.setIncidence(None,0,0)
    plotData = system.scanSpectrum(wlRange)
    pl.plot(plotData[1], plotData[2], label = str(plotData[0]))
    pl.legend(fontsize = 9)
    pl.xlabel("Wavelength /nm")
    pl.ylabel("Reflectivity")
