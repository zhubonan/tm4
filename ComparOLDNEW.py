# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 14:38:34 2016
A script to reproduce the data in Digital Color in Cellulose Nanocrystal Films
Ahu Gümrah Dumanli,*,†,‡ Hanne M. van der Kooij,† Gen Kamita,† Erwin Reisner,‡ Jeremy J. Baumberg,† Ullrich Steiner,†,§ and Silvia Vignolini*,‡
†Cavendish
@author: Bonan
"""

import numpy as np
import matplotlib.pyplot as pl
from multilayer import Uniaxial_Material, H_Layers
import multilayerNew as mul
#Set refractie index

mtl = Uniaxial_Material(1.60,1.55)
layer = H_Layers(mtl, 300, 30, 300)
rrr = []
#Iterate through frequencies
wlRange = np.linspace(200,1000,100)
#%%
for wavelength in wlRange:
    layer.set_incidence([0,0,1], wavelength)
    layer.doit()
    rrr.append(layer.prop.RCRR)
pl.cla()
pl.plot(wlRange, rrr, label = "Old")
#%%new classes
cellulose = mul.UniaxialMaterial(1.6,1.55)
helix= mul.HeliCoidalStructure(cellulose, 150,30,300)
helix2 = mul.HeliCoidalStructure(cellulose, 200,30,300)
air = mul.HomogeneousNondispersiveMaterial(1)
front = mul.IsotropicHalfSpace(air)
back = mul.IsotropicHalfSpace(air)
system = mul.OptSystem()
system.setHalfSpaces(front,back)
system.setStructure([helix,helix2])
#%%
#Iterate through frequencies
rrr = []
for wavelength in wlRange:
    system.setIncidence(wavelength)
    system.updateStructurePartialTransfer()
    system.getTransferMatrix()
    rrr.append(system.prop.RC[0,0])
pl.plot(wlRange, rrr, label = "New")
pl.legend()
    