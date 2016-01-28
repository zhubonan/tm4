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
cellulose = mul.UniaxialMaterial(1.586,1.524)
helix1= mul.HeliCoidalStructure(cellulose, 150,30,600)
helix2 = mul.HeliCoidalStructure(cellulose, 200,30,600)
air = mul.HomogeneousNondispersiveMaterial(1)
front = mul.IsotropicHalfSpace(air)
back = mul.IsotropicHalfSpace(air)
system = mul.OptSystem()
system.setHalfSpaces(front,back)
#%%Test selfplotting

