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
#Set refractie index

mtl = Uniaxial_Material(1.60,1.55)
layer = H_Layers(mtl, 300, 30, 300)
rrr = []
#Iterate through frequencies
for wavelength in np.linspace(200,1000,100):
    layer.set_incidence([0,0,1], wavelength)
    layer.doit()
    rrr.append(layer.prop.RCRR)
pl.cla()
pl.plot(np.linspace(200,1000,100), rrr)