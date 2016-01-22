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
mean = 1.555
birfringence = 0.062
no = mean - birfringence/2
ne = mean + birfringence/2
#Initialise the system
mtl = Uniaxial_Material(ne,no)
layer = H_Layers(mtl, 150, 10, 1400)
rrr = []
#Iterate through frequencies
for wavelength in np.linspace(400,750,100):
    layer.set_incidence([0,0,1], wavelength)
    layer.doit()
    rrr.append(layer.prop.RCRR)
pl.cla()
pl.plot(np.linspace(400,750,100), rrr)