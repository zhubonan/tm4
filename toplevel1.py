# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 22:45:56 2015

@author: Bonan
"""

import numpy as np
import matplotlib.pyplot as plt
from class_def import Layers, Material

#%% Initialise the instance object
a = [[200,300,1000], [1,1.2,1.5]]
b = [[200,300,1000], [1,1.3,1.5]]
m = U_Material(a,b)
l = Layers(m, 200, 10, 1000)
#%% calculation routine
lam = []
tss = []
tsp = []

for wavelength in np.linspace(400,900, 100):
    l.set_incidence([1,0,1], wavelength)
    l.doit()
    tss.append(l.coeff_modulus["tss"])
    tsp.append(l.coeff_modulus["tsp"])
    lam.append(wavelength)
#%%
plt.plot(lam,tss)
plt.plot(lam,tsp)
    