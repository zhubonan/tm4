# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 22:45:56 2015
Test codes
@author: Bonan
"""

import numpy as np
import matplotlib.pyplot as plt
from multilayer import H_Layers,U_Material

#%% Initialise the instance object
a = [[200,300,1000], [2,2,2]]
b = [[200,300,1000], [2.1,2.1,2.1]]
m = U_Material(a,b)
l = H_Layers(m, 600,50,2000)
incident = [1,1,1]
#%% calculation routine
lam = []
rrr = []
rll = []

for wavelength in np.linspace(400,900, 100):
    l.set_incidence(incident, wavelength)
    l.doit()
    rrr.append(l.coeff_modulus_LR["rpp"])
    rll.append(l.coeff_modulus_LR["rss"])
    lam.append(wavelength)
#%%
plt.plot(lam,rrr,label = "R-R reflectivity")
plt.plot(lam,rll, label = "L-L reflectivity")
plt.xlabel("Wavelength /nm")
plt.ylabel("Reflectivity")
plt.legend(loc = 2)
plt.title(str(incident)+ "  incidence")
    