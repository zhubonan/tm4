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
a = [[500,300,2000], [1.7**2,1.7**2,1.7**2]]
b = [[500,300,2000], [1.5**2,1.5**2,1.5**2]]
m = U_Material(a,b)
l = H_Layers(m, 650,13, 325*5)
incident = [0,0,1]
#%% calculation routine
lam = []
rrr = []
rll = []

for wavelength in np.linspace(650,1600, 100):
    l.set_incidence(incident, wavelength)
    l.doit()
    rrr.append(l.prop.RCRR)
    rll.append(l.prop.RCLL)
    lam.append(wavelength)
#%%
plt.figure(dpi = 100)
# for ploting the other results
#plt.plot(lbda*1e9, R_RR, label='R_RR_Berreman')
#plt.plot(lbda*1e9, R_th, label='R_theo')
plt.plot(lam,rrr,label = "R-R Yeh")
plt.plot(lam,rll, label = "L-L Yeh")
plt.xlabel("Wavelength /nm")
plt.ylabel("Reflectivity")
plt.legend(loc='center right')
plt.title(str(incident)+ "  incidence")

