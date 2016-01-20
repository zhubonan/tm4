# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 22:45:56 2015
Test codes
@author: Bonan
"""

import numpy as np
import matplotlib.pyplot as plt
from multilayer import H_Layers,Uniaxial_Material

#%% Initialise the instance object
a = [[500,300,2000], [1.7,1.7,1.7]]
b = [[500,300,2000], [1.5,1.5,1.5]]
m = Uniaxial_Material(a,b)
l = H_Layers(m, 650,13, 325*5)
incident = [0,0,1]
#%% calculation routine
def calc():
    global lam,rrr,rll
    lam = []
    rrr = []
    rll = []
    for wavelength in np.linspace(650,1600, 100):
        l.set_incidence(incident, wavelength)
        l.doit()
        rrr.append(l.prop.RCRR)
        rll.append(l.prop.RCLL)
        lam.append(wavelength)
    return 

#%%
# for ploting the other results
plt.plot(lbda*1e9, R_RR, label='R_RR_Berreman')
plt.plot(lbda*1e9, R_th, label='R_theo')
calc()
plt.plot(lam,rll,label = "L-L Yeh 50dvision")
l = H_Layers(m, 650,26, 325*5)
calc()
plt.plot(lam,rll,label = "L-L Yeh 25dvision")
l = H_Layers(m, 650,6, 325*5)
calc()
plt.plot(lam,rll,label = "L-L Yeh ~100 dvision")
plt.xlabel("Wavelength /nm")
plt.ylabel("Reflectivity")
plt.legend(loc='center right')
plt.title(str(incident)+ "  incidence")

