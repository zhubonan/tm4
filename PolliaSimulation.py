# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:04:21 2016

@author: Bonan
"""

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
a = [[100,1000, 2000], [1.6,1.6, 1.6]]
b = [[100,1000, 2000], [1.55,1.55, 1.55]]
m = Uniaxial_Material(a,b)
l = H_Layers(m, 250,25, 2000)
incident = [1,0,1]
#%% calculation routine
lam = []
rrr = []
rll = []

for wavelength in np.linspace(300,700, 100):
    l.set_incidence(incident, wavelength)
    l.doit()
    rrr.append(l.prop.RCRR)
    rll.append(l.prop.RCLL)
    lam.append(wavelength)
#%%
plt.figure(dpi = 100)
# for ploting the other results
#plt.plot(lbda*1e9, R_RR, label='R_RR_Berreman')
#plt.plot(lbda*1e9, R_th, label='R_theo').
plt.plot(lam,rrr,label = "R-R")
plt.plot(lam,rll, label = "L-L")
plt.xlabel("Wavelength /nm")
plt.ylabel("Reflectivity")
plt.legend(loc='center right')
plt.title(str(incident)+ "  incidence")

