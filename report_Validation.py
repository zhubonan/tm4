# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 19:12:43 2016

@author: Bonan
"""

import numpy as np
import simClasses as sc
import matplotlib.pyplot as pl
pl.rcParams['figure.figsize'] = (4,3)
import old.multilayer as scOld
from numpy import sqrt, pi, exp
#%%
wlRange = np.linspace(400,700,100)
ne, no = 1.6, 1.5
p = 180
h= 1800
div = 30
#%% Compute results
testMaterial = sc.UniaxialMaterial(ne,no)
air = sc.HomogeneousNondispersiveMaterial(1)
glass =sc.HomogeneousNondispersiveMaterial(1.6)
airHalf =   sc.IsotropicHalfSpace(air)
glassHalf = sc.IsotropicHalfSpace(glass)
matchedHalf = sc.IsotropicHalfSpace(sc.HomogeneousNondispersiveMaterial((no + ne)/2))
s = sc.OptSystem()
s.setHalfSpaces(matchedHalf, matchedHalf)
heli = sc.HeliCoidalStructure(testMaterial, p, h, div)
s.setStructure([heli])
def computeB():
    return s.scanSpectrum(wlRange,1)[1]
#%% Yeh's Method
l = scOld.H_Layers(scOld.Uniaxial_Material(ne,no),p,div,h)
def computeY():
    yRes = []
    for wl in wlRange:
        l.set_incidence([0,0,1], wl, (ne+no)/2,(ne+no)/2)
        l.doit()
        yRes.append(l.prop.RCRR)
    return yRes
#%% Analytical calculation for the power reflection coefficient
def computeT():
    q = pi/p
    alpha = q/(2*pi/wlRange)
    k0 = 2*pi/wlRange
    epsilon = (no**2+ne**2)/2
    delta = (no**2-ne**2)/2
    n2 = sqrt((alpha**2 + epsilon \
         - sqrt(4*epsilon*alpha**2+delta**2)).astype(complex))
    w = 1j*(ne**2-n2**2-alpha**2)/(2*alpha*n2) # not k0/c     
    return abs((w**2+1)*(1-exp(-2j*k0*n2*h)) \
               / (2*w*(1+exp(-2j*k0*n2*h)) \
               - 1j*(w**2-1)*(1-exp(-2j*k0*n2*h))))**2
if __name__ == "__main__":
    #%% Compute
    bRes = computeB()
    yRes = computeY()
    tRes = computeT()
    #%% Plotting
    pl.rcParams['axes.labelsize'] = 9
    pl.rcParams['figure.titlesize'] = 10
    pl.figure(figsize = (5,5*2/3))
    #pl.rcParams['a']
    pl.plot(wlRange, tRes, '-',label = "Theotical")
    pl.plot(wlRange, bRes, '+',label = "Berreman")
    pl.plot(wlRange, yRes, 'x',label = "Teh")
    pl.title("Comparing three models")
    pl.xlabel('Wavelength /nm')
    pl.ylabel('Reflectance L-L')
    pl.legend(loc = 2, fontsize = "small")
    pl.tight_layout()
    pl.savefig("../PartIIIReport/fig/validation.pdf")
