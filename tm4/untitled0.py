# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 20:14:12 2016
Testing the numberical stability of TM4 package
Note the transfer matrix formulation used here is by Yeh but should produce same
result as Berreman's
@author: Bonan Zhu
"""
import tm4.simClasses as sim
import numpy as np
import matplotlib.pyplot as pl
for stack in [100,200,400,600,4000,4300]:
    m = sim.UniaxialMaterial(1.416,1.424) # material that makes the helix
    sub = sim.HomogeneousNondispersiveMaterial(1.42) # material of the two half spaces
    h = sim.HeliCoidalStructure(m, 200,stack*200,d = 30) # setup the helicoidal structure
    front,back = sim.IsotropicHalfSpace(sub), sim.IsotropicHalfSpace(sub)
    s = sim.OptSystem()
    s.setHalfSpaces(front, back)
    s.setStructure([h])
    wlRange = np.linspace(540,580,201)
    res = s.scanSpectrum(wlRange)
    pl.plot(res[0],res[1],label = 'Stack:' + str(stack))
pl.xlim(540,580)
pl.legend(loc = 2)
