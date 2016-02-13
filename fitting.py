# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 20:07:53 2016
A script trys to fit the spectrum
@author: Bonan
"""

import simClasses as sim
import scipy.optimize as so 
import PitchProfile as pit
import matplotlib.pyplot as pl
import numpy as np
pl.rcParams['figure.figsize'] = (8,6)
pl.rcParams['savefig.dpi'] = 100
#%%
#define materials
cellulose = sim.UniaxialMaterial(1.586,1.524) # this might not be right
air = sim.HomogeneousNondispersiveMaterial(1)
# Set the glass as a practical backend
glass = sim.HomogeneousNondispersiveMaterial(1.55)
front = sim.IsotropicHalfSpace(air)
airback = sim.IsotropicHalfSpace(air)
glassback = sim.IsotropicHalfSpace(glass)
system = sim.OptSystem()
#%%
helix = sim.AnyHeliCoidalStructure(cellulose, 100, 1000)
helix.setPitchProfile(pit.PolyPitchProfile([0,1,0], 172, 1000))
system.setHalfSpaces(front, glassback)
system.setStructure([helix])
system.setIncidence(500,0,0)
#%% Plot spectrum
def merit(target, s):
    """
    A merit function to be minimised
    system: the initalised system to be used. Will try to call:
    s.setIncidence(i,self.Theta,self.Phi)
    s.updateStructurePartialTransfer()
    s.getTransferMatrix()
    
    Target: a 3xN array. First row is the wavelength, second is the measured value
    . The third row is for the tolerance
    """
    target = np.array(target)
    cV = np.array(s.scanSpectrum(target[0])[1]) # calculate the model values
    diff = cV - target[1]
    return np.sum((diff/target[2])**2)
#%% Test using [-0.002,0.0001,0] base 150 thickness 1000
target = np.array([[  4.0061224e+02,   4.2469388e+02,   4.5877551e+02,
          5.14285714e+02,   5.7693878e+02,   5.9102041e+02,
          6.1510204e+02,   6.55918367e+02],
       [  1.19886818e-02,   5.08965969e-02,   9.49555603e-02,
          1.08123569e-01,   8.82413602e-02,   5.35027923e-02,
          2.36891293e-02,   7.83742888e-03],
       [  1.00000000e+00,   1.00000000e+00,   1.00000000e+00,
          1.00000000e+00,   1.00000000e+00,   1.00000000e+00,
          1.00000000e+00,   1.00000000e+00]])
# Plot the target profile 
target[0] = target[0] + 30
pl.scatter(target[0], target[1])
#%%
def minFunc(coeff):
    helix.pitchProfile.setCoeff(np.append(coeff, [0]))
    return merit(target, system)
    
optResult = so.fmin_powell(minFunc, [0,0,0,0,0], full_output = 1)
#%%Plot the function
scanResult = system.scanSpectrum(np.linspace(400,700,50))
pl.plot(scanResult[0],scanResult[1]) 