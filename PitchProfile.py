# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:32:28 2016
Function to compute the angluar profile of a generalised helix structure
@author: Bonan
"""

import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import matplotlib.pyplot as pl
class PitchProfile():
    """A class represent the pitch profile of a generalised helix structure"""
    
    def getAngles():
        raise NotImplementedError
        
        
class AnyPitchProfile(PitchProfile):
    """A class represent the pitch profile of a generalised helix structure"""
    def __init__(self, pitchList, handness = 'left'):
        """Initialise by passing pitch profile and total thickness"""
        self.pitchInterp = interp1d(pitchList[0], pitchList[1], 
        bounds_error = False, fill_value = 0)
        self.handness = handness
        
    def _angleFunction(self, z):
        """The integrand of the angle(z) function"""
        if self.handness == 'left':
            return -np.pi / self.pitchInterp(z)
        elif self.handness == 'right':
            return np.pi / self.pitchInterp(z)
        else: raise RuntimeError('Need to be either left or right handded')

    def getAngles(self, zList):
        """Get angles of a series of z positions"""
        # May need to change the integrate method in order to get best performance
        zAngles = [sp.integrate.fixed_quad(self._angleFunction, 0, z)[0] for z in zList]
        return zAngles

class StepPitchProfile(PitchProfile):
    """Describe a pritch profile using a series steps"""
    def __init__(self, stepLocation, stepSize, basePitch, totalThickness):
        """Initialise the profile
        
        * stepList: 2*N list of [[z],[step]]. A adding a step moves segments on the two sides by -step,+step
        
        * basePitch: Basis pitch value
        
        * totalThickness: Total thickness of the profile
        
        """
        self.stepLocation = np.array(stepLocation)
        self.stepSize = np.array(stepSize)
        self.basePitch = basePitch
        self.totalThickness = totalThickness
        self.calcDeltaList()
        
    def calcDeltaList(self):
        """Calculate the change from basis for each segments"""
        delta = []
        lastDelta = 0
        lastStep = 0
        for i in range(self.stepSize.size):
            lastDelta += lastStep*2 - self.stepSize[i]
            delta.append(lastDelta)
            lastStep = self.stepSize[i]
        delta.append(lastDelta + lastStep*2) #added the delta of the final segment
        self.delta = np.array(delta)
    def getPitch(self, z):
        
        #Check if z lie on any point of steps
        return self.delta[np.searchsorted(self.stepLocation, z)]+self.basePitch
        
    def plotPitch(self, n= 50):
        
        depth = np.linspace(0,self.totalThickness, n)
        pl.plot(depth, self.getPitch(depth))
        
if __name__ == '__main__':
    
    p = StepPitchProfile([500],[10],500,1000)
    p.plotPitch()