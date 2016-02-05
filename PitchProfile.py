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

pi = np.pi
#%%
class PitchProfile():
    """An abstract basis class for pitch profile"""
    
    def getAngles(): 
        raise NotImplementedError
    
    def getPitch():
        raise NotImplementedError
        
    def plotPitch(self, n= 50):
        """Plot the pitch against depth"""
        depth = np.linspace(0,self.totalThickness, n)
        pl.plot(depth, self.getPitch(depth))
        
    def plotAngles(self, n = 50):
        """Plot the angles agasint depth"""
        depth = np.linspace(0,self.totalThickness, n)
        pl.plot(depth, self.getAngles(depth)) 
        
 #%%        
class AnyPitchProfile(PitchProfile):
    """A class represent the pitch profile of a generalised helix structure"""
    def __init__(self, pitchList, totalThickness ,handness = 'left'):
        """Initialise by passing pitch profile and total thickness"""
        self.pitchCallalbe = interp1d(pitchList[0], pitchList[1], 
        bounds_error = False, fill_value = 0)
        self.handness = handness
        self.totalThickness
        
    def _angleFunction(self, z):
        """The integrand of the angle(z) function"""
        if self.handness == 'left':
            return -np.pi / self.pitchCallalbe(z)
        elif self.handness == 'right':
            return np.pi / self.pitchCallalbe(z)
        else: raise RuntimeError('Need to be either left or right handded')

    def getAngles(self, zList):
        """Get angles of a series of z positions"""
        # May need to change the integrate method in order to get best performance
        zAngles = [sp.integrate.fixed_quad(self._angleFunction, 0, z)[0] for z in zList]
        return zAngles
        
    def getPitch(self, z):
        """Get the pitch at depth z"""
        return self.pitchCallalbe(z)
        

#%%
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
        
        
    def getAngle(self, zList):
        """Calculate the angle of rotation"""
        #Consider the complete segments' rotations
        if type(zList) == int:
            zList = [zList]
            angles = []
        for z in zList:    
            pitches = self.delta[np.arange(np.searchsorted(self.stepLocation,z))] + self.basePitch
            depthes = self.delta[np.arange(np.searchsorted(self.stepLocation,z))]
            # get segmentLength
            segment = []
            lastdepth = 0
            for depth in depthes:
                segment.append(depthes - lastdepth)
                lastdepth= depth
            # set the correct length for the last segment
            segment = np.array(segment)
            segment[-1] = z - depthes[-2]
            angles.append(np.sum(segment / pitches * np.pi * -1))
        return angles
#%%
class PolyPitchProfile(AnyPitchProfile):
    """A class that gives an quadratic profile of the pitch"""
    integOrder = 5  #default order of fixed_quad integration
    def __init__(self, polycoeff, basePitch, totalThickness, handness = 'left'):
        """Initialise the class by passing coefficient of polynomial, basePitch and
        total thickness of the profile
        
        polycoeff: a list contain the coefficient of the variation polynomial. The
        variation polynmial is added to the inverse of the pitch. poly = delta(1/pitch)
        
        basePitch: basis pitch we use
        
        totalThickness: totalThickness of the nematic layer
        """
        polycoeff = np.array(polycoeff)
        self.pitchPoly = sp.poly1d(polycoeff/totalThickness)
        self.basePitch = basePitch
        self.totalThickness = totalThickness
        self.handness = handness
    def getPitch(self, zList):
        """Get the pitch at a list of z location. """
        return self.pitchPoly(zList) + self.basePitch

    def getAngles(self, zList):
        if self.handness == 'left':
            sign = -1
        elif self.handness == 'right':
            sign= 1
        else:
            raise RuntimeError('handness is either left or right')
        angles = []
        pitchIntegrand = lambda x: 1/(self.pitchPoly(x) + self.basePitch) * sign
        for z in zList:
            # We use fast fixed qudrature method here
            angles.append(sp.integrate.fixed_quad(pitchIntegrand, 0, z, 
                                                  n = self.integOrder)[0] * pi)
        return angles
        
    def changeCoeff(self, newCoeff):
        """Change the coefficient of the variation polynomial
        note that VariedPitch = (1/basePitch + varPoly)**-1
        """
        self.pitchPoly = sp.poly1d(newCoeff)
        
    def setIntegOrder(self, order):
        self.integOrder = order
        
if __name__ == '__main__':
    
    #p = StepPitchProfile([250,500,750],[10,10,0],500,1000)
    #p.plotPitch(200)
    #pl.ylim(400,600)
    #pl.xlabel('Depth /nm')
    #pl.ylabel('Pitch /nm')
    
    p = PolyPitchProfile([0.1,0.1,0], 500, 1000)
    #p.plotPitch()
    p.plotAngles()