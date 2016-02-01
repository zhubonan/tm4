# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:32:28 2016
Function to compute the angluar profile of a generalised helix structure
@author: Bonan
"""

import numpy as np
import scipy as sp
from scipy.interpolate import interp1d

class PitchProfile():
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

if __name__ == '__main__':
    
    p = PitchProfile([[0,1000],[150,150]])
    print(p.getAngles(np.linspace(0,1000,50)))