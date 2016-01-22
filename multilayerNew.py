# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 22:51:47 2016
New class file
@author: Bonan
"""
from scipy.interpolate import interp1d
import scipy as sp
import numpy as np
import mathfuncNew as mfc
#%%
# Propagator for a homogeneous slab of material...
# need to stop using deprecated scipy functions
class Propagator():
    """
     A propagator class for easy access 
    """
    
    def __init__(self, method = 'Pade'):
        """the value of 'method':
        "linear" -> first order approximation of exp()
        "Padé"   -> Padé approximation of exp()
        "Taylor" -> Taylor development of exp()
        "eig"    -> calculation with eigenvalue decomposition
        The default method is Pade. It is the default method in scipy libray for 
        matrix exponential
        """
        self.method = 'Pade'
        
    def __call__(self, Delta, h, k0, q=None):
        """
        'Delta' : Delta matrix of the homogeneous material
        'h' : thickness of the homogeneous slab
        'k0' : wave vector in vacuum, k0 = ω/c
        'q' : order of the approximation method, if useful
        Returns a ndarry of the propagator for the division
        """
        if   self.method == "linear":    return np.identity(4) + 1j * h * k0 * Delta
        elif self.method == "Pade":      return sp.linalg.expm(1j * h * k0 * Delta, q)
        elif self.method == "Taylor":    return sp.linalg.expm3(1j * h * k0 * Delta, q+1)
        elif self.method == "eig":       return sp.linalg.expm2(1j * h * k0 * Delta)
            
#%% Defining Material Classes
class Material:
    """
    An class object that represent any type of anisotropic material
    """
    
    def __init__(self, a, b, c, kind = 'quadratic'):
        """
        Input some known points of the refractive index with respect to wavelength.
        known_e and known_o should be 2 x N array of known values of ne and no.
        """
        if type(a) != list and type(a) != np.ndarray:
            self.fa = lambda x: a
            self.fb = lambda x: b
            self.fc = lambda x: c
            return
            
        self.a , self.b, self.c = np.asarray(a), np.asarray(b), np.asarray(c)
        self.kind = kind
        self._build_interp()
    
    def _build_interp(self):
        # assigne function for interpolated
        self.fa = interp1d(self.a[0], self.a[1], self.kind)
        self.fb = interp1d(self.b[0], self.b[1], self.kind)
        self.fc = interp1d(self.c[0], self.c[1], self.kind)
        
    def e_diag(self, wl):
        """
        Return the calcuated dielectric tensor at given wave length.
        Optical axis is along the x direction by default
        """
        # construct the dielectric constant tensor
        e = np.diag([self.fa(wl), self.fb(wl), self.fc(wl)]
        )**2
        return e

    def __call__(self, wavelength):
        return self.e_diag(wavelength)
        
class UniaxialMaterial(Material):
    """
    An sub-class representing an uniaxial material(a!=b=c)
    """
    def __init__(self, e, o, kind = 'quadratic'):
        Material.__init__(self, e,o,o,kind)

class HomogeneousNondispersive_Material(Material):
    
    def __init__(self,n):
        self.n = n
    def e_diag(self, wavelength):
        return np.diag([self.n,self.n,self.n])**2

class HomogeneousDispersiveMaterial(Material):
    
    def __init__(self,n_array):
        Material.__init__(self, n_array,n_array,n_array, kind= 'quadratic')

#%% Defining Structure class

class Structure():
    """
    A superclass represent a type of structure
    """
    material = None
    def __init__(self):
        raise NotImplementedError
        
    def constructEpsilon():
        """A method to construct the relative dielectric tensor"""
        raise NotImplementedError
        
    def constructDelta():
        """A method to build the delta matrix"""
        raise NotImplementedError
        
    def constructPropogator():
        """A method a build the propagator matrices"""
        raise NotImplementedError
        
    def constructPartialTransfer():
        """A method to build the partial transfer matrix for the layer"""
        raise NotImplementedError
    def setPropagater(self, method = 'eig'):
        return
#        self.propagator = 
        
class HeliCoidalStructure(Structure):
    """A structure class represent the helicoidal structure"""
    
    _hasRemainder = False
    wl = None #dont specify wavelength initially
    Kx = None # reduced wave vector Kx = kx/k0
    propagtor = Propagator() #default propagator
    
    def __init__(self, material, pitch, d, t, handness = 'left'):
        """
        Initialise the structure by passing material pitch, division per 
        pitchand total thickness. Handness is left by default
        """
        self.divisionThickness = pitch / d
        self.setHandness(handness) 
        self.anglesRepeating = np.linspace(0, np.pi * self._handness, d, 
                                           endpoint = False)
        self.nOfReapting = int(t/d)
        remain = np.remainder(t, d)
        remainDivisions = int(remain/self.divisionThickness)
        #discard the overflowing bits
        remain = remainDivisions * self.divisionThickness
        if remain != 0:
            self._hasRemainder = True
            self.anglesRemaining = np.linspace(0, remain/pitch*np.pi*self._handness,
                                               remainDivisions, endpoint = False)
        self.material = material
    def setMethod(self, name = 'eig'):
        """set the propagator"""
        if name == 'eig':
            self.propagtor = mfc.eig
        
    def setHandness(self, H):
        
        if H == 'left':
            self._handness = -1
        elif H == 'right':
            self._handness = 1
        else: raise RuntimeError('Handness need to be either left or right')
        
    def constructEpsilon(self):
        """Build the epsilon for each layer"""
        self.epsilonRepeating = [np.dot(mfc.rotZ(theta),self.material(self.wl)) for 
        theta in self.anglesRepeating]
        self.epsilonRemaining = [np.dot(mfc.rotZ(theta),self.material(self.wl)) for 
        theta in self.anglesRemaining]
        
    def constructDelta(self):
        """Build the Delta matrix in Berreman's formulation"""
        self.deltaRepeating = [mfc.buildDeltaMatrix(e,self.Kx) for e in self.epsilonRepeating]
        self.deltaRepeating = [mfc.buildDeltaMatrix(e,self.Kx) for e in self.epsilonRepeating]
        
    def constructPropagtion(self):
        """Build the propagation matrix"""
        #self.propRepeating = [self.propagtor()]
        
        