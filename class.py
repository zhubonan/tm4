# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 17:02:54 2015
class definitions
@author: Bonan
"""
from scipy.interpolate import interp1d
import numpy as np
from mathfunc import *
class material:
    """
    An class object that represent a type of material
    For now the material is an uniaxial birefringent. Optical axis is aligned with
    x direction in the lab frame
    """
    
    def __init__(self, known_e, known_o, kind = 'quadratic'):
        """
        Input some known points of the referactive index with respect to wavelength.
        known_e and known_o should be 2 x N array of known values of ne and no.
        """
        self._known_e = np.asarray(known_e)
        self._known_o = np.asarray(known_o)
        self._build_interp()
    
    def _build_interp(self):
        # assigne function for interpolated
        self.fe = interp1d(self._known_e[0], self._known_e[1], kind)
        self.fo = interp1d(self._known_o[0], self._known_o[1], kind)
    
    def e_diag(self, wavelength):
        """
        Return the calcuated dielectric tensor at given wave length.
        Optical axis is along the x direction by default
        """
        e = np.diag([self.fe(wavelength), self.fo(wavelength), self.fo(wavelength)])
        return e

        
class layers():
    """
    An class object represent anisotropic multilayer filter
    """
    def __init__(self, material, pitch, layer_thickness, total_thickness):
        self._material = material
        self._pitch = pitch
        self._layer_thickness = layer_thickness
        self._total_thickness = total_thickness
        
    def update_e(self, wavelength):
        """
        calculate the dielectric matrix for all layers
        """
        self.e = construct_epsilon(self._material(wavelength), self._pitch, 
                                   self._layer_thickness, self._total_thickness)
    def update_p(self, )    
        
        
        
        
        