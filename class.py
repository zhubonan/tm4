# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 17:02:54 2015
class definitions
@author: Bonan
"""
from scipy.interpolate import interp1d
import scipy.constants as sc
import numpy as np
from mathfunc import *
#define constants

class Material:
    """
    An class object that represent a type of material
    For now the material is an uniaxial birefringent. Optical axis is aligned with
    x direction in the lab frame
    """
    
    def __init__(self, a, b, c, kind = 'quadratic'):
        """
        Input some known points of the relative dielectric constants with respect to wavelength.
        known_e and known_o should be 2 x N array of known values of ne and no.
        """
        self.a , self.b, self.c = np.asarray(a), np.asarray(b), np.asarray(c)
        self.kind = kind
        self._build_interp()
    
    def _build_interp(self):
        # assigne function for interpolated
        self.fa = interp1d(self.a[0], self.a[1], self.kind)
        self.fb = interp1d(self.b[0], self.b[1], self.kind)
        self.fc = interp1d(self.c[0], self.c[1], self.kind)
        
    def e_diag(self, wavelength):
        """
        Return the calcuated dielectric tensor at given wave length.
        Optical axis is along the x direction by default
        """
        # construct the dielectric constant tensor
        e = np.diag([self.fa(wavelength), self.fb(wavelength), self.fc(wavelength)]
        ) * sc.epsilon_0
        return e

    def __call__(self, wavelength):
        return self.e_diag(wavelength)
        
class U_Material(Material):
    """
    An sub-class representing an uniaxial material(a!=b=c)
    """
    def __init__(self, e, o, kind = 'quadratic'):
        Material.__init__(self, e,o,o,kind)
    

class Layers():
    """
    An class object represent anisotropic multilayer filter
    """
    def __init__(self, material, pitch, layer_thickness, total_thickness):
        self.material = material
        self._pitch = pitch
        self._layer_thickness = layer_thickness
        self._total_thickness = total_thickness
        
    def update_e(self):
        """
        calculate the dielectric matrix for all layers
        """
        self.e = construct_epsilon(self.material(self.wavelength), self._pitch, 
                                   self._layer_thickness, self._total_thickness)
        self.N = self.e.shape[0]
    
    def set_incidence(self, direction, wavelength):
        """
        set propagation constants a and b based on incident light direction
        """
        self.wavelength = wavelength
        self.omega = sc.c / wavelength * 2 * np.pi
        k = direction / np.sqrt(np.dot(direction, direction)) * 2 * np.pi / self.wavelength
        self.a = k[0]
        self.b = k[1]
        self._p_incident_p_polarised = normalise(np.cross([self.b, -self.a, 0], k))
        self._p_incident_s_polarised = normalise(np.array([self.b, -self.a, 0]))
        
    def update_D(self):
        self.k = [calc_k(e, self.a, self.b, self.omega) for e in self.e]
        self.p = [calc_p(self.e[i], self.k[i], self.omega) for i in range(self.N)]
        self.q = [calc_q(self.k[i], self.p[i], self.omega) for i in range(self.N)]
        self.D = np.asarray([calc_D(self.p[i], self.q[i]) for i in range(self.N)])
   
    def update_P(self):
        self.k = np.asarray(self.k)
        self.P = [np.diag(np.exp(1j* self._layer_thickness * self.k[i,:,2])) for i in range(self.N)]
    ###Writhe the transfer matrix constructor
    def calc_T_medium(self):
        T_medium = np.diag([1,1,1,1])
        for i in range(self.N):
            T_medium = T_medium.dot(self.D[i].dot(self.P[i].dot(np.linalg.inv(self.D[i]))))
        self.T_medium = T_medium
        
if __name__ == '__main__':
    # self-testing codes
    a = [[200e-9,300e-9,500e-9], [1,1.2,1.5]]
    b = [[200e-9,300e-9,500e-9], [1,1.3,1.5]]
    c = [[200e-9,300e-9,600e-9], [1,1.5,1.6]]
    m = U_Material(a,b)
    l = Layers(m, 109, 30, 1000)
    l.set_incidence([0,0,1], 450e-9)
    l.update_e()
    l.update_D()
    l.update_P()
    l.calc_T_medium()