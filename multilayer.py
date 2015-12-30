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
from numpy.linalg import inv

class Material:
    """
    An class object that represent any type of anisotropic material
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
        )
        return e

    def __call__(self, wavelength):
        return self.e_diag(wavelength)
        
class U_Material(Material):
    """
    An sub-class representing an uniaxial material(a!=b=c)
    """
    def __init__(self, e, o, kind = 'quadratic'):
        Material.__init__(self, e,o,o,kind)
    

class H_Seg():
    """
    An class object represent a stack of helicoidal arranged layers
    """
    def __init__(self, material, pitch, layer_thickness, total_thickness):
        """Initialise a multilayer system. The length scale"""
        self.material = material
        self._pitch = pitch
        self._layer_thickness = layer_thickness
        self._total_thickness = total_thickness
    def update_e(self):
        """
        calculate the relative dielectric matrix for all layers
        """
        self.e = construct_epsilon(self.material(self.wavelength), self._pitch, 
                                   self._layer_thickness, self._total_thickness)
        self.N = self.e.shape[0]
    
    def set_incidence(self, direction, wavelength):
        """
        Set propagation constants a and b based on incident light direction
        The problem is scaled such that the wavelength is one
        """
        self.wavelength = wavelength
        self.scale_factor = 2* np.pi /wavelength
        self.omega = sc.c / wavelength * 2 * np.pi
        # k is the incident wave vector here
        k = direction / np.sqrt(np.dot(direction, direction))
        self.a = k[0]
        self.b = k[1]
        # set incident light polarisation directions
        self._incident_p = incident_p(k)
        # assign 4 k vectors for four polarisations
        k0 = np.array([k, [k[0],k[1],-k[2]],k, [k[0],k[1],-k[2]]])
        # construct the dynamical matrix of the incident media
        self.D0 = calc_D(self._incident_p, calc_q(k0, self._incident_p))
        
    def update_D(self):
        self.k = [calc_k(e, self.a, self.b) for e in self.e]
        self.p = [calc_p(self.e[i], self.k[i]) for i in range(self.N)]
        self.q = [calc_q(self.k[i], self.p[i]) for i in range(self.N)]
        self.D = np.asarray([calc_D(self.p[i], self.q[i]) for i in range(self.N)])
        # add dynamic matrix of the incident/exiting medium to the end of the stack
        # no need for now as the class represent an segment only
        #self.D = np.append([self.D0], self.D, axis = 0)
        #self.D = np.append(self.D, [self.D0], axis = 0)
   
    def update_P(self):
        self.k = np.asarray(self.k)
        self.P = [np.diag(np.exp(1j* self._layer_thickness * self.k[i,:,2] * self.scale_factor)) for i in range(self.N)]

        
    def update_T(self):
        """
        Calcualte effective transfer matrix for each internal interface T_eff = D[N]D[N+1]P[N+1]
        """
        self.T_eff = [np.linalg.solve(self.D[i-1],self.D[i].dot(self.P[i])) for i in range(1,self.N)]
        # multiply terms of D[0]P[0] and D[N-1]-1 to the product        
        self.T_eff_total = np.dot(self.D[0],self.P[0]).dot(stack_dot(self.T_eff).dot(np.linalg.inv(self.D[-1])))
        
    def doit(self):
         self.update_e()
         self.update_D()
         self.update_P()
         self.update_T()

         
class H_Layers():
    """
    Multilayer structure with helicoidal arrangement
    Calculate dynamic and propagation matrix for one segment only
    """
    def __init__(self, material, pitch, layer_thickness, total_thickness):
        """
        Initialise self and sub-segments
        """
        self.material = material
        self._pitch = pitch
        # descard any incomplete layers
        total_thickness = total_thickness - total_thickness%layer_thickness
        self._layer_thickness = layer_thickness
        self._total_thickness = total_thickness
        # consider the case where there is no reminder part
        if total_thickness % pitch != 0 :
            self._flag_has_reminder = True
            self._reminder = H_Seg(material, pitch, layer_thickness, total_thickness % pitch)
        else: self._flag_has_reminder = False
        
        self._repeat = H_Seg(material, pitch, layer_thickness, pitch)
        self._N_of_units = int(total_thickness/pitch)
                  
    def set_incidence(self, direction, wavelength):
        """
        Set incidence wave parameters and pass on to sub-layers
        """
        H_Seg.set_incidence(self, direction, wavelength)
        self._repeat.set_incidence(direction, wavelength)
        if self._flag_has_reminder:
            self._reminder.set_incidence(direction, wavelength)
    
    def calc_T(self):
        #Calculate the total transfer matrix for the bulk material
        self._repeat.doit()
        if self._flag_has_reminder:
            self._reminder.doit()
        # calculate the layer part of the total transfer matrix
            T_layers = np.linalg.matrix_power(self._repeat.T_eff_total, 
                       self._N_of_units).dot(self._reminder.T_eff_total)
        else:
            T_layers = np.linalg.matrix_power(self._repeat.T_eff_total, 
                       self._N_of_units)
        # Now add dynamic matrix of the incident and exiting medium
        # Assume to be the vacuum for now
        self.T_total = np.linalg.solve(self.D0, T_layers.dot(self.D0))
        self.coeff = calc_coeff(self.T_total)
        self.coeff_modulus = self.coeff.copy()
        
    def LR_basis(self):
        """
        Caculate T_total and coefficients for LR polarised light
        """
        a = 1/math.sqrt(2)
        b = - 1j / math.sqrt(2)
        M = np.array([[a,0,a,0],[0,a,0,a],[b,0,-b,0],[0,b,0,-b]])
        N = np.array([[a,0,a,0],[0,0,0,0],[b,0,-b,0],[0,0,0,0]])
        self.T_total_LR = np.linalg.solve(M, self.T_total.dot(N))
        self.coeff_LR = calc_coeff(self.T_total_LR)
        self.coeff_modulus_LR = self.coeff_LR.copy()
        for i in self.coeff_modulus_LR:
            self.coeff_modulus_LR[i] = np.abs(self.coeff_modulus_LR[i])**2
    def doit(self):
        self.calc_T()
        self.LR_basis()
        
if __name__ == '__main__':
    # self-testing codes
    a = [[200,300,500], [1,1.2,1.5]]
    b = [[200,300,500], [1.1,1.3,1.6]]
    c = [[200,300,600], [1,1.5,1.6]]
    m = U_Material(a,b)
    l = H_Layers(m, 100, 10, 5009)
    l.set_incidence([0,0,1], 450)
    l.doit()
    