# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 17:02:54 2015
class definitions
@author: Bonan
"""
from scipy.interpolate import interp1d
import scipy.constants as sc
import numpy as np
from old.mathFuncLast import *
from abc import ABCMeta, abstractmethod

class Optical_Properties:
 
    def __init__(self, T_overall):
        """
        A conventnt class to store/access optical properties.
        Initialise the class by pass the overall transfer matrix(including incident
        and reflected light).
        """
        self.c_matrices = calc_coupling_matrices(T_overall)
        self.update_R()
        self.update_T()
        
    def update_R(self):
        """
        Update reflectivities
        """
        R_C = np.abs(self.c_matrices["Circular_r"])**2
        R_P = np.abs(self.c_matrices["Plane_r"])**2
        self.RCRR = R_C[1,1]
        self.RCRL = R_C[0,1]
        self.RCLR = R_C[1,0]
        self.RCLL = R_C[0,0]
        self.RPss = R_P[0,0]
        self.RPsp = R_P[1,0]
        self.RPps = R_P[0,1]
        self.RPpp = R_P[1,1]
        
    def update_T(self):
        """
        Update transmisivities
        """
        T_C = np.abs(self.c_matrices["Circular_r"])**2
        T_P = np.abs(self.c_matrices["Plane_r"])**2
        self.TCLL = T_C[0,0]
        self.TCLR = T_C[1,0]
        self.TCRR = T_C[1,1]
        self.TCRL = T_C[0,1]
        self.TPss = T_P[0,0]
        self.TPsp = T_P[1,0]
        self.TPps = T_P[0,1]
        self.TPpp = T_P[1,1]
        
    def calc_rot_of_pol(self, theta):
        """
        Calculate the rotation of polarisation for ellipsometery measurement
        simulation
        """
        # Disable the function for now
        return
        r_p = self.c_matrices["Plane_r"]
        inc_y = np.sin(theta/180*np.pi)
        inc_x = np.sqrt(1-y**2)
        ref_y,ref_x = r_p.dot(np.array([inc_y,inc_x]))
        return np.arctan2(ref_y, ref_x)
        
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
        
    def e_diag(self, wavelength):
        """
        Return the calcuated dielectric tensor at given wave length.
        Optical axis is along the x direction by default
        """
        # construct the dielectric constant tensor
        e = np.diag([self.fa(wavelength), self.fb(wavelength), self.fc(wavelength)]
        )**2
        return e

    def __call__(self, wavelength):
        return self.e_diag(wavelength)
        
class Uniaxial_Material(Material):
    """
    An sub-class representing an uniaxial material(a!=b=c)
    """
    def __init__(self, e, o, kind = 'quadratic'):
        Material.__init__(self, e,o,o,kind)

class Homogeneous_Nondispersive_Material(Material):
    
    def __init__(self,n):
        self.n = n
    def e_diag(self, wavelength):
        return np.diag([self.n,self.n,self.n])**2

class Homogeneous_Dispersive_Material(Material):
    
    def __init__(self,n_array):
        Material.__init__(self, n_array,n_array,n_array, kind= 'quadratic')
    
class Seg(metaclass = ABCMeta):
    """
    An abstract class represent a stack of birefringent layers
    """
    def __init__(self, material, **structure):
        self.material = material
        # store stucture_parameters
        self.structure_paras = structure
            
    def set_incidence(self, direction, wavelength):
        """
        Set propagation constants a and b based on incident light direction
        The problem is scaled such that the wavelength is one
        Universal to any structure
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
        # assign 4 k vectors for two polarisations
        k0 = np.array([k, [k[0],k[1],-k[2]],k, [k[0],k[1],-k[2]]])
        # construct the dynamical matrix of the incident media
        self.D0 = calc_D(self._incident_p, calc_q(k0, self._incident_p))
    
    @abstractmethod    
    def update_thickness(self):
        """
        Update the thickness information
        store information in self._thickness as an array. Structure specific
        """
        pass
    
    @abstractmethod
    def update_e(self):
        """
        Update the dielectric matrices. Structure specific
        """
        pass
    
    def update_D(self):
        self.k = [calc_k(e, self.a, self.b) for e in self.e]
        self.p = [calc_p(self.e[i], self.k[i]) for i in range(self.N)]
        self.q = [calc_q(self.k[i], self.p[i]) for i in range(self.N)]
        self.D = np.asarray([calc_D(self.p[i], self.q[i]) for i in range(self.N)])
        
    def update_P(self):
        self.k = np.asarray(self.k)
        self.P = [np.diag(np.exp(1j* self._thickness[i] * self.k[i,:,2] * self.scale_factor)) for i in range(self.N)]

    def update_T(self):
        """
        Calcualte effective transfer matrix for each internal interface T_eff = D[N]D[N+1]P[N+1]
        Then multiply remaining term to calculate the effective transfer matrix 
        of the segment. Incident and exiting mediums are not taken into consideration
        """
        self.T_eff = [np.linalg.solve(self.D[i-1],self.D[i].dot(self.P[i])) for i in range(1,self.N)]
        # multiply terms of D[0]P[0] and D[N-1]-1 to the product        
        self.T_eff_total = np.dot(self.D[0],self.P[0]).dot(stack_dot(self.T_eff).dot(np.linalg.inv(self.D[-1])))
    
    def doit(self):
        """
        Do as much calculation as possible
        """
        self.update_e()
        self.update_thickness()
        self.update_D()
        self.update_P()
        self.update_T()
         
         
class H_Seg(Seg):
    """
    An class object represent a stack of helicoidal arranged layers
    """
    def __init__(self, material, pitch, layer_thickness, total_thickness):
        
        self.material = material
        self.structure_paras = [pitch, layer_thickness, total_thickness]
        
    def update_e(self):
        """
        calculate the relative dielectric matrix for all layers
        """
        self.e = construct_epsilon_heli(self.material(self.wavelength), *self.structure_paras)
        # store the number of stacks
        self.N = self.e.shape[0]

    def update_thickness(self):
        """
        Build the array to represent the thickness data
        """
        self._thickness = np.full(self.e.shape[0], self.structure_paras[1])
        
    def update_T(self):
        Seg.update_T(self)
        self.T_total = np.linalg.solve(self.D0, self.T_eff_total.dot(self.D0))
        self.prop = Optical_Properties(self.T_total)
        
class H_Layers():
    """
    Multilayer structure with helicoidal arrangement.
    To speed up the calculation and increase precision, the stack is spilited into
    repeating units and remainder. The effective transfer matrix of the repeating unit is raised
    to power N where N is the number of repeats, instead of calculate transfer matrices
    for each interface.
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
        self.prop = Optical_Properties(self.T_total)

    def doit(self):
        self.calc_T()
        
if __name__ == '__main__':
    # self-testing codes
    a = [[200,300,500], [1.6,1.6,1.6]]
    b = [[200,300,500], [1.1,1.5,1.5]]
    c = [[200,300,600], [1,1.5,1.5]]
    m = Uniaxial_Material(a,b)
    l = H_Layers(m, 180, 6, 180)
    l.set_incidence([0,0,1], 500)
    l.doit()
    T1 = l._repeat.T_eff_total
    