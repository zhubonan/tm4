# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 22:51:47 2016
New class file
@author: Bonan
"""
import scipy as sp
import numpy as np
import tm4.mathFunc as mfc
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d
from multiprocessing import Pool
from functools import partial
from os import cpu_count
#%%
class customCall:
    """A convenient class to overlording calls to a object. It is used in Material
    Class in order to allow multiprocessing compatibility"""
    def __init__(self,value):
        self.v= value
        
    def __call__(self,*args):
        return self.v
        
class OpticalProperties:
    # Transformation matrix from the (s,p) basis to the (L,R) basis...
    C = 1 / np.sqrt(2) * np.array([[1, 1], [1j, -1j]])
    D = 1 / np.sqrt(2) * np.array([[1, 1], [-1j, 1j]])  #Transition from circular to plane for reflected wave
    TM = np.array([C,C])
    invC = np.array(sp.linalg.inv(C))
    invD = np.array(sp.linalg.inv(D))
    invTM = np.array([invD, invC])    
    def __init__(self, TOverall):
        """
        A conventnt class to store/access optical properties.
        Initialise the class by pass the overall transfer matrix(including incident
        and reflected light).
        """
        self.J = self.getJones(TOverall)
        self.Jc = self.circularJones()
        self.RP = self.J[0].conj()*self.J[0]
        self.RC = self.Jc[0].conj()*self.Jc[0]
        
    def getJones(self, T):
        """Returns the Jones matrices with linear polarisation basis
        J_ri is the Jones matrix in reflexion : [[r_pp, r_ps],
                                                 [r_sp, r_ss]]

        J_ti is the Jones matrix in transmission : [[t_pp, t_ps],
                                                    [t_sp, t_ss]]
        basis: [p, s]
        """
        # Extract element from the transfer matrix to form the Jones matrix
        J_it = T[2::-2,2::-2]
        J_ti = sp.linalg.inv(J_it)
        J_rt = T[3::-2,2::-2]
        
        # Then we have J_ri = J_rt * J_ti
        J_ri = np.dot(J_rt, J_ti)
        return (J_ri, J_ti)        

    def circularJones(self):
        """Returns the Jones matrices for circular polarization basis
        
        The Jones matrices for circular polarization are Jr^c = D⁻¹ Jr C  and
        Jt^c = C⁻¹ Jt C.
        Jones matrice for circular polarisation is in the form of:
                                         [[r_LL, r_LR],
                                          [r_LL, r_RR]]
        Returns : array of the same shape.
        """
        J = self.J
        Jc_ri = np.linalg.solve(self.D, J[0].dot(self.C))
        Jc_ti = np.linalg.solve(self.C, J[1].dot(self.C))
        return (Jc_ri, Jc_ti)     
    
    def applyAnalyser(self, angle, i = 0):
        """Return the Intensity after applying analyser. Assuming the incident 
        light is unpolarised.
        
        i: set to 0 for reflection, 1 for transimission. default is 0
        """
        Jp = mfc.polariserJ(angle)
        v = Jp.dot(self.J[i].dot([1,1]))
        return np.linalg.norm(v)**2

#%%
# Propagator for a homogeneous slab of material...
# need to stop using deprecated scipy functions
class Propagator():
    """
     A propagator class for easy access 
    """
    
    def __init__(self, method = 'eig', inv = True):
        """the value of 'method':
        "linear" -> first order approximation of exp()
        "Padé"   -> Padé approximation of exp()
        "Taylor" -> Taylor development of exp()
        "eig"    -> calculation with eigenvalue decomposition
        The default method is Pade. It is the default method in scipy libray for 
        matrix exponential
        """
        self.method = method
        if inv == True:
            self._i = -1
        else:
            self._i = 1
        self.q = 1
    def setMethod(self, method):
        self.method = method
    def setQ(self, q):
        self.q = q
    def __call__(self, Delta, h, k0,q):
        """
        'Delta' : Delta matrix of the homogeneous material
        'h' : thickness of the homogeneous slab
        'k0' : wave vector in vacuum, k0 = ω/c
        'q' : order of the approximation method, if useful
        Returns a ndarry of the propagator for the division
        """
        #Ignore the passed value 
        q = self.q
        if   self.method == "linear":    return np.identity(4) + 1j * h * k0 * Delta * self._i
        elif self.method == "Pade":      return sp.linalg.expm(1j * h * k0 * Delta * self._i, q)
        elif self.method == "Taylor":    return sp.linalg.expm3(1j * h * k0 * Delta * self._i, q + 1)
        elif self.method == "eig":       return sp.linalg.expm2(1j * h * k0 * Delta* self._i)
            
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
            self.fa = customCall(a)
            self.fb = customCall(b)
            self.fc = customCall(c)
            return
            
        self.a , self.b, self.c = np.asarray(a), np.asarray(b), np.asarray(c)
        self.kind = kind
        self._build_interp()
    
    def _build_interp(self):
        # assigne function for interpolated
        self.fa = interp1d(self.a[0], self.a[1], self.kind)
        self.fb = interp1d(self.b[0], self.b[1], self.kind)
        self.fc = interp1d(self.c[0], self.c[1], self.kind)
        
    def getTensor(self, wl):
        """
        Return the calcuated dielectric tensor at given wave length.
        Optical axis is along the x direction by default
        """
        # construct the dielectric constant tensor
        e = np.diag([self.fa(wl), self.fb(wl), self.fc(wl)]
        )**2
        return e

        
class UniaxialMaterial(Material):
    """
    An sub-class representing an uniaxial material(a!=b=c)
    """
    def __init__(self, e, o, kind = 'quadratic'):
        Material.__init__(self, e,o,o,kind)

class HomogeneousNondispersiveMaterial(Material):
    
    def __init__(self,n):
        self.n = n
    def getTensor(self, wavelength):
        return np.diag([self.n,self.n,self.n])**2
        
    def getRefractiveIndex(self, wavelength):
        return self.n
        
class HomogeneousDispersiveMaterial(Material):
    
    def __init__(self,n_array):
        Material.__init__(self, n_array,n_array,n_array, kind= 'quadratic')
    
    def getRefractiveIndex(self, wavelength):
        return self.fa(wavelength)
        
        
#%%HalfSpace class used in the front and back layer

class HalfSpace:
    """Homogeneous half-space with arbitrary permittivity.
    
    A HalfSpace must provide this method:
    * getTransitionMatrix(k0, Kx) : return transition matrix
    
    """

    material = None     # Material object

    def __init__(self, material=None):
        """Create a homogeneous half-space of the given material."""
        self.setMaterial(material)

    def setMaterial(self, material):
        """Defines the material for this half-space"""
        self.material = material

    def getTransitionMatrix(self, Kx, wl):
        """Returns transition matrix L.
        
        'Kx' : reduced wavenumber in the x direction, Kx = kx/k0
        'k0' : wavenumber in vacuum, k0 = ω/c

        Sort eigenvectors of the Delta matrix according to propagation
        direction first, then according to $y$ component. 

        Returns eigenvectors ordered like (s+,s-,p+,p-)
        """
        epsilon = self.material.getTensor(wl)
        Delta = mfc.buildDeltaMatrix(epsilon, Kx)
        q, Psi = sp.linalg.eig(Delta)

        # Sort according to z propagation direction, highest Re(q) first
        i = np.argsort(-np.real(q))
        q, Psi = q[i], Psi[:,i]     # Result should be (+,+,-,-)
        # For each direction, sort according to Ey component, highest Ey first
        i1 = np.argsort(-np.abs(Psi[1,:2]))
        i2 = 2 + np.argsort(-np.abs(Psi[1,2:]))
        i = np.hstack((i1,i2))   # Result should be (s+,p+,s-,p-)
        # Reorder
        i[[1,2]] = i[[2,1]]
        q, Psi = q[i], Psi[:,i]     # Result should be(s+,s-,p+,p-)

        # Adjust Ey in ℝ⁺ for 's', and Ex in ℝ⁺ for 'p'
        E = np.hstack((Psi[1,:2], Psi[0,2:]))
        nE = np.abs(E)
        c = np.ones_like(E)
        i = (nE != 0.0)
        c[i] = E[i]/nE[i]
        Psi = Psi * c
        # Normalize so that Ey = c1 + c2, analog to Ey = Eis + Ers
        # For an isotropic half-space, this should return the same matrix 
        # as IsotropicHalfSpace
        c = Psi[1,0] + Psi[1,1]
        if abs(c) == 0:
            c = 1.
        Psi = 2 * Psi / c
        return Psi

class IsotropicHalfSpace(HalfSpace):
    """Homogeneous Isotropic HalfSpace.
 
    * Provides transition matrix L and the inverse.
    
      Can be equally used for the front half-space (Theta = Thetai) or for the back 
      half-space (Theta = Thetat).

    * Provides relations between angle Theta and reduced wave vector Kx.

      'Theta' is the angle of the plane wave traveling to the right (angle measured
      with respect to z axis and oriented by y). '-Theta' is the angle of the wave 
      traveling to the left.
    """
    
    def getKxFromTheta(self, Theta, wl):
        """Returns the value of Kx.
        
        'Phi' : angle of the wave traveling to the right (radians)
        'k0' : wavenumber in vacuum

        kx = n k0 sin(Theta) : Real and constant throughout the structure. 
                           If n ∈ ℂ, then Theta ∈ ℂ
        Kx = kx/k0 = n sin(Theta) : Reduced wavenumber.
        """
        n = self.material.getRefractiveIndex(wl)
        Kx = n * np.sin(Theta)
        return Kx

    def getKzFromKx(self, Kx, wl):
        """Returns the value of Kz in the half-space, function of Kx
        
        'Kx' : Reduced wavenumber,      Kx = kx/k0 = n sin(Theta)
        'k0' : wavenumber in vacuum,    kx = n k0 sin(Theta)

        Returns : reduced wave number Kz = kz/k0
        """
        # Not vectorized. Could be? 
        # Test type(Kz2)
        n = self.material.getRefractiveIndex(wl)
        Kz2 = n**2 - Kx**2
        return np.sqrt(complex(Kz2))

    def getThetaFromKx(self, Kx, wl):
        """Returns the value of angle Phi according to the value of Kx.
        
        'Kx' : Reduced wavenumber,      Kx = kx/k0 = n sin(Theta)
        'wkl' : wavelength,    kx = n k0 sin(Theta)

        Returns : angle Theta in radians.
        """
        # May be vectorized when I have time?
        n = self.material.getRefractiveIndex(wl)
        sin_Phi = Kx/n
        if abs(sin_Phi) > 1:
            sin_Phi = complex(sin_Phi)
        Phi = np.arcsin(sin_Phi)
        return Phi

    def getTransitionMatrix(self, Kx, wl, inv=False):
        """Returns transition matrix L.
        
        'Kx' : Reduced wavenumber
        'k0' : wavenumber in vacuum
        'inv' : if True, returns inverse transition matrix L^-1

        Returns : transition matrix L
        """
        n = self.material.getRefractiveIndex(wl)
        
        sin_Theta = Kx/n
        if abs(sin_Theta) > 1:
            sin_Theta = complex(sin_Theta)
        cos_Theta = np.sqrt(1 - sin_Theta**2)
        if inv:
            return 0.5 * np.array( 
                    [[ 0        , 1, -1/(n*cos_Theta),  0   ],
                     [ 0        , 1,  1/(n*cos_Theta),  0   ],
                     [ 1/cos_Theta, 0,  0            ,  1/n ],
                     [ 1/cos_Theta, 0,  0            , -1/n ]])
        else:
            return np.array( 
                    [[0         , 0        , cos_Theta, cos_Theta],
                     [1         , 1        , 0      , 0      ],
                     [-n*cos_Theta, n*cos_Theta, 0      , 0      ],
                     [0         , 0        , n      , -n     ]])
    
#%% Defining Structure class

class Structure():
    """
    A superclass represent a type of structure
    """
    propagtor = Propagator(method = 'eig')

    def __init__(self):
        self.optParas = {}
        self.phyParas = {}
        
    def constructEpsilon():
        """A method to construct the relative dielectric tensor"""
        raise NotImplementedError
        
    def constructDelta():
        """A method to build the delta matrix"""
        raise NotImplementedError
        
    def constructPropogator():
        """A method a build the propagator matrices"""
        raise NotImplementedError
        
    def getPartialTransfer():
        """A method to build the partial transfer matrix for the layer"""
        raise NotImplementedError
#        self.propagator = 
    def getInfo():
        """A method to get info of the structure"""
        raise NotImplementedError
        
    def setWl(self, wl):
        self.optParas['wl'] = wl

    def setThickness(self, t):
        self.phyParas['t'] = t
    
    def setPhi(self,Phi):
        self.optParas['Phi'] = Phi
        
    def setOptParas(self, wl, Kx, Phi = 0):
        """Set up the optical parameters
        
        wl: wavelength of incident light
        
        Theta: incident angle
        
        Phi: azimuthal angle
        """
        self.optParas ={'wl': wl, 'Kx': Kx, 'Phi': Phi}
        
class HomogeneousStructure(Structure):
    """
    A class for a homogeneous layer. Can be used as basic calculation element for
    complex structures.
    """
    sType = 'homo'
    
    def __init__(self, material, t, aor = 0):
        super(HomogeneousStructure,self).__init__()
        self.setPhyParas(material, t, aor)
        self.info = {"Type":"Homogeneous", "TotalThickness": self.phyParas['t']}
        
    def setPhyParas(self, material, t, aor = 0):
        
        self.phyParas.update({'m':material, 't':t, 'aor': aor})
        
    def constructDelta(self):
        o = self.optParas
        wl, Phi, Kx= o['wl'], o['Phi'], o['Kx']
        e =self.phyParas['m'].getTensor(wl)
        e = mfc.rotedEpsilon(e, self.phyParas['aor']-Phi)
        self.delta = mfc.buildDeltaMatrix(e, Kx)
    
    def getPartialTransfer(self):
        self.constructDelta()
        self.partialTransfer = self.propagtor(self.delta, self.phyParas['t'],2*np.pi/self.optParas['wl'], q = None)
        return self.partialTransfer
        
    def getInfo(self):
        """Get infomation of the structure"""
        return {"Type":"Homegeneous", "TotalThickness": self.t}
        
#%%        


class Helix(Structure):
    """An class represent a CNC helix that is sliced into a stack of homogeneous
    layers
    """
    sType = 'helix'

    def __init__(self, *argv, **phyParas):
        """
        input arguments:
            
            material: material of the structure
            
            pitch: pitch of the helix
            
            d: number of division of the helix. Default is 30
            
            aor: intrinsic angle of rotation of the helix
        """
        super(Helix, self).__init__()
        self.phyParas.update(phyParas)

    def setPhyParas(self, material, pitch, t, d, handness, aor = 0):
        """Set the physical parameters of the class"""
        self.phyParas = {'Description':'Helix', 'd':d , 
        'sliceThickness': t / d, 'p': pitch, 't': t,
        'm': material, 'aor': aor, 'handness': handness}
        
    def getInfo(self):
        """Get infomation of the structure"""
        return self.phyParas
        
    def setPitch(self, pitch):
        """set the pitch of the helix"""
        self.phyParas.update({'p':pitch})
        
    def setThickness(self, t):
        """set the thicknes of the helix"""
        self.phyParas['t'] = t
        
        
    def _checkFastCondition(self):
        """
        Check the assumption is valid.
        Requirement:
    
        Material is uniaxial, normal incidence        
        """
        diag = np.diag(self.phyParas['m'].getTensor(self.optParas['wl']))
        if not (diag[0] != diag[1] and diag[1] == diag[2]):
            return False
        if self.optParas['Kx'] != 0:
            return False
        return True       
            
###### CORE algorithium : calculating the partial transfer matrix ###### 
    def getAngles(self):
        """Get the angle of roatation for each layer. Return an 1d array of the angles
        to be rotated for each layer. These are the physics angles to be rotated"""
        p = self.phyParas
        endAngle = p['t'] / p['p'] * p['handness'] * np.pi
        return np.linspace(0, endAngle, p['d'], endpoint = False) + p['aor']   
    
    def getSliceThickness(self):
        """
        Return the slice thickness based on current setting of total thickness
        and number of divisions
        """
        return self.phyParas['t'] / self.phyParas['d']
        
    def getPartialTransfer(self, q = None):
        """
        Build the partial transfer matrix, need to input wl and q
        """
        # If the thickness is zero, return the identity as the partial transfer matrix
        p = self.phyParas
        o = self.optParas
        if p['t'] == 0:
            return np.identity(4)
        # Check we can use the Fast Route
        
        if self._checkFastCondition == True:
            raise RuntimeWarning('Can use the HelixFast for faster calculation')
        # Normal Calculation routine
        # First we spawn a list of HomogeneousStructures
        sliceT = self.getSliceThickness()
        layerList = [HomogeneousStructure(p['m'],sliceT) for i in range(p['d'])]
        # Then we set the .Phi of each layer.
        PhiList = -(self.getAngles() - o['Phi']) # Note the reqired .Phi is the opposite of the angles of the phyiscs rotation
        PMatrices = []
        for layer, phi in zip(layerList, PhiList):
            layer.setOptParas(o['wl'], o['Kx'], phi)
            PMatrices.append(layer.getPartialTransfer())

        self.P = PMatrices
        #Take dot product for all 4x4 slices the first axis 
        return  mfc.stackDot(PMatrices)
            
        
class HelixCustom(Helix):
    """
    A class that allow custom orientation of each layer in a chiral medium
    """
    
    def __init__(self, **phyParas):
        super(HelixCustom, self).__init__(**phyParas)
        self.angleList = []
    
    def getAngles(self):
        """
        Get a list of angles. This function should also update the thickness
        """
        if self.angleList == []:
            raise RuntimeError('angles not calculated')
        return self.angleList
        
    def calcStandardAngles(self):
        """
        Calculated the standard angle list and store it in self.angleList
        """
        self.angleList = Helix.getAngles(self)
        
        
class HelixFast(Helix):
    """
    A class that uses alternative Yeh calculation routine when the incident is 
    normal and each layer is a uniaxial slab
    """
    def __init__(self, **phyParas):
        super(HelixFast, self).__init__(**phyParas)
        
    #Calculation routine
    def getPQVectors(self):
        """
        Return a list of polarisation vectors
        
        return: (p, q) where p and q are 3xN array of 3x1 vectors
        """
        angles = self.getAngles()
        p = mfc.vectorFromTheta(angles)
        q = mfc.rotZ(np.pi/2).dot(p)
        return p , q
        
    def getPandD(self):
        """
        Calculate the P and D matrices for each layers
        """
        wl = self.optParas['wl']
        d = self.phyParas['d']
        e = self.phyParas['m'].getTensor(wl)
        e_o = e[1,1] # Ordinary relative dielectric constant
        e_e = e[0,0] # Extra-ordinary
        # Wave vector for two modes note k_z = 1 we assumed during scaling
        k_o, k_e = np.sqrt(e_o), np.sqrt(e_e) 
        p, q = self.getPQVectors()
        # Initialise storage space for P and D matrix
        P, D = np.zeros((4,4,d), dtype = np.complex), np.zeros((4,4,d), dtype = np.float)
        # We take the order of k, p ,q pair to be: k_e, -k_e, k_o, -k_o
        # Note there is a angle of rotion of pi needed since p = k cross p
        r = mfc.rotZ(np.pi/2)
        p1, p2, p3, p4 = p, p, q, q  # The sign of polarisation vectors are arbitary
        # For q vectors need to put minus sign due to negative kz
        q1, q2, q3, q4 = q * k_e, -q * k_e, -p * k_o, p * k_o # q -> -p is rotation of pi/2
        # Assign values to the D matrices        
        D[:3:2,: ,:] = np.array([p1[:2], p2[:2], p3[:2], p4[:2]]).swapaxes(0,1)
        D[-1:0:-2] =  np.array([q1[:2], q2[:2], q3[:2], q4[:2]]).swapaxes(0,1)
        # Assign values to the P matrices
        s = self.getSliceThickness() # get the slice thickness
        k_0 = 2 * np.pi / wl # The magnitude of k vector in vaccum
        P = np.array([[np.exp( -1j * s * k_0 * k_e), 0 ,0 ,0],
              [0, np.exp( 1j * s * k_0 * k_e), 0, 0],
              [0, 0, np.exp( -1j * s * k_0 * k_o), 0],
              [0, 0, 0, np.exp( 1j * s * k_0 * k_o)]])
        # We now have P and D ready
        return P, D.transpose((2,0,1))
        
    def getPartialTransfer(self, q = None):
        """
        Get the partial transfer matrix with basis Ex, Ey, Hx, Hy
        """
        if self._checkFastCondition() == False:
            raise RuntimeError('Condition for fast calcution is not satisfied')
        P, D = self.getPandD()
        DInv = np.linalg.inv(D)
        # Calcuate the transition to basis for partial waves with k_e, -k_e, k_o, -k_o
        D0, DEnd = D[0], D[-1]
        # Rows: D with px, qy, py, qx, Tr need : px, py, qx, qy
        Tr0 = np.array([D0[0], D0[2], D0[3], D0[1]])
        TrEnd = np.array([DEnd[0], DEnd[2], DEnd[3], DEnd[1]])
        # Get the parital transfer matrix
        # Now begin the loop to calculate the partial transfer matrix
        # T for each layer is D[n] @ P[n] @ inv(D[n])
        # @ is the matrix multiplication but for backward compatibility still use
        # .dot sytex        
        n = self.phyParas['d']
        for i in range(n):
            if i == 0:
                # Here
                T = P.dot(DInv[0])
                continue
            if i == n-1:
                # Here
                T = T.dot(D[i]).dot(P)
                continue
            T = T.dot(D[i]).dot(P).dot(DInv[i])
        self.T_eff = T #For debug
        return Tr0.dot(T).dot(np.linalg.inv(TrEnd))
        # Change basis to be compatible with the rest of the code
#%%        
class HeliCoidalStructure(Helix):
    """
    A class that speed up the calculation by dividing itself into a repeating 
    unit and a remainder unit. The transfer matrix of the repeating unit can be 
    raised to power to get the effective transfer matrix.
    """
    sType = 'helix'
    #wl = None #dont specify wavelength initially

    def __init__(self, material, pitch, t, d= 30, handness = 'left', aor = 0):
        """
        Constructor for HeliCoidalStucture class
        
        t: Total thickness of the structure
        
        d: divisions per pitch/remainder unit
        """
        # Set handness of the helix
        if handness == 'left':
           h = -1
        elif handness == 'right':
            h = 1
        else: raise RuntimeError('Handness need to be either left or right')
        Helix.__init__(self)
        self.setPhyParas(material, pitch, t, d, h, aor = 0)
        self.phyParas['sliceThickness'] = pitch /d 
        
    def getPartialTransfer(self, q = None, updateEpsilon = True):
        """
        Build the partial transfer matrix, need to input wl and q
        """
        p = self.phyParas
        o = self.optParas
        r = np.remainder(p['t'], p['p'])
        if self._checkFastCondition() == True:
            unit, remainder = HelixFast(), HelixFast()
        else:
            unit, remainder = Helix(), Helix()
        # Need to use a copy for the sub helixs
        unit.setPhyParas(p['m'], p['p'], p['p'], p['d'], p['handness'], p['aor'])
        remainder.setPhyParas(p['m'],p['p'], r, p['d'], p['handness'], p['aor'])
        # Copy properties
        unit.optParas, remainder.optParas = o, o
        self.unit, self.remainder = unit, remainder
        unitT = unit.getPartialTransfer(None)
        remainderT = remainder.getPartialTransfer(None)
        n = int(p['t']/p['p'])
        return np.linalg.matrix_power(unitT,n).dot(remainderT)
        
#%% OptSystem Class definition

class OptSystem():
    """
    A class that represent a full set up. Contains multiple structures and front/
    back half-space. Optical properties should be extracted from a full system
    """
    # Attributes
    # wl: wavelength, Theta: Incident angle, Phi:Azimuthal angle
    Phi = 0
    Theta = 0
    wl = None
    def setStructure(self, strucList):
        """Set the Stucture of the System"""
        self.structures = strucList
        
    def setFrontHalfSpcae(self, frontHalf):
        
        self.frontHalfSpace = frontHalf
        
    def setBackHalfSpace(self, backHalf):
        
        self.backHalfSpace = backHalf
        
    def setHalfSpaces(self, front, back):
        
        self.setFrontHalfSpcae(front)
        self.setBackHalfSpace(back)
        
    def setIncidence(self, wl, Theta =0, Phi = 0):
        """Set the incidence conditions
        
        wl: wavelength of the incident light
        
        Theta: incident angle of the incidentl light
        
        Phi: azimuthal angle of the incident light"""
        self.wl = wl
        # Calculate Kx from the property of front half-space, Kx is then conserved
        # throughout the whole calculation
        self.Kx = self.frontHalfSpace.getKxFromTheta(Theta, self.wl)
        self.Theta = Theta
        self.Phi = Phi
        
    def updateStructurePartialTransfer(self):
        """Calculate the partial transfer matrices for all structures and update
        the overall partial transfer matrix
        This is reqired before calculating overall transfer matrix
        """
        overallPartial = np.identity(4)
        for s in self.structures:
            s.setOptParas(self.wl, self.Kx, self.Phi)
            # Apply matrix products
            overallPartial = overallPartial.dot(s.getPartialTransfer())
        self.overallPartial = overallPartial  
        
    def getTransferMatrix(self):
        """Calculate the overall transfer matrix of the optical system and update
        the prop property of the optSystem instance
        
        Return: overall transfer matrix (4 by 4)"""
        frontinv = self.frontHalfSpace.getTransitionMatrix(self.Kx,self.wl, inv = True)
        back = self.backHalfSpace.getTransitionMatrix(self.Kx, self.wl)
        self.transferMatrix = frontinv.dot(self.overallPartial.dot(back))
        self.prop = OpticalProperties(self.transferMatrix)
        return self.transferMatrix
    
    def getSubStructureInfo(self):
        """Print out the information about structures in the optsystem instance"""
        index = 0
        for s in self.structures:
            index += 1
            print("Layer " + str(index), s.getInfo())
            
    def setPitch(self, pitchList):
        """Change pitch of the structure in one go"""
        if type(pitchList) == int:
            pitchList = [pitchList]
            
        for i in range(len(pitchList)):
            self.structures[i].setPitch(pitchList[i])
            
    def setThickness(self, tList):
        """Set the thickness of all structures
        tList: a list containing the thicknesses of the structures in structures
        """
        if type(tList) == int:
            tList = [tList]
            
        for i in range(len(tList)):
            self.structures[i].setThickness(tList[i])
    
    def calcWl(self, wl, coupling = 'LL'):
            """A function to be called for calculating at a certain wl"""
            self.setIncidence(wl, self.Theta, self.Phi)
            self.updateStructurePartialTransfer()
            self.getTransferMatrix()
            if coupling == 'LL':
            # take real part only. Imag has some residual 
                return self.prop.RC[0,0].real
            elif coupling == 'RR':
                return self.prop.RC[1,1].real
            elif coupling == 'SS':
                return self.prop.RP[1,1].real
            elif coupling == 'PP':
                return self.prop.RP[0,0].real
            elif coupling == 'full':
                return self.prop
            else:
                raise RuntimeError('Output not specified')
    
    def scanSpectrum(self, wlList,  coreNum = 1, giveInfo = False, coupling = 'LL'):
        """Cacluate the respecon at the given wavelengths. 
       
        giveinfo: boolen, determine if return information about the strcuture
        
        useProp: if set to True the result will be a list of OptProperties. Userful
        if want to use the full information from calculation. e.g. Check the effect
        of using analyser
        """
        result = []
        calcWl = partial(self.calcWl, coupling = coupling)
        # Initialise multiprocessing
        if coreNum == 'auto':
            coreNum = cpu_count() - 1
        # Want to terminate the processes if anything goes wrong        
        if coreNum == 1:
            result = list(map(calcWl, wlList))
        else:
            with Pool(processes=coreNum) as pool:
                result = pool.map(calcWl, wlList)
        if giveInfo:
            return wlList, result, self.getSubStructureInfo()
        else: return wlList,result
      
#%%Some staff for convenience
air = HomogeneousNondispersiveMaterial(1)
airHalfSpace = IsotropicHalfSpace(air)
glass = HomogeneousNondispersiveMaterial(1.55)
glassHalfSpace = IsotropicHalfSpace(glass)
CNC = UniaxialMaterial(1.586,1.524)

if __name__ == "__main__":
#%%
    # For chekcing if two methods are consistent
    """
    from preset import *
    from time import clock    
    repeat = 5
    wlRange = np.linspace(400,800,200)
    cStart = clock()
    h1.phyParas['t'] = 150
    for i in range(repeat):
        res = s.scanSpectrum(wlRange,2)
        #pl.plot(res[0],res[1])
    timeUsed = (clock() - cStart) / repeat
    print(timeUsed)
    """
    #Test for angled incidence
    """
    from preset import *
    def run():
        s.setIncidence(500, 0,0)
        res1 = s.scanSpectrum(wlRange, 1)
        s.setIncidence(500,np.pi/3,0)
        res2 = s.scanSpectrum(wlRange, 1)
        pl.plot(res1[0],res1[1],label = 'Normal Incidence')
        pl.plot(res2[0],res2[1],label = '45 degree')
        pl.legend()
        return
    """
    # Testing the fast calcution routine
    testMaterial = UniaxialMaterial(1.6,1.5)        
    h = HelixFast(m = testMaterial, t = 1800, p = 180, d = 300, aor = 0, handness = -1)
    h.setOptParas(500, 0)
    P, D = h.getPandD()
    res = h.getPartialTransfer()
    T = h.T_eff
    from preset import *
    h1.setPhyParas(testMaterial, 180, 1800, 30, aor = 0, handness = -1)
    pl.plot(*s.scanSpectrum(wlRange,1))    