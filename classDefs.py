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
import matplotlib.pyplot as pl
#%%
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
        """Returns the Jones matrices.
        
        Returns : tuple (T_ri, T_ti)

        T_ri is the Jones matrix in reflexion : [[r_pp, r_ps],
                                                 [r_sp, r_ss]]

        T_ti is the Jones matrix in transmission : [[t_pp, t_ps],
                                                    [t_sp, t_ss]]
       
        Naming convention (Fujiwara, p. 220):
        't_ps' : transmitted p component for a s incident wave
        't_ss' : transmitted s component for a s incident wave
        ...

        Note: If all materials are isotropic, r_ps = r_sp = t_sp = t_ps = 0

        See also: 
        * extractCoefficient() to extract the desired coefficients.
        * circularJones() for circular polarization basis
        """
        # Extraction of T_it out of T. "2::-2" means integers {2,0}.
        T_it = T[2::-2,2::-2]
        T_ti = sp.linalg.inv(T_it)
        
        # Extraction of T_rt out of T. "3::-2" means integers {3,1}.
        T_rt = T[3::-2,2::-2]
        
        # Then we have T_ri = T_rt * T_ti
        T_ri = np.dot(T_rt, T_ti)
        return (T_ri, T_ti)        

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
    

#%%
# Propagator for a homogeneous slab of material...
# need to stop using deprecated scipy functions
class Propagator():
    """
     A propagator class for easy access 
    """
    
    def __init__(self, method = 'Pade', inv = True):
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
    def __call__(self, Delta, h, k0, q=None):
        """
        'Delta' : Delta matrix of the homogeneous material
        'h' : thickness of the homogeneous slab
        'k0' : wave vector in vacuum, k0 = ω/c
        'q' : order of the approximation method, if useful
        Returns a ndarry of the propagator for the division
        """
        if   self.method == "linear":    return np.identity(4) + 1j * h * k0 * Delta * self._i
        elif self.method == "Pade":      return sp.linalg.expm(1j * h * k0 * Delta * self._i, q)
        elif self.method == "Taylor":    return sp.linalg.expm3(1j * h * k0 * Delta * self._i, q+1)
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
    
    def get_Kx_from_Theta(self, Theta, wl):
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

    def get_Kz_from_Kx(self, Kx, wl):
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

    def get_Theta_from_Kx(self, Kx, wl):
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
        
    def getPartialTransfer():
        """A method to build the partial transfer matrix for the layer"""
        raise NotImplementedError
#        self.propagator = 
    def setKx(self, Kx):
        
        self.Kx = Kx
    def setWl(self, wl):
        self.wl = wl

        
class HybridStructure(Structure):
    """
    A structure representing a combination of other structures
    """
class HomogeneousStructure(Structure):

    propagtor = Propagator() #default propagator
    Phi = 0 #Set the angle phy to be zero in the start
    Kx = None
    
    def __init__(self, thickness, material = None):
        self.material = material
        self.thickness = thickness
    def constructDelta(self):
        self.delta = mfc.buildDeltaMatrix(self.material.getTensor(self.wl), self.Kx)
    
    def getPartialTransfer(self):
        self.constructDelta()
        self.partialTransfer = self.propagtor(self.delta, self.thickness, 2*np.pi/ self.wl)
        return self.partialTransfer
        
class HeliCoidalStructure(Structure):
    """A structure class represent the helicoidal structure"""
    
    _hasRemainder = False
    #wl = None #dont specify wavelength initially
    propagtor = Propagator(method = 'eig') #default propagator
    Phi = 0 #Set the angle phy to be zero in the start
    Kx = None
    
    def __init__(self, material, pitch, d, t, handness = 'left'):
        """
        Initialise the structure by passing material pitch, division per 
        pitchand total thickness. Handness is left by default
        """
        self.divisionThickness = pitch / d
        # Set handness of the helix
        if handness == 'left':
            self._handness = -1
        elif handness == 'right':
            self._handness = 1
        else: raise RuntimeError('Handness need to be either left or right')
        # Calculate the angles
        self.anglesRepeating = np.linspace(0, np.pi * self._handness, d, 
                                           endpoint = False)
        self.nOfReapting = int(t/pitch)
        remain = np.remainder(t, pitch)
        remainDivisions = int(remain/self.divisionThickness)
        #discard the overflowing parts
        remain = remainDivisions * self.divisionThickness
        if remain != 0:
            self._hasRemainder = True
            self.anglesRemaining = np.linspace(0, remain/pitch*np.pi*self._handness,
                                               remainDivisions, endpoint = False)
        self.material = material
        self.info = {"Type":"HeliCodidal", "Pitch":pitch, "DivisionPerPitch": d ,\
        "Handness":handness, "TotalThickness": remain + self.nOfReapting * pitch}

    def constructEpsilon(self,wl = None):
        """Build the epsilon for each layer"""
        self.wl = wl
        epsilon = self.material.getTensor(wl)
        self.epsilonRepeating = [mfc.rotedEpsilon(epsilon, theta) for theta in self.anglesRepeating]
        if self._hasRemainder:
            self.epsilonRemaining = [mfc.rotedEpsilon(epsilon, theta) for theta in self.anglesRemaining]
        
    def constructDelta(self):
        """Build the Delta matrix in Berreman's formulation"""
        self.deltaRepeating = [mfc.buildDeltaMatrix(e,self.Kx) for e in self.epsilonRepeating]
        if self._hasRemainder: self.deltaRemaining = [mfc.buildDeltaMatrix(e,self.Kx) for e in self.epsilonRemaining]
        
    def getPartialTransfer(self, q = None, updateEpsilon = True):
        """
        Build the partial transfer matrix, need to input wl and q
        """
        #Constructed needed matrices

        #Update the dielectric tensor stack. Can be supressed if needed e.g. at
        # a different psy angle
        if updateEpsilon:    
            self.constructEpsilon(self.wl)
        self.constructDelta()
        #Get propagation matrices, require devision thickness the wl has the same
        #unit. This effectively scales the problem
        d= self.divisionThickness
        k0 = 2*np.pi/self.wl
        PMatricesRepeating = [self.propagtor(delta, d, k0, q) for delta
                              in self.deltaRepeating]
        self.P = PMatricesRepeating
        
        if self._hasRemainder:
            PMatricesRemaining = [self.propagtor(delta, d, k0, q) for delta in 
                                  self.deltaRemaining]
        #Calculate the overall partial transfer matrix
        TRepeat = mfc.stackDot(PMatricesRepeating)
        self.PT = TRepeat
        if self._hasRemainder:
            TRemain = mfc.stackDot(PMatricesRemaining)
        else: TRemain = np.identity(4) 
        self.partialTransfer = np.linalg.matrix_power(TRepeat, self.nOfReapting).dot(TRemain)
        self.partialTransferParameters = {"wavelength":self.wl, "Kx":self.Kx, "phy": self.Phi}
        return self.partialTransfer.copy()
    

#%% OptSystem Class definition

class OptSystem():
    """
    A class that represent a full set up. Contains multiple structures and front/
    back half-space. Optical properties should be extracted from a full system
    """
    # Attributes
    # wl: wavelength, Theta: Incident angle, Phi:Azimuthal angle
    wl, Theta, Phi= None, None, None

    
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
        """Set the incidence conditions"""
        self.wl = wl
        # Calculate Kx from the property of front half-space, Kx is then conserved
        # throughout the whole calculation
        self.Kx = self.frontHalfSpace.get_Kx_from_Theta(Theta, self.wl)
        self.Theta = Theta
        self.Phi = Phi
        
    def updateStructurePartialTransfer(self):
        
        overallPartial = np.identity(4)
        for s in self.structures:
            s.setKx(self.Kx)
            s.setWl(self.wl)
            s.getPartialTransfer()
            # Apply matrix products
            overallPartial = overallPartial.dot(s.getPartialTransfer())
        self.overallPartial = overallPartial  
    def getTransferMatrix(self):
        
        frontinv = self.frontHalfSpace.getTransitionMatrix(self.Kx,self.wl, inv = True)
        back = self.backHalfSpace.getTransitionMatrix(self.Kx, self.wl)
        self.transferMatrix = frontinv.dot(self.overallPartial.dot(back))
        self.prop = OpticalProperties(self.transferMatrix)
        return self.transferMatrix
    
    def getSubStructureInfo(self):
        index = 0
        for s in self.structures:
            index += 1
            print("Layer " + str(index), s.info)
            
    def scanSpectrum(self, wlList, keyword = "L-L"):
        """Cacluate the respecon at the given wavelengths"""
        result = []
        for i in wlList:
            self.setIncidence(i,self.Theta,self.Phi)
            self.updateStructurePartialTransfer()
            self.getTransferMatrix()
            result.append(self.prop.RC[0,0])
        intel =["P:"+str(s.info["Pitch"]) + " T:" +str(s.info["TotalThickness"])
                for s in self.structures]
        return intel, wlList, result
        
#%%Some staff for convenience
air = HomogeneousNondispersiveMaterial(1)
airHalfSpace = IsotropicHalfSpace(air)

if __name__ == "__main__":
#%%    
    
    air = HomogeneousNondispersiveMaterial(1)
    glass = HomogeneousNondispersiveMaterial(1.5)
    celloluse = UniaxialMaterial(1.60,1.55)
    helix = HeliCoidalStructure(celloluse, 300, 30, 300)
    helixDouble = HeliCoidalStructure(celloluse, 100,10,200)
    slab = HomogeneousStructure(250 , glass)
    setup = OptSystem()
    front = IsotropicHalfSpace(air)
    back = IsotropicHalfSpace(air)
    setup.setHalfSpaces(front,back)
    #%%
    setup.setStructure([helix])
    result = []
    wavelength = np.linspace(200,1000,100)
    for wl in wavelength:
        setup.setIncidence(wl)
        setup.updateStructurePartialTransfer()
        setup.getTransferMatrix()
        result.append(setup.prop.RC[0,0])
    pl.plot(wavelength,result)
    