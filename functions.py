# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:42:02 2015
4 by 4 matrix formulation for birefringent thin layers
This file contains functions for mathematical manipulations
@author: Bonan
"""
import numpy as np
import math
def normalise(vector):
    return vector / math.sqrt(np.dot(vector, vector))
    
def v_p(epsilon, k, omega, mu = 1):
    """
    Constructe the polarisation vection with given wave vectors
    
    espilon: 3 by 3 dielectric tensor
    
    omega: frequency of wave
    
    k: wave vector (aplha, beta, gamma)
    
    mu: relative permitivity, default is one for non magnetic material
    
    return a 1-d numpy array representing the polarisation vector( of electric 
    field)
    """
    
    # use simplier names
    e, w, u, a, b, g = epsilon, omega, mu, k[0], k[1], k[2]
    
    # calculate each component of p 
    p_1 = (w**2 * u * e[1,1] - a**2 - g**2) * (w**2 * u * e[2,2] - a**2 - b**2) \
          - (w**2 * u * e[1,2]+ b * g)**2
    p_2 = (w**2 * u * e[1,2] + b * g) * (w**2 * u * e[2,0] + a * g) - \
          (w**2 * u * e[0,1] + a * b) * (w**2 * u * e[2,2] - a **2 - b**2)
    p_3 = (w**2 * u * e[0,1] + a * b) * (w**2 * u * e[1,2] + b * g) - \
          (w**2 * u * e[0,2] + a * g) * (w*82 * u * e[1,1] - a**2 - g**2)
    
    # normalise polarisation vector
    p = np.array([p_1,p_2,p_3])
    p = p / np.dot(p,p)**0.5      
    return p
    
def v_q(p, k, omega, mu = 1):
    """
    Construct the magnetic field polarisation vector from electric field's
    
    p: electric field polarisation vector
    
    k: wave vector
    
    omega: frequency of wave
    
    return a 1-d numpy array representing the magnetic field polarisation vector
    """
    
    # not q is not normalised here
    return 3e8 / omega / mu * np.cross(k, p) 
