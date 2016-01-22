# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 23:03:47 2016

@author: Bonan
"""


import numpy as np


def normalise(vector):
    """Return a normalised vector"""
    return vector / np.linalg.norm(vector)
    
def rotZ(theta):
    """
    Return rotation matrix that rotates with repect to z axis with theta degress
    """
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],[np.sin(theta), np.cos(theta), 0]
                    ,[0, 0, 1]])
    return rot
    
def stack_dot(array):
    """
    Take the matrix product along 1st axis
    """
    product = np.identity(len(array[0]))
    for i in array:
        product = product.dot(i)
    return product
    
def buildDeltaMatrix(eps, Kx):
    """Returns Delta matrix for given relative permitivity and reduced wave number.
    
    'Kx' : reduced wave number, Kx = kx/k0
    'eps' : relative permitivity tensor
    
    Returns : Delta 4x4 matrix, generator of infinitesimal translations
    """
    return np.array(
        [[-Kx * eps[2,0] / eps[2,2], -Kx * eps[2,1] / eps[2,2], 
          0, 1 - Kx**2 / eps[2,2]],
         [0, 0, -1, 0],
         [eps[1,2] * eps[2,0] / eps[2,2] - eps[1,0],
          Kx**2 - eps[1,1] + eps[1,2] * eps[2,1] / eps[2,2],
          0, Kx * eps[1,2] / eps[2,2]],
         [eps[0,0] - eps[0,2] * eps[2,0] / eps[2,2],
          eps[0,1] - eps[0,2] * eps[2,1] / eps[2,2],
          0, -Kx * eps[0,2] / eps[2,2]]])