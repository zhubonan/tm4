# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 12:09:32 2015
Some useful math functions
@author: Bonan
"""

import numpy as np
import math as m


def rot_z(theta):
    """
    Return rotation matrix that rotates with repect to z axis with theta degress
    """
    rot = np.array([[m.cos(theta), -m.sin(theta), 0],[m.sin(theta), m.cos(theta), 0]
                    ,[0, 0, 1]])
    return rot


def rot_angles(pitch, layer_t, thickness):
    """
    get a array containing the anlges base on the pitch and thickness given
    """
    # get the number of layers 
    n_l = m.modf(thickness/layer_t)
    return np.linspace(0,2*np.pi*thickness/pitch,n_l[1])
    
def construct_epsilon(epsilon, pitch, layer_t, thickness):
    angles = rot_angles(pitch, layer_t, thickness)
    return np.array([rot_z(i).dot(epsilon.dot(rot_z(-i))) for i in angles])
    