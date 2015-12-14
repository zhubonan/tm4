# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 12:09:32 2015
Some useful math functions
@author: Bonan
"""

import numpy as np
import math as m

def construct_epsilon(epsilon_diag, pitch, layer_t, thickness):
    """
    construct the dielectric matrices of all layers
    return a N*3*3 array where N is the number of layers
    """
    
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
        
        
    angles = rot_angles(pitch, layer_t, thickness)
    return np.array([rot_z(i).dot(epsilon_diag.dot(rot_z(-i))) for i in angles])
    
def calc_c(e, a, b , omega, u = 1.25e-6): # Check units
    """
    calculate the z components of 4 partial waves in medium
    
    e: dielectric tensor
    
    a,b: components of wavevector in direction of x and y direction
    
    omega: fequency
    
    return a list containting 4 roots for the z components of the partial waves
    """
    # assign names 
    x = e* u * omega**2
    x11, x12, x13 = x[0]
    x21, x22, x23 = x[1] 
    x31, x32, x33 = x[2]
    # calculate the coeffciency based on symbolic expression
    coef4 = b**2
    coef3 = b*x23 + b*x32
    coef2 = a**2*b**2 + a**2*x33 + a*b*x12 + a*b*x21 - b**2*x11 + b**2*x22 + \
    b**2*x33 - x22*x33 + x23*x32
    coef1 = a**3*x13 + a**3*x31 + a**2*b*x23 + a**2*b*x32 + 2*a*b**2*x13 + \
    2*a*b**2*x31 + a*x12*x23 - a*x13*x22 + a*x21*x32 - a*x22*x31 + b**3*x23 + \
    b**3*x32 - b*x11*x23 - b*x11*x32 + b*x12*x31 + b*x13*x21
    coef0 = a**4*x11 + a**3*b*x12 + a**3*b*x21 - a**2*b**4 + 2*a**2*b**2*x11 + \
    a**2*b**2*x22 - a**2*x11*x22 - a**2*x11*x33 + a**2*x12*x21 + a**2*x13*x31 + \
    a*b**3*x12 + a*b**3*x21 - a*b*x12*x33 + a*b*x13*x32 - a*b*x21*x33 + \
    a*b*x23*x31 - b**6 + b**4*x11 + b**4*x22 + b**4*x33 - b**2*x11*x22 - \
    b**2*x11*x33 + b**2*x12*x21 + b**2*x13*x31 - b**2*x22*x33 + b**2*x23*x32 \
    + x11*x22*x33 - x11*x23*x32 - x12*x21*x33 + x12*x23*x31 + x13*x21*x32 - x13*x22*x31
    # calculate the roots of the quartic equation
    c = np.roots([coef4, coef3, coef2, coef1, coef0])
    if len(c) == 2:
        return np.append(c,c)
    return c

def calc_k(e , a, b, omega, u = 1.25e-6):
    
    """
    A wrapper to calcualte k vector 
    """
    c = calc_c(e, a , b, omega, u)
    return np.array([[a, b, c[0]],[a, b , c[1]],[a,b,c[2]],[a,b,c[3]]])
    
def calc_p(e, k,  omega, u = 1.25e-6):
    """
    Calculate the polarisation vector based on the calculated wavevector and frequency
    equation(9.7-5)
    
    e: dielectric tensor
    
    k: 4x3 array of 4 k vectors
    """
    x = e* u * omega**2
    p = []
    for i in k:
        a = i[0]
        b = i[1]
        c = i[2]
        x11, x12, x13 = x[0]
        x21, x22, x23 = x[1] 
        x31, x32, x33 = x[2]
        p1 = (x22 - a**2 - c**2) * (x33 - a**2 - b**2) - (x23 + b * a )**2
        p2 = (x23 + b * c) * (x31 + a *c) - (x12 + a *b) * (x33 - a**2 - b **2)
        p3 = (x12 + a *b ) * (x23 + b *c) - (x13 + a * c) * (x22 - a**2 - c **2)
        # normalised the vector
        p_temp = np.array([p1,p2,p3])
        p.append(p_temp / np.sqrt(np.dot(p_temp, p_temp)))

    return np.array(p)
    
def calc_q(k , p, omega, u = 1.25e-6):
      
    return 3e8 / omega / u * np.cross(k ,p)

def calc_D(p,q):
    
    return np.array([p[:,0],q[:,1], p[:,1], q[:,0]])

def calc_P(k,t):
    return np.diag(np.exp(1j*t*k[:,2]))

def construct_D(e, a, b, omega, u = 1.25e-6):
    """
    construct the dynamic matrix for one layer with know dielectric tensor
    """
    k = calc_k(e, a, b, omega, u)
    p = calc_p(e , k , omega, u)
    q = calc_q(k , p, omega, u)
    return calc_D(p,q)
    
if __name__ == "__main__":
    
    e = np.diag([2,1,1]) * 8.85e-12
    a = 0
    b = 0
    omega = 500e12*6
    #print(calc_c(e,a,b,omega))
    k = calc_k(e,a,b,omega)
    p = calc_p(e, k, omega)
    q = calc_q(k , p, omega)
    D = calc_D(p,q)
    