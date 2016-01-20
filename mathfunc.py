# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 12:09:32 2015
@author: Bonan
"""

import numpy as np


def normalise(vector):
    """Return a normalised vector"""
    return vector / np.linalg.norm(vector)
    
def construct_epsilon_heli(epsilon_diag, pitch, layer_t, thickness):
    """
    construct the dielectric matrices of all layers
    return a N*3*3 array where N is the number of layers
    """
        
        
    angles = rot_angles(pitch, layer_t, thickness)
    return np.array([rot_z(i).dot(epsilon_diag.dot(rot_z(-i))) for i in angles])
    
def rot_z(theta):
    """
    Return rotation matrix that rotates with repect to z axis with theta degress
    """
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],[np.sin(theta), np.cos(theta), 0]
                    ,[0, 0, 1]])
    return rot


def rot_angles(pitch, layer_t, thickness):
    """
    get a array containing the anlges base on the pitch and thickness given
    """
    # get the number of layers 
    n_l = np.modf(thickness/layer_t)
    return np.linspace(0,2*np.pi*thickness/pitch,n_l[1], endpoint = False)

    
def calc_c(e, a, b , u = 1): # Check units
    """
    calculate the z components of 4 partial waves in medium
    
    e: dielectric tensor
    
    a,b: components of wavevector in direction of x and y direction
        
    return a list containting 4 roots for the z components of the partial waves
    """
    # assign names 
    x = e * u
    x11, x12, x13 = x[0]
    x21, x22, x23 = x[1] 
    x31, x32, x33 = x[2]
    # calculate the coeffciency based on symbolic expression
    coef4 = x33
    
    coef3 = a*x13 + a*x31 + b*x23 + b*x32
    
    coef2 = a**2*x11 + a**2*x33 + a*b*x12 + a*b*x21 + b**2*x22 + b**2*x33 - \
    x11*x33 + x13*x31 - x22*x33 + x23*x32
    
    coef1 = a**3*x13 + a**3*x31 + a**2*b*x23 + a**2*b*x32 + a*b**2*x13 + \
    a*b**2*x31 + a*x12*x23 - a*x13*x22 + a*x21*x32 - a*x22*x31 + b**3*x23 \
    + b**3*x32 - b*x11*x23 - b*x11*x32 + b*x12*x31 + b*x13*x21
    
    coef0 = a**4*x11 + a**3*b*x12 + a**3*b*x21 + a**2*b**2*x11 + a**2*b**2*x22 \
    - a**2*x11*x22 - a**2*x11*x33 + a**2*x12*x21 + a**2*x13*x31 + a*b**3*x12 + \
    a*b**3*x21 - a*b*x12*x33 + a*b*x13*x32 - a*b*x21*x33 + a*b*x23*x31 + \
    b**4*x22 - b**2*x11*x22 + b**2*x12*x21 - b**2*x22*x33 + b**2*x23*x32 +  \
    x11*x22*x33 - x11*x23*x32 - x12*x21*x33 + x12*x23*x31 + x13*x21*x32 - \
    x13*x22*x31
    
    # calculate the roots of the quartic equation
    c = np.roots([coef4, coef3, coef2, coef1, coef0])
    if len(c) == 2:
        return np.append(c,c)
    return c

def calc_k(e , a, b, u = 1):
    
    """
    A wrapper to calcualte k vector 
    """
    c = calc_c(e, a , b, u)
    return np.array([[a, b, c[0]],[a, b , c[1]],[a,b,c[2]],[a,b,c[3]]])
    
def calc_p(e, k,  u = 1): #something is wrong with this function. Not giving
#correct directions
    """
    Calculate the polarisation vector based on the calculated wavevector and frequency
    equation(9.7-5)
    
    e: dielectric tensor
    
    k: 4x3 array of 4 k vectors
    """
    x = e * u
    p = []
    x11, x12, x13 = x[0]
    x21, x22, x23 = x[1] 
    x31, x32, x33 = x[2]
    for i in k:
        a = i[0]
        b = i[1]
        c = i[2]
        coeff_m = np.array([[x11 - b**2 - c**2, x12 + a * b, x13 + a *c],
                            [x21 + a * b, x22 - a**2 - c **2, x23 + b *c ],
                            [x31 + a * c, x32 + b * c, x33 - a**2 - b**2]])

        # The function seems to return the normalised null spcae vector
        p.append(null(coeff_m))

    return np.array(p)

    
def calc_q(k , p,  u = 1):
    """
    calcualte the direction of q vector based on the k and q vectors given
    
    k: an 4x3 array of 4 k vectors
    
    p: an 4x3 array of 4 p vectors
    
    return a 4x3 array of 4 q vectors
    
    use a special unit for the magnetic field such that c/2pi/mu_0 = 1
    
    note these vectors are not normlised
    """
    
      
    return  np.cross(k ,p)  / u 

def calc_D(p,q):
    
    return np.array([p[:,0],q[:,1], p[:,1], q[:,0]])

def calc_P(k,t):
    return np.diag(np.exp(1j*t*k[:,2]))

def construct_D(e, a, b, omega, u = 1):
    """
    construct the dynamic matrix for one layer with know dielectric tensor
    """
    k = calc_k(e, a, b, omega, u)
    p = calc_p(e , k , omega, u)
    q = calc_q(k , p, omega, u)
    return calc_D(p,q)
    


def null(A, eps=1e-15):
    """
    Return the null vector of matrix A, usefull for calculate the p vector with
    known k vector
    """
    u, s, vh = np.linalg.svd(A)
    null_mask = (s <= eps)
    # relax the threshold if no null singularity is identified
    if null_mask.any() == False:
        return null(A, eps*10)
    
    null_space = np.compress(null_mask, vh, axis=0).flatten()     
    return np.transpose(null_space)

def incident_p(k):
    """
    Calculate the 4 polarisation vectors based on the incident wave wave vector.
    Assuming the medium is isotropic and polarastion is splited into s and p
    k is a 3-vector
    return a array of ps+,ps-,pp+,pp-
    """
    # For normal incidence, fix the direction of polarisation
    # p is aligned with x and s is aligned with y
    # note the p vector is reversed for reflected wave otherwise p,s,k don't form
    # a right hand set
    if k[0] == 0 and k[1] == 0:
        return np.array([[0,1,0],[0,1,0],[1,0,0],[-1,0,0]])
    # calculate the polarisation vectors see lab book for defined geometries
    # Note the normal vector is [0,0,-1] as the incident wave is traving in postive
    # z direction
    si = normalise(np.cross(k, [0,0,-1]))
    sr = si
    pi = normalise(np.cross(si, k))
    pr = normalise(np.cross(sr, [k[0], k[1], -k[2]]))
    # return a 4x3 array of the four polarisation vectors
    return np.array([si,sr,pi,pr])
    
def stack_dot(array):
    """
    Take the matrix product along 1st axis
    """
    product = np.identity(len(array[0]))
    for i in array:
        product = product.dot(i)
    return product

def calc_coeff(T):
    """
    Given the transfer matrix calculate the transmission and reflection coefficients
    Not currently in use
    """
    deno = (T[0,0] * T[2,2] - T[0,2] * T[2,0])
    rss = (T[1,0] * T[2,2] - T[1,2] * T[2,0])/deno
    rsp = (T[3,0] * T[2,2] - T[3,2] * T[2,0])/deno
    rps = (T[0,0] * T[1,2] - T[1,0] * T[0,2])/deno
    rpp = (T[0,0] * T[3,2] - T[3,0] * T[0,2])/deno
    tss = T[2,2]/deno
    tsp = -T[2,0]/deno
    tps = -T[0,2]/deno
    tpp = T[0,0]/deno
    return {"rss":rss, "rsp":rsp, "rps":rps, "rpp":rpp, "tss":tss, "tsp":tsp,
            "tps":tps, "tpp":tpp}
            
def calc_coupling_matrices(T):
    """
    Calculate the coupling matrix between reflected/transmitted light and incident light
    T is the overall transfer matrix of the system. Return a dictionary of coupling matrice
    Indice are always in the order of s,p or L,R
    Note p direction is aligned with x and s is aligned with y in the frame that 
    wave is traveling in the z direction. Refer to geometry guideline in the lab
    book. 
    """
    # Build the coupling matrice between transmitted light and incident/reflected light
    T_ti = np.array([[T[0,0], T[0,2]], [T[2,0], T[2,2]]])
    T_tr = np.array([[T[1,0], T[1,2]], [T[3,0], T[3,2]]])
    # Connect reflected light to incident light using the coupling to transmitted light
    T_ir = np.linalg.solve(T_ti, T_tr)
    T_it = np.linalg.inv(T_ti)
    # Switching to circular polarisation
    # Coupling matrix between planar and circular polarsiation T_cp * [L,R] = [S,P]
    T_cp = np.array([[1j, -1j],[1, 1]])*np.sqrt(1/2)
    T_ir_c = np.linalg.solve(T_cp, T_ir.dot(T_cp))
    T_it_c = np.linalg.solve(T_cp, T_it.dot(T_cp))
    
    coupling_matrices = {"Plane_r":T_ir, "Plane_t":T_it, "Circular_r":T_ir_c, 
                          "Circular_t":T_it_c}
    return coupling_matrices
   
if __name__ == "__main__":
    
    e = np.diag([2,1,1])
    a = 0.5
    b = 0.8
    #print(calc_c(e,a,b,omega))
    k = calc_k(e,a,b)
    p = calc_p(e, k)
    q = calc_q(k , p)
    D = calc_D(p,q)
    