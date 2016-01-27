# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:48:11 2016
Experimentally: fallback to the function in Yeh's book to calculate P with given k
@author: Bonan
"""
import numpy as np
def calc_p_alternative(e, kx, ky, kz):
    """
    e is the relative dielectric tensor 3*3
    kz is a 4 array containing 4 kz values
    return a 4*3 array contatining 4 p vectors in the same order a k given
    due to scaling there is no frequency dependence here
    """
    # extract values from the dielectric tensor
    exx, exy, exz = e[0]
    eyx, eyy, eyz = e[1]
    ezx, ezy, ezz = e[2]
    p = []
    for kz in kz:
        p1 = (eyy - kx**2 - kz**2)*(ezz - kx**2 - ky**2) - (eyz + ky*kz)**2
        p2 = (eyz + ky * kz)*(ezx + kx* kz)-(exy + kx * ky)*(ezz - kx**2 - ky**2)
        p3 = (exy + ky * kx)*(eyz + ky * kz) - (exz + kx *kz) *(eyy - kx**2 - kz**2)
        p.append([p1,p2,p3])
    return np.array(p)