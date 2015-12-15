# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 16:17:21 2015
Sympy script to calculate the polynomial for getting vertical component of k in a layer
@author: Bonan
"""

import sympy as sp

def coeff_kc():
    
    sp.var("x11,x12,x13,x21,x22,x23,x31,x32,x33")
    sp.var("a,b,c")
    M = sp.Matrix([[x11 - b**2 - c**2, x12 + a * b, x13 + a * c ], 
                   [x21 + a*b, x22 - a**2 - c**2, x23 + b*c], 
                   [x31 + a *c , x32 + b * c, x33 - a**2 - b **2]])
    
    p = sp.poly(sp.det(M), c)
    print(p.all_coeffs())

    return
    """
    results:
    c**4y term: x33
    c**3 term: a*x13 + a*x31 + b*x23 + b*x32
    c**2 term: a**2*x11 + a**2*x33 + a*b*x12 + a*b*x21 + b**2*x22 + b**2*x33 - x11*x33 + x13*x31 - x22*x33 + x23*x32
    c**1 term: a**3*x13 + a**3*x31 + a**2*b*x23 + a**2*b*x32 + a*b**2*x13 + a*b**2*x31 + a*x12*x23 - a*x13*x22 + a*x21*x32 - a*x22*x31 + b**3*x23 + b**3*x32 - b*x11*x23 - b*x11*x32 + b*x12*x31 + b*x13*x21
    constant term: a**4*x11 + a**3*b*x12 + a**3*b*x21 + a**2*b**2*x11 + a**2*b**2*x22 - a**2*x11*x22 - a**2*x11*x33 + a**2*x12*x21 + a**2*x13*x31 + a*b**3*x12 + a*b**3*x21 - a*b*x12*x33 + a*b*x13*x32 - a*b*x21*x33 + a*b*x23*x31 + b**4*x22 - b**2*x11*x22 + b**2*x12*x21 - b**2*x22*x33 + b**2*x23*x32 + x11*x22*x33 - x11*x23*x32 - x12*x21*x33 + x12*x23*x31 + x13*x21*x32 - x13*x22*x31
    """
"""    
sp.var("p1x_,p2x_,p3x_,p4x_,p1y_,p2y_,p3y_,p4y_")
sp.var("q1x_,q2x_,q3x_,q4x_,q1y_,q2y_,q3y_,q4y_")
sp.var("p1x_ex,p2x_ex,p3x_ex,p4x_ex, p1y_ex, p2y_ex, p3y_ex, p4y_ex")
sp.var("q1x_ex,q2x_ex,q3x_ex,q4x_ex, q1y_ex, q2y_ex, q3y_ex, q4y_ex")
sp.var("A1,A2,A3,A4")
"""

# p1x_ex = p1 dot x * A1 * exp(i gamma1 * t)
def transform_matrix(p_,q_, p , q , gamma , thickness):
    
    p1x_, p2x_, p3x_, p4x_ = p_[:,0]
    p1y_, p2y_, p3y_, p4y_ = p_[:,1]
    q1x_, q2x_, q3x_, q4x_ = q_[:,0]
    q1y_, q2y_, q3y_, q4y_ = q_[:,1]
    p1x, p2x, p3x, p4x = p[:,0] * np.exp(1j* gamma * thickness)
    p1y, p2y, p3y, p4y = p[:,1] * np.exp(1j* gamma * thickness)
    q1x, q2x, q3x, q4x = q[:,0] * np.exp(1j* gamma * thickness)
    q1y, q2y, q3y, q4y = q[:,1] * np.exp(1j* gamma * thickness)
    sp.var("A1,A2,A3,A4")
    M = sp.Matrix([[p1x_,p2x_,p3x_,p4x_, p1x*A1 + p2x*A2 + p3x*A3 + p4x*A4],
               [p1y_,p2y_,p3y_,p4y_, p1y*A1+ p2y * A2+p3y*A3 + p4y*A4],
               [q1x_,q2x_,q3x_,q4x_, q1x*A1 + q2x*A2 + q3x*A3 + q4x*A4],
               [q1y_,q2y_,q3y_,q4y_, q1y*A1+ q2y*A2+q3y*A3 + q4y*A4]])
    sp.var("A1_,A2_,A3_,A4_")
    global result
    result = sp.solve_linear_system_LU(M, (A1_,A2_,A3_,A4_))
    return

transform_matrix(l.p[0],l.q[0],l.p[1], l.q[1], l.k[1][:,2], 10e-9)