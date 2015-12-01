# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 16:17:21 2015
Sympy script to calculate the polynomial for getting vertical component of k in a layer
@author: Bonan
"""

import sympy as sp
sp.var("x11,x12,x13,x21,x22,x23,x31,x32,x33")
sp.var("a,b,c")
M = sp.Matrix([[x11 - b**2 - c**2, x12 + a * b, x13 + a * c ], 
               [x21 + a*b, x22 -det a**2 - b**2, x23 + b*c], 
               [x31 + a *c , x32 + b * c, x33 - a**2 - b **2]])

p = sp.poly(sp.det(M), c)
print(p.all_coeffs())


"""
results:
c**4 term: b**2
c**3 term: b*x23 + b*x32
c**2 term: a**2*b**2 + a**2*x33 + a*b*x12 + a*b*x21 - b**2*x11 + b**2*x22 + b**2*x33 - x22*x33 + x23*x32
c**1 term: a**3*x13 + a**3*x31 + a**2*b*x23 + a**2*b*x32 + 2*a*b**2*x13 + 2*a*b**2*x31 + a*x12*x23 - a*x13*x22 + a*x21*x32 - a*x22*x31 + b**3*x23 + b**3*x32 - b*x11*x23 - b*x11*x32 + b*x12*x31 + b*x13*x21
constant term: a**4*x11 + a**3*b*x12 + a**3*b*x21 - a**2*b**4 + 2*a**2*b**2*x11 + a**2*b**2*x22 - a**2*x11*x22 - a**2*x11*x33 + a**2*x12*x21 + a**2*x13*x31 + a*b**3*x12 + a*b**3*x21 - a*b*x12*x33 + a*b*x13*x32 - a*b*x21*x33 + a*b*x23*x31 - b**6 + b**4*x11 + b**4*x22 + b**4*x33 - b**2*x11*x22 - b**2*x11*x33 + b**2*x12*x21 + b**2*x13*x31 - b**2*x22*x33 + b**2*x23*x32 + x11*x22*x33 - x11*x23*x32 - x12*x21*x33 + x12*x23*x31 + x13*x21*x32 - x13*x22*x31
"""