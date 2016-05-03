# -*- coding: utf-8 -*-
"""
Created on Sun May  1 17:28:37 2016
Figure for the plot of the reflection band
@author: Bonan
"""
import numpy as np
import matplotlib.pyplot as pl
#%% Define the polynomial
ne = 1.6
no = 1.5
nBar = np.average([ne, no])
epsilonBar = (ne**2 + no**2) / 2
delta = (ne**2 - no**2) /2
def expr(beta):
    alpha = beta * nBar
    return alpha**2 + epsilonBar - (4 * alpha**2*epsilonBar + delta**2)**0.5
beta = np.linspace(0.9,1.1)
pl.figure(figsize = (5,5*2/3))
pl.plot(beta, expr(beta))
pl.plot(beta, np.zeros(50))
pl.title('Existence of reflection band')
pl.xlabel(r'Normalised wavelength $\lambda^\prime = \lambda / \bar{n}P$')
pl.ylabel(r'$m^2$',fontsize = 'medium')
pl.axvspan(0.968, 1.032, alpha = 0.2, ls='solid', lw = 1,color = 'b' )
pl.axhline(0,0.968, 1.032)
pl.annotate('Reflection Band', (0.9725,0.01))
pl.tight_layout()
#pl.savefig('../PartIIIReport/fig/reflectionBand.pdf')