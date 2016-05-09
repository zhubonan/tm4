# -*- coding: utf-8 -*-
"""
Created on Tue May  3 21:07:27 2016
Polt a vector field
For showing the rotation of the cholesteric director
@author: Bonan
"""
import numpy as np
import matplotlib.pyplot as pl

# Meshgrid of xyz
x, y = np.meshgrid([0], np.linspace(-450,270,30))

# Define uvw
P = 180
u = np.cos(np.pi * y / P)
v = 0 * y
# Plot
fig = pl.figure(figsize = (5,5*2/3))
ax = fig.add_subplot(111)
ax.quiver(x,y,u,v)
ax.set_xlim(-1,1)
ax.set_ylim(-450,270)
pl.axhspan(0,180,alpha = 0.2, color = 'r')
pl.axhspan(0,-360,alpha = 0.4, color = 'b')
pl.text(-0.75,80,'Apparant pitch \n 180nm')
pl.text(-0.75,-200,'Real pitch \n 360nm')
pl.xticks([])
pl.yticks([200,0,-200,-400])
pl.ylabel('Depth /nm')
fig.tight_layout() #Needed for ensuring labels are inside the saved image
#pl.savefig('../PartIIIReport/fig/LCPitch.pdf')