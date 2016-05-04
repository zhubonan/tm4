# -*- coding: utf-8 -*-
"""
Created on Wed May  4 12:39:18 2016
Plot the change of reflectance with the number of division
@author: Bonan
"""
import simClasses as sim
import matplotlib.pyplot as pl
import numpy as np
from preset import CNC
from multiprocessing import Pool
#Define the range of division we want
div = [2,3,10,20,30,35]
airHalf = sim.airHalfSpace
glassHalf = sim.glassHalfSpace
s = sim.OptSystem()
s.setHalfSpaces(airHalf, glassHalf)
wlRange = np.linspace(400,800,200)

def calcDivLL(div):
    h1 = sim.HeliCoidalStructure(CNC, 180, 1800, d = div)
    s.setStructure([h1])
    return s.scanSpectrum(wlRange,1)

def calcDivRR(div):
    h1 = sim.HeliCoidalStructure(CNC, 180, 1800, d = div)
    s.setStructure([h1])
    return s.scanSpectrum(wlRange,1,coupling='RR')
    
if __name__ == '__main__':
    from report_Paras import figPath
    with Pool(processes = 4) as pool:
        resLL = pool.map(calcDivLL, div)
        resRR = pool.map(calcDivRR,div)
        #%%
        pl.rc('font', size = 9)
        pl.rc('')
        pl.figure(figsize = (3.3,5))
        pl.subplot(211)
        for i, (wl, spec) in enumerate(resLL):
            pl.plot(wl,spec, label = 'd= ' + '{0:.0f}'.format(div[i]))
        pl.ylim(0,0.3)
        pl.ylabel('Reflectance LL')
        pl.legend(fontsize = 8)
        pl.title('Reflectance and number of division')
        pl.subplot(212)
        for i, (wl, spec) in enumerate(resRR):
            pl.plot(wl,spec, label = 'd= ' + '{0:.0f}'.format(div[i]))
        pl.ylim(0,0.3)
        pl.legend(fontsize = 8)
        pl.xlabel('Wavelength /nm')
        pl.ylabel('Reflectance RR')
        pl.tight_layout()
    #pl.savefig(figPath+'RvsDiv.pdf')