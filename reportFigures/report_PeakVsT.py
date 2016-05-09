# -*- coding: utf-8 -*-
"""
Created on Thu May  5 21:19:15 2016

@author: Bonan
"""

from preset import *
from matTools import specData

import scipy as sp
from scipy.stats import linregress
wlRange = np.linspace(450,650,100)
l = 5000
def calc(p):
    s.setThickness([l])
    # Define the range of pitch
    s.setPitch([p])
    return s.scanSpectrum(wlRange,1)[1]

if __name__ == '__main__':
    pList = np.linspace(150,200,50)
    from multiprocessing import Pool
    with Pool(processes = 7) as pool:
        res = pool.map(calc, pList)
    resArray = np.array(res).T
    spec = specData(wlRange,resArray)
    pData = spec.getPeaks()
    #%% Plotting
    pl.rc('font', size = 8)
    pl.rc('figure', figsize = (3.3,2.8))    
    pl.figure()
    pl.plot(pData[0],pData[1], 'x')
    pl.xlabel('Wavelength /nm')
    pl.xlim(wlRange[0],wlRange[-1])
    pl.ylabel('Peak Reflectance')
    pl.ylim(0,1)
    #pl.savefig('L' + str(l) + '.png')
    fitresult = linregress(pData)
    fiteq = lambda x: fitresult.slope * x + fitresult.intercept
    pl.plot(wlRange, fiteq(wlRange),'--')
    equation = 'y = ' + '{0:.2e}'.format(fitresult.slope) + ' x + ' + '{0:.2f}'.format(fitresult.intercept) + \
    '  R = ' + '{0:.3f}'.format(fitresult.rvalue)
    pl.annotate(equation, (pData[0][0], pData[1][0] + 0.05))