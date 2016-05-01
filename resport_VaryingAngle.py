# -*- coding: utf-8 -*-
"""
Created on Sun May  1 10:55:31 2016
Investigate the reflectance with varying angle of incidence
@author: Bonan
"""
import numpy as np
from multiprocessing import Pool
from matTools import specData
from preset import *

#%% Range of Angle
angleList = np.array([0,10,30,50,55,60,65,70,75,78,80,82,85])
angleList = angleList/180*np.pi
res = []
def getResultForAngle(angle):
    s.setIncidence(None, angle, 0)
    s.setPitch([180])
    print("Agnle " + str(angle) + "Finished", flush = True)
    return s.scanSpectrum(wlRange,1)[1]
    
if __name__ == '__main__':
    with Pool(processes = 4) as pool:
        res = pool.map(getResultForAngle, angleList) # Select the 2nd element of the result array
    resArray  = np.array(res).swapaxes(0,1)
    resData = specData(wlRange, resArray)
    #%%Plotting data
    pl.figure(figsize = (6,4))
    for i,angle in list(enumerate(angleList)):
        pl.plot(wlRange, resArray[:,i], label = 'theta = ' + '{0:.0f}'.format(angle/np.pi * 180))
    pl.title('Reflectance at various incident angle')
    pl.xlabel('Wavelength /nm')
    pl.ylabel('Reflectance')
    #%% Plot trace of the peaks
    pl.figure(figsize = (5,5*2/3))
    pl.plot(angleList,resData.getPeaks()[0],'o')
    pl.plot(angleList, np.cos(angleList) * 460)