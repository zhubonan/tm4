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
angleList = np.array(np.arange(0,90))
angleListR = angleList/180*np.pi
res = []
def getResultForAngle(angle):
    s.setIncidence(None, angle, 0)
    s.setPitch([180])
    s.setThickness([1800])
    print("Agnle " + str(angle) + "Finished", flush = True)
    return s.scanSpectrum(wlRange,1)[1]
    
if __name__ == '__main__':
    with Pool(processes = 4) as pool:
        res = pool.map(getResultForAngle, angleListR) # Select the 2nd element of the result array
    resArray  = np.array(res).swapaxes(0,1)
    resData = specData(wlRange, resArray)
    #np.save('angleDependence0_89',resArray)
    #%%Plotting data
    from report_Paras import figPath
    def plotWl():
        pl.figure(figsize = (3.3,2.8))
        pl.rc('font',size = 8)
        for i,angle in list(enumerate(angleListR))[:20]:
            pl.plot(wlRange, resArray[:,i], label = 'theta = ' + '{0:.0f}'.format(angle/np.pi * 180))
        pl.title('Reflectance at various incident angle')
        pl.xlabel('Wavelength /nm')
        pl.ylabel('Reflectance')
        pl.xlim(400,800)
    #plotWl()
    #pl.legend()
    #%% Plot trace of the peaks
    pl.figure(figsize = (3.3, 2.8))
    pl.plot(angleList,resData.getPeaks()[0],'-x')
    pl.title('Trace of peak')
    pl.ylabel('Wavelength /nm')
    pl.xlabel('Incidenct angle /degree')
    pl.tight_layout()
    #pl.savefig(figPath + 'TraceOfPeakVsAngles.pdf')
    #pl.plot(angleList, np.cos(angleList) * 460)
    #%% Plot the 2d representation
    pl.figure(figsize = (3.3,2.8))    
    pl.imshow(np.flipud(resArray.T), aspect = 'auto', extent = [400,800, 0,89],
              interpolation = 'none')
    pl.colorbar()
    pl.title('Spectrum vs incident angle')
    pl.ylabel('Incident angle /degree')
    pl.xlabel('Wavelength /nm')
    #pl.savefig(figPath + 'SpectrumVsIncidentAngle.pdf')