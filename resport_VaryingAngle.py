# -*- coding: utf-8 -*-
"""
Created on Sun May  1 10:55:31 2016
Investigate the reflectance with varying angle of incidence
@author: Bonan
"""
import numpy as np
from matTools import specData
from preset import *
import matplotlib.pyplot as pl
pl.rc('font', size = 8)
#%% Range of Angle
wlRange = np.linspace(300,900,1000)
angleList = np.array(np.arange(0,90))
angleListR = angleList/180*np.pi
res = []
def getResultForAngle(angle):
    s.setIncidence(None, angle, 0)
    s.setPitch([200])
    s.setThickness([2000])
    print("Agnle " + str(angle) + "Finished", flush = True)
    return s.scanSpectrum(wlRange,1)[1]
    
if __name__ == '__main__':
    # Actual computation of the data    
        
    import matplotlib.gridspec as gridspec
    from multiprocessing import Pool
    """    
    with Pool(processes = 4) as pool:
        res = pool.map(getResultForAngle, angleListR) # Select the 2nd element of the result array
    resArray  = np.array(res).swapaxes(0,1)    
    np.save('angleDependence0_89',resArray)
    """
    #%%Plotting data
    resArray = np.load('angleDependence0_89.npy') #Using saved data
    resData = specData(wlRange, resArray)
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
    pl.legend()
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
    gs = gridspec.GridSpec(1, 2, width_ratios=[10,1])
    fig = pl.figure(figsize = (3.3,2.8))
    ax1 = pl.subplot(gs[0])
    ax2 = pl.subplot(gs[1])
    plot1 = ax1.imshow(np.flipud(resArray.T), aspect = 'auto', extent = [wlRange[0],wlRange[-1], 0,89],
              interpolation = 'none')
    pl.colorbar(plot1, ax = ax1)
    ax1.set_title('Spectrum vs incident angle')
    ax1.set_ylabel('Incident angle /degree')
    ax1.set_xlabel('Wavelength /nm')
    ax1.set_xticks(np.arange(wlRange[0],wlRange[-1]+1,100))
    from colourTools import specToRGB
    RGB = np.empty((resArray.shape[1],3))
    
    for i,spec in enumerate(resArray.T):
        RGB[i] = specToRGB([wlRange, spec])
    RGB = np.array([RGB])
    ax2.imshow(np.flipud(RGB.swapaxes(0,1)),aspect = 0.2)
    ax2.set_title('RGB')
    pl.sca(ax2)    
    pl.xticks([])
    pl.yticks([])
    pl.tight_layout(pad = 0)
    pl.savefig(figPath + 'SpectrumVsIncidentAngle.pdf')
    """
    pl.plot(range(90),RGB[:,0],color = 'r')
    pl.plot(range(90),RGB[:,1],color = 'g')
    pl.plot(range(90),RGB[:,2],color = 'b')
    """