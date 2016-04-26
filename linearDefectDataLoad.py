# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:56:37 2016

@author: Bonan
"""
import numpy as np
import matplotlib.pyplot as pl
pl.rcParams['image.cmap'] = 'viridis'
wlRange = np.linspace(400,800,100)
#%% Calculate the average
def avgSpec(wlRange,spec):
    """Calculate the average of the spectrum stack. Averaging is along the 1st axis
    """
    avgData = np.average(spec, axis = 0)
    rng = np.linspace(wlRange[0], wlRange[-1], len(avgData))
    pl.plot(rng, avgData)
#%%
    
if __name__ == "__main__":
    path = "C:/Users/Bonan/OneDrive/Documents/PartIII/Project/linearDefectResults/"
    filename = "20160405-103908AligenTrue.npy"
    data = np.load(path + filename)
    pl.imshow(data)
    pl.figure()
    avgSpec(wlRange, data)
    #%% Compare with single spectrum
    from preset import s
    s.setThickness([2000])
    s.setPitch([180])
    res0 = s.scanSpectrum(wlRange)
    pl.plot(wlRange,res0[1], label = "180 helix")
    s.setPitch([150])
    s.setThickness([3000])
    res1 = s.scanSpectrum(wlRange)
    pl.plot(wlRange,res1[1],label = "150 helix")
    pl.legend()
    pl.xlabel("Wavelength /nm")
    pl.ylabel("Reflectivity(L-L)")
    pl.title("3 Layer defect averaged vs Single helix")
    #%%