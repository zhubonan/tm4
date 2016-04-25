# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:56:37 2016

@author: Bonan
"""
import numpy as np
import matplotlib.pyplot as pl
pl.rcParams['image.cmap'] = 'viridis'


#%% Calculate the average
def avgSpec(spec):
    """Calculate the average of the spectrum stack. Averaging is along the 1st axis
    """
    avgData = np.average(spec, axis = 0)
    pl.plot(np.arange(200), avgData)
#%%
    
if __name__ == "__main__":
    path = "C:/Users/Bonan/OneDrive/Documents/PartIII/Project/linearDefectResults/"
    filename = "20160404-202720AligenFalse.npy"
    data = np.load(path + filename)
    pl.imshow(data)
    pl.figure()
    avgSpec(data)