# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 21:27:20 2016
2D Scan of spectrum
@author: Bonan
"""
from matTools import scanData, specData
import numpy as np
from scipy.signal import find_peaks_cwt
from colourTools import specToRGB
path = r'C:\Users\Bonan\OneDrive\Documents\PartIII\Project\20160314\\' 
def getFilename(filename):
    return path + filename
    
print(getFilename('data'))
mat = scanData(getFilename('scanLite'))
spec = mat.getCurrentSpecData()
spec.cropSelf((400,700))
#%%
def convertToRGB(wl, spec):
    s, n = spec.shape # n is the number of spectrum
    RGB = np.zeros((n,3))
    for i in range(n):
        RGB[i,:] = specToRGB([wl, spec[:,i]])
    l = int(np.sqrt(n))
    RGBArray = RGB.reshape((l,l,3))
    return RGBArray

def customGreyScale(wl,spec, wlRange):
    """Convert the 2D Scan spec data to grey scale with given averaging range"""
    specObj = specData(wl,spec)
    specObj.cropSelf(wlRange)
    s, n = specObj.currentSpec.shape
    g = np.sum(specObj.spec, axis = 0)
    l = np.sqrt(n)
    gArray = g.reshape((l,l))
    return gArray
    
def findPeakCWT(wl, spec, widthRange,**args):
    pIndex = find_peaks_cwt(spec, np.arange(1,widthRange), **args)
    pValues = spec[pIndex]
    pWl = wl[pIndex]
    mxValue = np.max(pValues)
    return pWl[pValues == mxValue], mxValue