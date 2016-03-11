# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 18:02:50 2016
matlab for loading data from spectrometer
@author: Bonan
"""
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks_cwt
import numpy as np
import matplotlib.pyplot as pl
pi = np.pi


        
class specDataStack():
    """Analysis spectrum from a stack of spec data"""

    def __init__(self, wl, spec):
        if len(spec.shape) == 1:
            # convect to rowVector
            spec = spec[:,np.newaxis]
        
        self.n = spec.shape[1]
        self.wl = wl
        self.spec = spec
        self.range = [wl[0], wl[-1]]
        self.currentSpec = self.spec
        self.currentWl = self.wl
        self.currentRange = self.range
    
    @staticmethod
    def _cropIndices(wl, wlRange):
        """return the indices use for cropping spec data"""
        upper = np.where(wlRange[0] < wl)
        lower = np.where(wlRange[1] > wl)
        return np.intersect1d(lower, upper)    

    @staticmethod
    def _FTSpectrum(wl, spec, ns,  paddling, k = 'linear'):
        """Calculate the FFTed resampledSpectrum in 1/nm base so the axis should now
        be in nm. Can use this to see the effect of thin film interference and therby
        find the thickness of the film. The x axis is n*d(p)
        
        return an specData object
        """
        #Then we make the main to be zero
        meanValue = np.mean(spec)
        spec = spec - meanValue #So now the mean is zero
        # ZeroPaddling
        paddled = np.append(spec, np.zeros(paddling))
        fspec = fft(paddled)
        n= paddled.size
        # We want to scale the x axis such it is labed as nd. Note the ffted frequency
        # is in the unit of hertz.
        fX = fftfreq(len(paddled),abs( wl[1] - wl[0]))
        fY = np.abs(fspec)
        if n%2 == 0:
            t = n/2
        else:
            t = (n+1)/2
        fX = fX[0:t]
        fY = fY[0:t]
        return (fX,fY)

    @staticmethod    
    def _resampledSpectrum(wl, spec, ns, k):
        """Return a reSampled spectrum with 1/nm as the x axis, y axis is still
        the reflectrivity(as before)
        
        return an specData object"""
        #Resampling
        interp = interp1d(wl, spec, kind = k )
        #Get range
        newX = np.linspace(wl.min(),wl.max(), ns)
        newY = interp(newX)
        return (newX, newY)
        
    def findPeakIndex(self, widthRange = (1,100), **args):
        """Find the peak on the active data"""
        mxIndex = find_peaks_cwt(self.currentSpec, np.arange(*widthRange), **args)
        return mxIndex
        
        
    def getCropped(self,wlRange):
        """Return an specData object that is cropped with given range"""
        intersect = self._cropIndices(self.wl, wlRange)
        return specDataStack(self.wl[intersect, :],self.spec[intersect, :])
        
    def plot(self, showPeaks = False, peakWidthRange = (1,100), **argv):
        """Plot  the current data"""
        pl.figure()
        pl.plot(self.currentWl, self.currentSpec)
        if showPeaks:        
            peakId = self.findPeakIndex(peakWidthRange,**argv)
            peakCoord =  np.array([self.currentWl[peakId] , self.currentSpec[peakId]]).T
            for p in peakCoord:
                x = p[0]
                y = p[1]
                pl.annotate("{0:.2f}".format(x), (x,y))    
                
    def cropSelf(self,wlRange):
        """Crop the data. This will act one the oject itself, changing the current
        wl,spec and range"""
        intersect  = self._cropIndices(self.wl, wlRange)
        self.currentWl = self.wl[intersect]
        self.currentSpec= self.spec[intersect,:]
        self.currentRange = wlRange    

    
    def getResampledSpectrum(self, ns = 1000, k = 'linear'):
        """Return the Resampled specDataStack object with x axis using 1/wl. Respect
        the cropping on self"""
        thisSpec = self.currentSpec[:,0]
        resSpecData = self._resampledSpectrum(1/self.currentWl, thisSpec, ns, k)
        resSpec = (resSpecData[1])[:,np.newaxis] #resize to column vector
        for i in range(1,self.n):
            thisSpec = self.currentSpec[:,i]   
            thisResampled = self._resampledSpectrum(1/self.currentWl, thisSpec, ns, k)[1]
            resSpec = np.append(resSpec, thisResampled[:,np.newaxis], axis = 1)
        return specDataStack(resSpecData[0], resSpec)
        
    def append(self, other):
        self.spec = np.append(self.spec, other.spec, axis = 1)

            
    def __add__(self, other):
        newSpec = np.append(self.spec,other.spec, axis = 1)
        output = specDataStack(self.wl, newSpec)
        if self.currentRange == other.currentRange:
            output.cropSelf(self.currentRange)
        return output
        
class scanData:
    """An alternative way to load scan.mat into python"""
    def __init__(self, filename):
        self.scan = loadmat(filename, squeeze_me = True, struct_as_record = False,
                            appendmat = True)['scan']
        #scan is an 1D numpy object array contatining scipy mat_struct objects
        #This is preffered as it let data access easier since we don't care about saving
        self.noOfData = self.scan.size
        self.selectStruct(0)
        self.wl = self.currentStruct.wl

    def selectStruct(self, index):
        """Select the current struct from"""
        self.currentStruct = self.scan[index]
    
    def getIndexViaKeyword(self, keyword):
        """Get the indices of the entries with the keyword"""
        desc = self.currentStruct.desc.astype(np.str) #cast the array to string array
        return np.where(np.char.find(desc,keyword) != -1)
        
    def getCurrentSpecData(self):
        return specDataStack(self.wl, self.currentStruct.spec)
        
    def plotCurrentSpec(self, indices = None, keyword = None):
        if keyword != None:
            indicesViaKeyword = self.getIndexViaKeyword(keyword)
            if indices != None:
                indicesToUse = np.intersect1d(indicesViaKeyword, indices)
            else:
                indicesToUse = indicesViaKeyword
        elif indices != None:
            indicesToUse = indices
        else:
            thisSpecData = specDataStack(self.wl, self.currentStruct.spec)
            thisSpecData.plot()
            return
        thisSpecData = specDataStack(self.wl, self.currentStruct.spec[:,indicesToUse])
        thisSpecData.plot() 
    
if __name__ == '__main__':
    scan = scanData('scan')
    scan.selectStruct(2)
    s = scan.getCurrentSpecData()
    s.cropSelf((400,700))
    s2 = s + s 
    sr = s2.getResampledSpectrum()
    sr.plot()