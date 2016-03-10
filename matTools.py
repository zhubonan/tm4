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

class specData:
    """A class for storing/anaylsing/manipulating a single spectrum"""
    def __init__(self, wl, spec):
        self.wl = wl
        self.spec = spec
        self.range = [wl[0], wl[-1]]
        self.currentSpec = self.spec
        self.currentWl = self.wl
        self.currentRange = self.range
        
    def crop(self,wlRange):
        """Crop the data. This will act one the oject itself, changing the current
        wl,spec and range"""
        wl = self.wl
        spec = self.spec
        upper = np.where(wlRange[0] < wl)
        lower = np.where(wlRange[1] > wl)
        intersect = np.intersect1d(lower, upper)
        self.currentWl = wl[intersect]
        self.currentSpec= spec[intersect]
        self.currentRange = wlRange
        
    def findPeakIndex(self, widthRange = (1,100), **args):
        """Find the peak on the active data"""
        mxIndex = find_peaks_cwt(self.currentSpec, np.arange(*widthRange), **args)
        return mxIndex
    
    def getCropped(self,wlRange):
        """Return an specData object that is cropped with given range"""
        wl = self.wl
        spec = self.spec
        upper = np.where(wlRange[0] < wl)
        lower = np.where(wlRange[1] > wl)
        intersect = np.intersect1d(lower, upper)
        return specData(wl[intersect],spec[intersect])
        
    def getResampledSpectrum(self, ns = 10000, k = 'linear'):
        """Return a reSampled spectrum with 1/nm as the x axis, y axis is still
        the reflectrivity(as before)
        
        return an specData object"""
        spec = self.currentSpec
        #Change x coorindates to 1/nm
        wl = 1 / self.currentWl
        #Resampling
        interp = interp1d(wl, spec, kind = k )
        #Get range
        newX = np.linspace(wl.min(),wl.max(), ns)
        newY = interp(newX)
        return specData(newX, newY)
        
    def getCurrentSpectrum(self):
        """Return the current spectrum with a tuple (wl, spec)"""
        return [self.currentWl, self.currentSpec]
        
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
        
    def getFTSpectrum(self, ns = 10000, paddling = 100000, k = 'linear'):
        """Calculate the FFTed resampledSpectrum in 1/nm base so the axis should now
        be in nm. Can use this to see the effect of thin film interference and therby
        find the thickness of the film. The x axis is n*d(p)
        
        return an specData object
        """
        #Fist load the relvant data
        resampled = self.getResampledSpectrum(ns)
        wl,spec = resampled.wl, resampled.spec
        #Then we make the main to be zero
        meanValue = np.mean(spec)
        spec = spec - meanValue #So now the mean is zero
        # ZeroPaddling
        paddled = np.append(spec, np.zeros(paddling))
        fspec = fft(paddled)
        n= paddled.size
        # We want to scale the x axis such it is labed as nd. Note the ffted frequency
        # is in the unit of hertz.
        fX = fftfreq(len(paddled), wl[1] - wl[0])/2
        fY = np.abs(fspec)
        if n%2 == 0:
            t = n/2
        else:
            t = (n+1)/2
        fX = fX[0:t]
        fY = fY[0:t]
        return specData(fX,fY)


class scanData2:
    """An alternative way to load scan.mat into python"""
    def __init__(self, filename):
        self.scan = loadmat(filename, squeeze_me = True, struct_as_record = False,
                            appendmat = True)['scan']
        #scan is an 1D numpy object array contatining scipy mat_struct objects
        #This is preffered as it let data access easier since we don't care about saving
        self.noOfData = self.scan.size
        self.selectStruct(0)
        self.selectWl(0)

    def selectStruct(self, index):
        """Select the current struct from"""
        self.currentStruct = self.scan[index]
    
    def selectWl(self, index):
        """return the wl given the index of the element want to use"""
        wl = self.scan[index].wl
        if len(wl.shape) == 2:
            self.wl = wl[:,0]
        else:
            self.wl = wl
    def selectSpec(self,index,index2 = -1):
        """Select a single spec data"""
        if index2 != -1:
            self.currentSpec = (self.scan[index].spec)[:,index2]
        else:
            self.currentSpec = (self.scan[index].sepc)
        
    def plotCurrentStructSpec(self):
        spec = self.currentStruct.spec
        h = pl.figure()        
        pl.plot(self.wl, spec)
        return h
        
    def plotCurrentSpec(self):
        pl.plot(self.wl, self.currentSpec)
        pl.ylim(0,1)
        pl.xlim(400,800)
 

class scanData:
    """A class for easy control of scan.mat data"""
    def __init__(self, filename):
        """Loaded the .mat file"""
        self.data = loadmat(filename, squeeze_me = True, appendmat = True )['scan']
        self.noOfData = self.data['wl'].shape[0]
        
    @staticmethod
    def cropSpectrum(spectrum, wlRange):
        """Return the cropped spectrum given spectrum range, droping the noise"""
        wl = spectrum[0]
        spec = spectrum[1]
        upper = np.where(wlRange[0] < wl)
        lower = np.where(wlRange[1] > wl)
        intersect = np.intersect1d(lower, upper)
        return np.vstack((wl[intersect], spec[intersect]))
        
    def listAllDescription(self):
        count = 0
        for i in self.data['desc']:
            print(i, '   index:' ,count)
            count += 1
            
    def getSpectrumViaIndex(self, index, field = None):
        
        if field != None:
            spec = self.data['spec'][index][:, field]
            wl = self.data['wl'][index][:, field]
        else:
            spec = self.data['spec'][index]
            wl = self.data['wl'][index]
        #Each data['sepc'] contine a number of 1d arrays each for a recording
        return wl, spec
        
    def _getSpectrumViaKeyWrod(self, keyword):
        result = []
        for i in range(len(self.data['desc'])):
            if self.data['desc'][i].find(keyword) != -1:
                result.append(self.getSpectrumViaIndex(i))
        return np.array(result)
                
    def getDescViaIndex(self, index, field = None):
        if field != None:
            return self.data['desc'][index][field]
        else:
            return self.data['desc'][index]
        
    def getcroppedSpectrum(self, index, field = None, wlRange = [300,800] ):
        """Return the cropped spectrum given spectrum range, droping the noise"""
        if field != None:
            return self.cropSpectrum(self.getSpectrumViaIndex(index, field), wlRange)
        else:
            return self.cropSpectrum(self.getSpectrumViaIndex(index), wlRange)
        
    def getResampledSpectrum(self, index,  field = None, wlRange = [300,800], ns = 10000, k = 'linear'):
        """Return a reSampled spectrum with 1/nm as the x axis, y axis is still
        the reflectrivity(as before)"""
        spectrum = self.getcroppedSpectrum(index,field, wlRange)
        #Change x coorindates to 1/nm
        spectrum[0] = 1/spectrum[0]
        #Resampling
        interp = interp1d(spectrum[0], spectrum[1], kind = k )
        #Get range
        mx = spectrum[0].max()
        mn = spectrum[0].min()
        newX = np.linspace(mn,mx,ns)
        newY = interp(newX)
        
        return np.vstack((newX, newY))
    
    def getFTSpectrum(self, index, field = None, wlRange = [300,800],  ns = 10000, paddling = 100000, k = 'linear'):
        """Calculate the FFTed resampledSpectrum in 1/nm base so the axis should now
        be in nm. Can use this to see the effect of thin film interference and therby
        find the thickness of the film. The x axis is 
        """
        #Fist load the relvant data
        spectrum = self.getResampledSpectrum(index, field, wlRange, ns, k)
        #Then we make the main to be zero
        meanValue = np.mean(spectrum[1])
        spectrum[1] = spectrum[1] - meanValue #So now the mean is zero
        # ZeroPaddling
        paddled = np.append(spectrum[1], np.zeros(paddling))
        fspec = fft(paddled)
        n= paddled.size
        # We want to scale the x axis such it is labed as nd. Note the ffted frequency
        # is in the unit of hertz.
        fX = fftfreq(len(paddled), spectrum[0,1] - spectrum[0,0])/2
        fY = np.abs(fspec)
        if n%2 == 0:
            t = n/2
        else:
            t = (n+1)/2
        fX = fX[0:t]
        fY = fY[0:t]
        return fX,fY
               
if __name__ == '__main__':
    import matplotlib.pyplot as pl
    scan = scanData2('scan')
    scan.selectSpec(2,0)
    spectrum = specData(scan.wl,scan.currentSpec)
    spectrum.crop((400,700))
    spectrum.getFTSpectrum().plot()
    pl.xlim(0,10000)