# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 18:02:50 2016
matlab for loading data from spectrometer
@author: Bonan
"""
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.fftpack import fft, fftfreq
import numpy as np

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
        for i in self.data['desc']:
            print(i)
            
    def getSpectrumViaIndex(self, index):
        
        spec = self.data['spec'][index]
        #Each data['sepc'] contine a number of 1d arrays each for a recording
        wl = self.data['wl'][index]
        return wl, spec
        
    def getSpectrumViaKeyWrod(self, keyword):
        result = []
        for i in range(len(self.data['desc'])):
            if self.data['desc'][i].find(keyword) != -1:
                result.append(self.getSpectrumViaIndex(i))
        return np.array(result)
                
    def getDescViaIndex(self, index):
        return self.data['desc'][index]
        
    def getcroppedSpectrum(self, index, wlRange):
        """Return the cropped spectrum given spectrum range, droping the noise"""
        return self.cropSpectrum(self.getSpectrumViaIndex(index), wlRange)
        
    def getResampledFSpectrum(self, index, wlRange, ns = 2000, k = 'linear'):
        """Return a reSampled spectrum with 1/nm as the x axis, y axis is still
        the reflectrivity(as before)"""
        spectrum = self.getcroppedSpectrum(index,wlRange)
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
    
    def getFTSpectrum(self, index, wlRange):
        """Calculate the FFTed resampledSpectrum in 1/nm base so the axis should now
        be in nm. Can use this to see the effect of thin film interference and therby
        find the thickness of the film
        """
        #Fist load the relvant data
        spectrum = self.getResampledFSpectrum(index, wlRange)
        #Then we make the main to be zero
        meanValue = np.mean(spectrum[1])
        spectrum[1] = spectrum[1] - meanValue #So now the mean is zero
        n = np.size(spectrum[1])
        fspec = fft(spectrum[1])/2/np.pi
        fX = fftfreq(np.size(spectrum[0]), spectrum[0,1] - spectrum[0,0])
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
    data = scanData('scan')
    resampled = data.getcroppedSpectrum(2,[370,870])
    #resampled = data.cropSpectrum(resampled, [0.0012,0.0014])
    fted = data.getFTSpectrum(2,[700,800])
    #pl.plot(fted[0], fted[1])
    #pl.xlim(0,10000)
    pl.plot(resampled[0],resampled[1])
    pl.show()