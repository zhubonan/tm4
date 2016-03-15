"""
Created on Sun Feb 14 18:02:50 2016
matlab for loading data from spectrometer
@author: Bonan
"""
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks_cwt, savgol_filter
from colourTools import specToRGB
import numpy as np
import matplotlib.pyplot as pl
pi = np.pi


        
class specData():
    """Analysis spectrum from a stack of spec data"""

    def __init__(self, wl, spec, desc = None):
        if len(spec.shape) == 1:
            # convect to rowVector
            spec = spec[:,np.newaxis]
        
        self.n = spec.shape[1]
        self.wl = wl
        self.spec = spec
        self.range = [wl[0], wl[-1]]
        self.spec = self.spec
        self.desc = desc
        
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
        mxIndex = find_peaks_cwt(self.spec, np.arange(*widthRange), **args)
        return mxIndex
        
    def getIndexViaKeyword(self, keyword):
        """Get the indices of the entries with the keyword"""
        desc = self.desc.astype(np.str) #cast the array to string array
        return np.where(np.char.find(desc,keyword) != -1)
        
    def crop(self,wlRange):
        """Return an specData object that is cropped with given range"""
        intersect = self._cropIndices(self.wl, wlRange)
        return specData(self.wl[intersect],self.spec[intersect, :])
        
    def plot(self, indices, showPeaks = False, peakWidthRange = (10,30), **argv):
        """Plot  the current data"""
        pl.figure()
        spec = self.spec[:,indices]
        n = spec.shape[1]
        for i in range(n):
            pl.plot(self.wl, spec[:,i])
            if showPeaks:        
                peakId = find_peaks_cwt(spec[:,i], np.arange(*peakWidthRange), **argv)
                peakCoord =  np.array([self.wl[peakId] , spec[:,i][peakId]]).T
                for p in peakCoord:
                    x = p[0]
                    y = p[1]
                    pl.annotate("{0:.2f}".format(x), (x,y))
                
    def getSelected(self, indices = 'none', keyword = 'none'):
        if keyword != 'none':
            indicesViaKeyword = self.getIndexViaKeyword(keyword)
            if indices != None:
                indicesToUse = np.intersect1d(indicesViaKeyword, indices)
            else:
                indicesToUse = indicesViaKeyword
        elif indices != 'none':
            indicesToUse = indices
        else:
            raise RuntimeError("Need to have at least one argument")
        out = specData(self.wl, self.spec[:,indicesToUse])
        return out
        
    def cropSelf(self,wlRange):
        """Crop the data. This will act one the oject itself"""
        intersect  = self._cropIndices(self.wl, wlRange)
        self.wl = self.wl[intersect]
        self.spec= self.spec[intersect,:]
        self.range = wlRange

    
    def getFilteredSpectrum(self, window_length = 5, polyorder = 2, **args):
        s, n = self.spec.shape
        r = np.zeros((s, n))
        for i in range(n):
            r[:,i] = savgol_filter(self.spec[:,i], window_length, polyorder, **args)
        return specData(self.wl, r, self.desc)
        
    def testFilter(self, window_length, polyorder, n = 0, **sav_args):
        """Test parameters"""
        origin = self.spec[:,n]
        processed = savgol_filter(origin, window_length, polyorder, **sav_args)
        pl.figure()
        pl.subplot(211)        
        pl.plot(self.wl, origin)
        pl.title('Origin')
        pl.subplot(212)
        pl.plot(self.wl, processed)
        
    def getResampledSpectrum(self, ns = 1000, k = 'linear'):
        """Return the Resampled specData object with x axis using 1/wl. Respect
        the cropping on self"""
        thisSpec = self.spec[:,0]
        resSpecData = self._resampledSpectrum(1/self.wl, thisSpec, ns, k)
        resSpec = (resSpecData[1])[:,np.newaxis] #resize to column vector
        for i in range(1,self.n):
            thisSpec = self.spec[:,i]   
            thisResampled = self._resampledSpectrum(1/self.wl, thisSpec, ns, k)[1]
            resSpec = np.append(resSpec, thisResampled[:,np.newaxis], axis = 1)
        return specData(resSpecData[0], resSpec)
    
    def applyFunc2D(self, func, nOutput, *args):
        """Apply a function to each spectrum measurement and return an (n,n,X) array
        where X is the output of the function"""
        s, n = self.spec.shape
        g = np.zeros((n, nOutput))
        for i in range(n):
            g[i,:] = func(self.spec[:,i], *args)
        l = np.sqrt(n)
        gArray = g.reshape((l,l,nOutput))
        #Try to squeeze the array
        try:
            gArray = np.squeeze(gArray, 2)
        except: pass
        return gArray
        
    def get2DColouredImage(self, show = False):
        """Show image by converting spectrum to RGB"""
        s, n = self.spec.shape # n is the number of spectrum
        RGB = np.zeros((n,3))
        for i in range(n):
            RGB[i,:] = specToRGB([self.wl, self.spec[:,i]])
        l = int(np.sqrt(n))
        RGBArray = RGB.reshape((l,l,3))
        if show:
            pl.figure()
            pl.imshow(RGBArray)
            pl.title('RGBImage')
        return RGBArray
    
    def getPeaks(self, dim = 1 ,show = False):
        """Return the peak wavelength and peak height in a tuple"""
        s, n = self.spec.shape
        l = np.sqrt(int(n))
        peakH = np.max(self.spec, axis = 0)
        peakWl = np.zeros(n)
        for i in range(n):
            peakWl[i] = self.wl[self.spec[:,i] == peakH[i]]
        if dim == 2:
            peakH = peakH.reshape((l,l))
            peakWl = peakWl.reshape((l,l))
            if show == True:
                return disp2DPeakData(peakWl, peakH)
        return peakWl, peakH
        
    def get2DGreyScaleImage(self, show = False):
        """Get a 2D grey scale image by summing spec data"""
        s, n = self.spec.shape
        g = np.sum(self.spec, axis = 0)
        l = np.sqrt(n)
        gArray = g.reshape((l,l))
        if show:
            f = pl.figure()
            pl.imshow(gArray, cmap = pl.get_cmap('gray'))
            pl.colorbar()
            wlr = "{0:.0f}".format(self.range[0]) + " to " + "{0:.0f}".format(self.range[1]) 
            pl.title('Grey scale image using range ' + wlr)
            pl.show()
            return f
        return gArray
        
    def __add__(self, other):
        newSpec = np.append(self.spec,other.spec, axis = 1)
        output = specData(self.wl, newSpec)
        return output
#%%
def disp2DPeakData(pWlMat, pHMat):
      """
      Display the data of the peaks. pWlMat and pHMat are matrices of peak 
      wavelength and peak height"""
      pWlF = pl.figure()
      pl.title('Peak wavelength')
      pl.imshow(pWlMat)
      pl.colorbar()
      pHMatF = pl.figure()
      pl.title('Peak hieght')
      pl.imshow(pHMat)
      pl.colorbar()
      pl.show()
      return pWlF, pHMatF

#%%
class scanData:
    """An alternative way to load scan.mat into python"""
    def __init__(self, filename):
        self.scan = loadmat(filename, squeeze_me = False, struct_as_record = False,
                            appendmat = True)['scan']
        #scan is an 1D numpy object array contatining scipy mat_struct objects
        #This is preffered as it let data access easier since we don't care about saving
        self.noOfData = self.scan.size
        self.selectStruct(0)
        self.wl = self.scan[0,0].wl[:,0]

    def selectStruct(self, index):
        """Select the current struct from"""
        self.currentStruct = self.scan[index,0]
        
    def getCurrentSpecData(self):
        return specData(self.wl, self.currentStruct.spec, self.currentStruct.desc)
    
    @staticmethod
    def _convertDescField(desc):
        """Convert the desc field. If desc is a string then make it np array
        if desc is an array then squeeze the data structure"""
        if desc.size != 1:
            descList = list(map(lambda x: desc[x][0][0], range(desc.size)))
            out = np.array(descList)
        else:
            out = desc
        return out
        
    def getCombinedSpecData(self,indices = 'all'):
        """Combine selected spec data, and return an specData object"""
        if indices == 'all':
            indices = np.arange(self.scan.size)
        specList = map(lambda x: self.scan[x,0].spec, indices)
        specOut = np.concatenate(list(specList),axis = 1)
        descList = map(lambda x: self._convertDescField(self.scan[x,0].desc), indices)
        descOut = np.concatenate(list(descList))
        return specData(self.wl , specOut, descOut)
        
if __name__ == '__main__':
    scan = scanData('scan')
    comb = scan.getCombinedSpecData()
    #scan.selectStruct(2)
    #s = scan.getspecData()
    #s.cropSelf((400,700))
    #ss = s.getSelected(indices = np.arange(0,2,116))