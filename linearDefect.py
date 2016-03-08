"""
Python script for plotting and simulate 1-D defect structure. e.g. get sepctrum 
along a line
"""
import simClasses as sim
import matplotlib.pyplot as pl
import numpy as np
import copy as cp
from colourTools import specToRGB
from multiprocessing import Pool
from time import clock
pl.rcParams['figure.figsize'] = (8,6)
pl.rcParams['savefig.dpi'] = 100
#%%
#define materials
CNC = sim.UniaxialMaterial(1.586,1.524) # this might not be right
air = sim.HomogeneousNondispersiveMaterial(1)
cellulose = sim.HomogeneousNondispersiveMaterial(1.55)
# Set the glass as a practical backend
glass = sim.HomogeneousNondispersiveMaterial(1.55)
front = sim.IsotropicHalfSpace(air)
airhalf= sim.IsotropicHalfSpace(air)
glasshalf = sim.IsotropicHalfSpace(glass)
s = sim.OptSystem()
s.setHalfSpaces(airhalf,glasshalf)
heli = sim.HeliCoidalStructure
h1 = heli(CNC, 150, 1000)
spacer = sim.HomogeneousStructure(CNC, 10)
airSpacer = sim.HomogeneousStructure(air,10)
glassSpacer = sim.HomogeneousStructure(glass,10)
s.setStructure([h1])
wlRange = np.linspace(350,850,50)
#%% Functions
def plotSpectrum(OptSys, wlRange):
    result = OptSys.scanSpectrum(wlRange)
    pl.plot(result[0],result[1], color = specToRGB(result))
    return
    
def applyConfig(optSys, configArray, templateList):
    """Apply configuration to an Optsys object. configList should be an Nx2 
    ndarray object. Where N is the number of layers. Two columns are thickness and 
    pitch repectively.Example:
        [[1000 150]
         [1500 200]
         [2000 150]]
    
    The tempatelist will contain object of structures whose copy will be used.
    to added to the structure list of the optSys
    """
    # first check if the number of layers agree
    if len(configArray) != len(templateList):
        raise RuntimeError("Length of input don't agree")
    slist = []
    for i in range(len(configArray)):
        struc = cp.copy(templateList[i])
        if type(struc) == sim.HeliCoidalStructure:
            struc.setPitch(configArray[i][0])
        struc.setThickness(configArray[i][1])
        slist.append(struc)
    optSys.setStructure(slist)
    return

class CrossSection():
    """A class represent a CrossSection of the film"""
    def __init__(self, optSys,depth, length, nol):
        """Initialise a CrossSection object"""
        self.d = depth # depth, we assume the flim is flat
        self.l = length # Length of the cross-section
        self.n = nol # Number of layers
        self.s = optSys # The intantiated OptSystem object
        self.interface = [0] * (nol-1) # Preallocate the 
        self.wlList = None
        
    def setLayerTemplate(self, template):
        """Set up the type of layers using the a list of existing Strucuture objects.
        The length of list must match the number of layer. The Structure need to
        be instantiated(e.g. cannot directly pass a class)
        """
        if len(template) != self.n:
            raise RuntimeError("The tempate must match the total number of layers")
        self.tmp = template
        
    def setInterfaceFunction(self, func, i):
        """Set the function of the ith interface. The function must be callable 
        with a single argument"""
        self.interface[i] = func
        
    def calcPixelConfig(self, div):
        """This will calcuate the thickness of the layers under each pixel
        
        d: divisions of the cross-section
        """
        self.p = np.linspace(0, self.l, div) #location of the points
        # Calculate the location of interfaces
        if (np.array(self.interface) == 0).any():
            raise RuntimeError('At least one of the intface is undefined')
        
        h = np.empty((div, self.n - 1)) # array of the height of each interface
        for i in range(len(self.interface)):
            vec = np.vectorize(self.interface[i])
            h[:,i] = vec(self.p)
        # Now we need to convert the interface height to layer thickness
        self.h = h # The array of height of the interfaces
        h1 = np.append(h,np.zeros((div,1)),axis = 1)
        h2 = np.append(np.repeat([[self.d]], div, axis = 0),h, axis = 1)
        self.t = h2 - h1 # t is the array of thickness for each layer a all points
        return self.t
        
    def showLayerStructure(self):
        x = np.repeat([self.p], self.n, axis = 0).T
        # Calculate the ys to be plot
        y = np.append(self.h, np.repeat([[self.d]], len(self.p), axis = 0), axis = 1)
        pl.plot(x,y,'o-')
        pl.ylim(0,self.d+100)
        pl.title('Cross-section')
        pl.xlabel('distance /nm')
        pl.ylabel('Depth /nm')
        pl.show()
        
    def getSpectrum(self, wlList = None, showResult = True):
        """Calculate the spectrum for each point with given list of wavelengths"""
        if wlList == None:
            wlList = self.wlList
        result = []
        n = len(self.t)                
        for i in range(n):
            print('Calculating point ' + str(i+1) +  ' out of ' + str(n) + '...', flush = True)
            self.s.setStructure(self.tmp)
            self.s.setThickness(self.t[i])
            result.append(self.s.scanSpectrum(wlList, giveInfo = False)[1]) #Select the spec data only
        result = np.array(result)
        if showResult: self.plotResult(wlList, result)
        return result
        
    def setWlList(self,wlList):
        self.wlList = wlList
        
    def getResultForPoint(self, pointIndex, wlList = None):
        """This method is to be called for getting the spectrum for a certain point""" 
        if wlList == None:
            wlList = self.wlList
        self.s.setStructure(self.tmp)
        self.s.setThickness(self.t[pointIndex])
        result = self.s.scanSpectrum(wlList, giveInfo = False)[1]
        print('Calculation of point ' + str(pointIndex) + ' finished', flush = True)
        return result
        
    @staticmethod
    def plotResult(wlList, result):
        r = result
        pl.figure(1)
        pl.subplot(211)
        pl.plot(wlList,r.T)
        pl.xlim((wlList[0], wlList[-1]))
        pl.xlabel('Wavelength /nm')
        pl.ylabel('Reflectivity(L-L)')
        pl.subplot(212)    
        pl.imshow(r, aspect = 'auto', interpolation = 'none')
        pl.legend()
        pl.show()
        return
def f1(x):
    return 4500 - 2 *x
    
def f2(x):
    return 500 + 2 *x
    
if __name__ == '__main__':
    wlRange = np.linspace(450,670,200)
    h1 = heli(CNC,165,1000)
    h2 = heli(CNC, 200 ,1000)
    h2.Phi = np.pi/4;
    
    tmp = [h2,h1, h2]
    #%% Set layer structure
    c = CrossSection(s, 5000,1000,3)
    c.setInterfaceFunction(f1,0)
    c.setInterfaceFunction(f2,1)
    c.calcPixelConfig(50)
    c.setLayerTemplate(tmp)
    #r = c.getSpectrum(wlRange)
    c.setWlList(wlRange)
    if 1:
        t = clock()
        pool = Pool(processes = 2)
        res = pool.map(c.getResultForPoint, range(50))
        print(clock()-t)
    c.plotResult(wlRange, np.array(res))