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
import time
import matTools as mt
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
        
    def getResultForPoint(self, pointIndex, wlList = None, align = True):
        """This method is to be called for getting the spectrum for a certain point""" 
        if wlList == None:
            wlList = self.wlList
        self.s.setStructure(self.tmp)
        self.s.setThickness(self.t[pointIndex])
        ## Align each helix
        if align == True:
            end = 0
            for i, helix in enumerate(self.s.structures):
                if i > 0: helix.phyParas['aor'] = end
                #print(i,end)
                # Calculate the end angle, not the intrinstic anlge of rotation need to be added
                end = - helix.phyParas['t'] / helix.phyParas['p'] * np.pi + helix.phyParas['aor']    
        #print(self.s.getSubStructureInfo(), flush = True)
        result = self.s.scanSpectrum(wlList, giveInfo = False)[1]
        print('Calculation of point ' + str(pointIndex) + ' finished', flush = True)
        return result
#%% Plotting    
def plotResult(wlList, result, title):
    """Plot the result from a series calculation of spectrum
    wlList is a 1D array of wavelengths.
    
    result: a NxW array where N is the number of calculations and W is the number
    of wavelengths calculated
    """
    r = result
    pl.figure()
    #pl.subplot(211)
    #pl.plot(wlList,r.T)
    #pl.xlim((wlList[0], wlList[-1]))
    #pl.xlabel('Wavelength /nm')
    #pl.ylabel('Reflectivity(L-L)')
    pl.title(title)
    #pl.subplot(212)    
    pl.imshow(r, cmap = 'viridis',aspect = 'auto', interpolation = 'none',
              extent = [wlRange[0], wlRange[-1],  1000, 0])
    pl.xlabel("Wavelength /nm")
    pl.ylabel("Distance /nm")
    pl.colorbar()
    return

def showLayerStructure(c):
    pl.figure()
    x = c.p.T
    # Calculate the ys to be plot
    y = np.append(c.h, np.repeat([[c.d]], len(c.p), axis = 0), axis = 1)
    pl.plot(y,x,'o-')
    pl.xlim(0,c.d+100)
    pl.ylim((c.l,0))
    pl.title('With pitch '+ str([x.phyParas['p'] for x in c.tmp])
    + " incident from right TA= " +  "{0:.2f}".format(c.tmp[1].tiltParas['tiltAngle']))
    pl.xlabel('Height from bottom /nm')
    pl.ylabel('Distance /nm')
    pl.show()
        
#%%
def f1(x):
    return 3000 + 2 * x
    
def f2(x):
    return 1000 + 2 * x
    
def getSaveName():
    return "results\\" + time.strftime("%Y%m%d-%H%M%S")
    
    
if __name__ == '__main__':
    pitchesList1 = [[200,180],[200,180],[200,160],[180,200], [170,200], [160,200]]
    pitchesList2 = [[210,180],[210,170],[210,160],[180,210], [170,210], [160,210]]
    pitchesList3 = [[200,180],[180,200],[210,180],[180,210]]
    pitchesList4 = [[180,180]]
    tiltList = [0] * 4 
    nop = 200 #Number of points to sample along the defect
    for pitches, tilt in zip(pitchesList3,tiltList):
        wlRange = np.linspace(400,800,200)
        h1 = heli(CNC,pitches[0],1000) #The thickness doesn't matter here
        h2 = heli(CNC, pitches[1] ,1000)
        h2.setTilt(tilt,[0,1,0])
        h3 = heli(CNC, pitches[0] ,1000)
        tmp = [h1,h2, h3]
        #%% Set layer structure
        c = CrossSection(s, 5000,1000,3)
        c.setInterfaceFunction(f1,0)
        c.setInterfaceFunction(f2,1)
        c.calcPixelConfig(nop)
        c.setLayerTemplate(tmp)
        c.setWlList(wlRange)
        res = []
        #for i in range(4):c.
       #     res.append(c.getResultForPoint(i))
        #c.getResultForPoint(0)
        #%% We reserve the choice of wether align or not here
        from functools import partial
        for alignment in [True]:
            getPoint = partial(c.getResultForPoint, align = alignment)
            if 1:
                t = clock()
                pool = Pool(processes = 7)
                res = pool.map(getPoint, range(nop))
                np.save(getSaveName()+ "Aligen" + str(alignment), res)
                print(clock()-t)
            plotResult(wlRange, np.array(res), title= 'Alignment ' + str(alignment))
            pl.savefig((getSaveName()+ "Aligen" + str(alignment)))
            #%% Save plot of layer structure
            showLayerStructure(c)
            pl.savefig(getSaveName()+"Structure")
            #%% Plotting the colour band figure
            resArray = np.array(res)
            spec = mt.specData(wlRange,resArray.T)
            RGB = spec.getRGBArray()
            pl.figure()
            pl.imshow(RGB.reshape((nop,1,3)),aspect='auto', extent=[0,1,1000,0])
            pl.ylabel("Distance /nm")
            pl.title("RGB colour from spectrum as different distance")
            pl.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # labels along the bottom edge are off
            pl.savefig(getSaveName()+ "Aligen" + str(alignment) + "CBand")
            
        #%% Close the pool
        pool.close()
        pool.join()
        pl.figure()