# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 16:55:52 2016

@author: Bonan
"""

"""
Python script for plotting and simulate 1-D defect structure. e.g. get sepctrum 
along a line
"""
import simClasses as sim
import matplotlib.pyplot as pl
import numpy as np
import copy as cp
from colourTools import specToRGB
from time import clock
import time
import matTools as mt
from multiprocessing import Pool
from functools import partial
import matplotlib.gridspec as gridspec
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
        with a single argument
        
        func -- function of the interface
           i -- index of the interface
        """
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
        
        
    def getSpectrum(self, wlList = None, align = True):
        """Calculate the spectrum for each point with given list of wavelengths"""
        if wlList == None:
            wlList = self.wlList
        result = []
        n = len(self.t)                
        for i in range(n):
            print('Calculating point ' + str(i+1) +  ' out of ' + str(n) + '...', flush = True)
            result.append(self.getResultForPoint(i, None, align, False))
        result = np.array(result)
        return result
        
    def setWlList(self,wlList):
        self.wlList = wlList
        
    def getResultForPoint(self, pointIndex, wlList = None, align = True, showProgress = True):
        """This method is to be called for getting the spectrum for a certain point""" 
        if type(wlList) == type(None):
            wlList = self.wlList
        _s = self.s
        _s.setStructure(self.tmp)
        _s.setThickness(self.t[pointIndex])
        ## Align each helix
        if align == True:
            _s.structures[1].phyParas['aor'] = s.structures[0].phyParas['t'] / \
            _s.structures[0].phyParas['p']
            _s.structures[2].phyParas['aor'] = s.structures[0].phyParas['t'] / \
            _s.structures[0].phyParas['p']
            """
            for i, helix in enumerate(self.s.structures):
                if i > 0: helix.phyParas['aor'] = end
                # Calculate the end angle, not the intrinstic anlge of rotation need to be added
                end = - helix.phyParas['t'] / helix.phyParas['p'] * np.pi + helix.phyParas['aor']
            """
        #print(self.s.getSubStructureInfo(), flush = True)
        result = self.s.scanSpectrum(wlList, coreNum = 1,giveInfo = False)[1]
        if showProgress == True:
            print('Calculation of point ' + str(pointIndex + 1) + ' finished', flush = True)
        return result
#%% Plotting    
def plotResult(wlList, result, title):
    """Plot the result from a series calculation of spectrum
    wlList is a 1D array of wavelengths.
    
    result: a NxW array where N is the number of calculations and W is the number
    of wavelengths calculated
    """
    r = result
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
    pl.ylabel("Distance /a.u.")
    pl.colorbar()
    return

def showLayerStructure(c):
    x = c.p.T
    # Calculate the ys to be plot
    y = np.append(c.h, np.repeat([[c.d]], len(c.p), axis = 0), axis = 1)
    pl.plot(y,x,'o-')
    pl.xlim(0,c.d+100)
    pl.ylim((c.l,0))
    pl.title('With pitch '+ str([x.phyParas['p'] for x in c.tmp])
    + " incident from right")
    pl.xlabel('Height from bottom /nm')
    pl.ylabel('Distance /a.u.')
        
#%%
def f1(x):
    return 4000 - 1.5*x
    
def f2(x):
    return 1000 + 1.5*x
    
def getSaveName():
    return "results\\" + time.strftime("%Y%m%d-%H%M%S")
    
    
if __name__ == '__main__':
    from report_Paras import figPath
    pitchesList1 = [[180,200,180]]
    nop = 200 #Number of points to sample along the defect
    for pitches in pitchesList1:
        wlRange = np.linspace(400,800,200)
        h1 = heli(CNC,pitches[0],1000) #The thickness doesn't matter here
        h2 = heli(CNC, pitches[1] ,1000)
        h3 = heli(CNC, pitches[2], 1000)
        tmp = [h1,h2,h3]
        #%% Set layer structure
        c = CrossSection(s, 5000,1000,3)
        c.setInterfaceFunction(f1,0)
        c.setInterfaceFunction(f2,1)
        c.calcPixelConfig(nop)
        c.setLayerTemplate(tmp)
        c.setWlList(wlRange)
        #for i in range(4):c.
       #     res.append(c.getResultForPoint(i))
        #c.getResultForPoint(0)
        t = clock()
        # %%Calculation
        with Pool(processes = 3) as pool:
            calc = partial(c.getResultForPoint, wlList = wlRange)
            res = pool.map(calc, list(range(nop)))
        # %% Plottign
        name = getSaveName()
        np.save(name, res)
        print('time elipsed',clock()-t)
        #%%Plotting the structure
        pl.rc('figure', figsize = (3.3,2.8))
        pl.rc('font', size = 8)
        pl.figure()
        x = c.p.T
        # Calculate the ys to be plot
        y = np.append(c.h, np.repeat([[c.d]], len(c.p), axis = 0), axis = 1)
        pl.plot(y, x,'.-')
        pl.plot(np.zeros(x.shape[0]), x, '.-')
        pl.xlim(0, c.d+2000)
        pl.ylim((c.l, 0))
        pl.annotate('', (5000,500), (7000,500),
                    arrowprops=dict(facecolor='black', headwidth = 10, width = 1,headlength = 10))
        pl.title('Pitch ='+ str([x.phyParas['p'] for x in c.tmp])
        + " Incident from right")
        pl.xlabel('Height from bottom /nm')
        pl.ylabel('Distance /a.u.')
        #pl.savefig(figPath+"LinearDefect3StructureWithMismatch.pdf")
        #%%Plotting Combined image £££££££
        fig = pl.figure()
        gs = gridspec.GridSpec(1, 2, width_ratios=[10,1])
        ax1 = pl.subplot(gs[0])
        ax2 = pl.subplot(gs[1])
        sPlot = ax1.imshow(res, cmap = 'viridis',aspect = 'auto', interpolation = 'none',
                  extent = [wlRange[0], wlRange[-1],  1000, 0])
        ax1.set_xlabel("Wavelength /nm")
        ax1.set_ylabel("Distance /a.u.")
        ax1.set_xticks(np.arange(400,801,100))
        pl.colorbar(sPlot, ax = ax1)
        resArray = np.array(res)
        spec = mt.specData(wlRange,resArray.T)
        RGB = spec.getRGBArray()
        ax2.imshow(RGB.reshape((nop,1,3)),aspect='auto', extent=[0,1,1000,0])
        ax2.set_title("RGB")
        ax2.set_xticks([])
        ax2.set_yticks([])
        pl.tight_layout(pad = 0)
        #pl.savefig(figPath+ "LinearDefect3SpectrumWithMismatch.pdf")
        ####