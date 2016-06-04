# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 17:14:22 2016

@author: Bonan

Python script for plotting and simulate 1-D defect structure. e.g. get sepctrum 
along a line
"""
import matplotlib.pyplot as pl
import numpy as np
import matplotlib.gridspec as gridspec
from time import clock
import tm4.simClasses as sim
import tm4.matTools as mt
from tm4.colourTools import specToRGB
from tm4.linearDefect import CrossSection
pi = np.pi
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
    
class CrossSection3L(CrossSection):
    """A class represent a CrossSection of the film"""

        
    def alignHelices(self, pointIndex):
        """Align the helices"""
        s = self.s.structures
        s[1].phyParas['aor'] = 2 * pi * s[0].phyParas['t'] / \
        s[0].phyParas['p'] + s[0].phyParas['aor'] + 5 * pi * pointIndex / len(self.t)
        #s[2].phyParas['aor'] = 2 *pi * s[0].phyParas['t'] / \
        #s[0].phyParas['p'] + s[0].phyParas['aor']
  
        
#%%
def f1(x):
    return 3300 #- .2*x
    
def f2(x):
    return 1500 # +.2*x
    
    
if __name__ == '__main__':
    pitchesList1 = [[180,180,180]]
    nop = 100 #Number of points to sample along the defect
    for pitches in pitchesList1:
        wlRange = np.linspace(400,800,200)
        h1 = heli(CNC,pitches[0],1000) #The thickness doesn't matter here
        h2 = heli(CNC, pitches[1] ,1000)
        h3 = heli(CNC, pitches[2], 1000)
        tmp = [h1,h2,h3]
        #%% Set layer structure
        c = CrossSection3L(s, 5000,1000,3)
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
        res = c.getSpectrum(wlRange, 4)
        # Performance evaluation
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
        #pl.title('Pitch ='+ str([x.phyParas['p'] for x in c.tmp])
        #+ " Incident from right")
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