"""
Python script for plotting and simulate 1-D defect structure. e.g. get sepctrum 
along a line
"""
import simClasses as sim
import matplotlib.pyplot as pl
import numpy as np
import copy as cp
from colourTools import specToRGB
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
h1 = heli(CNC,150,1000)
s.setStructure([h1])
wlRange = np.linspace(350,850,100)
print('Followings are added to the scope')
print('Materials: CNC, air, cellulosem, glass')
print('HalfSpace: airhalf, glasshalf')
print('OptSystem:s')
print('heli as HeliCoidalStructure')
print('wlRange as 400 to 800 nm')
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
        struc = cp.deepcopy(templateList[i])
        if type(struc) == sim.HeliCoidalStructure:
            struc.setPitch(configArray[i][0])
        struc.setThickness(configArray[i][1])
        slist.append(struc)
    optSys.setStructure(slist)
    return
    
#if __name__ == '__main__':
config = [[150,1500],[200,600]]
tmp = [h1,h1]
applyConfig(s,config,tmp)
plotSpectrum(s,wlRange)