# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:41:09 2016
Add randomness to the helix
@author: Bonan
"""
from preset import s,  CNC, wlRange
import numpy as np
import simClasses as sim
import matplotlib.pyplot as pl
testMaterial = sim.UniaxialMaterial(1.6,1.5)
hrandom = sim.HelixCustom(p = 180, t = 1800, d = 90, aor = 0, m = CNC, handness = -1)
# Inject random variations
def calc(dphi):
    repeat = 10
    res = []
    for i in range(repeat):
        hrandom.calcStandardAngles()
        randomArray = np.random.standard_normal(hrandom.angleList.size) * dphi/180 * np.pi
        hrandom.angleList = hrandom.angleList + randomArray
        s.setStructure([hrandom])
        res.append(s.scanSpectrum(wlRange,1)[1])
    return res

#%%
if __name__ == '__main__':
    dphiList= np.linspace(0,30,5)
    from multiprocessing import Pool
    with Pool(processes = 4) as pool:
        res = np.array(pool.map(calc, dphiList))
    #%%Processing the data
    avg = np.average(res,1)
    std = np.std(res,1)
    #%%Plottings
    pl.figure(figsize = (3.3,3))
    pl.rc('font', size = 8)
    for dphi, _avg, _std in zip(dphiList, avg, std):
        pl.errorbar(wlRange, _avg, _std, label = r'$\delta \theta = $' + str(dphi))
    pl.xlabel('Wavelength /nm')
    pl.ylabel('Reflectance (LL)')
    pl.legend()
    pl.xlim(500,650)
    pl.ylim(0,0.3)    
    pl.title('Randomised helix')
    pl.tight_layout(pad = 0.5)    
    from report_Paras import figPath
    #pl.savefig(figPath + 'RandomHelix.pdf')
    