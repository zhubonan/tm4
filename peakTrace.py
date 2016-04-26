# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 21:05:22 2016
Ojective: Trace the maximum of the calculated spectrum with changing pitch and 
compare with the theoratical values
@author: Bonan
"""
import numpy as np
import matplotlib.pyplot as pl
from multiprocessing import Pool
from matTools import specData
from preset import s
pl.rcParams['savefig.dpi'] = 200
wlRange = np.linspace(400,800,100)

def calcPitch(pitch):
    """Calculate reflectivity at a certain pitch"""
    s.setPitch([pitch])
    res = specData(*s.scanSpectrum(wlRange)[:2])
    # Debug message
    print("{0:.0f}".format(pitch) + " finished", flush = True)
    return (pitch,) + res.getPeaks()

def plotData(filename, size = [4, 3]):
    pl.figure(figsize = size)
    pl.plot(peakArray[0], peakArray[1],'-o')
    pl.xlabel("Pitch /nm")
    pl.ylabel("Wavelength /nm")
    pl.title("Peak wavelength against pitch")
    pl.savefig("./reportFigs/peakVsPitch")
    pl.show()
    
if __name__ == "__main__":
    pool = Pool(processes = 4)
    s.setThickness([2000])
    pList = np.linspace(150,250,10)
    peakList = pool.map(calcPitch, pList) #Peak list in order or [wl,ref]
    pool.close()
    pool.join()
    peakArray = np.array(peakList).swapaxes(0,1)
    np.save("./reportData/peakVsPitch.npy", peakArray)
    #%%Plotting
    plotData("./reportData/peakVsPitch.npy", [4.5,4.5])