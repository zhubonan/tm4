# -*- coding: utf-8 -*-
"""
Created on Fri May  6 12:24:39 2016

@author: Bonan
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May  6 10:07:21 2016
For spectral Imaging data
@author: Bonan
"""
import matplotlib.pyplot as pl
import numpy as np
import scipy as sp
import matTools as mt
from report_Paras import figPath
pl.rc('figure', figsize = (3.3,2.8))
pl.rc('font', size = 8)
#Define file path
scan0502slim = r'C:\Users\Bonan\OneDrive\Documents\PartIII\Project\20160502\scanslim.mat'
scan0427slim = r'C:\Users\Bonan\OneDrive\Documents\PartIII\Project\20160427\scanslim.mat'
scan0408slim = r'C:\Users\Bonan\OneDrive\Documents\PartIII\Project\20160408\scanslim.mat'
# %%Load Data
test = mt.scanData(scan0427slim)
#%% Preview
stemp = test.getCurrentSpecData(5)
stemp.cropSelf([450,750])
stemp.setSpacialShape((51,51)) ## Define the shape of data
stemp.testFilter(21,2)
stempcf = stemp.getFilteredSpectrum(21,2)
stempcf.spec = stempcf.spec #Minus a backgroud
stempcf.get2DColouredImage('auto',1)
#%% Plot sliced data
start , finish = (0,0),(50,50) #Deine positions
num = 80 # Number of interpolation
dist = np.linalg.norm(np.array(finish)-np.array(start))
zi = stempcf.plotLineSlice(start, finish, num, show = 0)
RGBImage = stempcf.get2DColouredImage()
pl.figure(figsize = (3.3,5))
pl.imshow(RGBImage)
arrow = dict(facecolor='black', headwidth = 10, width = 2,headlength = 10)
pl.annotate('', (finish[0],finish[1]), (start[0],start[1]), arrowprops = arrow)
pl.xlabel(r'$\mathsf{\mu m} $', fontsize = 9)
pl.ylabel(r'$\mathsf{\mu m} $', fontsize = 9)
pl.title('RGB Image')
pl.tight_layout(pad = 0.5)
#pl.savefig(figPath + 'SpectroImagingRGB.pdf')
pl.figure(figsize = (3.3,3))
pl.imshow(zi.T, aspect='auto', interpolation='none', extent = [stempcf.range[0],stempcf.range[1],dist,0])
pl.ylabel(r'$\mathsf{\mu m} $', fontsize = 9)
pl.xlabel('Wavelength nm')
pl.colorbar()
pl.title('Spectrum along the slice')
pl.tight_layout(pad = 0.5)
#pl.savefig(figPath + 'SpectroImagingSlice.pdf')
#%% Plot some stats data
peakData = stempcf.getPeaks()
pl.figure(figsize = (3.3,2.8))
pl.plot(peakData[0],peakData[1],'.')
pl.xlabel('Wavelength /nm')
pl.ylabel('Reflectance')
pl.plot([518,590],[0.90,0.55],'k--', linewidth = 2)
pl.title('Height vs Position')
arrow2 = dict(facecolor='black', headwidth = 5, width = 2,headlength = 10)
pl.annotate('Cut-off', [540,0.8],[560,0.9],arrowprops = arrow2)
pl.tight_layout(pad = 0.5)
#pl.savefig(figPath + 'SpectroImagingStats.pdf')
#%%Plot an example spectrum
import simClasses as sim
from preset import s,wlRange, CNC
nBar = 1.55
dn = 0.025
m = sim.UniaxialMaterial(nBar + dn, nBar-dn)
#m = CNC
s.structures[0].phyParas['m'] = m
pl.figure(figsize = (3.3,2.8))
pl.plot(stempcf.wl, (stempcf.spec[:, 15*51+15])/2-0.025, label = 'Measured')
s.setPitch([179.5])
s.setThickness([2000])
pl.plot(*s.scanSpectrum(wlRange), label = 'Fit-Best')
s.structures[0].phyParas['m'] = CNC
pl.plot(*s.scanSpectrum(wlRange), label = 'Fit-Default')
pl.xlim(450,750)
pl.ylim(0,0.35)
pl.xlabel('Wavelength /nm')
pl.ylabel('L-L Reflectance')
pl.title('Measured Spectrum vs Fit')
pl.legend()
pl.tight_layout(pad = 0.5)
#pl.savefig(figPath + 'SpectroSingleFit.pdf')
