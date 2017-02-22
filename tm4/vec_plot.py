#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:36:54 2017
Trys to detect blue/red shift
@author: bonan
"""
import scipy as sp
import matplotlib.pyplot as pl
import numpy as np
import tm4.matTools as mt
pl.rc('figure', figsize = (3.3,2.8))
pl.rc('font', size = 8)
#Define file path
scanPath = 'DataFoler/20160427/scanslim.mat'
# %%Load Data
test = mt.scanData(scanPath)
#%% Preview
stemp = test.getCurrentSpecData(5)
stemp.cropSelf([450,750])
stemp.setSpacialShape((51,51)) ## Define the shape of data
stemp.testFilter(21,2)
stempcf = stemp.getFilteredSpectrum(21,2)
stempcf.spec = stempcf.spec #Minus a backgroud
stempcf.get2DColouredImage('auto',1)
#%% Try to extract information about blue/red shift
scrop = stempcf.crop((525,650))  # Crop, choose only the interested range

class GradientArray:
    
    def __init__(self, sdata, start, finish, npoints):
        self.construct_data(sdata)
        self.spec_start = start
        self.spec_finish = finish
        self.npoints = npoints
        
    @staticmethod
    def metric(shift, spec1, spec2, wlrange, spec_start, spec_finish, npoints):
        """Function to be minised
        This is the sum of different between shifted spectrum given a shift"""
        points = np.linspace(spec_start, spec_finish, npoints)
        data = np.interp(points + shift, wlrange, spec1)
        ref = np.interp(points, wlrange, spec2)
        metric = np.sum(np.abs(data - ref))
        return metric
            
            
    def compute_shift(self, spec1, spec2):
        """Compute the shift of spec1 repect to spec2.
        spec1, spec2 : 1D array
        Return shift in nm. Plus for blue shift, minus for red shift"""
    
        res = sp.optimize.minimize_scalar(self.metric, method='Golden',
                                          args= (spec1, spec2, self.wlrange,
                                                 self.spec_start, 
                                                 self.spec_finish,
                                                 self.npoints))
        return res.x
        
    def calc_shift(self, row, col):
        """Calcaulte spectral shift of adject points"""
        store = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                if i == 1 and j == 1:
                    store[1,1] = 0
                    continue
                store[i,j] = self.compute_shift(self.spec[row, col], 
                            self.spec[row+i-1, col+j-1])
        return store
                
    def construct_data(self, sdata):
        self.spec =   sdata.spec.reshape((sdata.spec.shape[0], *sdata.spacialShape))
        self.spec = np.moveaxis(self.spec, 0, -1) # change orde to be row, col, value
        self.wlrange = sdata.wl

    def calc(self):
        shift_store = np.zeros((*stempcf.spacialShape, 3 ,3))
        for row in range(1, stempcf.spacialShape[0] - 1):
            for col in range(1, stempcf.spacialShape[1] - 1):
                shift_matrix = self.calc_shift(row, col)
                shift_store[row, col] = shift_matrix
        return shift_store


#%% Test
g = GradientArray(stempcf, 550, 650 , 100)
g.calc_shift(25,25)
g.calc()
