#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 21:36:54 2017
Trys to detect blue/red shift
@author: bonan
"""
import scipy as sp
import numpy as np
import tm4.matTools as mt

class GradientArray:
    sample_range = 1
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
        n = self.sample_range * 2 + 1
        r = self.sample_range
        store = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == r and j == r:
                    store[r,r] = 0
                    continue
                store[i,j] = self.compute_shift(self.spec[row, col], 
                            self.spec[row+i-r, col+j-r])
        return store
                
    def construct_data(self, sdata):
        self.spec =   sdata.spec.reshape((sdata.spec.shape[0], *sdata.spacialShape))
        self.spec = np.moveaxis(self.spec, 0, -1) # change orde to be row, col, value
        self.wlrange = sdata.wl

    def calc(self):
        n = self.sample_range * 2 + 1
        r = self.sample_range
        shift_store = np.zeros((*self.spec.shape[:2], n ,n))
        gradient_row = np.zeros(self.spec.shape[:2])
        gradient_col = np.zeros(self.spec.shape[:2])
        for row in range(r, self.spec.shape[0] - r):
            for col in range(r, self.spec.shape[1] - r):
                shift_matrix = self.calc_shift(row, col)
                shift_store[row, col] = shift_matrix
                grad = np.gradient(shift_matrix, r)
                gradient_row[row, col] = grad[0][r,r]
                gradient_col[row, col] = grad[1][r,r]
        self.shift_array = shift_store
        self.gradient_arrays = (gradient_row, gradient_col)
        return shift_store