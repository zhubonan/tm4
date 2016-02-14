# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 18:02:50 2016
matlab for loading data from spectrometer
@author: Bonan
"""
from scipy.io import loadmat
import numpy as np

class scanData:
    """A class for easy control of scan.mat data"""
    def __init__(self, filename):
        self.data = loadmat(filename, squeeze_me = True, appendmat = True )['scan']
        
    def listAllDescription(self):
        for i in self.data['desc']:
            print(i)
            
    def getSpectrumViaIndex(self, index):
        
        spec = self.data['spec'][index]
        wl = self.data['wl'][index]
        return wl, spec
        
    def getSpectrumViaKeyWrod(self, keyword):
        result = []
        for i in range(len(self.data['desc'])):
            if self.data['desc'][i].find(keyword) != -1:
                result.append(self.getSpectrumViaIndex(i))
        return np.array(result)
                
    def getDescViaIndex(self, index):
        return self.data['desc'][index]