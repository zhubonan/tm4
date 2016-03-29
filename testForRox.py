# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:55:12 2016

@author: Bonan
"""

from preset import *
s.Theta = 1
res = s.scanSpectrum(wlRange)
pl.plot(res[0],res[1])