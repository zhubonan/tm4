# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:20:56 2016

@author: Bonan
"""

from cProfile import Profile
from pstats import Stats
prof = Profile()
prof.disable()
# Import the preset module 
from tm4.preset import s, wlRange
prof.enable()
#%% Code to be profiled
s.scanSpectrum(wlRange,1)


#%%
prof.disable()  # don't profile the generation of stats
prof.dump_stats('mystats.stats')
with open('mystats_output.txt', 'wt') as output:
    stats = Stats('mystats.stats', stream=output)
    stats.sort_stats('cumulative', 'time')
    stats.print_stats()