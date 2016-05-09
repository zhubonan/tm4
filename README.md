# TM4
4-by-4 transfer matrix method for birefringent multilayer filters

This is a work in progress to calculate the optical response for layered anistropic(birefringent materials)

simClasses -- Definition of classes for simulation
reportFigs -- For plotting figures shown in reportFigsmat
matTools   -- Spectral imaging data processing and data loading
colourTools -- Colour calculation using CIE standard
linearDefect -- Calculation of spectrum of defect structures.
preset -- Import this for fast set up

Quick start:

#This will generate a spectrum using default setup
from preset import *
pl.plot(*s.scanSpectrum(wlRange,1))