# TM4

Introduction
===============

A python package for using 4-by-4 transfer matrix methods to simulate birefringent multilayer structures.

The main topic of this research is to study the reflectance spectrum of helicoidal structures which were found in cellulose nanocrystals.

Python support: python 3.5+

Usage
======

Set set your `PYTHONPATH` environmental variable so `tm4` can be imported

Modules
--------

simClasses -- Definition of classes for simulation
matTools   -- Spectral imaging data processing and data loading
colourTools -- Colour visulisation from spectrum
linearDefect -- Calculation of spectrum of defect structures.
preset -- Convenience module for fast setup of simulation conditions

Quick start
-------------

This will generate a spectrum using default setup

from tm4.preset import *
pl.plot(*s.scanSpectrum(wlRange,1))
pl.show()
