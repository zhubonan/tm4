"""
Convenience module can be used for setting up a simulation quickly
"""

import matplotlib.pyplot as pl
import numpy as np

import tm4.simClasses as sim
from tm4.colourTools import specToRGB

pl.rcParams['figure.figsize'] = (8, 6)
pl.rcParams['savefig.dpi'] = 100

# define materials
CNC = sim.UniaxialMaterial(1.586, 1.524)  # this might not be right
air = sim.HomogeneousNondispersiveMaterial(1)
cellulose = sim.HomogeneousNondispersiveMaterial(1.55)

# Define half spaces
glass = sim.HomogeneousNondispersiveMaterial(1.55)
front = sim.IsotropicHalfSpace(air)
airhalf = sim.IsotropicHalfSpace(air)
glasshalf = sim.IsotropicHalfSpace(glass)

# Initialise and model system
s = sim.OptSystem()
s.setHalfSpaces(airhalf, glasshalf)  # Assign half-spaces

# Define the helicoidal layers
heli = sim.HeliCoidalStructure
h1 = heli(CNC, 180, 1000)  # Define a CNC layer of pitch 180 and thickness 1000
s.setStructure([h1])
wlRange = np.linspace(400, 800, 100)  # Commonly used spectral range
