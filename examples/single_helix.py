"""

Example script for computing the reflectance spectrum for the following setup:

==============================================
      |
      |  air - refractive index = 1
      V
----------------------------------------------

    Cellulose nanocrystals with:
    pitch = 180 nm, refractive indices=1.586, 1.524
    thickness = 1000 nm

----------------------------------------------

    glass

==============================================

"""

import matplotlib.pyplot as pl
import numpy as np
import tm4.simClasses as sim

pl.rcParams['figure.figsize'] = (8, 6)
pl.rcParams['savefig.dpi'] = 100

# define materials
CNC = sim.UniaxialMaterial(1.586, 1.524)  # The cellulose nanocrystal has birefringence 1.586 and 1.524 

# Define half spaces
air = sim.HomogeneousNondispersiveMaterial(1)  # Optical index is one for air
glass = sim.HomogeneousNondispersiveMaterial(1.55)   # Assume the substrate a glass with refractive index 1.55
airhalf = sim.IsotropicHalfSpace(air)
glasshalf = sim.IsotropicHalfSpace(glass)

# Initialise and model system
s = sim.OptSystem()
s.setHalfSpaces(airhalf, glasshalf)  # Assign half-spaces

# Define the helicoidal layers
heli = sim.HeliCoidalStructure
h1 = heli(CNC, 180, 1000)  # Define a CNC layer of pitch 180 (full pitch 360) and thickness 1000
s.setStructure([h1])  # The system has a structure which is our helix

# Define the wavelenghs for scan
wlRange = np.linspace(400, 800, 100)  # Commonly used spectral range

fig, axs = pl.subplots(5, 1, sharex=True, sharey=True)
for coupling, ax in zip(['LL', 'RR', 'LR', 'RL', 'all'], axs):
    wl, reflect = s.scanSpectrum(wlRange, coupling=coupling)
    ax.plot(wl, reflect, label=coupling)
    ax.legend()
    np.savetxt("single-helix-data-{}.dat".format(coupling), np.stack([wl, reflect], axis=1))
    
ax.set_ylabel("Reflectance")
ax.set_xlabel("Wavelength (nm)")
axs[0].set_title("Reflectance of a single CNC helix")
fig.savefig("single_helix.png", dpi=150)