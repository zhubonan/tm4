# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 10:40:39 2016
Converting reflectance spectrum to a CIE coordinate
@author: Bonan
"""
import numpy as np
from scipy import interpolate
import os
# Adobe RGB (1998) D65 as reference white
# http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_RGB.html
_RGB_to_XYZ = np.array([
    [0.5767309, 0.1855540, 0.1881852],
    [0.2973769, 0.6273491, 0.0752741],
    [0.0270343, 0.0706872, 0.9911085], ])
_XYZ_to_RGB = np.array([
    [2.0413690, -0.5649464, -0.3446944],
    [-0.9692660, 1.8760108, 0.0415560],
    [0.0134474, -0.1183897, 1.0154096], ])
# Load the
_dirname = os.path.dirname(__file__)
fn1 = os.path.join(_dirname, 'CIE_1931_XYZ.txt')
fn2 = os.path.join(_dirname, 'CIE_A.txt')
fn3 = os.path.join(_dirname, 'CIE_D65.txt')
CIE_XYZ_table = np.loadtxt(fn1).T  # Transpose column into rows
CIE_A = np.loadtxt(fn2).T
CIE_D65 = np.loadtxt(fn3).T


def splineInterp(xNew, xRaw, yRaw):
    """
    Compute the spline interpolation(cubic) of the data
    """
    tck = interpolate.splrep(xRaw, yRaw)
    return interpolate.splev(xNew, tck, der=0, ext=1)


def specToXYZ(spec, SI='D65'):
    """
    Calculate the XYZ coordinate of the spectrum input.
    It interpolates the charts to every wavelength that was inputed.
    By default the input spectrum was first eveloped using a SPD function
    to simulation illumination.

    spec: input spectrum, 2*N ndarray, 1st row must be the wavelength

    return: (X,Y,Z)
    """
    wl = spec[0]  # the input must have the 1st element as the wavelength
    XYZ = CIE_XYZ_table
    if SI == 'D65':
        interpSI = splineInterp(wl, CIE_D65[0], CIE_D65[1])
    if SI == 'A':
        interpSI = splineInterp(wl, CIE_A[0], CIE_A[1])
    else:
        interpSI = np.ones(len(wl))
    interpX = splineInterp(wl, XYZ[0], XYZ[1])
    interpY = splineInterp(wl, XYZ[0], XYZ[2])
    interpZ = splineInterp(wl, XYZ[0], XYZ[3])
    interpXYZ = np.array([interpX, interpY, interpZ])
    X, Y, Z = np.sum(spec[1] * interpSI * interpXYZ, axis=1)
    return X, Y, Z


def specToxyz(spec, SI='D65'):
    """
    Transfer spectrum into normalised x,y,z coordinates

    Return:  (x, y, z)
    """
    X, Y, Z = specToXYZ(spec, SI)
    x = X / (X + Y + Z)
    y = Y / (X + Y + Z)
    z = 1 - x - y
    return x, y, z


def specToRGB(spec, SI='D65', scale_factor=1):
    """
    Convert the spectrum(reflectivity) into an RGB value

    Return: (R,G,B)
    """
    XYZArray = specToxyz(spec, SI)
    RGBArray = np.dot(_XYZ_to_RGB, XYZArray).clip(0, 1)
    RGBArray *= scale_factor
    return tuple(RGBArray.clip(0, 1))


if __name__ == '__main__':
    # Testing of the module
    import matplotlib.pyplot as pl
    wlRange = np.linspace(400, 800, 100)
    example = np.sin((wlRange - 400) * np.pi / 400)
    spec = np.array([wlRange, example])
    c = specToRGB(spec)
    pl.plot(spec[0], spec[1] / spec[1].max(),
            label='Example distribution', color=c)
    print(c)
    # Use the D65 as the light source
    spec = CIE_D65
    c = specToRGB(spec, SI='D65')
    print('Test using D65 illumination. Should give R=G=B')
    print(c)
    pl.plot(spec[0], spec[1] / spec[1].max(),
            label='D65 distribution', color=np.array(c))
    pl.title('Coloured Spectrum')
    pl.legend()
