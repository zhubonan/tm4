# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
from multiprocessing import Pool
from functools import partial


class CrossSection():
    """A base class represent a CrossSection of the film"""

    def __init__(self, optSys, depth, length, nol):
        """Initialise a CrossSection object"""
        self.d = depth  # depth, we assume the flim is flat
        self.l = length  # Length of the cross-section
        self.n = nol  # Number of layers
        self.s = optSys  # The intantiated OptSystem object
        self.interface = [0] * (nol - 1)  # Preallocate the
        self.wlList = None

    def setLayerTemplate(self, template):
        """Set up the type of layers using the a list of existing Strucuture objects.
        The length of list must match the number of layer.
        The Structure need to be instantiated (e.g. cannot directly pass
         a class)
        """
        if len(template) != self.n:
            raise RuntimeError(
                "The tempate must match the total number of layers")
        self.tmp = template

    def setInterfaceFunction(self, func, i):
        """Set the function of the ith interface.
        The function must be callable with a single argument"""
        self.interface[i] = func

    def calcPixelConfig(self, div):
        """This will calcuate the thickness of the layers under each pixel

        d: divisions of the cross-section
        """
        self.p = np.linspace(0, self.l, div)  # location of the points
        # Calculate the location of interfaces
        if (np.array(self.interface) == 0).any():
            raise RuntimeError('At least one of the intface is undefined')

        # array of the height of each interface
        h = np.empty((div, self.n - 1))
        for i in range(len(self.interface)):
            vec = np.vectorize(self.interface[i])
            h[:, i] = vec(self.p)
        # Now we need to convert the interface height to layer thickness
        self.h = h  # The array of height of the interfaces
        h1 = np.append(h, np.zeros((div, 1)), axis=1)
        h2 = np.append(np.repeat([[self.d]], div, axis=0), h, axis=1)
        # t is the array of thickness for each layer a all points
        self.t = h2 - h1
        return self.t

    def getSpectrum(self, wlList=None, process=1, align=True, coupling='LL'):
        """
        Calculate the spectrum for each point with given list of wavelengths"""
        if wlList is None:
            wlList = self.wlList
        calc = partial(self.getResultForPoint, wlList=wlList, coupling=coupling)
        # If we only use one process then no need to use Pool class
        if process == 1:
            result = np.array(list(map(calc, range(len(self.p)))))
        else:
            with Pool(processes=process) as pool:
                res = pool.map(calc, list(range(len(self.p))))
                result = np.array(res)
        self.result = result
        return result

    def setWlList(self, wlList):
        self.wlList = wlList

    def alignHelices(self, pointIndex):
        """
        Align the helices

        pointIndex: Index of the point to be switched to
        """
        raise NotImplementedError('alignment need to be implement in subclass')

    def getResultForPoint(self, pointIndex, wlList=None,
                          align=True, showProgress=True, coupling='LL'):
        """
        This method is to be called for getting the spectrum for a
        certain point
        """
        if wlList is None:
            wlList = self.wlList
        self.s.setStructure(self.tmp)
        self.s.setThickness(self.t[pointIndex])
        # Align each helix
        if align is True:
            self.alignHelices(pointIndex)
        # print(self.s.getSubStructureInfo(), flush = True)
        result = self.s.scanSpectrum(wlList, coreNum=1, giveInfo=False, coupling=coupling)[1]
        if showProgress is True:
            print('Calculation of point ' +
                  str(pointIndex + 1) + ' finished', flush=True)
        return result

    def showLayerStructure(self):
        x = self.p.T
        # Calculate the ys to be plot
        y = np.append(self.h, np.repeat(
            [[self.d]], len(self.p), axis=0), axis=1)
        pl.plot(y, x, 'o-')
        pl.xlim(0, self.d + 100)
        pl.ylim((self.l, 0))
        pl.title('With pitch ' + str([x.phyParas['p'] for x in self.tmp])
                 + " incident from right")
        pl.xlabel('Height from bottom /nm')
        pl.ylabel('Distance /a.u.')
