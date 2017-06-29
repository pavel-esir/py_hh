# coding: utf-8
"""
Created on Thu Oct 13 23:02:27 2016

@author: Pavel Esir
"""
from __future__ import division
import matplotlib.pylab as pl
import numpy as np
from numpy import random

def getDelays(N, seed=0, x0=20.0):
    random.seed(seed)

    x1 = x0*random.rand(N)
    x2 = x0*random.rand(N)
    d = (x1 - x2)**2

    x1 = x0*random.rand(N)
    x2 = x0*random.rand(N)
    d += (x1 - x2)**2

    x1 = x0*random.rand(N)
    x2 = x0*random.rand(N)
    d += (x1 - x2)**2
    return np.sqrt(d)
if __name__ == "__main__":
    x0 = 10
    d = getDelays(2000, x0=x0)
    pl.figure(1)
    pl.hist(d, bins=100, range=(0, np.sqrt(3)*x0), normed=True, histtype='step')
    
    pl.show()
