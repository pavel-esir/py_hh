# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 19:35:53 2017

@author: pavel
"""
import numpy as np
from numpy import random
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def getConns(Nneur, Ncon, seed=0, x0=20.0):
    random.seed(seed)
    x1 = x0*random.rand(Nneur)
    x2 = x0*random.rand(Nneur)
    x3 = x0*random.rand(Nneur)

    pre = random.randint(0, Nneur, Ncon).astype('uint32')
    post = random.randint(0, Nneur, Ncon).astype('uint32')

    preX1 = x1[pre]
    postX1 = x1[post]
    preX2 = x2[pre]
    postX2 = x2[post]
    preX3 = x3[pre]
    postX3 = x3[post]

    delay = np.sqrt((preX1 - postX1)**2 + (preX2 - postX2)**2 + (preX3 - postX3)**2)
    return pre, post, delay, (x1, x2, x3)

def plotConns(coord, pre, post):
    x1, x2, x3 = coord
    preX1 = x1[pre]
    postX1 = x1[post]
    preX2 = x2[pre]
    postX2 = x2[post]
    preX3 = x3[pre]
    postX3 = x3[post]

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection='3d')
    for i in xrange(len(pre)):
        ax.plot([preX1[i], postX1[i]], [preX2[i], postX2[i]], [preX3[i], postX3[i]], '-o', lw=0.1, color='k', ms=2., mec='b')    
    ax.set_xticks((0, 10, 20))
    ax.set_yticks((0, 10, 20))
    ax.set_zticks((0, 10, 20))
    
if __name__ == "__main__":
    pre, post, delay, coord = getConns(100, 2000, 0)
    
    plotConns(coord, pre, post)
