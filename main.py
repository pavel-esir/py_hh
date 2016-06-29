# coding: utf-8
'''
Created on 29 июня 2016 г.

@author: Pavel Esir
'''

from py_hh import *
import numpy as np
import matplotlib.pylab as pl

SimTime = 1000.
h = 0.02
Tsim = int(SimTime//h)
tframes = np.linspace(0, SimTime, Tsim)

Nneur = 1
tau_psc = 0.2
w_p = 1.
rate = 100

setCalcParams(Tsim, Nneur, 1, h)

Vm = np.zeros(Nneur) + 32.906693
Vrec = np.zeros((Tsim, Nneur))
ns = np.zeros(Nneur) + 0.574676
ms = np.zeros(Nneur) + 0.913177
hs = np.zeros(Nneur) + 0.223994

Ie = np.zeros(Nneur) + 5.27
y = np.zeros(Nneur)
Isyn = np.zeros(Nneur)
d_w_p = np.zeros(Nneur) + np.e*w_p/tau_psc

setNeurVars(Vm, Vrec, ns, ms, hs)
setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, 0)

simulate()

pl.figure()
pl.plot(tframes, Vrec[:, 0])
pl.xlabel('time, ms')
pl.ylabel('Membrane potential, mV')
pl.show()

