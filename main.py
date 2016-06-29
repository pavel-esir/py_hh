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

Nneur = 2
tau_psc = 0.2
w_p = 1.
rate = 100

setCalcParams(Tsim, Nneur, 1, 1, h)

Vm = np.zeros(Nneur) + 32.906693
Vrec = np.zeros((Tsim, Nneur))
ns = np.zeros(Nneur) + 0.574676
ms = np.zeros(Nneur) + 0.913177
hs = np.zeros(Nneur) + 0.223994

spike_times = np.zeros((int(SimTime/10), Nneur), dtype='uint32')
weight = np.array([1.])
delay = np.array([0], dtype='uint32')
pre = np.array([0], dtype='uint32')
post = np.array([1], dtype='uint32')

Ie = np.zeros(Nneur)
Ie[0] = 5.27
y = np.zeros(Nneur)
Isyn = np.zeros(Nneur)
d_w_p = np.zeros(Nneur) + np.e*w_p/tau_psc

setNeurVars(Vm, Vrec, ns, ms, hs)
setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, 0)

setSpikeTimes(spike_times)

setConns(weight, delay,  pre, post)

simulate()

pl.figure()
pl.plot(tframes, Vrec[:, 0])
pl.xlabel('time, ms')
pl.ylabel('Membrane potential, mV')
pl.show()

