# coding: utf-8
'''
Created on 29 июня 2016 г.

@author: Pavel Esir
'''

from py_hh import *
import numpy as np
import matplotlib.pylab as pl
np.random.seed(0)

SimTime = 1000.
h = 0.04
Tsim = int(SimTime/h)

Nneur = 100
tau_psc = 0.2
w_p = 1.9
rate = 100
Ncon = 2000
recInt = 10
setCalcParams(Tsim, Nneur, Ncon, recInt, h)

Vm = np.zeros(Nneur) + 32.906693
Vrec = np.zeros((int(Tsim/recInt), Nneur))
ns = np.zeros(Nneur) + 0.574676
ms = np.zeros(Nneur) + 0.913177
hs = np.zeros(Nneur) + 0.223994

spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
weight = np.zeros(Ncon) + ([1.3*np.e/tau_psc])
#delay = np.array([0./h], dtype='uint32')
#pre = np.array([0], dtype='uint32')
#post = np.array([1], dtype='uint32')
pre = np.random.randint(0, Nneur, Ncon).astype('uint32')
post = np.random.randint(0, Nneur, Ncon).astype('uint32')
delay = np.random.lognormal(1.8, 1/6.0, Ncon).astype('uint32')
#delay = np.zeros(Ncon, dtype='uint32')

Ie = np.zeros(Nneur) + 5.27
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
pl.plot(np.linspace(0, SimTime, int(Tsim/recInt)), Vrec[:, ::3])
pl.xlabel('time, ms')
pl.ylabel('Membrane potential, mV')
pl.show()

