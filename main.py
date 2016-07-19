# coding: utf-8
'''
Created on 29 июня 2016 г.

@author: Pavel Esir
'''

from py_hh import *
import numpy as np
import matplotlib.pylab as pl
np.random.seed(0)

SimTime = 100.0
h = 0.02
Tsim = int(SimTime/h)
recInt = 2

tau_psc = 0.2  # ms
w_p = 2.0      # Poisson noise weith, pA
w_n = 2.0      # connection weight, pA

rate = 200     # Poisson noise rate, Hz (shouldn't  be 0)
Nneur = 100
Ncon = int(Nneur*Nneur*0.2)
pre = np.random.randint(0, Nneur, Ncon).astype('uint32')
post = np.random.randint(0, Nneur, Ncon).astype('uint32')
delay = np.random.lognormal(1.8, 1/6.0, Ncon).astype('uint32') # delay arrays in time frames

#Nneur = 2
#Ncon = 1
#pre = np.array([0], dtype='uint32')
#post = np.array([1], dtype='uint32')
#delay = np.array([10./h], dtype='uint32') # delay arrays in time frames
#rate = 0.1     # Poisson noise rate, Hz (shouldn't  be 0)

# int(SimTime/10) тут 10 это период в мс максимального ожидаемого интервала между спайками
spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

weight = np.zeros(Ncon) + w_n*np.e/tau_psc

setCalcParams(Tsim, Nneur, Ncon, recInt, h)

Vm = np.zeros(Nneur) + 32.906693
Vrec = np.zeros((int(Tsim/recInt), Nneur))
ns = np.zeros(Nneur) + 0.574676
ms = np.zeros(Nneur) + 0.913177
hs = np.zeros(Nneur) + 0.223994


#Ie = np.zeros(Nneur) + 5.27

Ie = np.zeros(Nneur)
Ie[0] = 5.27

y = np.zeros(Nneur)
Isyn = np.zeros(Nneur)
d_w_p = np.zeros(Nneur) + np.e*w_p/tau_psc

setNeurVars(Vm, Vrec, ns, ms, hs)
setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, 0)

setSpikeTimes(spike_times, num_spikes_neur)

setConns(weight, delay,  pre, post)

simulate()

pl.figure()
pl.plot(np.linspace(0, SimTime, int(Tsim/recInt)), Vrec[:, :])
pl.xlabel('time, ms')
pl.ylabel('Membrane potential, mV')
pl.show()
