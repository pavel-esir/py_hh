# coding: utf-8
'''
Created on 29 июня 2016 г.

@author: Pavel Esir
'''
fltSize = 'float32'

import py_hh as phh
import numpy as np
from numpy import random
import matplotlib.pylab as pl
from distribute_delays import getDelays

random.seed(0)
psn_seed = 0

SimTime = 1000.
h = 0.02
Tsim = int(SimTime/h)
recInt = 4

tau_psc = 0.2   # ms
w_p = 0.0      # Poisson noise weith, pA
w_n = 0.2       # connection weight, pA
rate = 185.0    # Poisson noise rate, Hz (shouldn't  be 0)

Nneur = 100
Ncon = int(Nneur*Nneur*0.2)

pre = random.randint(0, Nneur, Ncon).astype('uint32')
post = random.randint(0, Nneur, Ncon).astype('uint32')
delay = (np.ones(Ncon)*1/h).astype('uint32')
#delay = (getDelays(Ncon)*0.0/h).astype('uint32')
weight = np.zeros(Ncon, dtype=fltSize) + w_n*np.e/tau_psc

spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

Vrec = np.zeros((int((Tsim  + recInt - 1)/recInt), Nneur), dtype=fltSize)

#v0, n0, m0, h0 = 32.906693, 0.574676, 0.913177, 0.223994
##v0, n0, m0, h0 = -60.8457, 0.3763, 0.0833, 0.4636
#Vm = np.zeros(Nneur, dtype=fltSize) + v0
#ns = np.zeros(Nneur, dtype=fltSize) + n0
#ms = np.zeros(Nneur, dtype=fltSize) + m0
#hs = np.zeros(Nneur, dtype=fltSize) + h0

phases = random.choice(len(np.load('../Vm_cycle.npy')), Nneur)
Vm = np.load('../Vm_cycle.npy')[phases]
ns = np.load('../n_cycle.npy')[phases]
ms = np.load('../m_cycle.npy')[phases]
hs = np.load('../h_cycle.npy')[phases]

Ie = 7.134 + 0.01*random.randn(Nneur).astype(fltSize)
Ie[:int(Nneur/2)] = 0.0
#Ie = np.zeros(Nneur, dtype=fltSize) + 7.134
y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)
d_w_p = np.zeros(Nneur, dtype=fltSize) + w_p*np.e/tau_psc

NnumSp = 1
wInc = 0.0
nums = np.zeros(Nneur, dtype='uint32') + 1
incTimes = np.zeros((NnumSp, Nneur), dtype='uint32')
incSpWeights = np.zeros((NnumSp, Nneur), dtype=fltSize) + wInc*np.e/tau_psc
#%%
phh.setCalcParams(Tsim, Nneur, Ncon, recInt, h)

phh.setIncomSpikes(incTimes, nums, incSpWeights, NnumSp)
phh.setNeurVars(Vm, Vrec, ns, ms, hs)
phh.setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, psn_seed)

phh.setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)

phh.setConns(weight, delay,  pre, post)

phh.simulate()
#%%
pl.figure()
pl.plot(np.linspace(0, SimTime/1000, int((Tsim + recInt - 1)/recInt)), Vrec[:, ::10])
pl.xlabel('Time, s')
pl.ylabel('Membrane potential, mV')
pl.show()
pl.xlim((0, SimTime/1000))
#%%
# combine all spike
#spikes = spike_times[:num_spikes_neur[0], 0]
#senders = np.array([0]*num_spikes_neur[0])
#for i, nsp in zip(xrange(1, Nneur), num_spikes_neur[1:]):
#    spikes = np.concatenate((spikes, spike_times[:nsp, i]))
#    senders = np.concatenate((senders, [i]*nsp))
#
##pl.figure()
##pl.plot(spikes*h, senders, '.')
#
#pl.figure()
#pl.hist(spikes*h/1000, bins=int(SimTime/20), histtype='step')
#pl.xlabel('Time, s')
#
#pl.show()
