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

SimTime = 1000000.
h = 0.02
Tsim = int(SimTime/h)
recInt = np.iinfo(np.int32).max
#recInt = 4

w_ps = np.arange(1.92, 1.971, 0.01)
nw = len(w_ps)

N = 200
Nneur = N*nw

tau_psc = 0.2   # ms
w_n = 1.0       # connection weight, pA
I0 = 5.27

rate = np.zeros(Nneur, dtype=fltSize) + 185.0    # Poisson noise rate, Hz (shouldn't  be 0)

#pcon = 0.3
#N = 100
#w_n = 1
Ncon = int(N*N*0.15)
#Ncon = 2
pre = np.zeros(Ncon*nw, dtype='uint32')
post = np.zeros(Ncon*nw, dtype='uint32')
delay = np.zeros(Ncon*nw, dtype='uint32')
d_w_p = np.zeros(Nneur, dtype=fltSize)

np.random.seed(0)
preTmp = random.randint(0, N, Ncon).astype('uint32')
postTmp = random.randint(0, N, Ncon).astype('uint32')
delaysTmp = (getDelays(Ncon)/h).astype('uint32')

#preTmp = np.array([0, 1], dtype='uint32')
#postTmp = np.array([1, 0], dtype='uint32')
#delaysTmp = (np.array([6.,6.])/h).astype('uint32')

for idx, w_p in enumerate(w_ps):
    pre[idx*Ncon:(idx + 1)*Ncon] = preTmp + idx*N
    post[idx*Ncon:(idx + 1)*Ncon] = postTmp + idx*N
    delay[idx*Ncon:(idx + 1)*Ncon] = delaysTmp
    d_w_p[idx*N:(idx + 1)*N] = w_p*np.e/tau_psc
#%%
#import csv
#d = getDelays(Ncon)
#with open('delays.csv', 'wb') as csvfile:
#    wrtr = csv.writer(csvfile, delimiter=' ')
#    for i in xrange(Ncon):    
#        wrtr.writerow([pre[i], post[i], delay[i]*h])
#%%
weight = np.zeros(Ncon*nw, dtype=fltSize) + w_n*np.e/tau_psc

spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

Vrec = np.zeros((int((Tsim  + recInt - 1)/recInt), Nneur), dtype=fltSize)

#v0, n0, m0, h0 = 32.906693, 0.574676, 0.913177, 0.223994
v0, n0, m0, h0 = -60.8457, 0.3763, 0.0833, 0.4636
Vm = np.zeros(Nneur, dtype=fltSize) + v0
ns = np.zeros(Nneur, dtype=fltSize) + n0
ms = np.zeros(Nneur, dtype=fltSize) + m0
hs = np.zeros(Nneur, dtype=fltSize) + h0

#phases = random.choice(len(np.load('../Vm_cycle.npy')), Nneur)
#Vm = np.load('../Vm_cycle.npy')[phases].astype(fltSize)
#ns = np.load('../n_cycle.npy')[phases].astype(fltSize)
#ms = np.load('../m_cycle.npy')[phases].astype(fltSize)
#hs = np.load('../h_cycle.npy')[phases].astype(fltSize)

#Ie = 7.134 + 0.01*random.randn(Nneur).astype(fltSize)
Ie = 5.27 + 0.0*random.randn(Nneur).astype(fltSize)
#Ie[:int(Nneur/2)] = 0.0
#Ie = np.zeros(Nneur, dtype=fltSize) + 7.134
y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)

NnumSp = 1
wInc = 0.0
nums = np.zeros(Nneur, dtype='uint32') + 1
incTimes = np.zeros((NnumSp, Nneur), dtype='uint32')
incSpWeights = np.zeros((NnumSp, Nneur), dtype=fltSize) + wInc*np.e/tau_psc
#%%
phh.setCalcParams(Tsim, Nneur, Ncon*nw, recInt, h)

phh.setIncomSpikes(incTimes, nums, incSpWeights, NnumSp)
phh.setNeurVars(Vm, Vrec, ns, ms, hs)
phh.setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, psn_seed)

phh.setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)

phh.setConns(weight, delay,  pre, post)

phh.simulate()
#%%
#pl.figure()
#pl.plot(np.linspace(0, SimTime/1000, int((Tsim + recInt - 1)/recInt)), Vrec[:, 1])
#pl.xlabel('Time, s')
#pl.ylabel('Membrane potential, mV')
#pl.show()
#pl.xlim((0, SimTime/1000))
#%%
#pl.hold(True)
#for idx, w_p in enumerate(xrange(Nneur)):
#    pl.plot(np.linspace(0, SimTime/1000, int((Tsim + recInt - 1)/recInt)), Vrec[:, idx], label=str(idx))
#pl.legend()
#%%
# combine all spike
import numpy.ma as ma
#
spikes = spike_times[:num_spikes_neur[0], 0]
senders = np.array([0]*num_spikes_neur[0])
for i, nsp in zip(xrange(1, Nneur), num_spikes_neur[1:]):
    spikes = np.concatenate((spikes, spike_times[:nsp, i]))
    senders = np.concatenate((senders, [i]*nsp))

(f, ax) = pl.subplots(nw, 1, sharex=True, sharey=True)
#(f2, ax2) = pl.subplots(nw, 1, sharex=True, sharey=True)
#for idx, (a, a2) in enumerate(zip(ax, ax2)):
for idx, a in enumerate(ax):
    mask = ma.masked_inside(senders, N*idx, N*(idx + 1) - 1)
    a.hist(spikes[mask.mask]*h/1000, bins=int(SimTime/20), histtype='step')
#    a.plot(spikes[mask.mask]*h/1000, senders[mask.mask] - idx*N, '.')
    a.set_title(w_ps[idx])
    a.set_xlim((0, SimTime/1000))
    
#    a2.plot(np.linspace(0, SimTime/1000, int((Tsim + recInt - 1)/recInt)), Vrec[:, N*idx:N*(idx + 1)])
#    a2.set_xlim((0, SimTime/1000))
pl.show()
