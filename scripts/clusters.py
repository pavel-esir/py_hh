# coding: utf-8
'''
Created on 29 июня 2016 г.

@author: Pavel Esir
'''
from __future__ import division

fltSize = 'float32'

import py_hh as phh
import numpy as np
from numpy import random
import matplotlib.pylab as pl
from matplotlib.gridspec import GridSpec
from distribute_delays import getDelays
import numpy.ma as ma
fs = 26.

psn_seed = 0
iv_seed = 3

SimTime = 30000
h = 0.02
Tcutoff = np.iinfo(np.int32).max
Tsim = int(SimTime/h)
recInt = np.iinfo(np.int32).max
#recInt = 2

w_ps = np.arange(1.85, 2.001, 0.01)
w_ps = [1.97]
#w_ps = [0.5]
nw = len(w_ps)

N = 100
Nneur = N*nw

tau_psc = 0.2   # ms.
w_n = 1.3       # connection weight, pA
I0 = 5.27

rate = np.zeros(Nneur, dtype=fltSize) + 185.0    # Poisson noise rate, Hz (shouldn't be 0)

pcon = 0.2
Ncon = int(N*N*pcon)

pre = np.zeros(Ncon*nw, dtype='uint32')
post = np.zeros(Ncon*nw, dtype='uint32')
delay = np.zeros(Ncon*nw, dtype='uint32')
d_w_p = np.zeros(Nneur, dtype=fltSize)

np.random.seed(0)

preTmp = random.randint(0, N, Ncon).astype('uint32')
postTmp = random.randint(0, N, Ncon).astype('uint32')
delaysTmp = (getDelays(Ncon, 0)/h).astype('uint32')

for idx, w_p in enumerate(w_ps):
    pre[idx*Ncon:(idx + 1)*Ncon] = preTmp + idx*N
    post[idx*Ncon:(idx + 1)*Ncon] = postTmp + idx*N
    delay[idx*Ncon:(idx + 1)*Ncon] = delaysTmp
    d_w_p[idx*N:(idx + 1)*N] = w_p*np.e/tau_psc
#%%
weight = np.zeros(Ncon*nw, dtype=fltSize) + w_n*np.e/tau_psc

spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

Vrec = np.zeros((int((Tsim  + recInt - 1)/recInt), Nneur), dtype=fltSize)

##v0, n0, m0, h0 = 32.906693, 0.574676, 0.913177, 0.223994
#v0, n0, m0, h0 = -60.8457, 0.3763, 0.0833, 0.4636
#Vm = np.zeros(Nneur, dtype=fltSize) + v0
#ns = np.zeros(Nneur, dtype=fltSize) + n0
#ms = np.zeros(Nneur, dtype=fltSize) + m0
#hs = np.zeros(Nneur, dtype=fltSize) + h0

np.random.seed(iv_seed)
Vmpert = 0.1*np.random.rand(Nneur).astype(fltSize)
phases = random.choice(len(np.load('../Vm_cycle.npy')), Nneur)
Vm = np.load('../Vm_cycle.npy')[phases].astype(fltSize) + Vmpert
ns = np.load('../n_cycle.npy')[phases].astype(fltSize)
ms = np.load('../m_cycle.npy')[phases].astype(fltSize)
hs = np.load('../h_cycle.npy')[phases].astype(fltSize)

Ie = I0 + 0.0*random.randn(Nneur).astype(fltSize)
y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)

NnumSp = 1
wInc = 0.0
nums = np.zeros(Nneur, dtype='uint32') + 1
incTimes = np.zeros((NnumSp, Nneur), dtype='uint32')
incSpWeights = np.zeros((NnumSp, Nneur), dtype=fltSize) + wInc*np.e/tau_psc
#%%
phh.setCalcParams(Tsim, Tcutoff, Nneur, Ncon*nw, recInt, h)

phh.setIncomSpikes(incTimes, nums, incSpWeights, NnumSp)
phh.setNeurVars(Vm, Vrec, ns, ms, hs)
phh.setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, psn_seed)

phh.setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)

phh.setConns(weight, delay,  pre, post)

phh.simulate()
# combine all spike
spikes = spike_times[:num_spikes_neur[0], 0]
senders = np.array([0]*num_spikes_neur[0])

for i, nsp in zip(xrange(1, Nneur), num_spikes_neur[1:]):
    spikes = np.concatenate((spikes, spike_times[:nsp, i]))
    senders = np.concatenate((senders, [i]*nsp))
#%%
#colors = ['C0', 'C1', 'C2', 'C3']
#recIdx = [9, 40, 92]
##recIdx = [9, 48, 94]
#Nrec = len(recIdx)
#
#fig = pl.figure(figsize=(8*1.5, 6*1.5))
#gs = GridSpec(Nrec + 1, 1, height_ratios=[2.5, 1, 1, 1])
#ax = [[]]*(Nrec)
#
#axSpks = fig.add_subplot(gs[0, 0])
#
#mask = ma.masked_inside(senders, 0, N - 1)
#
#for i in xrange(Nrec):
#    ax[i] = fig.add_subplot(gs[0 + Nrec - i, 0], sharex = axSpks)
#    ax[i].plot(np.linspace(0, SimTime/1000., int((Tsim + recInt - 1)/recInt)), Vrec[:, recIdx[i]], lw=1., color=colors[i])
#    tm = spike_times[:num_spikes_neur[recIdx[i]], recIdx[i]]*h/1000
##    axSpks.scatter(tm, [recIdx[i]]*len(tm), s=50, edgecolors=colors[i], facecolors=colors[i], zorder=10)
#    pl.setp(ax[i].get_xticklabels(), visible=False)
#
#axSpks.plot(spikes[mask.mask]*h/1000, senders[mask.mask] - idx*N, '.k', ms = 5)
#pl.setp(axSpks.get_xticklabels(), visible=False)
#
#pl.setp(ax[0].get_xticklabels(), visible=True)
##axSpks.set_xlim([0.0, 2.5])
#axSpks.set_xlim([2.25, 2.5])
##axSpks.set_xlim([2.5, 3.])
#axSpks.set_ylim([0, 100])
#
#axSpks.set_yticks([0, 50, 100])
#
#axSpks.set_ylabel('Neuron number', fontsize=fs)
#ax[1].set_ylabel('$V_m\ (mV)$', fontsize=fs)
#ax[0].set_xlabel('$time\ (s)$', fontsize=fs)
#pl.subplots_adjust(left = 0.10, right = 0.9, top = 0.95, bottom = 0.11, hspace = 0.19)
#pl.savefig('cluster_iv_{}.png'.format(iv_seed), dpi=256.0)
#pl.savefig('cluster_iv_{}_zoomed.png'.format(iv_seed), dpi=256.0)
#pl.savefig('nodelay_synchrony.png', dpi=256.0)
#%%
fig = pl.figure(figsize=(8*1.5, 6*1.5))

gs = GridSpec(2, 1, height_ratios=[3, 1])
axSpks = fig.add_subplot(gs[0])
axHist = fig.add_subplot(gs[1], sharex=axSpks)

mask = ma.masked_inside(senders, 0, N - 1)
axHist.hist(spikes[mask.mask]*h/1000, bins=int(SimTime/20.), color='C3', ls='-', histtype='step', lw=3.)

axSpks.plot(spikes[mask.mask]*h/1000, senders[mask.mask], marker='.', ls='', color='k', ms = 5)

axSpks.set_xlim([4.0, 25.1])
axSpks.set_xticks(np.arange(4., 25.1, 3.))
axSpks.set_ylim([0, 100])

axHist.set_yticks([50, 100])
axHist.set_ylim([0, 100])
axSpks.set_yticks([0, 50, 100])

axSpks.set_ylabel('Neuron number', fontsize=fs)
axHist.set_ylabel('Firing rate', fontsize=fs)
axHist.set_xlabel('$time\ (s)$', fontsize=fs)
pl.setp(axSpks.get_xticklabels(), visible=False)
pl.subplots_adjust(top = 0.97, bottom = 0.10, left = 0.10, right = 0.95, hspace = 0.08)
pl.savefig('cluster_iv_{}_with_noise_{}.png'.format(iv_seed, w_ps[0]), dpi=256.0)
