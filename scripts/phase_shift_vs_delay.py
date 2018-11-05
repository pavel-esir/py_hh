# coding: utf-8
"""
Created on Tue Oct  4 11:42:51 2016

@author: Pavel Esir
"""
from __future__ import division
import py_hh as phh
import numpy as np
from numpy import random
import matplotlib.pylab as pl
fltSize = 'float32'

h = 0.02
T = 20.28
Ntrans = 80
SimTime = T*Ntrans
Tsim = int(SimTime/h)
recInt = np.iinfo(np.uint32).max

I = 5.27
w = 5.
w_p = 0.00
tau_psc = 0.2

dDelay = 0.1
dphi = 0.1
delaysRange = ((np.arange(0, T, dDelay))/h).astype('uint32')
phisRange = (np.arange(0, T, dphi)/h).astype('uint32')

Ndelays = len(delaysRange)
Nphis = len(phisRange)
Nneur = Ndelays*Nphis*2
Ncon = Ndelays*Nphis*2

rate = np.zeros(Nneur, dtype=fltSize) + 50.0 # Hz
Ie = np.zeros(Nneur, dtype=fltSize) + I

delays = np.zeros(Ncon, dtype='uint32')
pre = np.zeros(Ncon, dtype='uint32')
post = np.zeros(Ncon, dtype='uint32')

pre = np.arange(0, Ncon, dtype='uint32')

v0, n0, m0, h0 = 32.906693, 0.574676, 0.913177, 0.223994
Vm = np.zeros(Nneur, dtype=fltSize) + v0
ns = np.zeros(Nneur, dtype=fltSize) + n0
ms = np.zeros(Nneur, dtype=fltSize) + m0
hs = np.zeros(Nneur, dtype=fltSize) + h0

for i, phi in enumerate(phisRange):
    post[i*2*Ndelays:(i + 1)*2*Ndelays] = np.arange((i + 1)*2*Ndelays - 1, i*2*Ndelays - 1, -1, dtype='uint32')

    delays[2*i*Ndelays:(2*i + 1)*Ndelays] = delaysRange
    delays[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = delaysRange[::-1]

    Vm[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../Vm_cycle.npy')[phi]
    ns[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../n_cycle.npy')[phi]
    ms[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../m_cycle.npy')[phi]
    hs[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../h_cycle.npy')[phi]
#%%
weight = np.zeros(Ncon, dtype=fltSize) + w*np.e/tau_psc

spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

Vrec = np.zeros((int((Tsim  + recInt - 1)/recInt), Nneur), dtype=fltSize)

incTimes = np.zeros((1, Nneur), dtype='uint32')
nums = np.zeros(Nneur, dtype='uint32')
incSpWeights = np.zeros((1, Nneur), dtype=fltSize)

y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)
d_w_p = np.zeros(Nneur, dtype=fltSize) + w_p*np.e/tau_psc

phh.setCalcParams(Tsim, np.iinfo(np.uint32).max, Nneur, Ncon, recInt, h)
phh.setIncomSpikes(incTimes, nums, incSpWeights, 1)
phh.setNeurVars(Vm, Vrec, ns, ms, hs)
phh.setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, 0)
phh.setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)
phh.setConns(weight, delays,  pre, post)
phh.simulate()
#%%
dphiSteady = np.zeros((Nphis, Ndelays)) + np.nan
periodSteady = np.zeros((Nphis, Ndelays)) + np.nan
Nmeans = 10
for i in xrange(Nphis):
    for j in xrange(Ndelays):
        idx1 = 2*i*Ndelays + j
        idx2 = (2*(i + 1))*Ndelays - j - 1
        spkTime1 = h*spike_times[num_spikes_neur[idx1] - 1, idx1]
        spkTime2 = h*spike_times[num_spikes_neur[idx2] - 1, idx2]
        if (spkTime1 > SimTime - 2*T) and (spkTime2 > SimTime - 2*T):
#            dphiSteady[i, j] = np.mean(h*spike_times[num_spikes_neur[idx1] - Nmeans:num_spikes_neur[idx1] - 1, idx1] -
#                               h*spike_times[num_spikes_neur[idx2] - Nmeans:num_spikes_neur[idx2] - 1, idx2])
#            periodSteady[i, j] = np.mean(h*spike_times[num_spikes_neur[idx1] - Nmeans:num_spikes_neur[idx1] - 1, idx1] -
#                               h*spike_times[num_spikes_neur[idx1] - 1 - Nmeans:num_spikes_neur[idx1] - 2, idx1])
            dphiSteady[i, j] = spkTime1 - spkTime2
            periodSteady[i, j] = spkTime1 - h*spike_times[num_spikes_neur[idx1] - 2, idx1]
arr = np.abs(np.min([np.abs(dphiSteady), periodSteady - np.abs(dphiSteady)], axis=0))

mArr = np.ma.masked_invalid(arr)
periodSteady = np.ma.masked_invalid(periodSteady)
#%%
pl.figure()
pl.pcolormesh(delaysRange*h, phisRange*h, mArr)
cb = pl.colorbar()
#cb.set_clim(vmin=0., vmax=21.0)
pl.xlim((0, T))
pl.ylim((0, T))
pl.xlabel("Delay[ms]")
pl.ylabel(r"$\Delta \varphi$ [ms]")

#%%
# Plot slice from 2D diagramm
# Dependence of steady state phase shift from initial dphi for fixed delay
#d = 1.0
#pl.figure()
#pl.plot(phisRange*h, mArr[:, int(d//0.1)])
#pl.xlabel(r'$\Delta \phi$')
##pl.plot(delaysRange*h, mArr[0.2//1.0, :] % T)
##pl.xlabel(r'$Delay$')
#pl.xlim((0, T))
#pl.show()
#%%
pl.figure()
Nbins = int(mArr.max()//0.04)
tmpArr = np.zeros(Nbins + 2)

for i, d in enumerate(delaysRange):
    hst, bins = np.histogram(mArr[:, i].compressed(), bins=Nbins, range=(0, mArr.max()))
    tmpArr[1:-1] = hst
    if not mArr[:, i].mask.all():
        indices = np.argmax(tmpArr)
        dphi = bins[indices - 1]
        pl.scatter([d*h], bins[indices - 1], color='C0')

    hst, bins = np.histogram(periodSteady[:, i].compressed(), bins=Nbins, range=(0, periodSteady.max()))
    tmpArr[1:-1] = hst
    if not mArr[:, i].mask.all():
        indices = np.argmax(hst)
        pl.scatter([d*h], bins[indices - 1], color='C2')
        pl.scatter([d*h], (dphi + d*h) % bins[indices - 1], color='C1')
pl.hlines(12, 0, T)
pl.hlines(15, 0, T)
#pl.xlim(0, T)
#pl.ylim(0, 25)
#pl.ylim(0, mArr.max()+1)
#%%
#phiIdx = int(8.0/dphi)
#(f, ax) = pl.subplots(Ndelays, 1, sharex=True)
#for i, a in enumerate(ax):
#    a.plot(np.linspace(0, SimTime, int((Tsim + recInt - 1)/recInt)), Vrec[:, phiIdx*2*Ndelays + i])
#    a.plot(np.linspace(0, SimTime, int((Tsim + recInt - 1)/recInt)), Vrec[:, phiIdx*2*Ndelays + i + Ndelays])
#    a.set_xlim((0, SimTime))
#
#pl.xlabel('time, ms')
pl.show()
