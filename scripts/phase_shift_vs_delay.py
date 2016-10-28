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
Ntrans = 20
SimTime = T*Ntrans
Tsim = int(SimTime/h)
recInt = 400

I = 5.27
w = 1.
w_p = 0.00
rate = 50.0 # Hz
tau_psc = 0.2

dDelay = 0.1
dphi = 0.1
delaysRange = (np.arange(0, T, dDelay)/h).astype('uint32')
phisRange = (np.arange(0, T, dphi)/h).astype('uint32')

Ndelays = len(delaysRange)
Nphis = len(phisRange)
Nneur = Ndelays*Nphis*2
Ncon = Ndelays*Nphis*2

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
    delays[i*2*Ndelays:(2*i+1)*Ndelays] = delaysRange
    delays[(2*i+1)*Ndelays:(i+1)*2*Ndelays] = delaysRange[::-1]

    Vm[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../Vm_cycle.npy')[phi]
    ns[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../n_cycle.npy')[phi]
    ms[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../m_cycle.npy')[phi]
    hs[(2*i + 1)*Ndelays:2*(i + 1)*Ndelays] = np.load('../h_cycle.npy')[phi]
#%%
weight = np.zeros(Ncon, dtype=fltSize) + w*np.e/tau_psc

spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

Vrec = np.zeros((int((Tsim  + recInt - 1)/recInt), Nneur), dtype=fltSize)

Ie = np.zeros(Nneur, dtype=fltSize) + I

incTimes = np.zeros((1, Nneur), dtype='uint32')
nums = np.zeros(Nneur, dtype='uint32')
incSpWeights = np.zeros((1, Nneur), dtype=fltSize)

y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)
d_w_p = np.zeros(Nneur, dtype=fltSize) + w_p*np.e/tau_psc

phh.setCalcParams(Tsim, Nneur, Ncon, recInt, h)
phh.setIncomSpikes(incTimes, nums, incSpWeights, 1)
phh.setNeurVars(Vm, Vrec, ns, ms, hs)
phh.setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, 0)
phh.setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)
phh.setConns(weight, delays,  pre, post)
phh.simulate()
#%%
dphiSteady = np.zeros((Nphis, Ndelays)) + np.nan
for i in xrange(Nphis):
    for j in xrange(Ndelays):
        idx1 = 2*i*Ndelays + j
        idx2 = (2*(i + 1))*Ndelays - j - 1
#        print(idx1, idx2)
        if (num_spikes_neur[idx1] > Ntrans//2) and (num_spikes_neur[idx2] > Ntrans//2):
            dphiSteady[i, j] = (h*spike_times[num_spikes_neur[idx1] - 1, idx1] -
                                   h*spike_times[num_spikes_neur[idx2] - 1, idx2])
arr = np.abs(np.min([np.abs(dphiSteady), np.abs(dphiSteady) - T], axis=0))
mArr = np.ma.array(arr, mask=(dphiSteady != dphiSteady))
#%%
#pl.figure()
#pl.hist(mArr[:, 10//dDelay].compressed() % T, bins=100)
#%%
pl.figure()
#pl.pcolormesh(delaysRange*h, phisRange*h, mArr)
#pl.pcolormesh(delaysRange*h, phisRange*h, np.ma.masked_greater(mArr, 800))
pl.pcolormesh(delaysRange*h, phisRange*h, np.ma.masked_outside(mArr, -T, T))
cb = pl.colorbar()
cb.set_clim(vmin=-20., vmax=21.0)
pl.xlim((0, T))
pl.ylim((0, T))
pl.xlabel("Delay[ms]")
pl.ylabel(r"$\Delta \varphi$ [ms]")
pl.show()
#%%
#d = 1.0
#pl.figure()
#pl.plot(phisRange*h, mArr[:, d//0.1])
#pl.xlabel(r'$\Delta \phi$')
##pl.plot(delaysRange*h, mArr[0.2//1.0, :] % T)
##pl.xlabel(r'$Delay$')
#pl.xlim((0, T))
#pl.show()
#%%
#pl.figure()
#pl.hist(mArr[:, 6.0//0.1].compressed(), bins=100, range=(0, T))
#pl.show()
#%%
from scipy.signal import argrelextrema
for i, d in enumerate(delaysRange):
    hst, bins= np.histogram(mArr[:, i].compressed(), bins=100, range=(0, T + 1))
    indices = argrelextrema(hst, np.greater)[0]
#    print(len(indices))
    pl.scatter([d*h]*len(indices), bins[indices])
#pl.plot(bins[:-1], hst)

#%%
#phiIdx = int(8.0/dphi)
#(f, ax) = pl.subplots(Ndelays, 1, sharex=True)
#for i, a in enumerate(ax):
#    a.plot(np.linspace(0, SimTime, int((Tsim + recInt - 1)/recInt)), Vrec[:, phiIdx*2*Ndelays + i])
#    a.plot(np.linspace(0, SimTime, int((Tsim + recInt - 1)/recInt)), Vrec[:, phiIdx*2*Ndelays + i + Ndelays])
#    a.set_xlim((0, SimTime))
#
#pl.xlabel('time, ms')
#pl.show()
