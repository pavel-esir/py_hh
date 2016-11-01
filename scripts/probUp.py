# coding: utf-8
'''
Created on 1 ноября 2016 г.

@author: Pavel Esir
'''
fltSize = 'float32'
import py_hh as phh
import numpy as np
from numpy import random
import matplotlib.pylab as pl
#random.seed(0)
psnSeed = 0
#psnSeed += Nneur
Ntrans = 4 # length of transient process in periods, integer
T = 20.28
SimTime = (Ntrans+0.5)*T
SimTime = 100000.0
h = 0.02
Tsim = int(SimTime/h)
recInt = 20000000

tau_psc = 0.2  # ms
w_p = 2.0      # Poisson noise weith, pA
w_n = 0.0      # connection weight, pA

Nneur = 10000
Ncon = 2
pre = np.array([0, 1], dtype='uint32')
post = np.array([1, 0], dtype='uint32')
delay = (np.ones(Ncon)*0.0/h).astype('uint32') # delay arrays in time frames
rate = 200.0     # Poisson noise rate, Hz (shouldn't  be 0)
weight = np.zeros(Ncon, dtype=fltSize) + w_n*np.e/tau_psc

# int(SimTime/10) тут 10 это период в мс максимального ожидаемого интервала между спайками
spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

Vrec = np.zeros((int((Tsim  + recInt - 1)/recInt), Nneur), dtype=fltSize)

#v0, n0, m0, h0 = 32.906693, 0.574676, 0.913177, 0.223994
v0, n0, m0, h0 = -61.52724457, 0.37197131, 0.07901067, 0.47200188
Vm = np.zeros(Nneur, dtype=fltSize) + v0
ns = np.zeros(Nneur, dtype=fltSize) + n0
ms = np.zeros(Nneur, dtype=fltSize) + m0
hs = np.zeros(Nneur, dtype=fltSize) + h0

NnumSp = 1
wInc = 0.0
nums = np.zeros(Nneur, dtype='uint32') + NnumSp
incTimes = np.zeros((NnumSp, Nneur), dtype='uint32')
incSpWeights = np.zeros((NnumSp, Nneur), dtype=fltSize) + wInc*np.e/tau_psc
Ie = np.zeros(Nneur, dtype=fltSize) + 5.27

y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)
d_w_p = np.zeros(Nneur, dtype=fltSize) + w_p*np.e/tau_psc
#%%
phh.setCalcParams(Tsim, Nneur, Ncon, recInt, h)
phh.setIncomSpikes(incTimes, nums, incSpWeights, NnumSp)
phh.setNeurVars(Vm, Vrec, ns, ms, hs)
phh.setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, psnSeed)
phh.setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)
phh.setConns(weight, delay,  pre, post)
phh.simulate()
#%%
lenMeanActProf = int(4*T/(recInt*h))
meanActProf = np.zeros(lenMeanActProf)
m = 0
n = 0
spt = np.zeros(Nneur)
for i in xrange(Nneur):
    if num_spikes_neur[i]:
        spt[n] = spike_times[0, i]
        spkTimeIdxRec = spike_times[0, i]//recInt
#        if spt[n]*h > 4*T:
#            meanActProf += Vrec[spkTimeIdxRec - lenMeanActProf:spkTimeIdxRec, i]
#            m += 1
        n += 1
spt = spt[:n]
#meanActProf /= m
#pl.figure(1)
#pl.plot(np.linspace(0, 4*T, lenMeanActProf), meanActProf/m)
#%%
#pl.figure()
#pl.plot(np.linspace(0, SimTime, int((Tsim + recInt - 1)/recInt)), Vrec[:, ::1])
#pl.xlabel('time, ms')
#pl.ylabel('Membrane potential, mV')
#pl.show()
#pl.xlim((0, SimTime))
#np.save('h_cycle.npy', Vrec[4054:5069, 0])
#%%
#pl.figure()
## combine all spike
#spikes = spike_times[:num_spikes_neur[0], 0]
#senders = np.array([0]*num_spikes_neur[0])
#for i, nsp in zip(xrange(1, Nneur), num_spikes_neur[1:]):
#    spikes = np.concatenate((spikes, spike_times[:nsp, i]))
#    senders = np.concatenate((senders, [i]*nsp))
#
#pl.plot(spikes*h, senders, '.')
####pl.xlim((0, 200))
#pl.show()

#pl.hist(spikes*h, bins=int(SimTime/20), histtype='step')
#pl.xlabel('time, ms')
#pl.ylabel('Membrane potential, mV')
#pl.show()
#
