# coding: utf-8
'''
Created on 29 июня 2016 г.

@author: Pavel Esir
'''
fltSize = 'float64'
from py_hh_cpu import *
import numpy as np
import matplotlib.pylab as pl
np.random.seed(0)

Ntrans = 2 # length of transient process in periods, integer
T = 20.275
SimTime = (Ntrans+0.5)*T
h = 0.02
Tsim = int(SimTime/h)
recInt = 20

tau_psc = 0.2  # ms
w_p = 0.0      # Poisson noise weith, pA
w_n = 0.0      # connection weight, pA

rate = 0.1     # Poisson noise rate, Hz (shouldn't  be 0)
#Nneur = 100
#Ncon = int(Nneur*Nneur*0.2)
Nneur = int(T/h)
Ncon = 1
pre = np.random.randint(0, Nneur, Ncon).astype('uint32')
post = np.random.randint(0, Nneur, Ncon).astype('uint32')
delay = np.array([4.0/h]*Ncon).astype('uint32') # delay arrays in time frames
#delay = np.random.lognormal(1.8, 1/6.0, Ncon).astype('uint32') # delay arrays in time frames

#Nneur = 2
#Ncon = 1
#pre = np.array([0], dtype='uint32')
#post = np.array([1], dtype='uint32')
#delay = np.array([10./h], dtype='uint32') # delay arrays in time frames
#rate = 0.1     # Poisson noise rate, Hz (shouldn't  be 0)

# int(SimTime/10) тут 10 это период в мс максимального ожидаемого интервала между спайками
spike_times = np.zeros((int(SimTime/10) + 2, Nneur), dtype='uint32')
num_spikes_neur = np.zeros(Nneur, dtype='uint32')

weight = np.zeros(Ncon, dtype=fltSize) + w_n*np.e/tau_psc

setCalcParams(Tsim, Nneur, Ncon, recInt, h)

Vm = np.zeros(Nneur, dtype=fltSize) + 32.906693
Vrec = np.zeros((int((Tsim  + recInt - 1)/recInt), Nneur), dtype=fltSize)
ns = np.zeros(Nneur, dtype=fltSize) + 0.574676
ms = np.zeros(Nneur, dtype=fltSize) + 0.913177
hs = np.zeros(Nneur, dtype=fltSize) + 0.223994

NnumSp = 1
wInc = 1.3
nums = np.zeros(Nneur, dtype='uint32')
incTimes = np.zeros((NnumSp, Nneur), dtype='uint32')
incSpWeights = np.zeros((NnumSp, Nneur), dtype=fltSize) + wInc*np.e/tau_psc
dts = np.arange(0, Nneur, dtype='uint32')
for i in xrange(NnumSp):
    incTimes[i, :] = (i + 1)*dts
incSpWeights[0, :] = wInc*np.e/tau_psc
nums[:] = 1
#inc_spikes = {1: np.array([40, 50])}
#for k, v in inc_spikes.iteritems():
#    incTimes[0:len(v), k] = np.array(v/h, dtype='uint32')
#    incSpWeights[0:len(v), k] = wInc
#    nums[k] = len(v)
setIncomSpikes(incTimes, nums, incSpWeights, NnumSp)
#%%
#Ie = np.zeros(Nneur) + 5.27

Ie = np.zeros(Nneur, dtype=fltSize)
Ie[:] = 5.27

y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)
d_w_p = np.zeros(Nneur, dtype=fltSize) + w_p*np.e/tau_psc

setNeurVars(Vm, Vrec, ns, ms, hs)
setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, 0)

setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)

setConns(weight, delay,  pre, post)
#%%
simulate()
#%%
#pl.figure()
#pl.plot(np.linspace(0, SimTime, int((Tsim + recInt - 1)/recInt)), Vrec[:, ::2])
#pl.xlabel('time, ms')
#pl.ylabel('Membrane potential, mV')
#pl.show()
#pl.xlim((0, SimTime))
#%%
#pl.figure()
## combine all spike
#spikes = spike_times[:num_spikes_neur[0], 0]
#for i, nsp in zip(xrange(1, Nneur), num_spikes_neur[1:]):
#    spikes = np.concatenate((spikes, spike_times[:nsp, i]))
#
#pl.hist(spikes*h, bins=int(SimTime/20), histtype='step')
#pl.xlabel('time, ms')
##pl.ylabel('Membrane potential, mV')
#pl.show()
##
#%%
pl.figure(1)
msc = np.zeros(np.shape(spike_times), dtype='bool')
lastSpikeTime = np.zeros(Nneur)
for i, n in enumerate(num_spikes_neur):
    if n != 0:
        lastSpikeTime[i] = spike_times[n - 1, i]
        if n < Ntrans:
            lastSpikeTime[i] = np.nan
    else:
        lastSpikeTime[i] = np.nan

pl.plot(dts*h, np.ma.array((Ntrans*T + h*dts) - h*lastSpikeTime, mask=(lastSpikeTime != lastSpikeTime)), '-')
pl.plot(arange(0, T, 0.01), arange(0, T, 0.01))
pl.ylim((0, T))
pl.xlim((0, T))
pl.show()
