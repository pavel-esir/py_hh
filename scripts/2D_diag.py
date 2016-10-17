# coding: utf-8
'''
Created on 29 июня 2016 г.

@author: Pavel Esir
'''
fltSize = 'float32'
from py_hh import *
import numpy as np
import matplotlib.pylab as pl
np.random.seed(0)

Ntrans = 6 # length of transient process in periods, integer
#T = 20.275
T = 20.28
SimTime = (Ntrans+0.5)*T
h = 0.01
Tsim = int(SimTime/h)
recInt = 1000

tau_psc = 0.2  # ms
w_p = 0.0      # Poisson noise weith, pA
w_n = 0.0      # connection weight, pA
wInc = 1.3

rate = 0.1     # Poisson noise rate, Hz (shouldn't  be 0)
dT = 0.04
Nneur = int(T/dT)**2
Ncon = 1
pre = np.random.randint(0, Nneur, Ncon).astype('uint32')
post = np.random.randint(0, Nneur, Ncon).astype('uint32')
delay = np.array([4.0/h]*Ncon).astype('uint32') # delay arrays in time frames


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

NnumSp = 2
nums = np.zeros(Nneur, dtype='uint32') + NnumSp
incTimes = np.zeros((NnumSp, Nneur), dtype='uint32')
incSpWeights = np.zeros((NnumSp, Nneur), dtype=fltSize) + wInc*np.e/tau_psc

dts = np.linspace(0, int(T/h), int(np.sqrt(Nneur)), endpoint=False, dtype='uint32')
t1, t2 = np.meshgrid(dts, dts)
incTimes[0, :] = t1.reshape(Nneur)
incTimes[1, :] = t2.reshape(Nneur)
incTimes = np.sort(incTimes, axis=0)

setIncomSpikes(incTimes, nums, incSpWeights, NnumSp)
#%%
Ie = np.zeros(Nneur, dtype=fltSize) + 5.27

y = np.zeros(Nneur, dtype=fltSize)
Isyn = np.zeros(Nneur, dtype=fltSize)
d_w_p = np.zeros(Nneur, dtype=fltSize) + w_p*np.e/tau_psc

setNeurVars(Vm, Vrec, ns, ms, hs)
setCurrents(Ie, y, Isyn, rate, tau_psc, d_w_p, 0)

setSpikeTimes(spike_times, num_spikes_neur, np.shape(spike_times)[0]*Nneur)

setConns(weight, delay,  pre, post)
#%%
simulate()
#np.save('spikes_T_{}_h_{}_dT_{}_Ntrans_{}'.format(T, h, dT, Ntrans), (spike_times))
#%%
#pl.figure()
#pl.plot(np.linspace(0, SimTime, int((Tsim + recInt - 1)/recInt)), Vrec[:, :100:10])
#pl.xlabel('time, ms')
#pl.ylabel('Membrane potential, mV')
#pl.show()
#pl.xlim((0, SimTime))
#%%
#pl.figure(1)
#lastSpikeTime = np.zeros(Nneur)
#for i, n in enumerate(num_spikes_neur):
#    if n < Ntrans:
#            lastSpikeTime[i] = np.nan
#    else:
#        lastSpikeTime[i] = spike_times[Ntrans - 1, i]*h
#
#lastSpikeTime = lastSpikeTime.reshape(np.sqrt(Nneur), np.sqrt(Nneur))
#mArr = np.ma.array(lastSpikeTime, mask=(lastSpikeTime != lastSpikeTime))
#
#pl.pcolormesh(dts*h, dts*h, mArr, cmap='jet')
#pl.colorbar()
#pl.xlim((0, T))
#pl.ylim((0, T))
#
#pl.plot(dts*h, np.ma.array((Ntrans*T + h*dts) - h*lastSpikeTime, mask=(lastSpikeTime != lastSpikeTime)), '-')
##pl.plot(arange(0, T, 0.01), arange(0, T, 0.01))
#pl.ylim((0, T))
#pl.xlim((0, T))
#pl.show()
