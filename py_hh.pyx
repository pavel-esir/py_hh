'''
Created on 27 июня 2016 г.

@author: Pavel Esir
'''

cimport numpy as np
import numpy as np
import cython

cdef extern from "hh_main.h": 
    void set_spike_times(unsigned int *spike_time, unsigned int* num_spike_neur, unsigned int sz)

    void set_conns(float *weight, unsigned int *delay, unsigned int *pre, unsigned int *post)
    
    void set_calc_params(unsigned int Tsim, unsigned int cutoff_ns_tm, unsigned int Nneur, unsigned int Ncon, unsigned int recInt, float h)
    
    void set_neur_vars(float *V_m, float *Vrec, float *n_ch, float *m_ch, float *h_ch)
    
    void set_currents(float *I_e, float *y, float *I_syn, float *rate, float tau_psc, float *d_w_p, unsigned int seed)
    
    void set_incom_spikes(unsigned int *times, unsigned int *nums, float *weights, unsigned int MaxNumIncom)
    
    void simulate_cpp()

def setSpikeTimes(np.ndarray[np.uint32_t, ndim=2] spike_time, np.ndarray[np.uint32_t, ndim=1] num_spike_neur, unsigned int sz):
    set_spike_times(<unsigned int*> spike_time.data, <unsigned int*> num_spike_neur.data, sz)
    
def setConns(np.ndarray[np.float32_t, ndim=1] weight, np.ndarray[np.uint32_t, ndim=1] delay, np.ndarray[np.uint32_t, ndim=1]  pre, np.ndarray[np.uint32_t, ndim=1] post):
    set_conns(<float*> weight.data, <unsigned int*> delay.data, <unsigned int*> pre.data, <unsigned int*> post.data)

def setCalcParams(unsigned int Tsim, unsigned int cutoff_ns_tm, unsigned int Nneur, unsigned int Ncon, unsigned int recInt, float h):
    set_calc_params(Tsim, cutoff_ns_tm, Nneur, Ncon, recInt, h)

def setNeurVars(np.ndarray[np.float32_t, ndim=1] V_m, np.ndarray[np.float32_t, ndim=2] Vrec, 
                np.ndarray[np.float32_t, ndim=1] n_ch, 
                np.ndarray[np.float32_t, ndim=1] m_ch, 
                np.ndarray[np.float32_t, ndim=1] h_ch):
    set_neur_vars(<float*> V_m.data, <float*> Vrec.data, 
                  <float*> n_ch.data, <float*> m_ch.data, <float*> h_ch.data)

def setCurrents(np.ndarray[np.float32_t, ndim=1] I_e, np.ndarray[np.float32_t, ndim=1] y, 
                np.ndarray[np.float32_t, ndim=1] I_syn, np.ndarray[np.float32_t, ndim=1] rate, 
                float tau_psc, np.ndarray[np.float32_t, ndim=1] d_w_p, unsigned int seed):
    set_currents(<float*> I_e.data, <float*> y.data, <float*> I_syn.data, <float*> rate.data, tau_psc, <float*> d_w_p.data, seed)

def setIncomSpikes(np.ndarray[np.uint32_t, ndim=2] times, np.ndarray[np.uint32_t, ndim=1] nums, np.ndarray[np.float32_t, ndim=2] weights, unsigned int maxnum):
    set_incom_spikes(<unsigned int*> times.data, <unsigned int*> nums.data, <float*> weights.data, maxnum)

def simulate():
    simulate_cpp()
