'''
Created on 27 июня 2016 г.

@author: Pavel Esir
'''

cimport numpy as np
import numpy as np
import cython

cdef extern from "hh_main.h": 
    void set_spike_times(unsigned int *spike_time)

    void set_conns(double *weight, unsigned int *delay, unsigned int *pre, unsigned int *post)
    
    void set_calc_params(unsigned int Tsim, unsigned int Nneur, unsigned int Ncon, unsigned int recInt, double h)
    
    void set_neur_vars(double *V_m, double *Vrec, double *n_ch, double *m_ch, double *h_ch)
    
    void set_currents(double *I_e, double *y, double *I_syn, double rate, double tau_psc, double *d_w_p, unsigned int seed)
    
    void simulate_cpp()

def setSpikeTimes(np.ndarray[np.uint32_t, ndim=2] spike_time):
    set_spike_times(<unsigned int*> spike_time.data)
    
def setConns(np.ndarray[np.float64_t, ndim=1] weight, np.ndarray[np.uint32_t, ndim=1] delay, np.ndarray[np.uint32_t, ndim=1]  pre, np.ndarray[np.uint32_t, ndim=1] post):
    set_conns(<double*> weight.data, <unsigned int*> delay.data, <unsigned int*> pre.data, <unsigned int*> post.data)

def setCalcParams(unsigned int Tsim, unsigned int Nneur, unsigned int Ncon, unsigned int recInt, double h):
    set_calc_params(Tsim, Nneur, Ncon, recInt, h)

def setNeurVars(np.ndarray[np.float64_t, ndim=1] V_m, np.ndarray[np.float64_t, ndim=2] Vrec, 
                np.ndarray[np.float64_t, ndim=1] n_ch, 
                np.ndarray[np.float64_t, ndim=1] m_ch, 
                np.ndarray[np.float64_t, ndim=1] h_ch):
    set_neur_vars(<double*> V_m.data, <double*> Vrec.data, 
                  <double*> n_ch.data, <double*> m_ch.data, <double*> h_ch.data)

def setCurrents(np.ndarray[np.float64_t, ndim=1] I_e, np.ndarray[np.float64_t, ndim=1] y, 
                np.ndarray[np.float64_t, ndim=1] I_syn, 
                double rate, double tau_psc, np.ndarray[np.float64_t, ndim=1] d_w_p, unsigned int seed):
    set_currents(<double*> I_e.data, <double*> y.data, <double*> I_syn.data, rate, tau_psc, <double*> d_w_p.data, seed)

def simulate():
    simulate_cpp()