/*
 * hh_main_gpu.cpp
 *
 *  Created on: 20 июля 2016 г.
 *      Author: Pavel Esir
 */

#include <cmath>
#include "hh_main_gpu.h"
#include <cstdio>

#define Cm_    1.0 //  inverse of membrane capacity, 1/pF
#define g_Na  120.0 // nS
#define g_K   36.0
#define g_L   0.3
#define E_K   -77.0
#define E_Na  55.0
#define E_L   -54.4
#define V_peak 25.0

#define NBlockSz 128
#define SBlockSz 512

__device__ float get_random(unsigned int *seed){
    // Park-Miller generator
    // return random number homogeneously distributed in interval [0:1)
    unsigned long a = 16807;
    unsigned long m = 2147483647;
    unsigned long x = (unsigned long) *seed;
    x = (a * x) % m;
    *seed = (unsigned int) x;
    return ((float) x)/m;
}

__device__
float hh_Vm(float V, float n_ch, float m_ch, float h_ch, float I, float h){
    return (-g_K*(V - E_K)*n_ch*n_ch*n_ch*n_ch - g_Na*(V - E_Na)*m_ch*m_ch*m_ch*h_ch - g_L*(V - E_L) + I)*h*Cm_;
}

__device__
float hh_n_ch(float V, float n_ch, float h){
    float temp = 1.0f - exp(-(V + 55.0f)*0.1f);
    if (temp != 0.0f){
        return (.01f*(1.0f - n_ch)*(V + 55.0f)/temp - 0.125f*n_ch*exp(-(V + 65.0f)*0.0125f))*h;
    } else {
//      printf("dividing to zero while calculating n! \n");
//      to understand why it'so, calculate the limit for v/(1 - exp(v/10)) then v tend to 0
        return (0.1f*(1.0f - n_ch)- 0.125f*n_ch*exp(-(V + 65.0f)*0.0125f))*h;
    }
}

__device__
float hh_m_ch(float V, float m_ch, float h){
    float temp = 1.0f - exp(-(V + 40.0f)*0.1f);
    if (temp != 0.0f){
        return (0.1f*(1.0f - m_ch)*(V + 40.0f)/temp - 4.0f*m_ch*exp(-(V + 65.0f)*0.055555556f))*h;
    } else {
//      printf("dividing to zero while calculating  m! \n");
        return ((1.0f - m_ch) - 4.0f*m_ch*exp(-(V + 65.0f)*0.055555556f))*h;
    }
}

__device__
float hh_h_ch(float V, float h_ch, float h){
    return (0.07f*(1.0f - h_ch)*exp(-(V + 65.0f)*0.05f) - h_ch/(1.0f + exp(-(V + 35.0f)*0.1f)))*h;
}

__global__
void integrate_synapses(unsigned int t, unsigned int Ncon, unsigned int Nneur, unsigned int *pre_nidx, unsigned int *post_nidx, float *weight,
                        float *y, unsigned int *delay, unsigned int *num_spike_syn, unsigned int *num_spike_neur, unsigned int *spike_time){
    unsigned int s = blockIdx.x*blockDim.x + threadIdx.x;
    if (s < Ncon){
        // if we processed less spikes than presynaptic neuron generated
        // we need to check whether the new spikes arrive at this moment of time
        if (num_spike_syn[s] < num_spike_neur[pre_nidx[s]]){
            if (spike_time[Nneur*num_spike_syn[s] + pre_nidx[s]] == t - delay[s]){
                atomicAdd(&y[post_nidx[s]], weight[s]);
                num_spike_syn[s]++;
            }
        }
    }
}

__global__
void integrate_neurons(unsigned int t, unsigned int Nneur, float h, float rate, unsigned int *psn_seed, unsigned int *psn_time,
                       float exp_psc, float exp_psc_half, float tau_cor, NeurVars nv, RecordVars rv,
                       unsigned int *num_spike_neur, unsigned int *spike_time, IncSpikes incSpikes){
    unsigned int n = blockIdx.x*blockDim.x + threadIdx.x;
    if (n < Nneur){
        float I_syn_half = (nv.y[n]*h*0.5f + nv.Isyn[n])*exp_psc_half;


        // if where is poisson impulse on neuron
        while (psn_time[n] == t){
            nv.y[n] += nv.weight_p[n];
            // after taking logarithm from uniformly distributed from 0 to 1
            // random number we get exponentially distributed random number
            // for Poisson process time interals between impulses are exponentially distributed
            // sign of right part is negative hence here is "-="
            psn_time[n] += (unsigned int) (-1000.0f*log(get_random(psn_seed + n))/(rate*h));
        }

        while (incSpikes.numProcessed[n] < incSpikes.nums[n] && incSpikes.times[Nneur*incSpikes.numProcessed[n] + n] == t){
            nv.y[n] += incSpikes.weights[Nneur*incSpikes.numProcessed[n] + n];
            incSpikes.numProcessed[n] += 1;
        }
        float V_mem, n_channel, m_channel, h_channel;
        float v1, v2, v3, v4;
        float n1, n2, n3, n4;
        float m1, m2, m3, m4;
        float h1, h2, h3, h4;
        float Inoise_;
        float ns1, ns2, ns3, ns4;

        float dNoise = 0.0f;
    //    float dNoise = sqrtf(2.0f*h*D[n])*curand_normal(&state[n]);

        V_mem = nv.V[n];
        n_channel = nv.n[n];
        m_channel = nv.m[n];
        h_channel = nv.h[n];
        Inoise_ = nv.Inoise[n];
        v1 = hh_Vm(nv.V[n], nv.n[n], nv.m[n], nv.h[n], nv.Isyn[n] + nv.Inoise[n] + nv.Ie[n], h);
        n1 = hh_n_ch(nv.V[n], nv.n[n], h);
        m1 = hh_m_ch(nv.V[n], nv.m[n], h);
        h1 = hh_h_ch(nv.V[n], nv.h[n], h);
        ns1 = (-nv.Inoise[n]*h + dNoise)/tau_cor;
        nv.V[n] = V_mem + v1/2.0f;
        nv.n[n] = n_channel + n1/2.0f;
        nv.m[n] = m_channel + m1/2.0f;
        nv.h[n] = h_channel + h1/2.0f;
        nv.Inoise[n] = Inoise_ + ns1/2.0f;

        v2 = hh_Vm(nv.V[n], nv.n[n], nv.m[n], nv.h[n], I_syn_half + nv.Inoise[n] + nv.Ie[n], h);
        n2 = hh_n_ch(nv.V[n], nv.n[n], h);
        m2 = hh_m_ch(nv.V[n], nv.m[n], h);
        h2 = hh_h_ch(nv.V[n], nv.h[n], h);
        ns2 = (-nv.Inoise[n]*h + dNoise)/tau_cor;
        nv.V[n] = V_mem + v2/2.0f;
        nv.n[n] = n_channel + n2/2.0f;
        nv.m[n] = m_channel + m2/2.0f;
        nv.h[n] = h_channel + h2/2.0f;
        nv.Inoise[n] = Inoise_ + ns2/2.0f;


        v3 = hh_Vm(nv.V[n], nv.n[n], nv.m[n], nv.h[n], I_syn_half + nv.Inoise[n] + nv.Ie[n], h);
        n3 = hh_n_ch(nv.V[n], nv.n[n], h);
        m3 = hh_m_ch(nv.V[n], nv.m[n], h);
        h3 = hh_h_ch(nv.V[n], nv.h[n], h);
        ns3 = (-nv.Inoise[n]*h + dNoise)/tau_cor;
        nv.V[n] = V_mem + v3;
        nv.n[n] = n_channel + n3;
        nv.m[n] = m_channel + m3;
        nv.h[n] = h_channel + h3;
        nv.Inoise[n] = Inoise_ + ns3;

        nv.Isyn[n]  = (nv.y[n]*h + nv.Isyn[n])*exp_psc;
        nv.y[n] *= exp_psc;

        v4 = hh_Vm(nv.V[n], nv.n[n], nv.m[n], nv.h[n], nv.Isyn[n] + nv.Inoise[n] + nv.Ie[n], h);
        n4 = hh_n_ch(nv.V[n], nv.n[n], h);
        m4 = hh_m_ch(nv.V[n], nv.m[n], h);
        h4 = hh_h_ch(nv.V[n], nv.h[n], h);
        ns4 = (-nv.Inoise[n]*h + dNoise)/tau_cor;

        nv.V[n] = V_mem + (v1 + 2.0f*(v2 + v3) + v4)/6.0f;
        nv.n[n] = n_channel + (n1 + 2.0f*(n2 + n3) + n4)/6.0f;
        nv.m[n] = m_channel + (m1 + 2.0f*(m2 + m3) + m4)/6.0f;
        nv.h[n] = h_channel + (h1 + 2.0f*(h2 + h3) + h4)/6.0f;
        nv.Inoise[n] = Inoise_ + (ns1 + 2.0f*(ns2 + ns3) + ns4)/6.0f;

        // checking if there's spike on neuron
        if (nv.V[n] > V_peak && V_mem > nv.V[n] && nv.V_last[n] <= V_mem){
            // second condition is necessary in the presence of noise
            if (num_spike_neur[n] == 0 || t - spike_time[Nneur*(num_spike_neur[n] - 1) + n] > 5.0f/h){
                spike_time[Nneur*num_spike_neur[n] + n] = t;
                num_spike_neur[n]++;
            }
        }
        nv.V_last[n] = V_mem;

        if (t % rv.interval == 0){
            rv.V[Nneur*t/rv.interval + n] = nv.V[n];
//            rv.V[Nneur*t/rv.interval + n] = nv.Isyn[n];
        }
    }
}

__global__
void init_noise(unsigned int seed, unsigned int Nneur, float h, float rate, unsigned int *psn_seed, unsigned int *psn_time){
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < Nneur){
        psn_seed[i] = 100000*(seed + i + 1);
        psn_time[i] = 1 + (unsigned int) (-1000.0*log(get_random(psn_seed + i))/(rate*h));
    }
}

__global__
void fillFloatArr(unsigned int size, float *arr, float val){
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < size){
        arr[i] = val;
    }
}

void init_noise_gpu(unsigned int seed, unsigned int Nneur, float h, float rate, unsigned int *psn_seed, unsigned int *psn_time){
    init_noise<<<Nneur/NBlockSz + 1, NBlockSz>>>(seed, Nneur, h, rate, psn_seed, psn_time);
}

void fillFloatArr_gpu(unsigned int size, float *arr, float val){
    fillFloatArr<<<(size + NBlockSz - 1)/NBlockSz, NBlockSz>>>(size, arr, val);
}

void integrate_neurons_gpu(unsigned int t, unsigned int Nneur, float h, float rate, unsigned int *psn_seed, unsigned int *psn_time,
        float exp_psc, float exp_psc_half, float tau_cor, NeurVars nv, RecordVars rv, unsigned int *num_spike_neur, unsigned int *spike_time, IncSpikes incSpikes){
    integrate_neurons<<<(Nneur + NBlockSz - 1)/NBlockSz, NBlockSz>>>(t, Nneur, h, rate, psn_seed, psn_time, exp_psc, exp_psc_half, tau_cor, nv, rv,
            num_spike_neur, spike_time, incSpikes);
}

void integrate_synapses_gpu(unsigned int t, unsigned int Ncon, unsigned int Nneur, unsigned int *pre_nidx, unsigned int *post_nidx, float *weight,
                            float *y, unsigned int *delay, unsigned int *num_spike_syn, unsigned int *num_spike_neur, unsigned int *spike_time){
    integrate_synapses<<<(Ncon + SBlockSz - 1)/SBlockSz, SBlockSz>>>(t, Ncon, Nneur, pre_nidx, post_nidx, weight, y, delay, num_spike_syn, num_spike_neur, spike_time);
}

