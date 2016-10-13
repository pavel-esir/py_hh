/*
 * hh_main.cpp
 *
 *  Created on: 27 июня 2016 г.
 *      Author: Pavel Esir
 */

#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cstring>
#include "common_declrations.h"
#include "hh_main.h"
#include "hh_main_gpu.h"
#include <cstdio>
#include <iostream>

#define CUDA_CHECK_RETURN(value) {                                    \
    cudaError_t _m_cudaStat = value;                                  \
    if (_m_cudaStat != cudaSuccess) {                                 \
        fprintf(stderr, "Error %s at line %d in file %s\n",           \
                cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__); \
        exit(1);                                                      \
    }}

namespace hh{

    unsigned int Tsim;
    unsigned int Nneur;
    unsigned int Ncon;
    unsigned int recInt;
    float h;

    float *V_m;
    float *V_m_last;
    float *Vrec;
    float *Vrec_out;
    float *n_ch, *m_ch, *h_ch;
    unsigned int *spike_time, *num_spike_neur, *num_spike_syn;
    unsigned int *spike_time_out, *num_spike_neur_out;
    unsigned int numSpikeSz;

    IncSpikes incSpikes;

    float *I_e;
    float *y, *I_syn;
    float rate;
    float tau_psc;
    float *d_w_p;
    unsigned int *psn_time, *psn_seed;

    // spikes with definite times which are fed into neurons separately from Poisson
    unsigned int *ext_spk_times;
    unsigned int *ext_spk_num;

    float *Inoise;
    float tau_cor = 0.2;
    float *D;

    float *weight;
    unsigned int *delay;
    unsigned int *pre_nidx;
    unsigned int *post_nidx;

    RecordVars rv;
    NeurVars nv;
}

void set_calc_params(unsigned int Tsim, unsigned int Nneur, unsigned int Ncon, unsigned int recInt, float h){
    hh::Tsim = Tsim;
    hh::Nneur = Nneur;
    hh::Ncon = Ncon;
    hh::recInt = recInt;
    hh::h = h;
    cudaMalloc((void**) &hh::V_m_last, Nneur*sizeof(float));
    cudaMalloc((void**) &hh::psn_time, Nneur*sizeof(unsigned int));
    cudaMalloc((void**) &hh::psn_seed, Nneur*sizeof(unsigned int));
    cudaMalloc((void**) &hh::Inoise, Nneur*sizeof(float));
    cudaMalloc((void**) &hh::num_spike_syn, Ncon*sizeof(unsigned int));

    cudaMemset(hh::Inoise, 0, Nneur*sizeof(float));
    fillFloatArr_gpu(Nneur, hh::V_m_last, 0.0f);
    cudaMemset(hh::num_spike_syn, 0, Ncon*sizeof(unsigned int));
}

void set_spike_times(unsigned int *spike_time, unsigned int *num_spike_neur, unsigned int sz){
    hh::numSpikeSz = sz;
    cudaMalloc((void**) &hh::spike_time, sz*sizeof(unsigned int));
    cudaMalloc((void**) &hh::num_spike_neur, sizeof(unsigned int)*hh::Nneur);

    cudaMemcpy(hh::spike_time, spike_time, sz*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(hh::num_spike_neur, num_spike_neur, sizeof(unsigned int)*hh::Nneur, cudaMemcpyHostToDevice);

    hh::spike_time_out = spike_time;
    hh::num_spike_neur_out = num_spike_neur;
}

void set_conns(float *weight, unsigned int *delay, unsigned int *pre, unsigned int *post){
    cudaMalloc((void**) &hh::weight, sizeof(float)*hh::Ncon);
    cudaMalloc((void**) &hh::delay, sizeof(unsigned int)*hh::Ncon);
    cudaMalloc((void**) &hh::pre_nidx, sizeof(unsigned int)*hh::Ncon);
    cudaMalloc((void**) &hh::post_nidx, sizeof(unsigned int)*hh::Ncon);

    cudaMemcpy(hh::weight, weight, sizeof(float)*hh::Ncon, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::delay, delay, sizeof(unsigned int)*hh::Ncon, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::pre_nidx, pre, sizeof(unsigned int)*hh::Ncon, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::post_nidx, post, sizeof(unsigned int)*hh::Ncon, cudaMemcpyHostToDevice);
}

void set_neur_vars(float *V_m, float *Vrec, float *n_ch, float *m_ch, float *h_ch){
    CUDA_CHECK_RETURN(cudaMalloc((void**) &(hh::V_m), sizeof(float)*hh::Nneur));
    cudaMalloc((void**) &hh::Vrec, sizeof(float)*hh::Nneur*(hh::Tsim + hh::recInt - 1)/ hh::recInt);
    cudaMalloc((void**) &hh::n_ch, sizeof(float)*hh::Nneur);
    cudaMalloc((void**) &hh::m_ch, sizeof(float)*hh::Nneur);
    cudaMalloc((void**) &hh::h_ch, sizeof(float)*hh::Nneur);

    CUDA_CHECK_RETURN(cudaMemcpy(hh::V_m, V_m, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice));
    cudaMemcpy(hh::Vrec, Vrec, sizeof(float)*hh::Nneur*(hh::Tsim + hh::recInt - 1) / hh::recInt, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::n_ch, n_ch, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::m_ch, m_ch, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::h_ch, h_ch, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice);
    hh::Vrec_out = Vrec;
}

void set_currents(float *I_e, float *y, float *I_syn, float rate, float tau_psc, float *d_w_p, unsigned int seed){
    cudaMalloc((void**) &hh::I_e, sizeof(float)*hh::Nneur);
    cudaMalloc((void**) &hh::y, sizeof(float)*hh::Nneur);
    cudaMalloc((void**) &hh::I_syn, sizeof(float)*hh::Nneur);
    cudaMalloc((void**) &hh::d_w_p, sizeof(float)*hh::Nneur);

    cudaMemcpy(hh::I_e, I_e, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::y, y, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::I_syn, I_syn, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice);
    cudaMemcpy(hh::d_w_p, d_w_p, sizeof(float)*hh::Nneur, cudaMemcpyHostToDevice);

    hh::rate = rate;

    hh::tau_psc = tau_psc;

    using namespace hh;
    init_noise_gpu(seed, Nneur, h, rate, psn_seed, psn_time);
}

void set_incom_spikes(unsigned int *times, unsigned int *nums, float* weights, unsigned int MaxNumIncom){
    CUDA_CHECK_RETURN(cudaMalloc((void**) &hh::incSpikes.times, MaxNumIncom*hh::Nneur*sizeof(unsigned int)));
    CUDA_CHECK_RETURN(cudaMalloc((void**) &hh::incSpikes.weights, MaxNumIncom*hh::Nneur*sizeof(float)));
    CUDA_CHECK_RETURN(cudaMalloc((void**) &hh::incSpikes.nums, hh::Nneur*sizeof(unsigned int)));
    CUDA_CHECK_RETURN(cudaMalloc((void**) &hh::incSpikes.numProcessed, hh::Nneur*sizeof(unsigned int)));
    hh::incSpikes.MAXSZ = MaxNumIncom;

    CUDA_CHECK_RETURN(cudaMemcpy(hh::incSpikes.times, times, MaxNumIncom*hh::Nneur*sizeof(unsigned int), cudaMemcpyHostToDevice));
    cudaMemcpy(hh::incSpikes.weights, weights, MaxNumIncom*hh::Nneur*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(hh::incSpikes.nums, nums, hh::Nneur*sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemset(hh::incSpikes.numProcessed, 0, hh::Nneur*sizeof(unsigned int));
}

using namespace hh;

void simulate_cpp(){
    nv.V = V_m;
    nv.V_last = V_m_last;
    nv.n = n_ch;
    nv.m = m_ch;
    nv.h = h_ch;
    nv.Ie = I_e;
    nv.Isyn = I_syn;
    nv.y = y;
    nv.Inoise = Inoise;
    nv.weight_p = d_w_p;

    rv.V = Vrec;
    rv.interval = recInt;
    float exp_psc =  exp(-h/tau_psc);
    float exp_psc_half =  exp(-h*0.5/tau_psc);
    
    for (unsigned int t = 0; t < Tsim; t++){
        if (t % 50000 == 0){
            printf("%.3f\n", t*h);
        }
        
        integrate_neurons_gpu(t, Nneur, h, rate, psn_seed, psn_time, exp_psc, exp_psc_half, tau_cor, nv, rv, num_spike_neur, spike_time, incSpikes);
        CUDA_CHECK_RETURN(cudaGetLastError());

        integrate_synapses_gpu(t, Ncon, Nneur, pre_nidx, post_nidx, weight, y, delay, num_spike_syn, num_spike_neur, spike_time);
        CUDA_CHECK_RETURN(cudaGetLastError());

    }

    cudaMemcpy(Vrec_out, Vrec, sizeof(float)*Nneur*(Tsim + recInt - 1) / recInt, cudaMemcpyDeviceToHost);
    cudaMemcpy(spike_time_out, spike_time, numSpikeSz*sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaMemcpy(num_spike_neur_out, num_spike_neur, sizeof(unsigned int)*Nneur, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
}
