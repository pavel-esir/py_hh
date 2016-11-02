/*
 * hh_main.h
 *
 *  Created on: 21 июля 2016 г.
 *      Author: Pavel Esir
 */

#ifndef HH_MAIN_GPU_H_
#define HH_MAIN_GPU_H_

#include "common_declrations.h"

void init_noise_gpu(unsigned int seed, unsigned int Nneur, float h, float *rate, unsigned int *psn_seed, unsigned int *psn_time);

void fillFloatArr_gpu(unsigned int size, float *arr, float val);

void integrate_neurons_gpu(unsigned int t, unsigned int Nneur, float h, float *rate, unsigned int *psn_seed, unsigned int *psn_time,
        float exp_psc, float exp_psc_half, float tau_cor, NeurVars nv, RecordVars rv, unsigned int *num_spike_neur, unsigned int *spike_time, IncSpikes incSpikes);

void integrate_synapses_gpu(unsigned int t, unsigned int Ncon, unsigned int Nneur, unsigned int *pre_nidx, unsigned int *post_nidx, float *weight,
                            float *y, unsigned int *delay, unsigned int *num_spike_syn, unsigned int *num_spike_neur, unsigned int *spike_time);

#endif /* HH_MAIN_GPU_H_ */
