/*
 * hh_main.h
 *
 *  Created on: 27 июня 2016 г.
 *      Author: Pavel Esir
 */

#ifndef HH_MAIN_H_
#define HH_MAIN_H_

void set_calc_params(unsigned int Tsim, unsigned int Nneur, unsigned int Ncon, unsigned int recInt, float h);

void set_neur_vars(float *V_m, float *Vrec, float *n_ch, float *m_ch, float *h_ch);

void set_currents(float *I_e, float *y, float *I_syn, float *rate, float tau_psc, float *d_w_p, unsigned int seed);

void set_incom_spikes(unsigned int *times, unsigned int *nums, float *weights, unsigned int MaxNumIncom);

void set_spike_times(unsigned int *spike_time, unsigned int *num_spike_neur, unsigned int sz);

void set_conns(float *weight, unsigned int *delay, unsigned int *pre, unsigned int *post);

void simulate_cpp();

#endif /* HH_MAIN_H_ */
