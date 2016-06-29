/*
 * hh_main.h
 *
 *  Created on: 27 июня 2016 г.
 *      Author: Pavel Esir
 */

#ifndef HH_MAIN_H_
#define HH_MAIN_H_


void init_noise(unsigned int seed);

void set_delays(unsigned int *delays);

void set_calc_params(unsigned int Tsim, unsigned int Nneur, unsigned int recInt, double h);

void set_neur_vars(double *V_m, double *Vrec, double *n_ch, double *m_ch, double *h_ch);

void set_currents(double *I_e, double *y, double *I_syn, double rate, double tau_psc, double *d_w_p, unsigned int seed);

void simulate_cpp();

#endif /* HH_MAIN_H_ */
