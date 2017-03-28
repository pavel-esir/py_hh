/*
 * hh_main.cpp
 *
 *  Created on: 27 июня 2016 г.
 *      Author: Pavel Esir
 */

#include <cmath>
#include "hh_main_cpu.h"
#include <cstdio>

#define Cm_    1.0 //  inverse of membrane capacity, 1/pF
#define g_Na  120.0 // nS
#define g_K   36.0
#define g_L   0.3
#define E_K   -77.0
#define E_Na  55.0
#define E_L   -54.4
#define V_peak 25.0

namespace hh{

    unsigned int Tsim;
    unsigned int Nneur;
    unsigned int Ncon;
    unsigned int recInt;
    double h;

    double *V_m;
    double *V_m_last;
    double *Vrec;
    double *n_ch, *m_ch, *h_ch;
    unsigned int *spike_time, *num_spike_neur, *num_spike_syn;

    double *I_e;
    double *y, *I_syn;
    double *rate;
    double tau_psc;
    double exp_psc;
    double exp_psc_half;
    double *d_w_p;
    unsigned int *psn_time, *psn_seed;

    // spikes with definite times which are fed into neurons separately from Poisson
    unsigned int *ext_spk_times;
    unsigned int *ext_spk_num;

    double *Inoise;
    double tau_cor = 0.2;
    double *D;


    double *weight;
    unsigned int *delay;
    unsigned int *pre_nidx;
    unsigned int *post_nidx;

    unsigned int MAXSZ; // [maximum(over all neurons) number of incoming spikes
    unsigned int *incSpTimes; // size of times is  (MaxNumIncom, Nneur)
    double *incSpWeights;
    unsigned int *numIncomSpikes;
    unsigned int *numProcessed;

    double get_random(unsigned int *seed){
        // Park-Miller generator
        // return random number homogeneously distributed in interval [0:1)
        unsigned long a = 16807;
        unsigned long m = 2147483647;
        unsigned long x = (unsigned long) *seed;
        x = (a * x) % m;
        *seed = (unsigned int) x;
        return ((double) x)/m;
    }

    double hh_Vm(double V, double n_ch, double m_ch, double h_ch, double I, double h){
        return (-g_K*(V - E_K)*n_ch*n_ch*n_ch*n_ch - g_Na*(V - E_Na)*m_ch*m_ch*m_ch*h_ch - g_L*(V - E_L) + I)*h*Cm_;
    }

    double hh_n_ch(double V, double n_ch, double h){
        double temp = 1.0 - exp(-(V + 55.0)*0.1);
        if (temp != 0.0){
            return (.01*(1.0 - n_ch)*(V + 55.0)/temp - 0.125*n_ch*exp(-(V + 65.0)*0.0125))*h;
        } else {
    //      printf("dividing to zero while calculating n! \n");
    //      to understand why it'so, calculate the limit for v/(1 - exp(v/10)) then v tend to 0
            return (0.1*(1.0 - n_ch)- 0.125*n_ch*exp(-(V + 65.0)*0.0125))*h;
        }
    }

    double hh_m_ch(double V, double m_ch, double h){
        double temp = 1.0 - exp(-(V + 40.0)*0.1);
        if (temp != 0.0){
            return (0.1*(1.0 - m_ch)*(V + 40.0)/temp - 4.0*m_ch*exp(-(V + 65.0)*0.055555556))*h;
        } else {
    //      printf("dividing to zero while calculating  m! \n");
            return ((1.0 - m_ch) - 4.0*m_ch*exp(-(V + 65.0)*0.055555556))*h;
        }
    }

    double hh_h_ch(double V, double h_ch, double h){
        return (0.07*(1.0 - h_ch)*exp(-(V + 65.0)*0.05) - h_ch/(1.0 + exp(-(V + 35.0)*0.1)))*h;
    }


    void integrate_synapses(unsigned int t, unsigned int s){
        // if we processed less spikes than presynaptic neuron generated
        // we need to check whether the new spikes arrive at this moment of time
        if (num_spike_syn[s] < num_spike_neur[pre_nidx[s]]){
            if (spike_time[Nneur*num_spike_syn[s] + pre_nidx[s]] == t - delay[s]){
                y[post_nidx[s]] += weight[s];
                num_spike_syn[s]++;
            }
        }
    }

//    while (incSpikes.nums[n] != 0 && incSpikes.times[incSpikes.MAXSZ*incSpikes.numProcessed[n] + n] == t){
//        nv.y[n] += incSpikes.weights[incSpikes.MAXSZ*incSpikes.numProcessed[n] + n];
//        if (incSpikes.numProcessed[n] < incSpikes.nums[n]){
//            incSpikes.numProcessed[n] += 1;
//        } else {
//            break;
//        }
//    }

    void integrate_neurons(unsigned int t, unsigned int n){
        double I_syn_half = (y[n]*h*0.5 + I_syn[n])*exp_psc_half;

        // if where is poisson impulse on neuron
        while (psn_time[n] == t){
            y[n] += d_w_p[n];
            // after taking logarithm from uniformly distributed from 0 to 1
            // random number we get exponentially distributed random number
            // for Poisson process time interals between impulses are exponentially distributed
            // sign of right part is negative hence here is "-="
            psn_time[n] += (unsigned int) (-1000.0*log(get_random(psn_seed + n))/(rate[n]*h));
        }

        while (numProcessed[n] < numIncomSpikes[n] && incSpTimes[Nneur*numProcessed[n] + n] == t){
            y[n] += incSpWeights[Nneur*numProcessed[n] + n];
            numProcessed[n] += 1;
        }

        double V_mem, n_channel, m_channel, h_channel;
        double v1, v2, v3, v4;
        double n1, n2, n3, n4;
        double m1, m2, m3, m4;
        double h1, h2, h3, h4;
        double Inoise_;
        double ns1, ns2, ns3, ns4;

        double dNoise = 0.0;
    //    double dNoise = sqrtf(2.0f*h*D[n])*curand_normal(&state[n]);

        V_mem = V_m[n];
        n_channel = n_ch[n];
        m_channel = m_ch[n];
        h_channel = h_ch[n];
        Inoise_ = Inoise[n];
        v1 = hh_Vm(V_m[n], n_ch[n], m_ch[n], h_ch[n], I_syn[n] + Inoise[n] + I_e[n], h);
        n1 = hh_n_ch(V_m[n], n_ch[n], h);
        m1 = hh_m_ch(V_m[n], m_ch[n], h);
        h1 = hh_h_ch(V_m[n], h_ch[n], h);
        ns1 = (-Inoise[n]*h + dNoise)/tau_cor;
        V_m[n] = V_mem + v1/2.0;
        n_ch[n] = n_channel + n1/2.0;
        m_ch[n] = m_channel + m1/2.0;
        h_ch[n] = h_channel + h1/2.0;
        Inoise[n] = Inoise_ + ns1/2.0;

        v2 = hh_Vm(V_m[n], n_ch[n], m_ch[n], h_ch[n], I_syn_half + Inoise[n] + I_e[n], h);
        n2 = hh_n_ch(V_m[n], n_ch[n], h);
        m2 = hh_m_ch(V_m[n], m_ch[n], h);
        h2 = hh_h_ch(V_m[n], h_ch[n], h);
        ns2 = (-Inoise[n]*h + dNoise)/tau_cor;
        V_m[n] = V_mem + v2/2.0;
        n_ch[n] = n_channel + n2/2.0;
        m_ch[n] = m_channel + m2/2.0;
        h_ch[n] = h_channel + h2/2.0;
        Inoise[n] = Inoise_ + ns2/2.0;


        v3 = hh_Vm(V_m[n], n_ch[n], m_ch[n], h_ch[n], I_syn_half + Inoise[n] + I_e[n], h);
        n3 = hh_n_ch(V_m[n], n_ch[n], h);
        m3 = hh_m_ch(V_m[n], m_ch[n], h);
        h3 = hh_h_ch(V_m[n], h_ch[n], h);
        ns3 = (-Inoise[n]*h + dNoise)/tau_cor;
        V_m[n] = V_mem + v3;
        n_ch[n] = n_channel + n3;
        m_ch[n] = m_channel + m3;
        h_ch[n] = h_channel + h3;
        Inoise[n] = Inoise_ + ns3;


        I_syn[n]  = (y[n]*h + I_syn[n])*exp_psc;
        y[n] *= exp_psc;

        v4 = hh_Vm(V_m[n], n_ch[n], m_ch[n], h_ch[n], I_syn[n] + Inoise[n] + I_e[n], h);
        n4 = hh_n_ch(V_m[n], n_ch[n], h);
        m4 = hh_m_ch(V_m[n], m_ch[n], h);
        h4 = hh_h_ch(V_m[n], h_ch[n], h);
        ns4 = (-Inoise[n]*h + dNoise)/tau_cor;

        V_m[n] = V_mem + (v1 + 2.0*(v2 + v3) + v4)/6.0;
        n_ch[n] = n_channel + (n1 + 2.0*(n2 + n3) + n4)/6.0;
        m_ch[n] = m_channel + (m1 + 2.0*(m2 + m3) + m4)/6.0;
        h_ch[n] = h_channel + (h1 + 2.0*(h2 + h3) + h4)/6.0;
        Inoise[n] = Inoise_ + (ns1 + 2.0*(ns2 + ns3) + ns4)/6.0;

        // checking if there's spike on neuron
        if (V_m[n] > V_peak && V_mem > V_m[n] && V_m_last[n] <= V_mem){
            // second condition is necessary in the presence of noise
            if (num_spike_neur[n] == 0 || t - spike_time[Nneur*(num_spike_neur[n] - 1) + n] > 5.0/h){
                spike_time[Nneur*num_spike_neur[n] + n] = t;
                num_spike_neur[n]++;
            }
        }
        V_m_last[n] = V_mem;

        if (t % recInt == 0){
            Vrec[Nneur*t/recInt + n] = V_m[n];
//            Vrec[Nneur*t/recInt + n] = I_syn[n];
        }
    }
}

void set_calc_params(unsigned int Tsim, unsigned int Nneur, unsigned int Ncon, unsigned int recInt, double h){
    hh::Tsim = Tsim;
    hh::Nneur = Nneur;
    hh::Ncon = Ncon;
    hh::recInt = recInt;
    hh::h = h;
    hh::V_m_last = new double[Nneur]();
    for (unsigned int i = 0; i < Nneur; i++){
        hh::V_m_last[i] = 100.0;
    }
    hh::psn_time = new unsigned int[Nneur];
    hh::psn_seed = new unsigned int[Nneur];
    hh::Inoise = new double[Nneur]();
//     hh::num_spike_neur = new unsigned int[Nneur]();
    hh::num_spike_syn = new unsigned int[Ncon]();
}

void set_incom_spikes(unsigned int *times, unsigned int *nums, double *weights, unsigned int MaxNumIncom){
    hh::MAXSZ = MaxNumIncom;
    hh::incSpTimes = times;
    hh::incSpWeights = weights;
    hh::numIncomSpikes = nums;
    hh::numProcessed = new unsigned int[hh::Nneur]();
}

void set_spike_times(unsigned int *spike_time, unsigned int *num_spike_neur, unsigned int sz){
    hh::spike_time = spike_time;
    hh::num_spike_neur = num_spike_neur;
}

void set_conns(double *weight, unsigned int *delay, unsigned int *pre, unsigned int *post){
    hh::weight = weight;
    hh::delay = delay;
    hh::pre_nidx = pre;
    hh::post_nidx = post;
}

void set_neur_vars(double *V_m, double *Vrec, double *n_ch, double *m_ch, double *h_ch){
    hh::V_m = V_m;
    hh::Vrec = Vrec;
    hh::n_ch = n_ch;
    hh::m_ch = m_ch;
    hh::h_ch = h_ch;
}

void set_currents(double *I_e, double *y, double *I_syn, double *rate, double tau_psc, double *d_w_p, unsigned int seed){
    hh::I_e = I_e;
    hh::y = y;
    hh::I_syn = I_syn;
    hh::rate = rate;
    hh::d_w_p = d_w_p;
    init_noise(seed);
    hh::tau_psc = tau_psc;
    hh::exp_psc = exp(-hh::h/tau_psc);
    hh::exp_psc_half = exp(-hh::h*0.5/tau_psc);
}


using namespace hh;

void simulate_cpp(){
    for (unsigned int t = 0; t < Tsim; t++){
        for (unsigned int i = 0; i < Nneur; i++){
            integrate_neurons(t, i);
        }
        for (unsigned int i = 0; i < Ncon; i++){
            integrate_synapses(t, i);
        }
    }
}

void init_noise(unsigned int seed){
    for (unsigned int i = 0; i< Nneur; i++){
        psn_seed[i] = 100000*(seed + i + 1);
        psn_time[i] = 1 + (unsigned int) (-1000.0*log(get_random(psn_seed + i))/(rate[i]*h));
    }
}
