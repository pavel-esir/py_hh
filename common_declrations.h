/*
 * common_declrations.h
 *
 *  Created on: Jul 21, 2016
 *      Author: pavel
 */

#ifndef COMMON_DECLRATIONS_H_
#define COMMON_DECLRATIONS_H_


struct NeurVars{
    float *V;
    float *V_last;
    float *n;
    float *m;
    float *h;
    float *Ie;
    float *Isyn;
    float *y;
    float *Inoise;
    float *weight_p;
};

struct RecordVars{
    float *V;
    unsigned int interval;
};

struct IncSpikes{
    unsigned int MAXSZ; // [maximum(over all neurons) number of incoming spikes
    unsigned int *times; // size of times is  (MaxNumIncom, Nneur)
    float *weights;
    unsigned int *nums;
    unsigned int *numProcessed;
};

#endif /* COMMON_DECLRATIONS_H_ */
