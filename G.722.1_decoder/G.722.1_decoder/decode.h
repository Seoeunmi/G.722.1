#ifndef DECODE_H
#define DECODE_H

#include<stdio.h>
#include<math.h>
#include<string.h>

#define N 320
#define pi 3.141592

void DECODING(float *mlt, float *Q_rms, int k[][20], short *rms_index, FILE *(*fin), short *differential_rms_index, int Category[][14], int vector_index[][10], short *size);

void RMS_INDEX_Q_RMS(short *rms_index, short *differential_rms_index, float *Q_rms);

void NOISE_FILL_7(int r, float *mlt, float *Q_rms, int *k[][20]);

void AMPLITUDE(int *i, short *rms_index, short *count, short *size, FILE *(*fin), short *differential_rms_index);


#endif