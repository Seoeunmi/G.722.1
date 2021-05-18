#ifndef MLT_H
#define MLT_H

#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 320
#define pi 3.141592

void MLT(short *input, float *old, float *out);

float W(int n);

void RMS(float *mlt, float *rms, float *Q_rms);

float Q_RMS(int r, float rms);

void RMS_value(float *mlt, float *rms);

void RMS_INDEX(float *rms, float *rms_index, float *Q_rms);

void DIFF(float *rms_index, float *diff_rms);

void HUFFMAN(int *bit, int *code, float *diff_rms);

void CATEGORY(float *bit, float *rms_index, int Category[][14]);

void SQVH(float *mlt, float *rms, float *Q_rms, int Category[][14], int *number_of_region_bits, int *bit, int *select_categorization, int vector_index[][10], int k[][20]);

int MLT_SQVH_BITCOUNT(int a, int b);

int MLT_SQVH_CODE(int a, int b);

void BITSTREAM(float *rms_index, int *number_of_region_bits, int *bit, int *code, int *select_categorization, int vector_index[][10], int Category[][14], int k[][20], float *mlt, FILE *(*fout2));

void IMLT(float *mlt, float *u_old, short *out, float *u);
#endif 