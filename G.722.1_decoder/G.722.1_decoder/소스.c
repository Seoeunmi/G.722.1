#define _crt_secure_no_warnings //fopen 컴파일러에러처리방지
#pragma warning(disable:4996) //fcloseall 컴파일러에러처리방지

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include "mlt.h"
#define N 320
#define pi 3.141592



int main()
{
	FILE *fin, *fout;
	short input,size =0, count = 0, buffer, save, select_categorization = 0;
	short rms_index[14] = { 0, }, differential_rms_index[14] = { 0, };
	int Category[16][14], vector_index[14][10];
	int i = 0, a = 0,n,r;
	int kmax[7] = { 13, 9, 6, 4, 3, 2, 1 };
	int vd[7] = { 2, 2, 2, 4, 4, 5, 5 };
	int vpr[7] = { 10, 10, 10, 5, 5, 4, 4 };
	float stepsize[7] = { pow(2, -1.5), pow(2, -1.0), pow(2, -0.5), pow(2, 0), pow(2, 0.5), pow(2, 1.0), pow(2, 1.5) };
	float deadzone_rounding[7] = { 0.3, 0.33, 0.36, 0.39, 0.42, 0.45, 0.5 };
	int k[14][20];
	float mlt[N] = { 0, }, Q_rms[14];
	float old[N] = { 0, }, u[N], u_old[N / 2] = { 0, };
	short y[N];

	fin = fopen("encoding.raw", "rb");
	fout = fopen("decoding1.raw", "wb");

	srand(time(NULL));

	while (1)
	{
		DECODING(mlt, Q_rms, k, rms_index, &(*fin), differential_rms_index, Category, vector_index, &size);
		IMLT(mlt, u_old, y, u);

		fwrite(y, 2, N, fout);
		if (size == 0) break;

		memmove(u_old, &u[160], sizeof(float)*(N / 2));
	}
	fcloseall();
	return 0;
}