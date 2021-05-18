#define _crt_secure_no_warnings //fopen 컴파일러에러처리방지
#pragma warning(disable:4996) //fcloseall 컴파일러에러처리방지

#include<stdio.h>
#include<math.h>
#include<string.h>
#include "mlt.h"
#define N 320
#define pi 3.141592

int main()
{
	FILE *fin, *fout, *fout2;
	short input[N] = { 0, }, y[N], size;
	float old[N] = { 0, }, mlt[N], u[N], u_old[N / 2] = { 0, }, rms[14], rms_index[14], differential_rms_index[14] = { 0, };
	float Q_rms[14], x;
	int Category[16][14], bit[14], code[14], number_of_region_bits[14];
	int select_categorization, k[14][20], k_correct[14][20];
	int vector_index[14][10];
	short bit_stream, buffer, correct;
	int vpr[7] = { 10, 10, 10, 5, 5, 4, 4 };
	int vd[7] = { 2, 2, 2, 4, 4, 5, 5 };
	float stepsize[7] = { pow(2, -1.5), pow(2, -1.0), pow(2, -0.5), pow(2, 0), pow(2, 0.5), pow(2, 1.0), pow(2, 1.5) };
	float deadzone_rounding[7] = { 0.3, 0.33, 0.36, 0.39, 0.42, 0.45, 0.5 };
	

	fin = fopen("input16k.raw", "rb");
	fout = fopen("mlt2.raw", "wb");
	fout2 = fopen("encoding.raw", "wb");

	for (int h = 0;; h++)
	{
		size = fread(input, 2, N, fin);
		if (size != N)
			for (int i = size; i < N; i++) input[i] = 0;

		MLT(input, old, mlt);
		RMS_value(mlt, rms);
		RMS_INDEX(rms, rms_index, Q_rms);
		DIFF(rms_index, differential_rms_index);
		HUFFMAN(bit, code, differential_rms_index);
		CATEGORY(bit, rms_index, Category);
		SQVH(mlt, rms, Q_rms, Category, number_of_region_bits, bit, &select_categorization, vector_index, k, k_correct);
		BITSTREAM(rms_index, number_of_region_bits, bit, code, &select_categorization,vector_index, Category,k,mlt, &(*fout2));
		IMLT(mlt, u_old, y, u);

		if (h != 0) fwrite(y, 2, N, fout);
		if (size == 0) break;

		memmove(u_old, &u[160], sizeof(float)*(N / 2));
		for (int i = 0; i < N; i++) old[i] = input[i];
	}
	fcloseall();
	return 0;
}