#include "mlt.h"
#include "Huffman.h"
#include "MLT_SQVH.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void MLT(short *input, float *old, float *out)
{
	float v[N];
	for (int n = 0; n < 160; n++)
	{
		v[n] = W(159 - n) * old[159 - n] + W(160 + n) * old[160 + n];
		v[n + 160] = W(319 - n) * input[n] - W(n) * input[319 - n];
	}

	for (int m = 0; m < N; m++)   // mlt
	{
		out[m] = 0;
		for (int n = 0; n < N; n++)
			out[m] += sqrt(2.0 / 320)*cos(pi / 320 * (n + 0.5)*(m + 0.5))*v[n];
	}
}

float W(int n)
{
	float w;
	w = sin(pi / 640 * (n + 0.5));
	return w;
}

void RMS(float *mlt, float *rms, float *Q_rms)
{
	float  mlt_rms[N - 40];

	for (int r = 0; r < 14; r++)
	{
		for (int i = 0; i < 20; i++)
		{
			mlt_rms[20 * r + i] = mlt[20 * r + i] / Q_rms[r];
		}
	}

	for (int r = 0; r < 14; r++)
	{
		for (int i = 0; i < 20; i++)
			mlt[20 * r + i] = mlt_rms[20 * r + i] * Q_rms[r];
	}
	for (int r = 14; r < 16; r++)
	{
		for (int i = 0; i < 20; i++)
			mlt[20 * r + i] = 0;
	}
}

void RMS_value(float *mlt, float *rms)
{
	float rms_in;
	for (int r = 0; r < 14; r++)
	{
		rms_in = 0;
		for (int i = 0; i < 20; i++)
		{
			rms_in += mlt[20 * r + i] * mlt[20 * r + i];
		}
		rms[r] = sqrt(rms_in / 20.0);
	}
}

void RMS_INDEX(float *rms, float *rms_index, float *Q_rms)
{
	for (int r = 0; r < 14; r++)
	{
		if (r == 0)
			for (int i = 1; i < 32; i++)
			{
				if (((pow(2.0, (i - 0.5 + 2) / 2.0)) <= rms[r]) && (rms[r] < pow(2.0, (i + 0.5 + 2) / 2.0)))
				{
					rms_index[r] = i;
					Q_rms[r] = pow(2, (i + 2) / 2.0);
				}
			}

		else
			for (int i = -8; i < 32; i++)
			{
				if (((pow(2.0, (i - 0.5 + 2) / 2.0)) <= rms[r]) && (rms[r] < pow(2.0, (i + 0.5 + 2) / 2.0)))
				{
					rms_index[r] = i;
					Q_rms[r] = pow(2, (i + 2) / 2.0);
				}
			}
	}
}

void DIFF(float *rms_index, float *diff_rms)
{
	int j;

	for (int r = 12; r >= 0; r--)
	{
		if (rms_index[r] < rms_index[r + 1] - 11)
			rms_index[r] = rms_index[r + 1] - 11;
	}
	for (int r = 1; r < 14; r++)
	{
		j = rms_index[r] - rms_index[r - 1];
		if (j < -12)
		{
			j = -12;
			rms_index[r] = rms_index[r - 1] + j;
		}
		diff_rms[r] = j;
	}
}

void HUFFMAN(int *bit, int *code, float *diff_rms)
{
	int j;
	for (int r = 0; r < 14; r++)
	{
		j = (int)(diff_rms[r]);
		bit[r] = differential_region_power_bits[r][j + 12];
		code[r] = differential_region_power_codes[r][j + 12];
	}
}

void CATEGORY(short amplitude_envelope_bits, short *rms_index, int Category[][14])
{
	/* STEP 1 & 2 */
	int initial_categorization[14], max_category[14], min_category[14], temp_category_balances[32];
	int offset = -32, delta = 32;
	int expected_bits, max_expected_bits, min_expected_bits, max_pointer, min_pointer, categorization_count;
	int number_of_available_bits, estimated_number_of_available_bits;
	int expected_bits_table[8] = { 52, 47, 43, 37, 29, 22, 16, 0 };
	int sum_bit = 0;
	int temp, temp_region;

	/* STEP 0 */
	number_of_available_bits = 480 - amplitude_envelope_bits - 4;
	
	if (number_of_available_bits > 320)
		estimated_number_of_available_bits = 320 + ((number_of_available_bits - 320)* 5.0 / 8);
	else estimated_number_of_available_bits = number_of_available_bits;

	while (delta > 0)
	{
		/* STEP 3 & 4 */
		for (int r = 0; r < 14; r++)
		{
			initial_categorization[r] = (offset + delta - rms_index[r]) / 2.0;
			if (initial_categorization[r] < 0) initial_categorization[r] = 0;
			if (initial_categorization[r] > 7) initial_categorization[r] = 7;
		}

		/* STEP 5 */
		expected_bits = 0;
		for (int r = 0; r < 14; r++)
		{
			expected_bits += expected_bits_table[initial_categorization[r]];
		}

		/* STEP 6 */
		if (expected_bits >= estimated_number_of_available_bits - 32)
			offset = offset + delta;

		/* STEP 7 */
		delta = delta / 2;
	} /* STEP 8*/

	/* STEP 9 & 10 */
	for (int r = 0; r < 14; r++)
	{
		initial_categorization[r] = (offset - rms_index[r]) / 2;
		if (initial_categorization[r] < 0) initial_categorization[r] = 0;
		if (initial_categorization[r] > 7) initial_categorization[r] = 7;
	}

	/* STEP 11 */
	expected_bits = 0;
	for (int r = 0; r < 14; r++)
	{
		expected_bits += expected_bits_table[initial_categorization[r]];
	}

	/* STEP 12 */
	for (int r = 0; r < 14; r++)
	{
		max_category[r] = initial_categorization[r];
		min_category[r] = initial_categorization[r];

		max_expected_bits = expected_bits;
		min_expected_bits = expected_bits;

		max_pointer = 16;
		min_pointer = 16;
		categorization_count = 1;
	}

	while (categorization_count < 16)
	{
		/* STEP 13 */
		if ((max_expected_bits + min_expected_bits) <= 2 * estimated_number_of_available_bits)
		{
			/* STEP 14 */
			for (int r = 0; r < 14; r++)
			{
				if (max_category[r] > 0)
				{
					temp = (offset - rms_index[r] - 2 * max_category[r]);
					temp_region = r;
					break;
				}
			}

			for (int r = temp_region + 1; r < 14; r++)
			{
				if (max_category[r] > 0)
				{
					if (temp > (offset - rms_index[r] - 2 * max_category[r]))
					{
						temp = (offset - rms_index[r] - 2 * max_category[r]);
						temp_region = r;
					}
				}
			}
			/* STEP 15 */
			max_pointer = max_pointer - 1;
			temp_category_balances[max_pointer] = temp_region;
			max_expected_bits = max_expected_bits - expected_bits_table[max_category[temp_region]];
			max_category[temp_region] = max_category[temp_region] - 1;
			max_expected_bits = max_expected_bits + expected_bits_table[max_category[temp_region]];
		}
		else
		{
			/* STEP 16 */
			for (int r = 0; r < 14; r++)
			{
				if (min_category[r] < 7)
				{
					temp = (offset - rms_index[r] - 2 * min_category[r]);
					temp_region = r;
					break;
				}
			}

			for (int r = temp_region + 1; r < 14; r++)
			{
				if (min_category[r] < 7)
				{
					if (temp <= (offset - rms_index[r] - 2 * min_category[r]))
					{
						temp = (offset - rms_index[r] - 2 * min_category[r]);
						temp_region = r;
					}
				}
			}
			/* STEP 17 */
			temp_category_balances[min_pointer] = temp_region;
			min_pointer = min_pointer + 1;
			min_expected_bits = min_expected_bits - expected_bits_table[min_category[temp_region]];
			min_category[temp_region] = min_category[temp_region] + 1;
			min_expected_bits = min_expected_bits + expected_bits_table[min_category[temp_region]];
		}
		categorization_count = categorization_count + 1;
	} /* STEP 18 */

	/* STEP 19 */
	for (int r = 0; r < 14; r++)
	{
		Category[0][r] = max_category[r];
	}

	for (int n = 1; n < 16; n++)
	{
		/* STEP 20 & 21 */
		for (int r = 0; r < 14; r++)
		{
			Category[n][r] = Category[n - 1][r];
		}

		/* STEP 22 */
		Category[n][temp_category_balances[max_pointer]] = (Category[n][temp_category_balances[max_pointer]] + 1);

		/* STEP 23 */
		max_pointer = max_pointer + 1;
	}
}

void SQVH(float *mlt, float *rms, float *Q_rms, int Category[][14], int *number_of_region_bits, int *bit, int *select_categorization, int vector_index[][10], int k[][20])
{
	float stepsize[7] = { pow(2, -1.5), pow(2, -1.0), pow(2, -0.5), pow(2, 0), pow(2, 0.5), pow(2, 1.0), pow(2, 1.5) };
	float deadzone_rounding[7] = { 0.3, 0.33, 0.36, 0.39, 0.42, 0.45, 0.5 };
	int kmax[7] = { 13, 9, 6, 4, 3, 2, 1 };
	float x;
	int vd[7] = { 2, 2, 2, 4, 4, 5, 5 };
	int vpr[7] = { 10, 10, 10, 5, 5, 4, 4 };

	int bit_num = 4 + 5;
	int mlt_num = 0;

	for (int m = 0; m < 16; m++)
	{
		for (int r = 0; r < 14; r++)
		{
			x = 1.0 / (stepsize[Category[m][r]] * Q_rms[r]);
			if (Category[m][r] == 7)
				for (int i = 0; i < 20; i++)
				{
					k[r][i] = 0;
				}

			for (int i = 0; i < 20; i++)
			{
				k[r][i] = MIN((int)(x * fabs(mlt[20 * r + i]) + deadzone_rounding[Category[m][r]]), kmax[Category[m][r]]);
			}
		}

		for (int r = 0; r < 14; r++)
		{
			for (int n = 0; n < vpr[Category[m][r]]; n++)
			{
				vector_index[r][n] = 0;
				for (int j = 0; j < vd[Category[m][r]]; j++)
				{
					vector_index[r][n] += k[r][n * vd[Category[m][r]] + j] * pow(kmax[Category[m][r]] + 1, (vd[Category[m][r]] - (j + 1)));
				}
			}
		}

		for (int r = 0; r < 14; r++)
		{
			number_of_region_bits[r] = 0;
			for (int i = 0; i < 20; i++)
			{
				if (k[r][i] != 0) number_of_region_bits[r] += 1;
			}

			for (int n = 0; n < vpr[Category[m][r]]; n++)
			{
				number_of_region_bits[r] += MLT_SQVH_BITCOUNT(Category[m][r], vector_index[r][n]);
			}
		}

		/*Select_Categorization*/
		bit_num = 9;
		for (int i = 1; i < 14; i++)
		{
			bit_num += bit[i];
		}

		mlt_num = 0;
		for (int r = 0; r < 14; r++)
		{
			mlt_num += number_of_region_bits[r];
		}

		if ((480 - bit_num) >= mlt_num)
		{
			*select_categorization = m;
			break;
		}
		else *select_categorization = 15;
	}
}

int MLT_SQVH_BITCOUNT(int a, int b)
{
	switch (a)
	{
	case 0:
		return mlt_sqvh_bitcount_category_0[b];
		break;
	case 1:
		return mlt_sqvh_bitcount_category_1[b];
		break;
	case 2:
		return mlt_sqvh_bitcount_category_2[b];
		break;
	case 3:
		return mlt_sqvh_bitcount_category_3[b];
		break;
	case 4:
		return mlt_sqvh_bitcount_category_4[b];
		break;
	case 5:
		return mlt_sqvh_bitcount_category_5[b];
		break;
	case 6:
		return mlt_sqvh_bitcount_category_6[b];
		break;
	default:
		break;
	}
}

int MLT_SQVH_CODE(int a, int b)
{
	switch (a)
	{
	case 0:
		return  mlt_sqvh_code_category_0[b];
		break;
	case 1:
		return  mlt_sqvh_code_category_1[b];
		break;
	case 2:
		return  mlt_sqvh_code_category_2[b];
		break;
	case 3:
		return  mlt_sqvh_code_category_3[b];
		break;
	case 4:
		return  mlt_sqvh_code_category_4[b];
		break;
	case 5:
		return  mlt_sqvh_code_category_5[b];
		break;
	case 6:
		return  mlt_sqvh_code_category_6[b];
		break;
	default:
		break;
	}
}


void BITSTREAM(float *rms_index, int *number_of_region_bits, int *bit, int *code, int select_categorization, int vector_index[][10], int Category[][14], int k[][20])
{
	FILE *fout2;
	short bit_stream = 0;
	short count = 0;
	short buffer, save;
	int  n = 1, vector_bit, k_bit_dec, k_bit_bin;

	int vpr[7] = { 10, 10, 10, 5, 5, 4, 4 };

	fout2 = fopen("encoding.raw", "wb");


	//amplitude envelope bits , 마지막 카운트 값 = 11 , bit_stream = 00000 00000111001
	for (int i = 0; i < 14; i++)
	{
		n = 1;
		if (i == 0)
		{
			bit_stream = rms_index[0];
			count = 5;
			bit_stream = (bit_stream << 1);
		}
		else
		{
			for (int j = 1; j <= bit[i]; j++)
			{
				buffer = (code[i] >> (bit[i] - j));
				n++;
				buffer = buffer & 1;
				bit_stream |= buffer;
				count++;
				if (count == 16)
				{
					fwrite(&bit_stream, sizeof(short), 1, fout2);
					count = 0;
					bit_stream = 0;
					for (int a = n; a <= bit[i]; a++)
					{
						buffer = (code[i] >> (bit[i] - a));
						buffer = buffer & 1;
						bit_stream |= buffer;
						count++;
						bit_stream = (bit_stream << 1);
					}
					break;
				}
				bit_stream = (bit_stream << 1);
			}
		}
	}

	// Categorization control bits 
	n = 1;
	for (int j = 1; j <= 4; j++)
	{
		buffer = (select_categorization >> 4 - j);
		n++;
		buffer = buffer & 1;
		bit_stream |= buffer;
		count++;

		if (count == 16)
		{
			fwrite(&bit_stream, sizeof(short), 1, fout2);
			count = 0;
			bit_stream = 0;
			for (int a = n; a <= 4; a++)
			{
				buffer = (select_categorization >> 4 - a);
				buffer = buffer & 1;
				bit_stream |= buffer;
				count++;
				bit_stream = (bit_stream << 1);
			}
			break;
		}
		bit_stream = (bit_stream << 1);
	}

	// MLT coefficients bits
	for (int r = 0; r < 14; r++)
	{
		for (int m = 0; m < vpr[Category[select_categorization][r]]; m++)
		{
			vector_bit = MLT_SQVH_BITCOUNT(Category[select_categorization][r], vector_index[r][m]);

			n = 1;
			for (int j = 1; j <= vector_bit; j++)
			{
				buffer = (vector_index[r][m] >> vector_bit - j);
				n++;
				buffer = buffer & 1;
				bit_stream |= buffer;
				count++;

				if (count == 16)
				{
					fwrite(&bit_stream, sizeof(short), 1, fout2);
					count = 0;
					bit_stream = 0;
					for (int a = n; a <= vector_bit; a++)
					{
						buffer = (vector_index[r][m] >> vector_bit - a);
						buffer = buffer & 1;
						bit_stream |= buffer;
						count++;
						bit_stream = (bit_stream << 1);
					}
					break;
				}
				bit_stream = (bit_stream << 1);
			}
		}

		/* k가 0이 아닐때의 부호 보내기*/
		k_bit_dec = 0;
		for (int i = 0; i < 20; i++)
		{
			if (k[r][i] != 0) k_bit_dec += 1;
		}
		k_bit_bin = K_BIT_NUM(k_bit_dec);

		n = 1;
		for (int j = 1; j <= k_bit_bin; j++)
		{
			buffer = (k_bit_dec >> k_bit_bin - j);
			n++;
			buffer = buffer & 1;
			bit_stream |= buffer;
			count++;

			if (count == 16)
			{
				fwrite(&bit_stream, sizeof(short), 1, fout2);
				count = 0;
				bit_stream = 0;
				for (int a = n; a <= k_bit_bin; a++)
				{
					buffer = (k_bit_dec >> k_bit_bin - a);
					buffer = buffer & 1;
					bit_stream |= buffer;
					count++;
					bit_stream = (bit_stream << 1);
				}
				break;
			}
			bit_stream = (bit_stream << 1);
		}
	}
}

int K_BIT_NUM(int num)
{
	char binary[10];
	char position = 0;
	int k_count = 0;

	while (1)
	{
		binary[position] = num % 2;
		num = num / 2;
		position++;
		k_count++;

		if (num == 0) break;
	}
	return k_count;
}

void IMLT(float *mlt, float *u_old, short *out, float *u)
{
	for (int n = 0; n < N; n++)   //  imlt
	{
		u[n] = 0;
		for (int m = 0; m < N; m++)
			u[n] += sqrt(2.0 / 320)*cos(pi / 320 * (m + 0.5)*(n + 0.5))*mlt[m];
	}
	for (int n = 0; n < 160; n++)
	{
		out[n] = W(n) * u[159 - n] + W(319 - n) * u_old[n];
		out[n + 160] = W(160 + n) * u[n] - W(159 - n) * u_old[159 - n];
	}
}