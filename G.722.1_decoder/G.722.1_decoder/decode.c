#include "decode.h"
#include "mlt.h"
#include "HUFFDECODE.h"

int MLT_DECODER_TREE(int a, int b, int c)
{
	switch (a)
	{
	case 0:
		return mlt_decoder_tree_category_0[b][c];
		break;
	case 1:
		return mlt_decoder_tree_category_1[b][c];
		break;
	case 2:
		return mlt_decoder_tree_category_2[b][c];
		break;
	case 3:
		return mlt_decoder_tree_category_3[b][c];
		break;
	case 4:
		return mlt_decoder_tree_category_4[b][c];
		break;
	case 5:
		return mlt_decoder_tree_category_5[b][c];
		break;
	case 6:
		return mlt_decoder_tree_category_6[b][c];
		break;
	default:
		break;
	}
}

void DECODING(float *mlt, float *Q_rms, int k[][20], short *rms_index, FILE *(*fin), short *differential_rms_index, int Category[][14], int vector_index[][10],short *size)
{
	short input, count = 0, buffer, save, select_categorization = 0, amplitude_envelope_bits;
	int i = 0, a = 0, n, r;
	int kmax[7] = { 13, 9, 6, 4, 3, 2, 1 };
	int vd[7] = { 2, 2, 2, 4, 4, 5, 5 };
	int vpr[7] = { 10, 10, 10, 5, 5, 4, 4 };
	float stepsize[7] = { pow(2, -1.5), pow(2, -1.0), pow(2, -0.5), pow(2, 0), pow(2, 0.5), pow(2, 1.0), pow(2, 1.5) };
	float deadzone_rounding[7] = { 0.3, 0.33, 0.36, 0.39, 0.42, 0.45, 0.5 };

	amplitude_envelope_bits = 0;
	count = 0;
	i = 0;
	rms_index[0] = 0;

	*size = fread(&input, sizeof(short), 1, fin);
	count++;

	// rms_index[0] 구하기
	for (i = 0; i < 5; i++)
	{
		rms_index[0] = (rms_index[0] << 1);
		buffer = input;
		buffer = (input << i);
		buffer = buffer & 32768;
		if (buffer == -32768)
			buffer = 1;
		rms_index[0] |= buffer;
	}
	amplitude_envelope_bits = 5;

	// differential_rms_index 구하기
	for (int r = 1; r < 14; r++)
	{
		if (i == 16)
		{
			fread(&input, sizeof(short), 1, fin);
			count++;
			i = 0;
		}

		buffer = input;
		buffer = (input << i);
		buffer = buffer & 32768;
		if (buffer == -32768)
			buffer = 1;
		i++;
		amplitude_envelope_bits++;
		save = differential_region_power_decoder_tree[r][0][buffer];

		while (1)
		{
			if (save >= 0)
			{
				if (i == 16)
				{
					fread(&input, sizeof(short), 1, fin);
					count++;
					i = 0;
				}
				buffer = input;
				buffer = (input << i);
				buffer = buffer & 32768;
				if (buffer == -32768)
					buffer = 1;
				i++;
				amplitude_envelope_bits++;
				save = differential_region_power_decoder_tree[r][save][buffer];
			}
			else
			{
				save = -save;
				save = save - 12;
				differential_rms_index[r] = save;
				break;
			}
		}
	}
	// rms_index , Q_rms 구하기
	RMS_INDEX_Q_RMS(rms_index, differential_rms_index, Q_rms);

	// Determining categorization
	CATEGORY(amplitude_envelope_bits, rms_index, Category);
	select_categorization = 0;
	for (int a = 0; a < 4; a++)
	{
		if (i == 16)
		{
			count++;
			fread(&input, sizeof(short), 1, fin);
			i = 0;
		}
		select_categorization = (select_categorization << 1);
		buffer = input;
		buffer = (input << i);
		buffer = buffer & 32768;
		if (buffer == -32768)
			buffer = 1;
		i++;
		select_categorization |= buffer;
	}

	// Decoding MLT coefficients
	for (r = 0; r < 14; r++)
	{
		if (Category[select_categorization][r] == 7)
		{
			for (int j = 0; j < 20; j++)
			{
				k[r][j] = 0;
				mlt[20 * r + j] = Q_rms[r] * 0.707107;
			}
			continue;
		}

		for (int n = 0; n < vpr[Category[select_categorization][r]]; n++)
		{
			// VECTOR INDEX 구하기
			if (i == 16)
			{
				count++;
				if (count == 31) goto INS;
				fread(&input, sizeof(short), 1, fin);
				i = 0;
			}

			buffer = input;
			buffer = (input << i);
			buffer = buffer & 32768;
			if (buffer == -32768)
				buffer = 1;
			i++;
			save = MLT_DECODER_TREE(Category[select_categorization][r], 0, buffer);

			while (1)
			{
				if (save > 0)
				{
					if (i == 16)
					{
						count++;
						if (count == 31) goto INS;
						fread(&input, sizeof(short), 1, fin);
						i = 0;
					}
					buffer = input;
					buffer = (input << i);
					buffer = buffer & 32768;
					if (buffer == -32768)
						buffer = 1;
					i++;
					save = MLT_DECODER_TREE(Category[select_categorization][r], save, buffer);
				}
				else
				{
					save = -save;
					vector_index[r][n] = save;
					break;
				}
			}

			// K 값 구하기
			for (int j = 0; j < vd[Category[select_categorization][r]]; j++)
			{
				a = (n + 1)*vd[Category[select_categorization][r]] - j - 1;

				if (save == 0)
					k[r][a] = 0;
				else
					k[r][a] = (int)(vector_index[r][n] / pow(kmax[Category[select_categorization][r]] + 1, j)) % (kmax[Category[select_categorization][r]] + 1);
			}

			for (int j = 0; j < vd[Category[select_categorization][r]]; j++)
			{
				a = n*vd[Category[select_categorization][r]] + j;

				// MLT 구하기
				if (k[r][a] == 0)
					mlt[20 * r + a] = 0;
				else
					mlt[20 * r + a] = (k[r][a]) * (stepsize[Category[select_categorization][r]] * Q_rms[r]);

				// MLT 부호
				if (mlt[20 * r + a] != 0)
				{
					if (i == 16)
					{
						count++;
						if (count == 31) goto INS;
						fread(&input, sizeof(short), 1, fin);
						i = 0;
					}
					buffer = input;
					buffer = (input << i);
					buffer = buffer & 32768;
					if (buffer == -32768)
						buffer = 1;
					i++;

					if (buffer == 1)
						mlt[20 * r + a] = mlt[20 * r + a];
					else
						mlt[20 * r + a] = -mlt[20 * r + a];

					// noise fill
					if ((Category[select_categorization][r] == 5) && (k[r][a] == 0))
					{
						mlt[20 * r + a] = Q_rms[r] * 0.176777;

						if ((rand() % 2) == 0)
							mlt[20 * r + a] = mlt[20 * r + a];
						else
							mlt[20 * r + a] = -mlt[20 * r + a];
					}
					if ((Category[select_categorization][r] == 6) && (k[r][a] == 0))
					{
						mlt[20 * r + a] = Q_rms[r] * 0.25;

						if ((rand() % 2) == 0)
							mlt[20 * r + a] = mlt[20 * r + a];
						else
							mlt[20 * r + a] = -mlt[20 * r + a];
					}
				}
			}
		}
	}

	while (count < 30)
	{
		fread(&input, sizeof(short), 1, fin);
		count++;
	}

INS: NOISE_FILL_7(r, mlt, Q_rms, k);
}

void RMS_INDEX_Q_RMS(short *rms_index, short *differential_rms_index, float *Q_rms)
{
	// rms_index 구하기
	for (int r = 1; r < 14; r++)
		rms_index[r] = rms_index[r - 1] + differential_rms_index[r];

	// Q_rms 구하기
	for (int j = 0; j < 14; j++)
	{
		Q_rms[j] = pow(2, (rms_index[j] + 2) / 2.0);
	}
}

void NOISE_FILL_7(int r, float *mlt, float *Q_rms, int *k[][20])
{
	while (r < 14)
	{
		for (int j = 0; j < 20; j++)
		{
			k[r][j] = 0;
			mlt[20 * r + j] = Q_rms[r] * 0.707107;
		}
		r++;
	}
}