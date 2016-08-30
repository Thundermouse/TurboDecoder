#include<stdio.h>
#include <stdlib.h>
#include<time.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vml.h> 
#include <windows.h>
#include <xmmintrin.h>
#include <immintrin.h> //AVX
#include "TurboCode.h"
#include "CYHTurbo.h"
#include "Tools.h"
#define SEED  7777777
extern int Mem;
int main()
{
	/*****************************����ռ�***************************************/
	/************��Ҫ�˹��趨�Ĳ���************/
	int Gate = 0;
	double SNR_H;
	double SNR_L;
	double SNR_Pace;
	int iterCNum = 1;
	int Cycle = 1;
	int	Len = 640;
	int			CodeLen = 8400;
	short		*msg;
	short		*dec;
	int			*msg32;
	short		*code;
	double		*llr_in;
	double		*llr_out, *llr_out1;
	double		*store_llr;
	double		*awgn;
	short		*llr1_sign, *llr2_sign;
	int			*interleave_table;
	int			*ratematch_table;
	double		ber;
	int			errors = 0;
	int			i, k;
	clock_t	turbodec_clocks, turbodec_begin, turbodec_end;
	FILE		*fp;
	VSLStreamStatePtr stream, strGauss;
	auto		status = 0;
	double		SNR, SNRLinear;
	/********************Ϊ����������ռ�************/
	msg = (short*)_mm_malloc(Len*sizeof(short), sizeof(__m256i));     //_aligned_malloc(a,b)Ϊ��b��������ڴ�
	msg32 = (int*)_mm_malloc(Len*sizeof(int), sizeof(__m256i));
	dec = (short*)_mm_malloc(Len*sizeof(short), sizeof(__m256i));
	code = (short*)_mm_malloc(CodeLen*sizeof(short), sizeof(__m256i));
	llr_in = (double*)_mm_malloc(CodeLen*sizeof(double), sizeof(__m256i));
	llr_out = (double*)_mm_malloc(CodeLen*sizeof(double), sizeof(__m256i));
	llr_out1 = (double*)_mm_malloc(CodeLen*sizeof(double), sizeof(__m256i));
	store_llr = (double*)_mm_malloc((Len + Mem)*sizeof(double), sizeof(__m256i));
	awgn = (double*)_mm_malloc(CodeLen*sizeof(double), sizeof(__m256i));

	llr1_sign = (short*)_mm_malloc(CodeLen*sizeof(short), sizeof(__m256i));
	llr2_sign = (short*)_mm_malloc(CodeLen*sizeof(short), sizeof(__m256i));

	//memset(msg, 0, Len*sizeof(short));

	// Initialize the stream
	status = vslNewStream(&stream, VSL_BRNG_SFMT19937, SEED);
	status = vslNewStream(&strGauss, VSL_BRNG_SFMT19937, SEED);
	interleave_table = (int*)_mm_malloc(Len*sizeof(int), sizeof(__m256i));
	ratematch_table = (int*)_mm_malloc(CodeLen*sizeof(int), sizeof(__m256i));
	//ratematch_table = (int*)_mm_malloc(3*(Len+4)*sizeof(int), sizeof(__m256i));
	while (TRUE)
	{
		printf("�Ƿ����Turbo����棿1/0\n");
		scanf("%d", &Gate);
		if (Gate == 0) break;
		/*****************************��һ���������*********************************/
		double SNR_H;
		double SNR_L;
		double SNR_Pace;
		int iterCNum = 1;
		int Cycle = 1;
		int	Len = 640;
		printf("���������������\n");
		scanf("%d", &iterCNum);
		printf("���������SNR(dB)�����SNR(dB)���м��ÿո�ֿ���\n");
		scanf("%lf%lf", &SNR_L, &SNR_H);
		printf("������SNR��������[0.2]��\n");
		scanf("%lf", &SNR_Pace);
		printf("������鳤�ȣ�\n");
		scanf("%d", &Len);
		printf("������ÿ��SNRѭ��������\n");
		scanf("%d", &Cycle);
		/******************��֯���趨*****************
		1.�����֯����RandomInterleavetableGenerator(Len, interleave_table);
		2.Ĭ�Ͻ�֯����fp = fopen("interleave_table.dat", "rb");
		              fread(interleave_table, sizeof(int), Len, fp);
		              fclose(fp);
		*/
		RandomInterleavetableGenerator(Len, interleave_table);
		/***************����ƥ���趨*****************/
		
		fp = fopen("ratematch_table.dat", "rb");
		fread(ratematch_table, sizeof(int), CodeLen, fp);
		fclose(fp);
		
		//RatematchtableGenerator(Len + 4, ratematch_table);
		/**********************�������ѭ��************************/
		for (SNR=SNR_L; SNR<=SNR_H;SNR=SNR+SNR_Pace)
		{ 
			printf("��������������������������������������������������������\n");
			turbodec_clocks = 0;
			turbodec_begin = 0;
			turbodec_end = 0;
			SNRLinear = pow(10, SNR / 10); //?
			errors = 0;
			for (k = 0; k < Cycle; k++)
			{
				// ��������������
				status = viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, Len, msg32, 0.5);

				for (i = 0; i < Len; i++)
				{
					msg[i] = (short)msg32[i];
				}

				turbo_encode_rm(msg, code, Len, CodeLen, interleave_table, ratematch_table);

				status = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, strGauss, CodeLen, awgn, 0, sqrt(SNRLinear));

				for (i = 0; i < CodeLen; i++)
				{
					llr_in[i] = 2 * code[i] - 1 + awgn[i];
					llr_in[i] = 2 * llr_in[i] / SNRLinear;
				}

				for (i = 0; i < (Len + Mem); i++)
				{
					store_llr[i] = 0;
				}
				;
				turbodec_begin = clock();
				//turbo_decode_rm(dec, llr_out1, llr_in, store_llr, Len, CodeLen, 4, interleave_table, ratematch_table);
				Decoder(dec, llr_out1, llr_in, Len, CodeLen, iterCNum, interleave_table, ratematch_table);
				//turbo_decode_rm_avx(dec, llr_out, llr_in, store_llr, Len, CodeLen, 4, interleave_table, ratematch_table);
				turbodec_end = clock();
				turbodec_clocks += (turbodec_end - turbodec_begin);
				for (i = 0; i < Len; i++)
				{
					if (msg[i] != dec[i])
					{
						errors++;
					}
				}
			}//Cycle
			ber = (double)errors / (Cycle*Len);
			printf("����������������������ǰSNR:%0.2f(dB)��������������������\n", SNR);
			printf("���������%d\n�������%0.4lf\n", errors, ber);
			printf("ִ��ʱ��%0.4lf��\n", (double)turbodec_clocks / CLOCKS_PER_SEC);
		}//SNR_H
		printf("��������������������������������������������������������\n");
		printf("|                      ����ִ�����                    |\n");
		printf("��������������������������������������������������������\n\n");
	}
	/*********************************��ʱ��������*****************************************/
	//RandomInterleavetableGenerator(Len, interleave_table);
	//printf("%s\n", IsTable(interleave_table, Len) ?  "��Table":"����Table");
	/*
	for (int L = 40; L < 10000; ++L)
	{
		int*table1 = (int*)malloc(L * sizeof(int));
		int*table2 = (int*)malloc(L * sizeof(int));
		SubBlockInterleavetableGenerator(L, table1, 0);
		SubBlockInterleavetableGenerator(L, table2, 1);
		printf("%d", CmpVector(table1, table2, L));
		free(table1);
		free(table2);
		if (L % 100 == 0) printf("\n");
	}
	*/     //��֯�����
	/*
	int K = 640;
	int*table = (int*)malloc(3 * (K+4)*sizeof(int));
	int*streamin = (int*)malloc(3 * (K + 4)*sizeof(int));
	int*streamout = (int*)malloc(3 * (K + 4)*sizeof(int));
	int*s = (int*)malloc((K + 4)*sizeof(int));
	int*p1 = (int*)malloc((K + 4)*sizeof(int));
	int*p2 = (int*)malloc((K + 4)*sizeof(int));
	for (int i = 0; i < K + 4; ++i)
	{
		s[i] = i;
		p1[i] = i;
		p2[i] = i;
	}
	RatematchtableGenerator(K + 4, table);
	Ratematching(s, p1, p2, table, streamin,K+4);
	Ratedematching(streamin, table, streamout, K + 4);
	for (int i = 0; i < 3 * (K + 4); ++i)
	{
		printf("%d ", streamout[i]);
		if (i % 40 == 0 && i != 0)printf("\n");
	}
	free(table);
	free(s);
	free(p1);
	free(p2);
	free(streamin);
	free(streamout);
	*/
	/**************************************************************************************/
	/********************�ͷſռ�*********************/
	status = vslDeleteStream(&stream);
	status = vslDeleteStream(&strGauss);
	_mm_free(interleave_table);
	_mm_free(ratematch_table);
	_mm_free(msg32);
	_mm_free(msg);
	_mm_free(code);
	_mm_free(dec);
	_mm_free(llr_in);
	_mm_free(llr_out);
	_mm_free(llr_out1);
	_mm_free(store_llr);
	_mm_free(awgn);
	_mm_free(llr1_sign);
	_mm_free(llr2_sign);
	return 0;
}