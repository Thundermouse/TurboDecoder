#ifndef _TURBOCODE_H_
#define _TURBOCODE_H_

void turbo_encode_rm(short* msg,
	short* code,
	int Len,
	int CodeLen,
	int* interleave_table,
	int* ratematch_table);

void turbo_decode_rm(short* dec,
	double* llr_out,
	double* llr_in,
	double* store_llr,
	int Len,
	int CodeLen,
	int iterCNum,
	int* interleave_table,
	int* ratematch_table);

void turbo_decode_rm_fix(short* dec,
	double* llr_out,
	double* llr_in,
	double* store_llr,
	int Len,
	int CodeLen,
	int iterCNum,
	int* interleave_table,
	int* ratematch_table);

void turbo_decode_rm_avx(short* dec,
	double* llr_out,
	double* llr_in,
	double* store_llr,
	int Len,
	int CodeLen,
	int iterCNum,
	int* interleave_table,
	int* ratematch_table);

void turbo_decode_rm_avx2(short* dec,
	double* llr_out,
	double* llr_in,
	double* store_llr,
	int Len,
	int CodeLen,
	int iterCNum,
	int* interleave_table,
	int* ratematch_table);

#endif