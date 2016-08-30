#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <xmmintrin.h>
#include <immintrin.h>
#include "turboCode.h"

const short OUTM[16] = { 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1 };
const int NSM[16] = { 0, 4, 4, 0, 5, 1, 1, 5, 2, 6, 6, 2, 7, 3, 3, 7 };
const int LSM[16] = { 0, 1, 3, 2, 4, 5, 7, 6, 1, 0, 2, 3, 5, 4, 6, 7 };
const short TAILM[8] = { 0, 1, 1, 0, 0, 1, 1, 0 };
const int Mem = 3;

void interleave_short(short* msg, int Len, int *interleave_table)
{
	short	*tmp;
	int		i;

	tmp = (short*)malloc(Len*sizeof(short));

	for (i = 0; i<Len; i++)
	{
		tmp[i] = msg[interleave_table[i]];
	}

	for (i = 0; i<Len; i++)
	{
		msg[i] = tmp[i];
	}

	free(tmp);
}

void deinterleave_short(short* msg, int Len, int *interleave_table)
{
	short	*tmp;
	int		i;

	tmp = (short*)malloc(Len*sizeof(short));

	for (i = 0; i<Len; i++)
	{
		tmp[interleave_table[i]] = msg[i];
	}

	for (i = 0; i<Len; i++)
	{
		msg[i] = tmp[i];
	}

	free(tmp);
}

void	rm_interleave(int CodeLen, short *input, short *output, int *inter_index)
{
	int i;

	for (i = 0; i<CodeLen; i++)
	{
		output[i] = input[inter_index[i]];
	}
}


void rsc_encode(short* msg, short* parity, int Len)
{
	int		i;
	int		state = 0;

	for (i = 0; i<Len - Mem; i++)
	{
		parity[i] = OUTM[state * 2 + msg[i]];
		state = NSM[state * 2 + msg[i]];
	}

	for (i = Len - Mem; i<Len; i++)
	{
		msg[i] = TAILM[state];
		parity[i] = OUTM[state * 2 + msg[i]];
		state = NSM[state * 2 + msg[i]];
	}
}

void turbo_encode_rm(	short* msg,
						short* code,
						int Len,
						int CodeLen,
						int* interleave_table,
						int* ratematch_table)
{
	short		*tmpmsg;
	short		*parity1;
	short		*parity2;
	short		*tmptailmsg;
	short		*tmptailp1;
	short		*tmptailp2;
	short		*code_tmp;

	int			i;

	tmpmsg = (short*)malloc((Len + 4)*sizeof(short));
	parity1 = (short*)malloc((Len + 4)*sizeof(short));
	parity2 = (short*)malloc((Len + 4)*sizeof(short));
	code_tmp = (short*)malloc(3 * (Len + 4)*sizeof(short));
	tmptailmsg = (short*)malloc(4 * sizeof(short));
	tmptailp1 = (short*)malloc(4 * sizeof(short));
	tmptailp2 = (short*)malloc(4 * sizeof(short));

	for (i = 0; i<Len; i++)
	{
		tmpmsg[i] = msg[i];
	}

	for (i = Len; i<Len + 4; i++)
	{
		tmpmsg[i] = 0;
	}

	rsc_encode(tmpmsg, parity1, Len + Mem);
	
	// TailBits
	// msg:		x0,p1,x0',p1'
	// parity1:	p0,x2,p0',x2'
	// parity2:	x1,p2,x1',p2'
	
	tmptailmsg[0] = tmpmsg[Len];
	tmptailp1[0] = parity1[Len];
	tmptailp2[0] = tmpmsg[Len + 1];
	tmptailmsg[1] = parity1[Len + 1];
	tmptailp1[1] = tmpmsg[Len + 2];
	tmptailp2[1] = parity1[Len + 2];

	interleave_short(tmpmsg, Len, interleave_table);

	rsc_encode(tmpmsg, parity2, Len + Mem);

	//transmitted bits for trellis termination encoder 2
	tmptailmsg[2] = tmpmsg[Len];
	tmptailp1[2] = parity2[Len];
	tmptailp2[2] = tmpmsg[Len + 1];
	tmptailmsg[3] = parity2[Len + 1];
	tmptailp1[3] = tmpmsg[Len + 2];
	tmptailp2[3] = parity2[Len + 2];

	deinterleave_short(tmpmsg, Len, interleave_table);

	// adding tailbits to code bits
	for (i = 0; i<4; i++)
	{
		tmpmsg[i + Len] = tmptailmsg[i];
		parity1[i + Len] = tmptailp1[i];
		parity2[i + Len] = tmptailp2[i];
	}
	//Ratematching(tmpmsg, parity1, parity2, ratematch_table, code, Len+4);
	
	for (i = 0; i<Len + 4; i++)
	{
		code_tmp[3 * i] = tmpmsg[i];
		code_tmp[3 * i + 1] = parity1[i];
		code_tmp[3 * i + 2] = parity2[i];
	}

	rm_interleave(CodeLen, code_tmp, code, ratematch_table);
	

	free(tmpmsg);
	free(parity1);
	free(parity2);
	free(code_tmp);
	free(tmptailmsg);
	free(tmptailp1);
	free(tmptailp2);
}

void interleave_double(double* msg, int Len, int *interleave_table)
{
	double	*tmp;
	int		i;

	tmp = (double*) malloc(Len*sizeof(double));

	for(i=0; i<Len; i++)
	{
		tmp[i] = msg[interleave_table[i]];
	}

	for(i=0; i<Len; i++)
	{
		msg[i] = tmp[i];
	}

	free(tmp);
}

void deinterleave_double(double* msg, int Len, int *interleave_table)
{
	double	*tmp;
	int		i;

	tmp = (double*) malloc(Len*sizeof(double));

	for(i=0; i<Len; i++)
	{
		tmp[interleave_table[i]] = msg[i];
	}

	for(i=0; i<Len; i++)
	{
		msg[i] = tmp[i];
	}

	free(tmp);
}

//rate dematching; rate dematch with table inter_index
void	rate_dematching(int Len, int CodeLen, double *input, double *output, int *inter_index)
{
	int i;
	int j;
	for(i=0; i<3*Len; i++)
	{
		output[i] = 0;
	}
	for(i=0; i<CodeLen; i++)
	{
		
		output[inter_index[i]] += input[i];
	}
}

void	rm_interleave_doub(int Len, int CodeLen, double *input, double *output, int *inter_index)
{
	int i;
	
	for(i=0; i<CodeLen; i++)
	{
		output[i] = input[inter_index[i]];
	}
}

double maxo(double t0, double t1)
{
	//if (t0 >= t1)
	//{
	//	return t0 + log(1+exp(-fabs(t0-t1)));
	//}
	//else
	//{
	//	return t1 + log(1+exp(-fabs(t0-t1)));
	//}
	if (t0 >= t1)
	{
		return t0;
	}
	else
	{
		return t1;
	}
}


double ACSO(double a0, double g0, double a1, double g1)
{
	double	t0, t1;

	t0 = a0 + g0;
	t1 = a1 + g1;
	
	return maxo(t0, t1);
}

void turbo_decode_rm(	short* dec, 
						double* llr_out, 
						double* llr_in, 
						double* store_llr, 
						int Len, 
						int CodeLen,
						int iterCNum, 
						int* interleave_table,
						int* ratematch_table)
{
	double		*llr_in_rm, *llr_out_rm, *llr_s1, *llr_s2, *llr_p1, *llr_p2, *llr_es1, *llr_es2, *llr_ep1, *llr_ep2;
	double		*alfa1, *beta1, *gs1, *gp1, *alfa2, *beta2, *gs2, *gp2;

	int			i, s, iter_idx;
	double		llr_s_tail[4], llr_p1_tail[4], llr_p2_tail[4];
	double		inf = 128;
	double		limit = 8;
	int			StateNum = 8;

	llr_in_rm	=	(double*) malloc(3*(Len+4)*sizeof(double));
	llr_out_rm	=	(double*) malloc(3*(Len+4)*sizeof(double));
	llr_s1		=	(double*) malloc((Len+4)*sizeof(double));
	llr_s2		=	(double*) malloc((Len+4)*sizeof(double));
	llr_p1		=	(double*) malloc((Len+4)*sizeof(double));
	llr_p2		=	(double*) malloc((Len+4)*sizeof(double));
	llr_es1		=	(double*) malloc((Len+Mem)*sizeof(double));
	llr_es2		=	(double*) malloc((Len+Mem)*sizeof(double));
	llr_ep1		=	(double*) malloc((Len+Mem)*sizeof(double));
	llr_ep2		=	(double*) malloc((Len+Mem)*sizeof(double));

	alfa1 = (double*) malloc((Len+Mem)*StateNum*sizeof(double));
	beta1 = (double*) malloc((Len+Mem)*StateNum*sizeof(double));
	alfa2 = (double*) malloc((Len+Mem)*StateNum*sizeof(double));
	beta2 = (double*) malloc((Len+Mem)*StateNum*sizeof(double));
	gs1 = (double*) malloc((Len+Mem)*sizeof(double));
	gp1 = (double*) malloc((Len+Mem)*sizeof(double));
	gs2 = (double*) malloc((Len+Mem)*sizeof(double));
	gp2 = (double*) malloc((Len+Mem)*sizeof(double));

	//rate match....
	rate_dematching(Len+4, CodeLen, llr_in, llr_in_rm, ratematch_table);

	for(i=0; i<Len+4; i++)
	{
		llr_s1[i] = llr_in_rm[3*i];
		llr_s2[i] = llr_in_rm[3*i];
		llr_p1[i] = llr_in_rm[3*i+1];
		llr_p2[i] = llr_in_rm[3*i+2];
	}

	//process tail bit
	for(i=0; i<4; i++)
	{
		llr_s_tail[i] = llr_s1[Len+i];
		llr_p1_tail[i] = llr_p1[Len+i];
		llr_p2_tail[i] = llr_p2[Len+i];
	}
	
	llr_s1[Len]		= llr_s_tail[0];
	llr_s1[Len+1]	= llr_p2_tail[0];
	llr_s1[Len+2]	= llr_p1_tail[1];

	llr_p1[Len]		= llr_p1_tail[0];
	llr_p1[Len+1]	= llr_s_tail[1];
	llr_p1[Len+2]	= llr_p2_tail[1];

	llr_s2[Len]		= llr_s_tail[2];
	llr_s2[Len+1]	= llr_p2_tail[2];
	llr_s2[Len+2]	= llr_p1_tail[3];

	llr_p2[Len]		= llr_p1_tail[2];
	llr_p2[Len+1]	= llr_s_tail[3];
	llr_p2[Len+2]	= llr_p2_tail[3];
	/*
	//接收的尾比特与输入的尾比特要重排，有一定的映射关系
	double LLR_copys1 = llr_s1[Len + 1];
	//LLR_s1[Len] = LLR_s1[Len];//llr_s_tail[0];
	llr_s1[Len + 1] = llr_p2[Len];//llr_p2_tail[0];   无影响
	llr_s1[Len + 2] = llr_p1[Len + 1];//llr_p1_tail[1]; 无影响

	double LLR_copyp1 = llr_p1[Len + 2];
	//LLR_p1[Len] = LLR_p1[Len];//llr_p1_tail[0];              
	llr_p1[Len + 1] = LLR_copys1;//LLR_s1[Len + 1];//llr_s_tail[1]    被覆盖 ;
	llr_p1[Len + 2] = llr_p2[Len + 1];//llr_p2_tail[1];       

	llr_s2[Len] = llr_s1[Len];//llr_s_tail[2];
	llr_s2[Len + 1] = llr_p2[Len + 2];//llr_p2_tail[2];        无影响
	llr_s2[Len + 2] = llr_p1[Len + 3];//llr_p1_tail[3];         无影响

	llr_p2[Len] = LLR_copyp1;//LLR_p1[Len + 2];//llr_p1_tail[2];    被覆盖         
	llr_p2[Len + 1] = llr_s1[Len + 3]; //llr_s_tail[3];          无影响
	llr_p2[Len + 2] = llr_p2[Len + 3]; //llr_p2_tail[3];         无影响
	*/
	interleave_double(llr_s2, Len, interleave_table);

	for(i=0; i<Len; i++)
	{
		llr_es2[i] = store_llr[i];
	}

	for(iter_idx=0; iter_idx<iterCNum; iter_idx++)
	{
		// decoder 1
		// initialize for alfa and beta 
		alfa1[0] = 0;
		alfa1[1] = -inf;
		alfa1[2] = -inf;
		alfa1[3] = -inf;
		alfa1[4] = -inf;
		alfa1[5] = -inf;
		alfa1[6] = -inf;
		alfa1[7] = -inf;

		beta1[(Len+Mem-1)*StateNum+0] = 0;
		beta1[(Len+Mem-1)*StateNum+1] = -inf;
		beta1[(Len+Mem-1)*StateNum+2] = -inf;
		beta1[(Len+Mem-1)*StateNum+3] = -inf;
		beta1[(Len+Mem-1)*StateNum+4] = -inf;
		beta1[(Len+Mem-1)*StateNum+5] = -inf;
		beta1[(Len+Mem-1)*StateNum+6] = -inf;
		beta1[(Len+Mem-1)*StateNum+7] = -inf;

		/* initialize for gs and gp */
		for(i=0; i<Len; i++)
		{
			gs1[i] = 0.5*llr_s1[i] + 0.5*llr_es2[i];
			gp1[i] = 0.5*llr_p1[i];
		}
		for(i=Len; i<Len+Mem; i++)
		{
			gs1[i] = 0.5*llr_s1[i];
			gp1[i] = 0.5*llr_p1[i];
		}

		/* compute the alfa */
		for(i=1; i<Len+Mem; i++)
		{
			alfa1[i*StateNum+0] = ACSO(alfa1[(i-1)*StateNum+0], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+1], +gs1[i-1]+gp1[i-1]);
			alfa1[i*StateNum+1] = ACSO(alfa1[(i-1)*StateNum+3], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+2], +gs1[i-1]-gp1[i-1]);
			alfa1[i*StateNum+2] = ACSO(alfa1[(i-1)*StateNum+4], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+5], +gs1[i-1]-gp1[i-1]);
			alfa1[i*StateNum+3] = ACSO(alfa1[(i-1)*StateNum+7], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+6], +gs1[i-1]+gp1[i-1]);
			alfa1[i*StateNum+4] = ACSO(alfa1[(i-1)*StateNum+1], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+0], +gs1[i-1]+gp1[i-1]);
			alfa1[i*StateNum+5] = ACSO(alfa1[(i-1)*StateNum+2], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+3], +gs1[i-1]-gp1[i-1]);
			alfa1[i*StateNum+6] = ACSO(alfa1[(i-1)*StateNum+5], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+4], +gs1[i-1]-gp1[i-1]);
			alfa1[i*StateNum+7] = ACSO(alfa1[(i-1)*StateNum+6], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+7], +gs1[i-1]+gp1[i-1]);
		}

		/* compute the beta */
		for(i=Len+Mem-2; i>=0; i--)
		{
			beta1[i*StateNum+0] = ACSO(beta1[(i+1)*StateNum+0], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+4], +gs1[i+1]+gp1[i+1]);
			beta1[i*StateNum+1] = ACSO(beta1[(i+1)*StateNum+4], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+0], +gs1[i+1]+gp1[i+1]);
			beta1[i*StateNum+2] = ACSO(beta1[(i+1)*StateNum+5], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+1], +gs1[i+1]-gp1[i+1]);
			beta1[i*StateNum+3] = ACSO(beta1[(i+1)*StateNum+1], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+5], +gs1[i+1]-gp1[i+1]);
			beta1[i*StateNum+4] = ACSO(beta1[(i+1)*StateNum+2], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+6], +gs1[i+1]-gp1[i+1]);
			beta1[i*StateNum+5] = ACSO(beta1[(i+1)*StateNum+6], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+2], +gs1[i+1]-gp1[i+1]);
			beta1[i*StateNum+6] = ACSO(beta1[(i+1)*StateNum+7], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+3], +gs1[i+1]+gp1[i+1]);
			beta1[i*StateNum+7] = ACSO(beta1[(i+1)*StateNum+3], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+7], +gs1[i+1]+gp1[i+1]);
		}

		/* compute the llres1 */
		for(i=0; i<Len+Mem; i++)
		{
			double nom,denom;
			double nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;		

			nom1	=	maxo(alfa1[i*StateNum+0]+gp1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gp1[i]+beta1[i*StateNum+0]);
			nom2	=	maxo(alfa1[i*StateNum+2]-gp1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]-gp1[i]+beta1[i*StateNum+5]);
			nom3	=	maxo(alfa1[i*StateNum+4]-gp1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]-gp1[i]+beta1[i*StateNum+2]);
			nom4	=	maxo(alfa1[i*StateNum+6]+gp1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gp1[i]+beta1[i*StateNum+7]);

			denom1	=	maxo(alfa1[i*StateNum+0]-gp1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gp1[i]+beta1[i*StateNum+4]);
			denom2	=	maxo(alfa1[i*StateNum+2]+gp1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]+gp1[i]+beta1[i*StateNum+1]);
			denom3	=	maxo(alfa1[i*StateNum+4]+gp1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]+gp1[i]+beta1[i*StateNum+6]);
			denom4	=	maxo(alfa1[i*StateNum+6]-gp1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gp1[i]+beta1[i*StateNum+3]);

			nom = maxo(maxo(nom1,nom2),maxo(nom3,nom4));
			denom = maxo(maxo(denom1,denom2),maxo(denom3,denom4));

			llr_es1[i] = nom - denom;
			
			if (llr_es1[i] > limit)
			{
				llr_es1[i] = limit;
			}

			if (llr_es1[i] < -limit)
			{
				llr_es1[i] = -limit;
			}
		}
		/* end of decoder1 */

		interleave_double(llr_es1, Len, interleave_table);

		/* decoder2 */
		/* initialize for alfa and beta */
		alfa2[0] = 0;
		alfa2[1] = -inf;
		alfa2[2] = -inf;
		alfa2[3] = -inf;
		alfa2[4] = -inf;
		alfa2[5] = -inf;
		alfa2[6] = -inf;
		alfa2[7] = -inf;

		beta2[(Len+Mem-1)*StateNum+0] = 0;
		beta2[(Len+Mem-1)*StateNum+1] = -inf;
		beta2[(Len+Mem-1)*StateNum+2] = -inf;
		beta2[(Len+Mem-1)*StateNum+3] = -inf;
		beta2[(Len+Mem-1)*StateNum+4] = -inf;
		beta2[(Len+Mem-1)*StateNum+5] = -inf;
		beta2[(Len+Mem-1)*StateNum+6] = -inf;
		beta2[(Len+Mem-1)*StateNum+7] = -inf;

		/* initialize for gp and gs */
		for(i=0; i<Len; i++)
		{
			gs2[i] = 0.5*llr_s2[i] + 0.5*llr_es1[i];
			gp2[i] = 0.5*llr_p2[i];
		}
		for(i=Len; i<Len+Mem; i++)
		{
			gs2[i] = 0.5*llr_s2[i];
			gp2[i] = 0.5*llr_p2[i];
		}

		/* compute the alfa2 */
		for(i=1; i<Len+Mem; i++)
		{
			alfa2[i*StateNum+0] = ACSO(alfa2[(i-1)*StateNum+0], -gs2[i-1]-gp2[i-1], alfa2[(i-1)*StateNum+1], +gs2[i-1]+gp2[i-1]);
			alfa2[i*StateNum+1] = ACSO(alfa2[(i-1)*StateNum+3], -gs2[i-1]+gp2[i-1], alfa2[(i-1)*StateNum+2], +gs2[i-1]-gp2[i-1]);
			alfa2[i*StateNum+2] = ACSO(alfa2[(i-1)*StateNum+4], -gs2[i-1]+gp2[i-1], alfa2[(i-1)*StateNum+5], +gs2[i-1]-gp2[i-1]);
			alfa2[i*StateNum+3] = ACSO(alfa2[(i-1)*StateNum+7], -gs2[i-1]-gp2[i-1], alfa2[(i-1)*StateNum+6], +gs2[i-1]+gp2[i-1]);
			alfa2[i*StateNum+4] = ACSO(alfa2[(i-1)*StateNum+1], -gs2[i-1]-gp2[i-1], alfa2[(i-1)*StateNum+0], +gs2[i-1]+gp2[i-1]);
			alfa2[i*StateNum+5] = ACSO(alfa2[(i-1)*StateNum+2], -gs2[i-1]+gp2[i-1], alfa2[(i-1)*StateNum+3], +gs2[i-1]-gp2[i-1]);
			alfa2[i*StateNum+6] = ACSO(alfa2[(i-1)*StateNum+5], -gs2[i-1]+gp2[i-1], alfa2[(i-1)*StateNum+4], +gs2[i-1]-gp2[i-1]);
			alfa2[i*StateNum+7] = ACSO(alfa2[(i-1)*StateNum+6], -gs2[i-1]-gp2[i-1], alfa2[(i-1)*StateNum+7], +gs2[i-1]+gp2[i-1]);
		}

		/* compute the beta */
		for(i=Len+Mem-2; i>=0; i--)
		{
			beta2[i*StateNum+0] = ACSO(beta2[(i+1)*StateNum+0], -gs2[i+1]-gp2[i+1], beta2[(i+1)*StateNum+4], +gs2[i+1]+gp2[i+1]);
			beta2[i*StateNum+1] = ACSO(beta2[(i+1)*StateNum+4], -gs2[i+1]-gp2[i+1], beta2[(i+1)*StateNum+0], +gs2[i+1]+gp2[i+1]);
			beta2[i*StateNum+2] = ACSO(beta2[(i+1)*StateNum+5], -gs2[i+1]+gp2[i+1], beta2[(i+1)*StateNum+1], +gs2[i+1]-gp2[i+1]);
			beta2[i*StateNum+3] = ACSO(beta2[(i+1)*StateNum+1], -gs2[i+1]+gp2[i+1], beta2[(i+1)*StateNum+5], +gs2[i+1]-gp2[i+1]);
			beta2[i*StateNum+4] = ACSO(beta2[(i+1)*StateNum+2], -gs2[i+1]+gp2[i+1], beta2[(i+1)*StateNum+6], +gs2[i+1]-gp2[i+1]);
			beta2[i*StateNum+5] = ACSO(beta2[(i+1)*StateNum+6], -gs2[i+1]+gp2[i+1], beta2[(i+1)*StateNum+2], +gs2[i+1]-gp2[i+1]);
			beta2[i*StateNum+6] = ACSO(beta2[(i+1)*StateNum+7], -gs2[i+1]-gp2[i+1], beta2[(i+1)*StateNum+3], +gs2[i+1]+gp2[i+1]);
			beta2[i*StateNum+7] = ACSO(beta2[(i+1)*StateNum+3], -gs2[i+1]-gp2[i+1], beta2[(i+1)*StateNum+7], +gs2[i+1]+gp2[i+1]);
		}

		/* compute the llres2 */
		for(i=0; i<Len+Mem; i++)
		{
			double nom,denom;
			double nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;		

			nom1	=	maxo(alfa2[i*StateNum+0]+gp2[i]+beta2[i*StateNum+4], alfa2[i*StateNum+1]+gp2[i]+beta2[i*StateNum+0]);
			nom2	=	maxo(alfa2[i*StateNum+2]-gp2[i]+beta2[i*StateNum+1], alfa2[i*StateNum+3]-gp2[i]+beta2[i*StateNum+5]);
			nom3	=	maxo(alfa2[i*StateNum+4]-gp2[i]+beta2[i*StateNum+6], alfa2[i*StateNum+5]-gp2[i]+beta2[i*StateNum+2]);
			nom4	=	maxo(alfa2[i*StateNum+6]+gp2[i]+beta2[i*StateNum+3], alfa2[i*StateNum+7]+gp2[i]+beta2[i*StateNum+7]);

			denom1	=	maxo(alfa2[i*StateNum+0]-gp2[i]+beta2[i*StateNum+0], alfa2[i*StateNum+1]-gp2[i]+beta2[i*StateNum+4]);
			denom2	=	maxo(alfa2[i*StateNum+2]+gp2[i]+beta2[i*StateNum+5], alfa2[i*StateNum+3]+gp2[i]+beta2[i*StateNum+1]);
			denom3	=	maxo(alfa2[i*StateNum+4]+gp2[i]+beta2[i*StateNum+2], alfa2[i*StateNum+5]+gp2[i]+beta2[i*StateNum+6]);
			denom4	=	maxo(alfa2[i*StateNum+6]-gp2[i]+beta2[i*StateNum+7], alfa2[i*StateNum+7]-gp2[i]+beta2[i*StateNum+3]);

			nom = maxo(maxo(nom1,nom2),maxo(nom3,nom4));
			denom = maxo(maxo(denom1,denom2),maxo(denom3,denom4));

			llr_es2[i] = nom - denom;
			
			if (llr_es2[i] > limit)
			{
				llr_es2[i] = limit;
			}

			if (llr_es2[i] < -limit)
			{
				llr_es2[i] = -limit;
			}
		}
		/* end of decoder2 */

		deinterleave_double(llr_es2, Len, interleave_table);
	}

	for(i=0; i<Len; i++)
	{
		store_llr[i] = llr_es2[i];	
	}

	deinterleave_double(llr_es1, Len, interleave_table);

	/* compute the llr_ep1 */
	for(i=0; i<Len+Mem; i++)
	{
		double nom,denom;
		double nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;	

		nom1	=	maxo(alfa1[i*StateNum+0]+gs1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gs1[i]+beta1[i*StateNum+0]);
		nom2	=	maxo(alfa1[i*StateNum+2]-gs1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]-gs1[i]+beta1[i*StateNum+1]);
		nom3	=	maxo(alfa1[i*StateNum+4]-gs1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]-gs1[i]+beta1[i*StateNum+6]);
		nom4	=	maxo(alfa1[i*StateNum+6]+gs1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gs1[i]+beta1[i*StateNum+7]);

		denom1	=	maxo(alfa1[i*StateNum+0]-gs1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gs1[i]+beta1[i*StateNum+4]);
		denom2	=	maxo(alfa1[i*StateNum+2]+gs1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]+gs1[i]+beta1[i*StateNum+5]);
		denom3	=	maxo(alfa1[i*StateNum+4]+gs1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]+gs1[i]+beta1[i*StateNum+2]);
		denom4	=	maxo(alfa1[i*StateNum+6]-gs1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gs1[i]+beta1[i*StateNum+3]);

		nom = maxo(maxo(nom1,nom2),maxo(nom3,nom4));
		denom = maxo(maxo(denom1,denom2),maxo(denom3,denom4));

		llr_ep1[i] = nom - denom;

		if (llr_ep1[i] > limit)
		{
			llr_ep1[i] = limit;
		}
		if (llr_ep1[i] < -limit)
		{
			llr_ep1[i] = -limit;
		}
	}

	/* compute the llr_ep2 */
	for(i=0; i<Len+Mem; i++)
	{
		double nom,denom;
		double nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;	

		nom1	=	maxo(alfa2[i*StateNum+0]+gs2[i]+beta2[i*StateNum+4], alfa2[i*StateNum+1]+gs2[i]+beta2[i*StateNum+0]);
		nom2	=	maxo(alfa2[i*StateNum+2]-gs2[i]+beta2[i*StateNum+5], alfa2[i*StateNum+3]-gs2[i]+beta2[i*StateNum+1]);
		nom3	=	maxo(alfa2[i*StateNum+4]-gs2[i]+beta2[i*StateNum+2], alfa2[i*StateNum+5]-gs2[i]+beta2[i*StateNum+6]);
		nom4	=	maxo(alfa2[i*StateNum+6]+gs2[i]+beta2[i*StateNum+3], alfa2[i*StateNum+7]+gs2[i]+beta2[i*StateNum+7]);

		denom1	=	maxo(alfa2[i*StateNum+0]-gs2[i]+beta2[i*StateNum+0], alfa2[i*StateNum+1]-gs2[i]+beta2[i*StateNum+4]);
		denom2	=	maxo(alfa2[i*StateNum+2]+gs2[i]+beta2[i*StateNum+1], alfa2[i*StateNum+3]+gs2[i]+beta2[i*StateNum+5]);
		denom3	=	maxo(alfa2[i*StateNum+4]+gs2[i]+beta2[i*StateNum+6], alfa2[i*StateNum+5]+gs2[i]+beta2[i*StateNum+2]);
		denom4	=	maxo(alfa2[i*StateNum+6]-gs2[i]+beta2[i*StateNum+7], alfa2[i*StateNum+7]-gs2[i]+beta2[i*StateNum+3]);

		nom = maxo(maxo(nom1,nom2),maxo(nom3,nom4));
		denom = maxo(maxo(denom1,denom2),maxo(denom3,denom4));

		llr_ep2[i] = nom - denom;

		if (llr_ep2[i] > limit)
		{
			llr_ep2[i] = limit;
		}

		if (llr_ep2[i] < -limit)
		{
			llr_ep2[i] = -limit;
		}
	}

	for(i=0; i<Len+Mem; i++)
	{
		llr_s1[i] = llr_s1[i] + llr_es1[i] + llr_es2[i];
		llr_p1[i] = llr_p1[i] + llr_ep1[i];
		llr_p2[i] = llr_p2[i] + llr_ep2[i];
	}	
	
	// Tailbiting
	llr_s_tail[0]	=	llr_s1[Len];
	llr_p1_tail[0]	=	llr_p1[Len];
	llr_p2_tail[0]	=	llr_s1[Len+1];

	llr_s_tail[1]	=	llr_p1[Len+1];
	llr_p1_tail[1]	=	llr_s1[Len+2];
	llr_p2_tail[1]	=	llr_p1[Len+2];

	llr_s_tail[2]	=	llr_s1[Len];
	llr_p1_tail[2]	=	llr_p2[Len];
	llr_p2_tail[2]	=	llr_s1[Len+1];

	llr_s_tail[3]	=	llr_p2[Len+1];
	llr_p1_tail[3]	=	llr_s1[Len+2];
	llr_p2_tail[3]	=	llr_p2[Len+2];

	for(i=0; i<4; i++)
	{
		llr_s1[Len+i] = llr_s_tail[i];
		llr_p1[Len+i] = llr_p1_tail[i];
		llr_p2[Len+i] = llr_p2_tail[i];
	}
	/* decision */
	for(i=0; i<Len; i++)
	{
		dec[i] = (llr_s1[i] > 0)? 1:0;
	}

	for(i=0; i<Len+4; i++)
	{
		llr_out_rm[3*i] = llr_s1[i];
		llr_out_rm[3*i+1] = llr_p1[i];
		llr_out_rm[3*i+2] = llr_p2[i];
	}

	for(i=0; i<3*(Len+4); i++)
	{
		llr_out_rm[i] = (llr_out_rm[i]>64)? 64:llr_out_rm[i];
		llr_out_rm[i] = (llr_out_rm[i]<-64)? -64:llr_out_rm[i];			
	}

	rm_interleave_doub(Len+4, CodeLen, llr_out_rm, llr_out, ratematch_table);

	free(llr_in_rm);
	free(llr_out_rm);
	free(llr_s1);
	free(llr_s2);
	free(llr_p1);
	free(llr_p2);
	free(llr_es1);
	free(llr_es2);
	free(llr_ep1);
	free(llr_ep2);
	free(alfa1);
	free(beta1);
	free(gs1);
	free(gp1);
	free(alfa2);
	free(beta2);
	free(gs2);
	free(gp2);
}



//rate dematching; rate dematch with table inter_index
void	rate_dematching_fix(int Len, int CodeLen, int *input, int *output, int *inter_index)
{
	int i;
	int j;

	for(i=0; i<3*Len; i++)
	{
		output[i] = 0;
	}
	for(i=0; i<CodeLen; i++)
	{
		
		output[inter_index[i]] += input[i];
	}
}

void	rm_interleave_fix(int Len, int CodeLen, int *input, int *output, int *inter_index)
{
	int i;
	
	for(i=0; i<CodeLen; i++)
	{
		output[i] = input[inter_index[i]];
	}
}

void interleave_fix(int* msg, int Len, int *interleave_table)
{
	int	*tmp;
	int		i;

	tmp = (int*) malloc(Len*sizeof(int));

	for(i=0; i<Len; i++)
	{
		tmp[i] = msg[interleave_table[i]];
	}

	for(i=0; i<Len; i++)
	{
		msg[i] = tmp[i];
	}

	free(tmp);
}

void deinterleave_fix(int* msg, int Len, int *interleave_table)
{
	int	*tmp;
	int		i;

	tmp = (int*) malloc(Len*sizeof(int));

	for(i=0; i<Len; i++)
	{
		tmp[interleave_table[i]] = msg[i];
	}

	for(i=0; i<Len; i++)
	{
		msg[i] = tmp[i];
	}

	free(tmp);
}

int maxo_fix(int t0, int t1)
{
	//if (t0 >= t1)
	//{
	//	return t0 + log(1+exp(-fabs(t0-t1)));
	//}
	//else
	//{
	//	return t1 + log(1+exp(-fabs(t0-t1)));
	//}
	if (t0 >= t1)
	{
		return t0;
	}
	else
	{
		return t1;
	}
}

int maxo_avx(short t0, short t1)
{
	//if (t0 >= t1)
	//{
	//	return t0 + log(1+exp(-fabs(t0-t1)));
	//}
	//else
	//{
	//	return t1 + log(1+exp(-fabs(t0-t1)));
	//}
	if (t0 >= t1)
	{
		return t0;
	}
	else
	{
		return t1;
	}
}

int ACSO_fix(int a0, int g0, int a1, int g1)
{
	int	t0, t1;

	t0 = a0 + g0;
	t1 = a1 + g1;
	
	return maxo_fix(t0, t1);
}

int	fix_GetInt(int a, int width)
{
	if ((width<0)||(width>32))
	{
		printf("The width is out of range\n");
		return 0;
	}

	
	return a&((1<<width)-1);

}

void turbo_decode_rm_fix(short* dec, 
						double* llr_out, 
						double* llr_in, 
						double* store_llr, 
						int Len, 
						int CodeLen,
						int iterCNum, 
						int* interleave_table,
						int* ratematch_table)
{
	int			*llr_in_fix, *llr_out_fix;
	int			*llr_in_rm, *llr_out_rm, *llr_s1, *llr_s2, *llr_p1, *llr_p2, *llr_es1, *llr_es2, *llr_ep1, *llr_ep2;
	int			*alfa1, *beta1, *gs1, *gp1, *alfa2, *beta2, *gs2, *gp2, *g0, *g1;

	int			i, s, iter_idx;
	int			llr_s_tail[4], llr_p1_tail[4], llr_p2_tail[4];
	int			inf = 4095;
	int			limit = 255;
	int			qbits = 4;
	int			bitwidth = 20;
	int			StateNum = 8;

	int			*tmp_llr_es1 = (int*) malloc((Len+Mem)*sizeof(int));

	llr_in_fix	=	(int*) malloc(CodeLen*sizeof(int));
	llr_out_fix	=	(int*) malloc(CodeLen*sizeof(int));

	llr_in_rm	=	(int*) malloc(3*(Len+4)*sizeof(int));
	llr_out_rm	=	(int*) malloc(3*(Len+4)*sizeof(int));
	llr_s1		=	(int*) malloc((Len+4)*sizeof(int));
	llr_s2		=	(int*) malloc((Len+4)*sizeof(int));
	llr_p1		=	(int*) malloc((Len+4)*sizeof(int));
	llr_p2		=	(int*) malloc((Len+4)*sizeof(int));
	llr_es1		=	(int*) malloc((Len+Mem)*sizeof(int));
	llr_es2		=	(int*) malloc((Len+Mem)*sizeof(int));
	llr_ep1		=	(int*) malloc((Len+Mem)*sizeof(int));
	llr_ep2		=	(int*) malloc((Len+Mem)*sizeof(int));

	alfa1 = (int*) malloc((Len+Mem)*StateNum*sizeof(int));
	beta1 = (int*) malloc((Len+Mem)*StateNum*sizeof(int));
	alfa2 = (int*) malloc((Len+Mem)*StateNum*sizeof(int));
	beta2 = (int*) malloc((Len+Mem)*StateNum*sizeof(int));
	gs1 = (int*) malloc((Len+Mem)*sizeof(int));
	gp1 = (int*) malloc((Len+Mem)*sizeof(int));
	gs2 = (int*) malloc((Len+Mem)*sizeof(int));
	gp2 = (int*) malloc((Len+Mem)*sizeof(int));
	g0 = (int*) malloc((Len+Mem)*sizeof(int));
	g1 = (int*) malloc((Len+Mem)*sizeof(int));

	for(i=0; i<CodeLen; i++)
	{
		llr_in_fix[i] = (int) (llr_in[i]*(1<<qbits));

		if (llr_in_fix[i] > limit)
		{
			llr_in_fix[i] = limit;
		}

		if (llr_in_fix[i] < -limit)
		{
			llr_in_fix[i] = -limit;
		}

	}

	//rate match....
	rate_dematching_fix(Len+4, CodeLen, llr_in_fix, llr_in_rm, ratematch_table);

	for(i=0; i<Len+4; i++)
	{
		llr_s1[i] = llr_in_rm[3*i];
		llr_s2[i] = llr_in_rm[3*i];
		llr_p1[i] = llr_in_rm[3*i+1];
		llr_p2[i] = llr_in_rm[3*i+2];
	}

	//process tail bit
	for(i=0; i<4; i++)
	{
		llr_s_tail[i] = llr_s1[Len+i];
		llr_p1_tail[i] = llr_p1[Len+i];
		llr_p2_tail[i] = llr_p2[Len+i];
	}
	llr_s1[Len]		= llr_s_tail[0];
	llr_s1[Len+1]	= llr_p2_tail[0];
	llr_s1[Len+2]	= llr_p1_tail[1];

	llr_p1[Len]		= llr_p1_tail[0];
	llr_p1[Len+1]	= llr_s_tail[1];
	llr_p1[Len+2]	= llr_p2_tail[1];

	llr_s2[Len]		= llr_s_tail[2];
	llr_s2[Len+1]	= llr_p2_tail[2];
	llr_s2[Len+2]	= llr_p1_tail[3];

	llr_p2[Len]		= llr_p1_tail[2];
	llr_p2[Len+1]	= llr_s_tail[3];
	llr_p2[Len+2]	= llr_p2_tail[3];

	interleave_fix(llr_s2, Len, interleave_table);

	for(i=0; i<Len; i++)
	{
		llr_es2[i] = (int)	(store_llr[i]*(1<<qbits));
	}

	for(iter_idx=0; iter_idx<iterCNum; iter_idx++)
	{
		// decoder 1
		// initialize for alfa and beta 
		alfa1[0] = 0;
		alfa1[1] = -inf;
		alfa1[2] = -inf;
		alfa1[3] = -inf;
		alfa1[4] = -inf;
		alfa1[5] = -inf;
		alfa1[6] = -inf;
		alfa1[7] = -inf;

		beta1[(Len+Mem-1)*StateNum+0] = 0;
		beta1[(Len+Mem-1)*StateNum+1] = -inf;
		beta1[(Len+Mem-1)*StateNum+2] = -inf;
		beta1[(Len+Mem-1)*StateNum+3] = -inf;
		beta1[(Len+Mem-1)*StateNum+4] = -inf;
		beta1[(Len+Mem-1)*StateNum+5] = -inf;
		beta1[(Len+Mem-1)*StateNum+6] = -inf;
		beta1[(Len+Mem-1)*StateNum+7] = -inf;

		/* initialize for gs and gp */
		for(i=0; i<Len; i++)
		{
			gs1[i] = (llr_s1[i]>>1) + (llr_es2[i]>>1);
			gp1[i] = (llr_p1[i]>>1);
			
			g0[i] = - gs1[i] -  gp1[i]; 
			g1[i] = - gs1[i] +  gp1[i]; 
		}
		for(i=Len; i<Len+Mem; i++)
		{
			gs1[i] = (llr_s1[i]>>1);
			gp1[i] = (llr_p1[i]>>1);

			g0[i] = - gs1[i] -  gp1[i]; 
			g1[i] = - gs1[i] +  gp1[i]; 
		}

		/* compute the alfa and beta */
		for(i=1; i<Len+Mem; i++)
		{
			int		j = Len+Mem-1-i;
			int		at0[8], at1[8], as0[8], as1[8];
			int		bt0[8], bt1[8], bs0[8], bs1[8];
			int		ct0[8], ct1[8], dt0[8], dt1[8];

			int		*pa = alfa1 + (i-1)*StateNum;
			int		*pb = beta1 + (j+1)*StateNum;

			int		*pta, *ptb;

			int		nom,denom;

			//if (i>=j-1)
			//{
			//	printf("\n");
			//}

			// add
			at0[0] = pa[0] + g0[i-1];
			at0[1] = pa[1] + g0[i-1];
			at0[2] = pa[2] + g1[i-1];
			at0[3] = pa[3] + g1[i-1];
			at0[4] = pa[4] + g1[i-1];
			at0[5] = pa[5] + g1[i-1];
			at0[6] = pa[6] + g0[i-1];
			at0[7] = pa[7] + g0[i-1];

			at1[0] = pa[0] - g0[i-1];
			at1[1] = pa[1] - g0[i-1];
			at1[2] = pa[2] - g1[i-1];
			at1[3] = pa[3] - g1[i-1];
			at1[4] = pa[4] - g1[i-1];
			at1[5] = pa[5] - g1[i-1];
			at1[6] = pa[6] - g0[i-1];
			at1[7] = pa[7] - g0[i-1];

			// shuffle
			as0[0] = at0[0];
			as0[1] = at0[3];
			as0[2] = at0[4];
			as0[3] = at0[7];
			as0[4] = at0[1];
			as0[5] = at0[2];
			as0[6] = at0[5];
			as0[7] = at0[6];

			as1[0] = at1[1];
			as1[1] = at1[2];
			as1[2] = at1[5];
			as1[3] = at1[6];
			as1[4] = at1[0];
			as1[5] = at1[3];
			as1[6] = at1[4];
			as1[7] = at1[7];

			// max
			pa = pa + StateNum;

			pa[0] = (as0[0]>as1[0])?as0[0]:as1[0];
			pa[1] = (as0[1]>as1[1])?as0[1]:as1[1];
			pa[2] = (as0[2]>as1[2])?as0[2]:as1[2];
			pa[3] = (as0[3]>as1[3])?as0[3]:as1[3];
			pa[4] = (as0[4]>as1[4])?as0[4]:as1[4];
			pa[5] = (as0[5]>as1[5])?as0[5]:as1[5];
			pa[6] = (as0[6]>as1[6])?as0[6]:as1[6];
			pa[7] = (as0[7]>as1[7])?as0[7]:as1[7];

			//alfa1[i*StateNum+0] = ACSO_fix(alfa1[(i-1)*StateNum+0], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+1], +gs1[i-1]+gp1[i-1]);
			//alfa1[i*StateNum+1] = ACSO_fix(alfa1[(i-1)*StateNum+3], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+2], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+2] = ACSO_fix(alfa1[(i-1)*StateNum+4], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+5], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+3] = ACSO_fix(alfa1[(i-1)*StateNum+7], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+6], +gs1[i-1]+gp1[i-1]);
			//alfa1[i*StateNum+4] = ACSO_fix(alfa1[(i-1)*StateNum+1], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+0], +gs1[i-1]+gp1[i-1]);
			//alfa1[i*StateNum+5] = ACSO_fix(alfa1[(i-1)*StateNum+2], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+3], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+6] = ACSO_fix(alfa1[(i-1)*StateNum+5], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+4], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+7] = ACSO_fix(alfa1[(i-1)*StateNum+6], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+7], +gs1[i-1]+gp1[i-1]);
		
	/*		fix_GetInt(alfa1[i*StateNum+0],bitwidth);
			fix_GetInt(alfa1[i*StateNum+1],bitwidth);
			fix_GetInt(alfa1[i*StateNum+2],bitwidth);
			fix_GetInt(alfa1[i*StateNum+3],bitwidth);
			fix_GetInt(alfa1[i*StateNum+4],bitwidth);
			fix_GetInt(alfa1[i*StateNum+5],bitwidth);
			fix_GetInt(alfa1[i*StateNum+6],bitwidth);
			fix_GetInt(alfa1[i*StateNum+7],bitwidth);*/

			// beta

			// add
			bt0[0] = pb[0] + g0[j+1];
			bt0[1] = pb[1] + g1[j+1];
			bt0[2] = pb[2] + g1[j+1];
			bt0[3] = pb[3] + g0[j+1];
			bt0[4] = pb[4] + g0[j+1];
			bt0[5] = pb[5] + g1[j+1];
			bt0[6] = pb[6] + g1[j+1];
			bt0[7] = pb[7] + g0[j+1];

			bt1[0] = pb[0] - g0[j+1];
			bt1[1] = pb[1] - g1[j+1];
			bt1[2] = pb[2] - g1[j+1];
			bt1[3] = pb[3] - g0[j+1];
			bt1[4] = pb[4] - g0[j+1];
			bt1[5] = pb[5] - g1[j+1];
			bt1[6] = pb[6] - g1[j+1];
			bt1[7] = pb[7] - g0[j+1];

			// shuffle
			bs0[0] = bt0[0];
			bs0[1] = bt0[4];
			bs0[2] = bt0[5];
			bs0[3] = bt0[1];
			bs0[4] = bt0[2];
			bs0[5] = bt0[6];
			bs0[6] = bt0[7];
			bs0[7] = bt0[3];

			bs1[0] = bt1[4];
			bs1[1] = bt1[0];
			bs1[2] = bt1[1];
			bs1[3] = bt1[5];
			bs1[4] = bt1[6];
			bs1[5] = bt1[2];
			bs1[6] = bt1[3];
			bs1[7] = bt1[7];

			// max
			pb = pb - StateNum;

			pb[0] = (bs0[0]>bs1[0])?bs0[0]:bs1[0];
			pb[1] = (bs0[1]>bs1[1])?bs0[1]:bs1[1];
			pb[2] = (bs0[2]>bs1[2])?bs0[2]:bs1[2];
			pb[3] = (bs0[3]>bs1[3])?bs0[3]:bs1[3];
			pb[4] = (bs0[4]>bs1[4])?bs0[4]:bs1[4];
			pb[5] = (bs0[5]>bs1[5])?bs0[5]:bs1[5];
			pb[6] = (bs0[6]>bs1[6])?bs0[6]:bs1[6];
			pb[7] = (bs0[7]>bs1[7])?bs0[7]:bs1[7];

	/*		fix_GetInt(beta1[j*StateNum+0],bitwidth);
			fix_GetInt(beta1[j*StateNum+1],bitwidth);
			fix_GetInt(beta1[j*StateNum+2],bitwidth);
			fix_GetInt(beta1[j*StateNum+3],bitwidth);
			fix_GetInt(beta1[j*StateNum+4],bitwidth);
			fix_GetInt(beta1[j*StateNum+5],bitwidth);
			fix_GetInt(beta1[j*StateNum+6],bitwidth);
			fix_GetInt(beta1[j*StateNum+7],bitwidth);*/

			// llres1
			ptb = beta1 + (i-1)*StateNum;

			ct0[0] = as0[0] + ptb[0];
			ct0[1] = as0[1] + ptb[1];
			ct0[2] = as0[2] + ptb[2];
			ct0[3] = as0[3] + ptb[3];
			ct0[4] = as0[4] + ptb[4];
			ct0[5] = as0[5] + ptb[5];
			ct0[6] = as0[6] + ptb[6];
			ct0[7] = as0[7] + ptb[7];

			ct1[0] = as1[0] + ptb[0];
			ct1[1] = as1[1] + ptb[1];
			ct1[2] = as1[2] + ptb[2];
			ct1[3] = as1[3] + ptb[3];
			ct1[4] = as1[4] + ptb[4];
			ct1[5] = as1[5] + ptb[5];
			ct1[6] = as1[6] + ptb[6];
			ct1[7] = as1[7] + ptb[7];

			nom = maxo_fix(maxo_fix(maxo_fix(ct1[0],ct1[1]),maxo_fix(ct1[2],ct1[3])),maxo_fix(maxo_fix(ct1[4],ct1[5]),maxo_fix(ct1[6],ct1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(ct0[0],ct0[1]),maxo_fix(ct0[2],ct0[3])),maxo_fix(maxo_fix(ct0[4],ct0[5]),maxo_fix(ct0[6],ct0[7])));

			llr_es1[i-1] = nom - denom - (gs1[i-1]<<1);
			//llr_es1[i-1] = nom - denom;

			llr_es1[i-1] = (llr_es1[i-1]>limit)? limit:llr_es1[i-1];
			llr_es1[i-1] = (llr_es1[i-1]<-limit)? -limit:llr_es1[i-1];

			pta = alfa1 + (j+1)*StateNum;

			dt0[0] = bs0[0] + pta[0];
			dt0[1] = bs0[1] + pta[1];
			dt0[2] = bs0[2] + pta[2];
			dt0[3] = bs0[3] + pta[3];
			dt0[4] = bs0[4] + pta[4];
			dt0[5] = bs0[5] + pta[5];
			dt0[6] = bs0[6] + pta[6];
			dt0[7] = bs0[7] + pta[7];

			dt1[0] = bs1[0] + pta[0];
			dt1[1] = bs1[1] + pta[1];
			dt1[2] = bs1[2] + pta[2];
			dt1[3] = bs1[3] + pta[3];
			dt1[4] = bs1[4] + pta[4];
			dt1[5] = bs1[5] + pta[5];
			dt1[6] = bs1[6] + pta[6];
			dt1[7] = bs1[7] + pta[7];

			nom = maxo_fix(maxo_fix(maxo_fix(dt1[0],dt1[1]),maxo_fix(dt1[2],dt1[3])),maxo_fix(maxo_fix(dt1[4],dt1[5]),maxo_fix(dt1[6],dt1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(dt0[0],dt0[1]),maxo_fix(dt0[2],dt0[3])),maxo_fix(maxo_fix(dt0[4],dt0[5]),maxo_fix(dt0[6],dt0[7])));
		
			llr_es1[j+1] = nom - denom - (gs1[j+1]<<1);
			//llr_es1[j+1] = nom - denom;
			
			llr_es1[j+1] = (llr_es1[j+1]>limit)? limit:llr_es1[j+1];
			llr_es1[j+1] = (llr_es1[j+1]<-limit)? -limit:llr_es1[j+1];
		}

		{
			int		*pta, *ptb;
			int		nom, denom;
			int		dt0[8], dt1[8];

			pta = alfa1;
			ptb = beta1;

			dt0[0] = pta[0]-gs1[0]-gp1[0]+ptb[0];
			dt0[1] = pta[1]-gs1[0]-gp1[0]+ptb[4];
			dt0[2] = pta[2]-gs1[0]+gp1[0]+ptb[5];
			dt0[3] = pta[3]-gs1[0]+gp1[0]+ptb[1];
			dt0[4] = pta[4]-gs1[0]+gp1[0]+ptb[2];
			dt0[5] = pta[5]-gs1[0]+gp1[0]+ptb[6];
			dt0[6] = pta[6]-gs1[0]-gp1[0]+ptb[7];
			dt0[7] = pta[7]-gs1[0]-gp1[0]+ptb[3];

			dt1[0] = pta[0]+gs1[0]+gp1[0]+ptb[4];
			dt1[1] = pta[1]+gs1[0]+gp1[0]+ptb[0];
			dt1[2] = pta[2]+gs1[0]-gp1[0]+ptb[1];
			dt1[3] = pta[3]+gs1[0]-gp1[0]+ptb[5];
			dt1[4] = pta[4]+gs1[0]-gp1[0]+ptb[6];
			dt1[5] = pta[5]+gs1[0]-gp1[0]+ptb[2];
			dt1[6] = pta[6]+gs1[0]+gp1[0]+ptb[3];
			dt1[7] = pta[7]+gs1[0]+gp1[0]+ptb[7];

			nom = maxo_fix(maxo_fix(maxo_fix(dt1[0],dt1[1]),maxo_fix(dt1[2],dt1[3])),maxo_fix(maxo_fix(dt1[4],dt1[5]),maxo_fix(dt1[6],dt1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(dt0[0],dt0[1]),maxo_fix(dt0[2],dt0[3])),maxo_fix(maxo_fix(dt0[4],dt0[5]),maxo_fix(dt0[6],dt0[7])));
		
			llr_es1[0] = nom - denom - (gs1[0]<<1);
			//llr_es1[0] = nom - denom;
			
			llr_es1[0] = (llr_es1[0]>limit)? limit:llr_es1[0];
			llr_es1[0] = (llr_es1[0]<-limit)? -limit:llr_es1[0];

			pta = alfa1 + (Len+Mem-1)*StateNum;
			ptb = beta1 + (Len+Mem-1)*StateNum;

			dt0[0] = pta[0]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[0];
			dt0[1] = pta[1]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[4];
			dt0[2] = pta[2]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[5];
			dt0[3] = pta[3]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[1];
			dt0[4] = pta[4]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[2];
			dt0[5] = pta[5]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[6];
			dt0[6] = pta[6]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[7];
			dt0[7] = pta[7]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[3];

			dt1[0] = pta[0]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[4];
			dt1[1] = pta[1]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[0];
			dt1[2] = pta[2]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[1];
			dt1[3] = pta[3]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[5];
			dt1[4] = pta[4]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[6];
			dt1[5] = pta[5]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[2];
			dt1[6] = pta[6]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[3];
			dt1[7] = pta[7]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[7];

			nom = maxo_fix(maxo_fix(maxo_fix(dt1[0],dt1[1]),maxo_fix(dt1[2],dt1[3])),maxo_fix(maxo_fix(dt1[4],dt1[5]),maxo_fix(dt1[6],dt1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(dt0[0],dt0[1]),maxo_fix(dt0[2],dt0[3])),maxo_fix(maxo_fix(dt0[4],dt0[5]),maxo_fix(dt0[6],dt0[7])));
		
			llr_es1[Len+Mem-1] = nom - denom - (gs1[Len+Mem-1]<<1);
			//llr_es1[Len+Mem-1] = nom - denom;
			
			llr_es1[Len+Mem-1] = (llr_es1[Len+Mem-1]>limit)? limit:llr_es1[Len+Mem-1];
			llr_es1[Len+Mem-1] = (llr_es1[Len+Mem-1]<-limit)? -limit:llr_es1[Len+Mem-1];
		}

		///* compute the beta */
		//for(i=Len+Mem-2; i>=0; i--)
		//{
		//	beta1[i*StateNum+0] = ACSO_fix(beta1[(i+1)*StateNum+0], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+4], +gs1[i+1]+gp1[i+1]);
		//	beta1[i*StateNum+1] = ACSO_fix(beta1[(i+1)*StateNum+4], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+0], +gs1[i+1]+gp1[i+1]);
		//	beta1[i*StateNum+2] = ACSO_fix(beta1[(i+1)*StateNum+5], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+1], +gs1[i+1]-gp1[i+1]);
		//	beta1[i*StateNum+3] = ACSO_fix(beta1[(i+1)*StateNum+1], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+5], +gs1[i+1]-gp1[i+1]);
		//	beta1[i*StateNum+4] = ACSO_fix(beta1[(i+1)*StateNum+2], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+6], +gs1[i+1]-gp1[i+1]);
		//	beta1[i*StateNum+5] = ACSO_fix(beta1[(i+1)*StateNum+6], -gs1[i+1]+gp1[i+1], beta1[(i+1)*StateNum+2], +gs1[i+1]-gp1[i+1]);
		//	beta1[i*StateNum+6] = ACSO_fix(beta1[(i+1)*StateNum+7], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+3], +gs1[i+1]+gp1[i+1]);
		//	beta1[i*StateNum+7] = ACSO_fix(beta1[(i+1)*StateNum+3], -gs1[i+1]-gp1[i+1], beta1[(i+1)*StateNum+7], +gs1[i+1]+gp1[i+1]);
		//
		//}

		/* compute the llres1 */
		for(i=0; i<Len+Mem; i++)
		{
			int nom,denom;
			int nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;		

			nom1	=	maxo_fix(alfa1[i*StateNum+0]+gp1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gp1[i]+beta1[i*StateNum+0]);
			nom2	=	maxo_fix(alfa1[i*StateNum+2]-gp1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]-gp1[i]+beta1[i*StateNum+5]);
			nom3	=	maxo_fix(alfa1[i*StateNum+4]-gp1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]-gp1[i]+beta1[i*StateNum+2]);
			nom4	=	maxo_fix(alfa1[i*StateNum+6]+gp1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gp1[i]+beta1[i*StateNum+7]);

			denom1	=	maxo_fix(alfa1[i*StateNum+0]-gp1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gp1[i]+beta1[i*StateNum+4]);
			denom2	=	maxo_fix(alfa1[i*StateNum+2]+gp1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]+gp1[i]+beta1[i*StateNum+1]);
			denom3	=	maxo_fix(alfa1[i*StateNum+4]+gp1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]+gp1[i]+beta1[i*StateNum+6]);
			denom4	=	maxo_fix(alfa1[i*StateNum+6]-gp1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gp1[i]+beta1[i*StateNum+3]);

			//nom1	=	maxo_fix(alfa1[i*StateNum+0]+gs1[i]+gp1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gs1[i]+gp1[i]+beta1[i*StateNum+0]);
			//nom2	=	maxo_fix(alfa1[i*StateNum+2]+gs1[i]-gp1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]+gs1[i]-gp1[i]+beta1[i*StateNum+5]);
			//nom3	=	maxo_fix(alfa1[i*StateNum+4]+gs1[i]-gp1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]+gs1[i]-gp1[i]+beta1[i*StateNum+2]);
			//nom4	=	maxo_fix(alfa1[i*StateNum+6]+gs1[i]+gp1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gs1[i]+gp1[i]+beta1[i*StateNum+7]);

			//denom1	=	maxo_fix(alfa1[i*StateNum+0]-gs1[i]-gp1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gs1[i]-gp1[i]+beta1[i*StateNum+4]);
			//denom2	=	maxo_fix(alfa1[i*StateNum+2]-gs1[i]+gp1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]-gs1[i]+gp1[i]+beta1[i*StateNum+1]);
			//denom3	=	maxo_fix(alfa1[i*StateNum+4]-gs1[i]+gp1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]-gs1[i]+gp1[i]+beta1[i*StateNum+6]);
			//denom4	=	maxo_fix(alfa1[i*StateNum+6]-gs1[i]-gp1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gs1[i]-gp1[i]+beta1[i*StateNum+3]);

			nom = maxo_fix(maxo_fix(nom1,nom2),maxo_fix(nom3,nom4));
			denom = maxo_fix(maxo_fix(denom1,denom2),maxo_fix(denom3,denom4));

			tmp_llr_es1[i] = nom - denom;
			
			if (tmp_llr_es1[i] > limit)
			{
				tmp_llr_es1[i] = limit;
			}

			if (tmp_llr_es1[i] < -limit)
			{
				tmp_llr_es1[i] = -limit;
			}
		}
		/* end of decoder1 */

		//{
		//	FILE *fp;

		//	fp = fopen("turbo_fix.txt","at");

		//	for(i=0; i<Len+Mem; i++)
		//	{
		//		fprintf(fp,"%8d, %8d, %8d, %8d, %8d, %8d\n", i, llr_es1[i], g0[i], g1[i], alfa1[i*StateNum+0], alfa1[i*StateNum+1]);
		//	}

		//	fclose(fp);
		//}

		interleave_fix(llr_es1, Len, interleave_table);

		/* decoder2 */
		/* initialize for alfa and beta */
		alfa2[0] = 0;
		alfa2[1] = -inf;
		alfa2[2] = -inf;
		alfa2[3] = -inf;
		alfa2[4] = -inf;
		alfa2[5] = -inf;
		alfa2[6] = -inf;
		alfa2[7] = -inf;

		beta2[(Len+Mem-1)*StateNum+0] = 0;
		beta2[(Len+Mem-1)*StateNum+1] = -inf;
		beta2[(Len+Mem-1)*StateNum+2] = -inf;
		beta2[(Len+Mem-1)*StateNum+3] = -inf;
		beta2[(Len+Mem-1)*StateNum+4] = -inf;
		beta2[(Len+Mem-1)*StateNum+5] = -inf;
		beta2[(Len+Mem-1)*StateNum+6] = -inf;
		beta2[(Len+Mem-1)*StateNum+7] = -inf;

		/* initialize for gp and gs */
		for(i=0; i<Len; i++)
		{
			gs2[i] = (llr_s2[i]>>1) + (llr_es1[i]>>1);
			gp2[i] = (llr_p2[i]>>1);

			g0[i] = - gs2[i] -  gp2[i]; 
			g1[i] = - gs2[i] +  gp2[i]; 
		}
		for(i=Len; i<Len+Mem; i++)
		{
			gs2[i] = (llr_s2[i]>>1);
			gp2[i] = (llr_p2[i]>>1);

			g0[i] = - gs2[i] -  gp2[i]; 
			g1[i] = - gs2[i] +  gp2[i]; 
		}

		/* compute the alfa2 */
		/* compute the alfa2, beta2 and LLR */
		for(i=1; i<Len+Mem; i++)
		{
			int		j = Len+Mem-1-i;
			int		at0[8], at1[8], as0[8], as1[8];
			int		bt0[8], bt1[8], bs0[8], bs1[8];
			int		ct0[8], ct1[8], dt0[8], dt1[8];

			int		*pa = alfa2 + (i-1)*StateNum;
			int		*pb = beta2 + (j+1)*StateNum;

			int		*pta, *ptb;

			int		nom,denom;

			//if (i>=j-1)
			//{
			//	printf("\n");
			//}

			// add
			at0[0] = pa[0] + g0[i-1];
			at0[1] = pa[1] + g0[i-1];
			at0[2] = pa[2] + g1[i-1];
			at0[3] = pa[3] + g1[i-1];
			at0[4] = pa[4] + g1[i-1];
			at0[5] = pa[5] + g1[i-1];
			at0[6] = pa[6] + g0[i-1];
			at0[7] = pa[7] + g0[i-1];

			at1[0] = pa[0] - g0[i-1];
			at1[1] = pa[1] - g0[i-1];
			at1[2] = pa[2] - g1[i-1];
			at1[3] = pa[3] - g1[i-1];
			at1[4] = pa[4] - g1[i-1];
			at1[5] = pa[5] - g1[i-1];
			at1[6] = pa[6] - g0[i-1];
			at1[7] = pa[7] - g0[i-1];

			// shuffle
			as0[0] = at0[0];
			as0[1] = at0[3];
			as0[2] = at0[4];
			as0[3] = at0[7];
			as0[4] = at0[1];
			as0[5] = at0[2];
			as0[6] = at0[5];
			as0[7] = at0[6];

			as1[0] = at1[1];
			as1[1] = at1[2];
			as1[2] = at1[5];
			as1[3] = at1[6];
			as1[4] = at1[0];
			as1[5] = at1[3];
			as1[6] = at1[4];
			as1[7] = at1[7];

			// max
			pa = pa + StateNum;

			pa[0] = (as0[0]>as1[0])?as0[0]:as1[0];
			pa[1] = (as0[1]>as1[1])?as0[1]:as1[1];
			pa[2] = (as0[2]>as1[2])?as0[2]:as1[2];
			pa[3] = (as0[3]>as1[3])?as0[3]:as1[3];
			pa[4] = (as0[4]>as1[4])?as0[4]:as1[4];
			pa[5] = (as0[5]>as1[5])?as0[5]:as1[5];
			pa[6] = (as0[6]>as1[6])?as0[6]:as1[6];
			pa[7] = (as0[7]>as1[7])?as0[7]:as1[7];
		
	/*		fix_GetInt(alfa1[i*StateNum+0],bitwidth);
			fix_GetInt(alfa1[i*StateNum+1],bitwidth);
			fix_GetInt(alfa1[i*StateNum+2],bitwidth);
			fix_GetInt(alfa1[i*StateNum+3],bitwidth);
			fix_GetInt(alfa1[i*StateNum+4],bitwidth);
			fix_GetInt(alfa1[i*StateNum+5],bitwidth);
			fix_GetInt(alfa1[i*StateNum+6],bitwidth);
			fix_GetInt(alfa1[i*StateNum+7],bitwidth);*/

			// beta

			// add
			bt0[0] = pb[0] + g0[j+1];
			bt0[1] = pb[1] + g1[j+1];
			bt0[2] = pb[2] + g1[j+1];
			bt0[3] = pb[3] + g0[j+1];
			bt0[4] = pb[4] + g0[j+1];
			bt0[5] = pb[5] + g1[j+1];
			bt0[6] = pb[6] + g1[j+1];
			bt0[7] = pb[7] + g0[j+1];

			bt1[0] = pb[0] - g0[j+1];
			bt1[1] = pb[1] - g1[j+1];
			bt1[2] = pb[2] - g1[j+1];
			bt1[3] = pb[3] - g0[j+1];
			bt1[4] = pb[4] - g0[j+1];
			bt1[5] = pb[5] - g1[j+1];
			bt1[6] = pb[6] - g1[j+1];
			bt1[7] = pb[7] - g0[j+1];

			// shuffle
			bs0[0] = bt0[0];
			bs0[1] = bt0[4];
			bs0[2] = bt0[5];
			bs0[3] = bt0[1];
			bs0[4] = bt0[2];
			bs0[5] = bt0[6];
			bs0[6] = bt0[7];
			bs0[7] = bt0[3];

			bs1[0] = bt1[4];
			bs1[1] = bt1[0];
			bs1[2] = bt1[1];
			bs1[3] = bt1[5];
			bs1[4] = bt1[6];
			bs1[5] = bt1[2];
			bs1[6] = bt1[3];
			bs1[7] = bt1[7];

			// max
			pb = pb - StateNum;

			pb[0] = (bs0[0]>bs1[0])?bs0[0]:bs1[0];
			pb[1] = (bs0[1]>bs1[1])?bs0[1]:bs1[1];
			pb[2] = (bs0[2]>bs1[2])?bs0[2]:bs1[2];
			pb[3] = (bs0[3]>bs1[3])?bs0[3]:bs1[3];
			pb[4] = (bs0[4]>bs1[4])?bs0[4]:bs1[4];
			pb[5] = (bs0[5]>bs1[5])?bs0[5]:bs1[5];
			pb[6] = (bs0[6]>bs1[6])?bs0[6]:bs1[6];
			pb[7] = (bs0[7]>bs1[7])?bs0[7]:bs1[7];

	/*		fix_GetInt(beta1[j*StateNum+0],bitwidth);
			fix_GetInt(beta1[j*StateNum+1],bitwidth);
			fix_GetInt(beta1[j*StateNum+2],bitwidth);
			fix_GetInt(beta1[j*StateNum+3],bitwidth);
			fix_GetInt(beta1[j*StateNum+4],bitwidth);
			fix_GetInt(beta1[j*StateNum+5],bitwidth);
			fix_GetInt(beta1[j*StateNum+6],bitwidth);
			fix_GetInt(beta1[j*StateNum+7],bitwidth);*/

			// llres1
			ptb = beta2 + (i-1)*StateNum;

			ct0[0] = as0[0] + ptb[0];
			ct0[1] = as0[1] + ptb[1];
			ct0[2] = as0[2] + ptb[2];
			ct0[3] = as0[3] + ptb[3];
			ct0[4] = as0[4] + ptb[4];
			ct0[5] = as0[5] + ptb[5];
			ct0[6] = as0[6] + ptb[6];
			ct0[7] = as0[7] + ptb[7];

			ct1[0] = as1[0] + ptb[0];
			ct1[1] = as1[1] + ptb[1];
			ct1[2] = as1[2] + ptb[2];
			ct1[3] = as1[3] + ptb[3];
			ct1[4] = as1[4] + ptb[4];
			ct1[5] = as1[5] + ptb[5];
			ct1[6] = as1[6] + ptb[6];
			ct1[7] = as1[7] + ptb[7];

			nom = maxo_fix(maxo_fix(maxo_fix(ct1[0],ct1[1]),maxo_fix(ct1[2],ct1[3])),maxo_fix(maxo_fix(ct1[4],ct1[5]),maxo_fix(ct1[6],ct1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(ct0[0],ct0[1]),maxo_fix(ct0[2],ct0[3])),maxo_fix(maxo_fix(ct0[4],ct0[5]),maxo_fix(ct0[6],ct0[7])));

			llr_es2[i-1] = nom - denom - (gs2[i-1]<<1);
			//llr_es1[i-1] = nom - denom;

			llr_es2[i-1] = (llr_es2[i-1]>limit)? limit:llr_es2[i-1];
			llr_es2[i-1] = (llr_es2[i-1]<-limit)? -limit:llr_es2[i-1];

			pta = alfa2 + (j+1)*StateNum;

			dt0[0] = bs0[0] + pta[0];
			dt0[1] = bs0[1] + pta[1];
			dt0[2] = bs0[2] + pta[2];
			dt0[3] = bs0[3] + pta[3];
			dt0[4] = bs0[4] + pta[4];
			dt0[5] = bs0[5] + pta[5];
			dt0[6] = bs0[6] + pta[6];
			dt0[7] = bs0[7] + pta[7];

			dt1[0] = bs1[0] + pta[0];
			dt1[1] = bs1[1] + pta[1];
			dt1[2] = bs1[2] + pta[2];
			dt1[3] = bs1[3] + pta[3];
			dt1[4] = bs1[4] + pta[4];
			dt1[5] = bs1[5] + pta[5];
			dt1[6] = bs1[6] + pta[6];
			dt1[7] = bs1[7] + pta[7];

			nom = maxo_fix(maxo_fix(maxo_fix(dt1[0],dt1[1]),maxo_fix(dt1[2],dt1[3])),maxo_fix(maxo_fix(dt1[4],dt1[5]),maxo_fix(dt1[6],dt1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(dt0[0],dt0[1]),maxo_fix(dt0[2],dt0[3])),maxo_fix(maxo_fix(dt0[4],dt0[5]),maxo_fix(dt0[6],dt0[7])));
		
			llr_es2[j+1] = nom - denom - (gs2[j+1]<<1);
			//llr_es1[j+1] = nom - denom;
			
			llr_es2[j+1] = (llr_es2[j+1]>limit)? limit:llr_es2[j+1];
			llr_es2[j+1] = (llr_es2[j+1]<-limit)? -limit:llr_es2[j+1];
		}

		{
			int		*pta, *ptb;
			int		nom, denom;
			int		dt0[8], dt1[8];

			pta = alfa2;
			ptb = beta2;

			dt0[0] = pta[0]-gs2[0]-gp2[0]+ptb[0];
			dt0[1] = pta[1]-gs2[0]-gp2[0]+ptb[4];
			dt0[2] = pta[2]-gs2[0]+gp2[0]+ptb[5];
			dt0[3] = pta[3]-gs2[0]+gp2[0]+ptb[1];
			dt0[4] = pta[4]-gs2[0]+gp2[0]+ptb[2];
			dt0[5] = pta[5]-gs2[0]+gp2[0]+ptb[6];
			dt0[6] = pta[6]-gs2[0]-gp2[0]+ptb[7];
			dt0[7] = pta[7]-gs2[0]-gp2[0]+ptb[3];

			dt1[0] = pta[0]+gs2[0]+gp2[0]+ptb[4];
			dt1[1] = pta[1]+gs2[0]+gp2[0]+ptb[0];
			dt1[2] = pta[2]+gs2[0]-gp2[0]+ptb[1];
			dt1[3] = pta[3]+gs2[0]-gp2[0]+ptb[5];
			dt1[4] = pta[4]+gs2[0]-gp2[0]+ptb[6];
			dt1[5] = pta[5]+gs2[0]-gp2[0]+ptb[2];
			dt1[6] = pta[6]+gs2[0]+gp2[0]+ptb[3];
			dt1[7] = pta[7]+gs2[0]+gp2[0]+ptb[7];

			nom = maxo_fix(maxo_fix(maxo_fix(dt1[0],dt1[1]),maxo_fix(dt1[2],dt1[3])),maxo_fix(maxo_fix(dt1[4],dt1[5]),maxo_fix(dt1[6],dt1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(dt0[0],dt0[1]),maxo_fix(dt0[2],dt0[3])),maxo_fix(maxo_fix(dt0[4],dt0[5]),maxo_fix(dt0[6],dt0[7])));
		
			llr_es2[0] = nom - denom - (gs2[0]<<1);
			//llr_es1[0] = nom - denom;
			
			llr_es2[0] = (llr_es2[0]>limit)? limit:llr_es2[0];
			llr_es2[0] = (llr_es2[0]<-limit)? -limit:llr_es2[0];

			pta = alfa2 + (Len+Mem-1)*StateNum;
			ptb = beta2 + (Len+Mem-1)*StateNum;

			dt0[0] = pta[0]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[0];
			dt0[1] = pta[1]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[4];
			dt0[2] = pta[2]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[5];
			dt0[3] = pta[3]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[1];
			dt0[4] = pta[4]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[2];
			dt0[5] = pta[5]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[6];
			dt0[6] = pta[6]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[7];
			dt0[7] = pta[7]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[3];

			dt1[0] = pta[0]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[4];
			dt1[1] = pta[1]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[0];
			dt1[2] = pta[2]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[1];
			dt1[3] = pta[3]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[5];
			dt1[4] = pta[4]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[6];
			dt1[5] = pta[5]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[2];
			dt1[6] = pta[6]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[3];
			dt1[7] = pta[7]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[7];

			nom = maxo_fix(maxo_fix(maxo_fix(dt1[0],dt1[1]),maxo_fix(dt1[2],dt1[3])),maxo_fix(maxo_fix(dt1[4],dt1[5]),maxo_fix(dt1[6],dt1[7])));
			denom = maxo_fix(maxo_fix(maxo_fix(dt0[0],dt0[1]),maxo_fix(dt0[2],dt0[3])),maxo_fix(maxo_fix(dt0[4],dt0[5]),maxo_fix(dt0[6],dt0[7])));
		
			llr_es2[Len+Mem-1] = nom - denom - (gs2[Len+Mem-1]<<1);
			//llr_es1[Len+Mem-1] = nom - denom;
			
			llr_es2[Len+Mem-1] = (llr_es2[Len+Mem-1]>limit)? limit:llr_es2[Len+Mem-1];
			llr_es2[Len+Mem-1] = (llr_es2[Len+Mem-1]<-limit)? -limit:llr_es2[Len+Mem-1];
		}
		/* end of decoder2 */

		//{
		//	FILE *fp;

		//	fp = fopen("turbo_fix2.txt","at");

		//	for(i=0; i<Len+Mem; i++)
		//	{
		//		fprintf(fp,"%8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d\n", i, llr_es2[i], g0[i], g1[i], 
		//			alfa2[i*StateNum+0], alfa2[i*StateNum+1],alfa2[i*StateNum+2], alfa2[i*StateNum+3],alfa2[i*StateNum+4], alfa2[i*StateNum+5],alfa2[i*StateNum+6], alfa2[i*StateNum+7],
		//			beta2[i*StateNum+0], beta2[i*StateNum+1],beta2[i*StateNum+2], beta2[i*StateNum+3],beta2[i*StateNum+4], beta2[i*StateNum+5],beta2[i*StateNum+6], beta2[i*StateNum+7]);
		//	}

		//	fclose(fp);
		//}

		deinterleave_fix(llr_es2, Len, interleave_table);
	}

	for(i=0; i<Len; i++)
	{
		store_llr[i] = (double) llr_es2[i]/(double) (1<<qbits);	
	}

	deinterleave_fix(llr_es1, Len, interleave_table);

	/* compute the llr_ep1 */
	for(i=0; i<Len+Mem; i++)
	{
		int nom,denom;
		int nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;	

		nom1	=	maxo_fix(alfa1[i*StateNum+0]+gs1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gs1[i]+beta1[i*StateNum+0]);
		nom2	=	maxo_fix(alfa1[i*StateNum+2]-gs1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]-gs1[i]+beta1[i*StateNum+1]);
		nom3	=	maxo_fix(alfa1[i*StateNum+4]-gs1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]-gs1[i]+beta1[i*StateNum+6]);
		nom4	=	maxo_fix(alfa1[i*StateNum+6]+gs1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gs1[i]+beta1[i*StateNum+7]);

		denom1	=	maxo_fix(alfa1[i*StateNum+0]-gs1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gs1[i]+beta1[i*StateNum+4]);
		denom2	=	maxo_fix(alfa1[i*StateNum+2]+gs1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]+gs1[i]+beta1[i*StateNum+5]);
		denom3	=	maxo_fix(alfa1[i*StateNum+4]+gs1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]+gs1[i]+beta1[i*StateNum+2]);
		denom4	=	maxo_fix(alfa1[i*StateNum+6]-gs1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gs1[i]+beta1[i*StateNum+3]);

		nom = maxo_fix(maxo_fix(nom1,nom2),maxo_fix(nom3,nom4));
		denom = maxo_fix(maxo_fix(denom1,denom2),maxo_fix(denom3,denom4));

		llr_ep1[i] = nom - denom;

		if (llr_ep1[i] > limit)
		{
			llr_ep1[i] = limit;
		}
		if (llr_ep1[i] < -limit)
		{
			llr_ep1[i] = -limit;
		}
	}

	/* compute the llr_ep2 */
	for(i=0; i<Len+Mem; i++)
	{
		int nom,denom;
		int nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;	

		nom1	=	maxo_fix(alfa2[i*StateNum+0]+gs2[i]+beta2[i*StateNum+4], alfa2[i*StateNum+1]+gs2[i]+beta2[i*StateNum+0]);
		nom2	=	maxo_fix(alfa2[i*StateNum+2]-gs2[i]+beta2[i*StateNum+5], alfa2[i*StateNum+3]-gs2[i]+beta2[i*StateNum+1]);
		nom3	=	maxo_fix(alfa2[i*StateNum+4]-gs2[i]+beta2[i*StateNum+2], alfa2[i*StateNum+5]-gs2[i]+beta2[i*StateNum+6]);
		nom4	=	maxo_fix(alfa2[i*StateNum+6]+gs2[i]+beta2[i*StateNum+3], alfa2[i*StateNum+7]+gs2[i]+beta2[i*StateNum+7]);

		denom1	=	maxo_fix(alfa2[i*StateNum+0]-gs2[i]+beta2[i*StateNum+0], alfa2[i*StateNum+1]-gs2[i]+beta2[i*StateNum+4]);
		denom2	=	maxo_fix(alfa2[i*StateNum+2]+gs2[i]+beta2[i*StateNum+1], alfa2[i*StateNum+3]+gs2[i]+beta2[i*StateNum+5]);
		denom3	=	maxo_fix(alfa2[i*StateNum+4]+gs2[i]+beta2[i*StateNum+6], alfa2[i*StateNum+5]+gs2[i]+beta2[i*StateNum+2]);
		denom4	=	maxo_fix(alfa2[i*StateNum+6]-gs2[i]+beta2[i*StateNum+7], alfa2[i*StateNum+7]-gs2[i]+beta2[i*StateNum+3]);

		nom = maxo_fix(maxo_fix(nom1,nom2),maxo_fix(nom3,nom4));
		denom = maxo_fix(maxo_fix(denom1,denom2),maxo_fix(denom3,denom4));

		llr_ep2[i] = nom - denom;

		if (llr_ep2[i] > limit)
		{
			llr_ep2[i] = limit;
		}

		if (llr_ep2[i] < -limit)
		{
			llr_ep2[i] = -limit;
		}
	}

	for(i=0; i<Len+Mem; i++)
	{
		llr_s1[i] = llr_s1[i] + llr_es1[i] + llr_es2[i];
		llr_p1[i] = llr_p1[i] + llr_ep1[i];
		llr_p2[i] = llr_p2[i] + llr_ep2[i];
	}	

	// Tailbiting
	llr_s_tail[0]	=	llr_s1[Len];
	llr_p1_tail[0]	=	llr_p1[Len];
	llr_p2_tail[0]	=	llr_s1[Len+1];

	llr_s_tail[1]	=	llr_p1[Len+1];
	llr_p1_tail[1]	=	llr_s1[Len+2];
	llr_p2_tail[1]	=	llr_p1[Len+2];

	llr_s_tail[2]	=	llr_s1[Len];
	llr_p1_tail[2]	=	llr_p2[Len];
	llr_p2_tail[2]	=	llr_s1[Len+1];

	llr_s_tail[3]	=	llr_p2[Len+1];
	llr_p1_tail[3]	=	llr_s1[Len+2];
	llr_p2_tail[3]	=	llr_p2[Len+2];

	for(i=0; i<4; i++)
	{
		llr_s1[Len+i] = llr_s_tail[i];
		llr_p1[Len+i] = llr_p1_tail[i];
		llr_p2[Len+i] = llr_p2_tail[i];
	}

	/* decision */
	for(i=0; i<Len; i++)
	{
		dec[i] = (llr_s1[i] > 0)? 1:0;
	}

	for(i=0; i<Len+4; i++)
	{
		llr_out_rm[3*i] = llr_s1[i];
		llr_out_rm[3*i+1] = llr_p1[i];
		llr_out_rm[3*i+2] = llr_p2[i];
	}

	for(i=0; i<3*(Len+4); i++)
	{
		llr_out_rm[i] = (llr_out_rm[i]>64)? 64:llr_out_rm[i];
		llr_out_rm[i] = (llr_out_rm[i]<-64)? -64:llr_out_rm[i];			
	}

	rm_interleave_fix(Len+4, CodeLen, llr_out_rm, llr_out_fix, ratematch_table);

	for(i=0; i<CodeLen; i++)
	{
		llr_out[i] = (double) llr_out_fix[i]/(double) (1<<qbits);
	}

	free(llr_in_rm);
	free(llr_out_rm);
	free(llr_s1);
	free(llr_s2);
	free(llr_p1);
	free(llr_p2);
	free(llr_es1);
	free(llr_es2);
	free(llr_ep1);
	free(llr_ep2);
	free(alfa1);
	free(beta1);
	free(gs1);
	free(gp1);
	free(alfa2);
	free(beta2);
	free(gs2);
	free(gp2);
	free(llr_in_fix);
	free(llr_out_fix);
	free(g0);
	free(g1);
	free(tmp_llr_es1);
}

void	rate_dematching_short(int Len, int CodeLen, short *input, short *output, int *inter_index)
{
	int i;
	int j;

	for(i=0; i<3*Len; i++)
	{
		output[i] = 0;
	}
	for(i=0; i<CodeLen; i++)
	{
		
		output[inter_index[i]] += input[i];
	}
}

//void	rate_dematching_short2(int Len, int CodeLen, short *input, short *output, int *inter_index)
//{
//	int			i, j;
//	int			*tmp_32 = (int*) _mm_malloc(3 * Len*sizeof(int), sizeof(__m256i));
//	__m256i		zeros_256i;
//	__m256i		*pout = (__m256i*) output;
//	__m256i		*p256i;
//
//	for (i=0; i < 3 * Len; i++)
//	{
//		tmp_32[i] = 0;
//	}
//
//	for (i = 0; i < CodeLen; i += 16)
//	{
//		__m256i vindex, outputa_256i, outputb_256i, output16_256i, input16_256i;
//
//		vindex = _mm256_load_si256((__m256i*)&inter_index[i]);
//
//		outputa_256i = _mm256_i32gather_epi32(tmp_32, vindex, 1);
//
//		vindex = _mm256_load_si256((__m256i*)&inter_index[i+8]);
//
//		outputb_256i = _mm256_i32gather_epi32(tmp_32, vindex, 1);
//
//		output16_256i = _mm256_packs_epi32(outputa_256i, outputb_256i);
//
//		output16_256i = _mm256_permute4x64_epi64(output16_256i, 0xd8);
//
//		input16_256i = _mm256_load_si256((__m256i*)&input[i]);
//
//		output16_256i = _mm256_add_epi16(output16_256i, input16_256i);
//
//		_mm256_store_si256((__m256i*)&output[i], output16_256i);
//	}
//
//	for (i; i<CodeLen; i++)
//	{
//		output[inter_index[i]] += input[i];
//	}
//
//	_mm_free(tmp_32);
//}

void	rm_interleave_short(int CodeLen, short *input, short *output, int *inter_index)
{
	int i;
	
	for(i=0; i<CodeLen; i++)
	{
		output[i] = input[inter_index[i]];
	}
}

void turbo_decode_rm_avx2(short* dec, 
						double* llr_out, 
						double* llr_in, 
						double* store_llr, 
						int Len, 
						int CodeLen,
						int iterCNum, 
						int* interleave_table,
						int* ratematch_table)
{
	short		*llr_in_fix, *llr_out_fix;
	short		*llr_in_rm, *llr_out_rm, *llr_s1, *llr_s2, *llr_p1, *llr_p2, *llr_es1, *llr_es2, *llr_ep1, *llr_ep2;
	short		*alfa1, *beta1, *gs1, *gp1, *alfa2, *beta2, *gs2, *gp2, *g0, *g1;

	int			i, s, iter_idx;
	short		llr_s_tail[4], llr_p1_tail[4], llr_p2_tail[4];
	int			inf = 4095;
	int			limit = 255;
	int			qbits = 4;
	int			bitwidth = 15;
	int			StateNum = 8;

	short		*tmp_llr_es1;

	tmp_llr_es1		=	(short*) malloc((Len+Mem)*sizeof(short));

	llr_in_fix	=	(short*) malloc(CodeLen*sizeof(short));
	llr_out_fix	=	(short*) malloc(CodeLen*sizeof(short));

	llr_in_rm	=	(short*) malloc(3*(Len+4)*sizeof(short));
	llr_out_rm	=	(short*) malloc(3*(Len+4)*sizeof(short));
	llr_s1		=	(short*) malloc((Len+4)*sizeof(short));
	llr_s2		=	(short*) malloc((Len+4)*sizeof(short));
	llr_p1		=	(short*) malloc((Len+4)*sizeof(short));
	llr_p2		=	(short*) malloc((Len+4)*sizeof(short));
	llr_es1		=	(short*) malloc((Len+Mem)*sizeof(short));
	llr_es2		=	(short*) malloc((Len+Mem)*sizeof(short));
	llr_ep1		=	(short*) malloc((Len+Mem)*sizeof(short));
	llr_ep2		=	(short*) malloc((Len+Mem)*sizeof(short));

	alfa1 = (short*) malloc((Len+Mem)*StateNum*sizeof(short));
	beta1 = (short*) malloc((Len+Mem)*StateNum*sizeof(short));
	alfa2 = (short*) malloc((Len+Mem)*StateNum*sizeof(short));
	beta2 = (short*) malloc((Len+Mem)*StateNum*sizeof(short));
	gs1 = (short*) malloc((Len+Mem)*sizeof(short));
	gp1 = (short*) malloc((Len+Mem)*sizeof(short));
	gs2 = (short*) malloc((Len+Mem)*sizeof(short));
	gp2 = (short*) malloc((Len+Mem)*sizeof(short));
	g0 = (short*) malloc((Len+Mem)*sizeof(short));
	g1 = (short*) malloc((Len+Mem)*sizeof(short));

	for(i=0; i<CodeLen; i++)
	{
		llr_in_fix[i] = (short) (llr_in[i]*(1<<qbits));

		if (llr_in_fix[i] > limit)
		{
			llr_in_fix[i] = limit;
		}

		if (llr_in_fix[i] < -limit)
		{
			llr_in_fix[i] = -limit;
		}

	}

	//rate match....
	rate_dematching_short(Len+4, CodeLen, llr_in_fix, llr_in_rm, ratematch_table);

	for(i=0; i<Len+4; i++)
	{
		llr_s1[i] = llr_in_rm[3*i];
		llr_s2[i] = llr_in_rm[3*i];
		llr_p1[i] = llr_in_rm[3*i+1];
		llr_p2[i] = llr_in_rm[3*i+2];
	}

	//process tail bit
	for(i=0; i<4; i++)
	{
		llr_s_tail[i] = llr_s1[Len+i];
		llr_p1_tail[i] = llr_p1[Len+i];
		llr_p2_tail[i] = llr_p2[Len+i];
	}
	llr_s1[Len]		= llr_s_tail[0];
	llr_s1[Len+1]	= llr_p2_tail[0];
	llr_s1[Len+2]	= llr_p1_tail[1];

	llr_p1[Len]		= llr_p1_tail[0];
	llr_p1[Len+1]	= llr_s_tail[1];
	llr_p1[Len+2]	= llr_p2_tail[1];

	llr_s2[Len]		= llr_s_tail[2];
	llr_s2[Len+1]	= llr_p2_tail[2];
	llr_s2[Len+2]	= llr_p1_tail[3];

	llr_p2[Len]		= llr_p1_tail[2];
	llr_p2[Len+1]	= llr_s_tail[3];
	llr_p2[Len+2]	= llr_p2_tail[3];

	interleave_short(llr_s2, Len, interleave_table);

	for(i=0; i<Len; i++)
	{
		llr_es2[i] = (short) (store_llr[i]*(1<<qbits));
	}

	for(iter_idx=0; iter_idx<iterCNum; iter_idx++)
	{
		// decoder 1
		// initialize for alfa and beta 
		alfa1[0] = 0;
		alfa1[1] = -inf;
		alfa1[2] = -inf;
		alfa1[3] = -inf;
		alfa1[4] = -inf;
		alfa1[5] = -inf;
		alfa1[6] = -inf;
		alfa1[7] = -inf;

		beta1[(Len+Mem-1)*StateNum+0] = 0;
		beta1[(Len+Mem-1)*StateNum+1] = -inf;
		beta1[(Len+Mem-1)*StateNum+2] = -inf;
		beta1[(Len+Mem-1)*StateNum+3] = -inf;
		beta1[(Len+Mem-1)*StateNum+4] = -inf;
		beta1[(Len+Mem-1)*StateNum+5] = -inf;
		beta1[(Len+Mem-1)*StateNum+6] = -inf;
		beta1[(Len+Mem-1)*StateNum+7] = -inf;

		/* initialize for gs and gp */
		for(i=0; i<Len; i++)
		{
			gs1[i] = (llr_s1[i]>>1) + (llr_es2[i]>>1);
			gp1[i] = (llr_p1[i]>>1);
			
			g0[i] = - gs1[i] -  gp1[i]; 
			g1[i] = - gs1[i] +  gp1[i]; 
		}
		for(i=Len; i<Len+Mem; i++)
		{
			gs1[i] = (llr_s1[i]>>1);
			gp1[i] = (llr_p1[i]>>1);

			g0[i] = - gs1[i] -  gp1[i]; 
			g1[i] = - gs1[i] +  gp1[i]; 
		}

		/* compute the alfa and beta */
		for(i=1; i<Len+Mem; i++)
		{
			int		j = Len+Mem-1-i;
			short		at0[8], at1[8], as0[8], as1[8];
			short		bt0[8], bt1[8], bs0[8], bs1[8];
			short		ct0[8], ct1[8], dt0[8], dt1[8];

			short		*pa = alfa1 + (i-1)*StateNum;
			short		*pb = beta1 + (j+1)*StateNum;

			short		*pta, *ptb;

			short		nom,denom;

			//if (i>=j-1)
			//{
			//	printf("\n");
			//}

			// add
			at0[0] = pa[0] + g0[i-1];
			at0[1] = pa[1] + g0[i-1];
			at0[2] = pa[2] + g1[i-1];
			at0[3] = pa[3] + g1[i-1];
			at0[4] = pa[4] + g1[i-1];
			at0[5] = pa[5] + g1[i-1];
			at0[6] = pa[6] + g0[i-1];
			at0[7] = pa[7] + g0[i-1];

			at1[0] = pa[0] - g0[i-1];
			at1[1] = pa[1] - g0[i-1];
			at1[2] = pa[2] - g1[i-1];
			at1[3] = pa[3] - g1[i-1];
			at1[4] = pa[4] - g1[i-1];
			at1[5] = pa[5] - g1[i-1];
			at1[6] = pa[6] - g0[i-1];
			at1[7] = pa[7] - g0[i-1];

			// shuffle
			as0[0] = at0[0];
			as0[1] = at0[3];
			as0[2] = at0[4];
			as0[3] = at0[7];
			as0[4] = at0[1];
			as0[5] = at0[2];
			as0[6] = at0[5];
			as0[7] = at0[6];

			as1[0] = at1[1];
			as1[1] = at1[2];
			as1[2] = at1[5];
			as1[3] = at1[6];
			as1[4] = at1[0];
			as1[5] = at1[3];
			as1[6] = at1[4];
			as1[7] = at1[7];

			// max
			pa = pa + StateNum;

			//update alfa
			pa[0] = (as0[0]>as1[0])?as0[0]:as1[0];
			pa[1] = (as0[1]>as1[1])?as0[1]:as1[1];
			pa[2] = (as0[2]>as1[2])?as0[2]:as1[2];
			pa[3] = (as0[3]>as1[3])?as0[3]:as1[3];
			pa[4] = (as0[4]>as1[4])?as0[4]:as1[4];
			pa[5] = (as0[5]>as1[5])?as0[5]:as1[5];
			pa[6] = (as0[6]>as1[6])?as0[6]:as1[6];
			pa[7] = (as0[7]>as1[7])?as0[7]:as1[7];

			//nomalize
			if (i % 16 == 0)
			{
				short	tmp = pa[0];
				
				pa[0] -= tmp;
				pa[1] -= tmp;
				pa[2] -= tmp;
				pa[3] -= tmp;
				pa[4] -= tmp;
				pa[5] -= tmp;
				pa[6] -= tmp;
				pa[7] -= tmp;
			}
			

			//alfa1[i*StateNum+0] = ACSO_fix(alfa1[(i-1)*StateNum+0], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+1], +gs1[i-1]+gp1[i-1]);
			//alfa1[i*StateNum+1] = ACSO_fix(alfa1[(i-1)*StateNum+3], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+2], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+2] = ACSO_fix(alfa1[(i-1)*StateNum+4], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+5], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+3] = ACSO_fix(alfa1[(i-1)*StateNum+7], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+6], +gs1[i-1]+gp1[i-1]);
			//alfa1[i*StateNum+4] = ACSO_fix(alfa1[(i-1)*StateNum+1], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+0], +gs1[i-1]+gp1[i-1]);
			//alfa1[i*StateNum+5] = ACSO_fix(alfa1[(i-1)*StateNum+2], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+3], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+6] = ACSO_fix(alfa1[(i-1)*StateNum+5], -gs1[i-1]+gp1[i-1], alfa1[(i-1)*StateNum+4], +gs1[i-1]-gp1[i-1]);
			//alfa1[i*StateNum+7] = ACSO_fix(alfa1[(i-1)*StateNum+6], -gs1[i-1]-gp1[i-1], alfa1[(i-1)*StateNum+7], +gs1[i-1]+gp1[i-1]);
		
			//fix_GetInt(alfa1[i*StateNum+0],bitwidth);
			//fix_GetInt(alfa1[i*StateNum+1],bitwidth);
			//fix_GetInt(alfa1[i*StateNum+2],bitwidth);
			//fix_GetInt(alfa1[i*StateNum+3],bitwidth);
			//fix_GetInt(alfa1[i*StateNum+4],bitwidth);
			//fix_GetInt(alfa1[i*StateNum+5],bitwidth);
			//fix_GetInt(alfa1[i*StateNum+6],bitwidth);
			//fix_GetInt(alfa1[i*StateNum+7],bitwidth);

			// beta

			// add
			bt0[0] = pb[0] + g0[j+1];
			bt0[1] = pb[1] + g1[j+1];
			bt0[2] = pb[2] + g1[j+1];
			bt0[3] = pb[3] + g0[j+1];
			bt0[4] = pb[4] + g0[j+1];
			bt0[5] = pb[5] + g1[j+1];
			bt0[6] = pb[6] + g1[j+1];
			bt0[7] = pb[7] + g0[j+1];

			bt1[0] = pb[0] - g0[j+1];
			bt1[1] = pb[1] - g1[j+1];
			bt1[2] = pb[2] - g1[j+1];
			bt1[3] = pb[3] - g0[j+1];
			bt1[4] = pb[4] - g0[j+1];
			bt1[5] = pb[5] - g1[j+1];
			bt1[6] = pb[6] - g1[j+1];
			bt1[7] = pb[7] - g0[j+1];

			// shuffle
			bs0[0] = bt0[0];
			bs0[1] = bt0[4];
			bs0[2] = bt0[5];
			bs0[3] = bt0[1];
			bs0[4] = bt0[2];
			bs0[5] = bt0[6];
			bs0[6] = bt0[7];
			bs0[7] = bt0[3];

			bs1[0] = bt1[4];
			bs1[1] = bt1[0];
			bs1[2] = bt1[1];
			bs1[3] = bt1[5];
			bs1[4] = bt1[6];
			bs1[5] = bt1[2];
			bs1[6] = bt1[3];
			bs1[7] = bt1[7];

			// max
			pb = pb - StateNum;

			pb[0] = (bs0[0]>bs1[0])?bs0[0]:bs1[0];
			pb[1] = (bs0[1]>bs1[1])?bs0[1]:bs1[1];
			pb[2] = (bs0[2]>bs1[2])?bs0[2]:bs1[2];
			pb[3] = (bs0[3]>bs1[3])?bs0[3]:bs1[3];
			pb[4] = (bs0[4]>bs1[4])?bs0[4]:bs1[4];
			pb[5] = (bs0[5]>bs1[5])?bs0[5]:bs1[5];
			pb[6] = (bs0[6]>bs1[6])?bs0[6]:bs1[6];
			pb[7] = (bs0[7]>bs1[7])?bs0[7]:bs1[7];

			//nomalize
			if (i % 16 == 0)
			{
				short	tmp = pb[0];
				
				pb[0] -= tmp;
				pb[1] -= tmp;
				pb[2] -= tmp;
				pb[3] -= tmp;
				pb[4] -= tmp;
				pb[5] -= tmp;
				pb[6] -= tmp;
				pb[7] -= tmp;
			}

	/*		fix_GetInt(beta1[j*StateNum+0],bitwidth);
			fix_GetInt(beta1[j*StateNum+1],bitwidth);
			fix_GetInt(beta1[j*StateNum+2],bitwidth);
			fix_GetInt(beta1[j*StateNum+3],bitwidth);
			fix_GetInt(beta1[j*StateNum+4],bitwidth);
			fix_GetInt(beta1[j*StateNum+5],bitwidth);
			fix_GetInt(beta1[j*StateNum+6],bitwidth);
			fix_GetInt(beta1[j*StateNum+7],bitwidth);*/

			// llres1
			ptb = beta1 + (i-1)*StateNum;

			ct0[0] = as0[0] + ptb[0];
			ct0[1] = as0[1] + ptb[1];
			ct0[2] = as0[2] + ptb[2];
			ct0[3] = as0[3] + ptb[3];
			ct0[4] = as0[4] + ptb[4];
			ct0[5] = as0[5] + ptb[5];
			ct0[6] = as0[6] + ptb[6];
			ct0[7] = as0[7] + ptb[7];

			ct1[0] = as1[0] + ptb[0];
			ct1[1] = as1[1] + ptb[1];
			ct1[2] = as1[2] + ptb[2];
			ct1[3] = as1[3] + ptb[3];
			ct1[4] = as1[4] + ptb[4];
			ct1[5] = as1[5] + ptb[5];
			ct1[6] = as1[6] + ptb[6];
			ct1[7] = as1[7] + ptb[7];

			nom = maxo_avx(maxo_avx(maxo_avx(ct1[0],ct1[1]),maxo_avx(ct1[2],ct1[3])),maxo_avx(maxo_avx(ct1[4],ct1[5]),maxo_avx(ct1[6],ct1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(ct0[0],ct0[1]),maxo_avx(ct0[2],ct0[3])),maxo_avx(maxo_avx(ct0[4],ct0[5]),maxo_avx(ct0[6],ct0[7])));

			llr_es1[i-1] = nom - denom - (gs1[i-1]<<1);
			//llr_es1[i-1] = nom - denom;

			llr_es1[i-1] = (llr_es1[i-1]>limit)? limit:llr_es1[i-1];
			llr_es1[i-1] = (llr_es1[i-1]<-limit)? -limit:llr_es1[i-1];

			pta = alfa1 + (j+1)*StateNum;

			dt0[0] = bs0[0] + pta[0];
			dt0[1] = bs0[1] + pta[1];
			dt0[2] = bs0[2] + pta[2];
			dt0[3] = bs0[3] + pta[3];
			dt0[4] = bs0[4] + pta[4];
			dt0[5] = bs0[5] + pta[5];
			dt0[6] = bs0[6] + pta[6];
			dt0[7] = bs0[7] + pta[7];

			dt1[0] = bs1[0] + pta[0];
			dt1[1] = bs1[1] + pta[1];
			dt1[2] = bs1[2] + pta[2];
			dt1[3] = bs1[3] + pta[3];
			dt1[4] = bs1[4] + pta[4];
			dt1[5] = bs1[5] + pta[5];
			dt1[6] = bs1[6] + pta[6];
			dt1[7] = bs1[7] + pta[7];

			nom = maxo_avx(maxo_avx(maxo_avx(dt1[0],dt1[1]),maxo_avx(dt1[2],dt1[3])),maxo_avx(maxo_avx(dt1[4],dt1[5]),maxo_avx(dt1[6],dt1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(dt0[0],dt0[1]),maxo_avx(dt0[2],dt0[3])),maxo_avx(maxo_avx(dt0[4],dt0[5]),maxo_avx(dt0[6],dt0[7])));
		
			llr_es1[j+1] = nom - denom - (gs1[j+1]<<1);
			//llr_es1[j+1] = nom - denom;
			
			llr_es1[j+1] = (llr_es1[j+1]>limit)? limit:llr_es1[j+1];
			llr_es1[j+1] = (llr_es1[j+1]<-limit)? -limit:llr_es1[j+1];
		}

		{
			short		*pta, *ptb;
			short		nom, denom;
			short		dt0[8], dt1[8];

			pta = alfa1;
			ptb = beta1;

			dt0[0] = pta[0]-gs1[0]-gp1[0]+ptb[0];
			dt0[1] = pta[1]-gs1[0]-gp1[0]+ptb[4];
			dt0[2] = pta[2]-gs1[0]+gp1[0]+ptb[5];
			dt0[3] = pta[3]-gs1[0]+gp1[0]+ptb[1];
			dt0[4] = pta[4]-gs1[0]+gp1[0]+ptb[2];
			dt0[5] = pta[5]-gs1[0]+gp1[0]+ptb[6];
			dt0[6] = pta[6]-gs1[0]-gp1[0]+ptb[7];
			dt0[7] = pta[7]-gs1[0]-gp1[0]+ptb[3];

			dt1[0] = pta[0]+gs1[0]+gp1[0]+ptb[4];
			dt1[1] = pta[1]+gs1[0]+gp1[0]+ptb[0];
			dt1[2] = pta[2]+gs1[0]-gp1[0]+ptb[1];
			dt1[3] = pta[3]+gs1[0]-gp1[0]+ptb[5];
			dt1[4] = pta[4]+gs1[0]-gp1[0]+ptb[6];
			dt1[5] = pta[5]+gs1[0]-gp1[0]+ptb[2];
			dt1[6] = pta[6]+gs1[0]+gp1[0]+ptb[3];
			dt1[7] = pta[7]+gs1[0]+gp1[0]+ptb[7];

			nom = maxo_avx(maxo_avx(maxo_avx(dt1[0],dt1[1]),maxo_avx(dt1[2],dt1[3])),maxo_avx(maxo_avx(dt1[4],dt1[5]),maxo_avx(dt1[6],dt1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(dt0[0],dt0[1]),maxo_avx(dt0[2],dt0[3])),maxo_avx(maxo_avx(dt0[4],dt0[5]),maxo_avx(dt0[6],dt0[7])));
		
			llr_es1[0] = nom - denom - (gs1[0]<<1);
			//llr_es1[0] = nom - denom;
			
			llr_es1[0] = (llr_es1[0]>limit)? limit:llr_es1[0];
			llr_es1[0] = (llr_es1[0]<-limit)? -limit:llr_es1[0];

			pta = alfa1 + (Len+Mem-1)*StateNum;
			ptb = beta1 + (Len+Mem-1)*StateNum;

			dt0[0] = pta[0]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[0];
			dt0[1] = pta[1]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[4];
			dt0[2] = pta[2]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[5];
			dt0[3] = pta[3]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[1];
			dt0[4] = pta[4]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[2];
			dt0[5] = pta[5]-gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[6];
			dt0[6] = pta[6]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[7];
			dt0[7] = pta[7]-gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[3];

			dt1[0] = pta[0]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[4];
			dt1[1] = pta[1]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[0];
			dt1[2] = pta[2]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[1];
			dt1[3] = pta[3]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[5];
			dt1[4] = pta[4]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[6];
			dt1[5] = pta[5]+gs1[Len+Mem-1]-gp1[Len+Mem-1]+ptb[2];
			dt1[6] = pta[6]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[3];
			dt1[7] = pta[7]+gs1[Len+Mem-1]+gp1[Len+Mem-1]+ptb[7];

			nom = maxo_avx(maxo_avx(maxo_avx(dt1[0],dt1[1]),maxo_avx(dt1[2],dt1[3])),maxo_avx(maxo_avx(dt1[4],dt1[5]),maxo_avx(dt1[6],dt1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(dt0[0],dt0[1]),maxo_avx(dt0[2],dt0[3])),maxo_avx(maxo_avx(dt0[4],dt0[5]),maxo_avx(dt0[6],dt0[7])));
		
			llr_es1[Len+Mem-1] = nom - denom - (gs1[Len+Mem-1]<<1);
			//llr_es1[Len+Mem-1] = nom - denom;
			
			llr_es1[Len+Mem-1] = (llr_es1[Len+Mem-1]>limit)? limit:llr_es1[Len+Mem-1];
			llr_es1[Len+Mem-1] = (llr_es1[Len+Mem-1]<-limit)? -limit:llr_es1[Len+Mem-1];
		}

		/* compute the llres1 */
		for(i=0; i<Len+Mem; i++)
		{
			short nom,denom;
			short nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;		

			nom1	=	maxo_fix(alfa1[i*StateNum+0]+gp1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gp1[i]+beta1[i*StateNum+0]);
			nom2	=	maxo_fix(alfa1[i*StateNum+2]-gp1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]-gp1[i]+beta1[i*StateNum+5]);
			nom3	=	maxo_fix(alfa1[i*StateNum+4]-gp1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]-gp1[i]+beta1[i*StateNum+2]);
			nom4	=	maxo_fix(alfa1[i*StateNum+6]+gp1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gp1[i]+beta1[i*StateNum+7]);

			denom1	=	maxo_fix(alfa1[i*StateNum+0]-gp1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gp1[i]+beta1[i*StateNum+4]);
			denom2	=	maxo_fix(alfa1[i*StateNum+2]+gp1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]+gp1[i]+beta1[i*StateNum+1]);
			denom3	=	maxo_fix(alfa1[i*StateNum+4]+gp1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]+gp1[i]+beta1[i*StateNum+6]);
			denom4	=	maxo_fix(alfa1[i*StateNum+6]-gp1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gp1[i]+beta1[i*StateNum+3]);

		/*	nom1	=	maxo_avx(alfa1[i*StateNum+0]+gs1[i]+gp1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gs1[i]+gp1[i]+beta1[i*StateNum+0]);
			nom2	=	maxo_avx(alfa1[i*StateNum+2]+gs1[i]-gp1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]+gs1[i]-gp1[i]+beta1[i*StateNum+5]);
			nom3	=	maxo_avx(alfa1[i*StateNum+4]+gs1[i]-gp1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]+gs1[i]-gp1[i]+beta1[i*StateNum+2]);
			nom4	=	maxo_avx(alfa1[i*StateNum+6]+gs1[i]+gp1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gs1[i]+gp1[i]+beta1[i*StateNum+7]);

			denom1	=	maxo_avx(alfa1[i*StateNum+0]-gs1[i]-gp1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gs1[i]-gp1[i]+beta1[i*StateNum+4]);
			denom2	=	maxo_avx(alfa1[i*StateNum+2]-gs1[i]+gp1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]-gs1[i]+gp1[i]+beta1[i*StateNum+1]);
			denom3	=	maxo_avx(alfa1[i*StateNum+4]-gs1[i]+gp1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]-gs1[i]+gp1[i]+beta1[i*StateNum+6]);
			denom4	=	maxo_avx(alfa1[i*StateNum+6]-gs1[i]-gp1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gs1[i]-gp1[i]+beta1[i*StateNum+3]);*/

			nom = maxo_avx(maxo_avx(nom1,nom2),maxo_avx(nom3,nom4));
			denom = maxo_avx(maxo_avx(denom1,denom2),maxo_avx(denom3,denom4));

			tmp_llr_es1[i] = nom - denom;
			
			if (tmp_llr_es1[i] > limit)
			{
				tmp_llr_es1[i] = limit;
			}

			if (tmp_llr_es1[i] < -limit)
			{
				tmp_llr_es1[i] = -limit;
			}
		}

		/*for(i=0; i<Len; i++)
		{
			printf("%d, %d, %d \n", i, llr_es1[i], tmp_llr_es1[i]);
		}*/

		//{
		//	FILE *fp;

		//	fp = fopen("turbo_avx.txt","at");

		//	for(i=0; i<Len+Mem; i++)
		//	{
		//		fprintf(fp,"%8d, %8d, %8d, %8d, %8d, %8d\n", i, llr_es1[i], g0[i], g1[i], alfa1[i*StateNum+0], alfa1[i*StateNum+1]);
		//	}

		//	fclose(fp);
		//}

		interleave_short(llr_es1, Len, interleave_table);

		/* decoder2 */
		/* initialize for alfa and beta */
		alfa2[0] = 0;
		alfa2[1] = -inf;
		alfa2[2] = -inf;
		alfa2[3] = -inf;
		alfa2[4] = -inf;
		alfa2[5] = -inf;
		alfa2[6] = -inf;
		alfa2[7] = -inf;

		beta2[(Len+Mem-1)*StateNum+0] = 0;
		beta2[(Len+Mem-1)*StateNum+1] = -inf;
		beta2[(Len+Mem-1)*StateNum+2] = -inf;
		beta2[(Len+Mem-1)*StateNum+3] = -inf;
		beta2[(Len+Mem-1)*StateNum+4] = -inf;
		beta2[(Len+Mem-1)*StateNum+5] = -inf;
		beta2[(Len+Mem-1)*StateNum+6] = -inf;
		beta2[(Len+Mem-1)*StateNum+7] = -inf;

		/* initialize for gp and gs */
		for(i=0; i<Len; i++)
		{
			gs2[i] = (llr_s2[i]>>1) + (llr_es1[i]>>1);
			gp2[i] = (llr_p2[i]>>1);

			g0[i] = - gs2[i] -  gp2[i]; 
			g1[i] = - gs2[i] +  gp2[i]; 
		}
		for(i=Len; i<Len+Mem; i++)
		{
			gs2[i] = (llr_s2[i]>>1);
			gp2[i] = (llr_p2[i]>>1);

			g0[i] = - gs2[i] -  gp2[i]; 
			g1[i] = - gs2[i] +  gp2[i]; 
		}

		/* compute the alfa2 */
		/* compute the alfa2, beta2 and LLR */
		for(i=1; i<Len+Mem; i++)
		{
			int			j = Len+Mem-1-i;
			short		at0[8], at1[8], as0[8], as1[8];
			short		bt0[8], bt1[8], bs0[8], bs1[8];
			short		ct0[8], ct1[8], dt0[8], dt1[8];

			short		*pa = alfa2 + (i-1)*StateNum;
			short		*pb = beta2 + (j+1)*StateNum;

			short		*pta, *ptb;

			short		nom,denom;

			//if (i>=j-1)
			//{
			//	printf("\n");
			//}

			// add
			at0[0] = pa[0] + g0[i-1];
			at0[1] = pa[1] + g0[i-1];
			at0[2] = pa[2] + g1[i-1];
			at0[3] = pa[3] + g1[i-1];
			at0[4] = pa[4] + g1[i-1];
			at0[5] = pa[5] + g1[i-1];
			at0[6] = pa[6] + g0[i-1];
			at0[7] = pa[7] + g0[i-1];

			at1[0] = pa[0] - g0[i-1];
			at1[1] = pa[1] - g0[i-1];
			at1[2] = pa[2] - g1[i-1];
			at1[3] = pa[3] - g1[i-1];
			at1[4] = pa[4] - g1[i-1];
			at1[5] = pa[5] - g1[i-1];
			at1[6] = pa[6] - g0[i-1];
			at1[7] = pa[7] - g0[i-1];

			// shuffle
			as0[0] = at0[0];
			as0[1] = at0[3];
			as0[2] = at0[4];
			as0[3] = at0[7];
			as0[4] = at0[1];
			as0[5] = at0[2];
			as0[6] = at0[5];
			as0[7] = at0[6];

			as1[0] = at1[1];
			as1[1] = at1[2];
			as1[2] = at1[5];
			as1[3] = at1[6];
			as1[4] = at1[0];
			as1[5] = at1[3];
			as1[6] = at1[4];
			as1[7] = at1[7];

			// max
			pa = pa + StateNum;

			pa[0] = (as0[0]>as1[0])?as0[0]:as1[0];
			pa[1] = (as0[1]>as1[1])?as0[1]:as1[1];
			pa[2] = (as0[2]>as1[2])?as0[2]:as1[2];
			pa[3] = (as0[3]>as1[3])?as0[3]:as1[3];
			pa[4] = (as0[4]>as1[4])?as0[4]:as1[4];
			pa[5] = (as0[5]>as1[5])?as0[5]:as1[5];
			pa[6] = (as0[6]>as1[6])?as0[6]:as1[6];
			pa[7] = (as0[7]>as1[7])?as0[7]:as1[7];

			//nomalize
			if (i % 16 == 0)
			{
				short	tmp = pa[0];
				
				pa[0] -= tmp;
				pa[1] -= tmp;
				pa[2] -= tmp;
				pa[3] -= tmp;
				pa[4] -= tmp;
				pa[5] -= tmp;
				pa[6] -= tmp;
				pa[7] -= tmp;
			}

			//fix_GetInt(alfa2[i*StateNum+0],bitwidth);
			//fix_GetInt(alfa2[i*StateNum+1],bitwidth);
			//fix_GetInt(alfa2[i*StateNum+2],bitwidth);
			//fix_GetInt(alfa2[i*StateNum+3],bitwidth);
			//fix_GetInt(alfa2[i*StateNum+4],bitwidth);
			//fix_GetInt(alfa2[i*StateNum+5],bitwidth);
			//fix_GetInt(alfa2[i*StateNum+6],bitwidth);
			//fix_GetInt(alfa2[i*StateNum+7],bitwidth);

			// beta

			// add
			bt0[0] = pb[0] + g0[j+1];
			bt0[1] = pb[1] + g1[j+1];
			bt0[2] = pb[2] + g1[j+1];
			bt0[3] = pb[3] + g0[j+1];
			bt0[4] = pb[4] + g0[j+1];
			bt0[5] = pb[5] + g1[j+1];
			bt0[6] = pb[6] + g1[j+1];
			bt0[7] = pb[7] + g0[j+1];

			bt1[0] = pb[0] - g0[j+1];
			bt1[1] = pb[1] - g1[j+1];
			bt1[2] = pb[2] - g1[j+1];
			bt1[3] = pb[3] - g0[j+1];
			bt1[4] = pb[4] - g0[j+1];
			bt1[5] = pb[5] - g1[j+1];
			bt1[6] = pb[6] - g1[j+1];
			bt1[7] = pb[7] - g0[j+1];

			// shuffle
			bs0[0] = bt0[0];
			bs0[1] = bt0[4];
			bs0[2] = bt0[5];
			bs0[3] = bt0[1];
			bs0[4] = bt0[2];
			bs0[5] = bt0[6];
			bs0[6] = bt0[7];
			bs0[7] = bt0[3];

			bs1[0] = bt1[4];
			bs1[1] = bt1[0];
			bs1[2] = bt1[1];
			bs1[3] = bt1[5];
			bs1[4] = bt1[6];
			bs1[5] = bt1[2];
			bs1[6] = bt1[3];
			bs1[7] = bt1[7];

			// max
			pb = pb - StateNum;

			pb[0] = (bs0[0]>bs1[0])?bs0[0]:bs1[0];
			pb[1] = (bs0[1]>bs1[1])?bs0[1]:bs1[1];
			pb[2] = (bs0[2]>bs1[2])?bs0[2]:bs1[2];
			pb[3] = (bs0[3]>bs1[3])?bs0[3]:bs1[3];
			pb[4] = (bs0[4]>bs1[4])?bs0[4]:bs1[4];
			pb[5] = (bs0[5]>bs1[5])?bs0[5]:bs1[5];
			pb[6] = (bs0[6]>bs1[6])?bs0[6]:bs1[6];
			pb[7] = (bs0[7]>bs1[7])?bs0[7]:bs1[7];

			//nomalize
			if (i % 16 == 0)
			{
				short	tmp = pb[0];
				
				pb[0] -= tmp;
				pb[1] -= tmp;
				pb[2] -= tmp;
				pb[3] -= tmp;
				pb[4] -= tmp;
				pb[5] -= tmp;
				pb[6] -= tmp;
				pb[7] -= tmp;
			}

			//fix_GetInt(beta2[j*StateNum+0],bitwidth);
			//fix_GetInt(beta2[j*StateNum+1],bitwidth);
			//fix_GetInt(beta2[j*StateNum+2],bitwidth);
			//fix_GetInt(beta2[j*StateNum+3],bitwidth);
			//fix_GetInt(beta2[j*StateNum+4],bitwidth);
			//fix_GetInt(beta2[j*StateNum+5],bitwidth);
			//fix_GetInt(beta2[j*StateNum+6],bitwidth);
			//fix_GetInt(beta2[j*StateNum+7],bitwidth);

			// llres2
			ptb = beta2 + (i-1)*StateNum;

			ct0[0] = as0[0] + ptb[0];
			ct0[1] = as0[1] + ptb[1];
			ct0[2] = as0[2] + ptb[2];
			ct0[3] = as0[3] + ptb[3];
			ct0[4] = as0[4] + ptb[4];
			ct0[5] = as0[5] + ptb[5];
			ct0[6] = as0[6] + ptb[6];
			ct0[7] = as0[7] + ptb[7];

			ct1[0] = as1[0] + ptb[0];
			ct1[1] = as1[1] + ptb[1];
			ct1[2] = as1[2] + ptb[2];
			ct1[3] = as1[3] + ptb[3];
			ct1[4] = as1[4] + ptb[4];
			ct1[5] = as1[5] + ptb[5];
			ct1[6] = as1[6] + ptb[6];
			ct1[7] = as1[7] + ptb[7];

			nom = maxo_avx(maxo_avx(maxo_avx(ct1[0],ct1[1]),maxo_avx(ct1[2],ct1[3])),maxo_avx(maxo_avx(ct1[4],ct1[5]),maxo_avx(ct1[6],ct1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(ct0[0],ct0[1]),maxo_avx(ct0[2],ct0[3])),maxo_avx(maxo_avx(ct0[4],ct0[5]),maxo_avx(ct0[6],ct0[7])));

			llr_es2[i-1] = nom - denom - (gs2[i-1]<<1);
			//llr_es2[i-1] = nom - denom;

			llr_es2[i-1] = (llr_es2[i-1]>limit)? limit:llr_es2[i-1];
			llr_es2[i-1] = (llr_es2[i-1]<-limit)? -limit:llr_es2[i-1];

			pta = alfa2 + (j+1)*StateNum;

			dt0[0] = bs0[0] + pta[0];
			dt0[1] = bs0[1] + pta[1];
			dt0[2] = bs0[2] + pta[2];
			dt0[3] = bs0[3] + pta[3];
			dt0[4] = bs0[4] + pta[4];
			dt0[5] = bs0[5] + pta[5];
			dt0[6] = bs0[6] + pta[6];
			dt0[7] = bs0[7] + pta[7];

			dt1[0] = bs1[0] + pta[0];
			dt1[1] = bs1[1] + pta[1];
			dt1[2] = bs1[2] + pta[2];
			dt1[3] = bs1[3] + pta[3];
			dt1[4] = bs1[4] + pta[4];
			dt1[5] = bs1[5] + pta[5];
			dt1[6] = bs1[6] + pta[6];
			dt1[7] = bs1[7] + pta[7];

			nom = maxo_avx(maxo_avx(maxo_avx(dt1[0],dt1[1]),maxo_avx(dt1[2],dt1[3])),maxo_avx(maxo_avx(dt1[4],dt1[5]),maxo_avx(dt1[6],dt1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(dt0[0],dt0[1]),maxo_avx(dt0[2],dt0[3])),maxo_avx(maxo_avx(dt0[4],dt0[5]),maxo_avx(dt0[6],dt0[7])));
		
			llr_es2[j+1] = nom - denom - (gs2[j+1]<<1);
			//llr_es1[j+1] = nom - denom;
			
			llr_es2[j+1] = (llr_es2[j+1]>limit)? limit:llr_es2[j+1];
			llr_es2[j+1] = (llr_es2[j+1]<-limit)? -limit:llr_es2[j+1];
		}

		{
			short		*pta, *ptb;
			short		nom, denom;
			short		dt0[8], dt1[8];

			pta = alfa2;
			ptb = beta2;

			dt0[0] = pta[0]-gs2[0]-gp2[0]+ptb[0];
			dt0[1] = pta[1]-gs2[0]-gp2[0]+ptb[4];
			dt0[2] = pta[2]-gs2[0]+gp2[0]+ptb[5];
			dt0[3] = pta[3]-gs2[0]+gp2[0]+ptb[1];
			dt0[4] = pta[4]-gs2[0]+gp2[0]+ptb[2];
			dt0[5] = pta[5]-gs2[0]+gp2[0]+ptb[6];
			dt0[6] = pta[6]-gs2[0]-gp2[0]+ptb[7];
			dt0[7] = pta[7]-gs2[0]-gp2[0]+ptb[3];

			dt1[0] = pta[0]+gs2[0]+gp2[0]+ptb[4];
			dt1[1] = pta[1]+gs2[0]+gp2[0]+ptb[0];
			dt1[2] = pta[2]+gs2[0]-gp2[0]+ptb[1];
			dt1[3] = pta[3]+gs2[0]-gp2[0]+ptb[5];
			dt1[4] = pta[4]+gs2[0]-gp2[0]+ptb[6];
			dt1[5] = pta[5]+gs2[0]-gp2[0]+ptb[2];
			dt1[6] = pta[6]+gs2[0]+gp2[0]+ptb[3];
			dt1[7] = pta[7]+gs2[0]+gp2[0]+ptb[7];

			nom = maxo_avx(maxo_avx(maxo_avx(dt1[0],dt1[1]),maxo_avx(dt1[2],dt1[3])),maxo_avx(maxo_avx(dt1[4],dt1[5]),maxo_avx(dt1[6],dt1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(dt0[0],dt0[1]),maxo_avx(dt0[2],dt0[3])),maxo_avx(maxo_avx(dt0[4],dt0[5]),maxo_avx(dt0[6],dt0[7])));
		
			llr_es2[0] = nom - denom - (gs2[0]<<1);
			//llr_es1[0] = nom - denom;
			
			llr_es2[0] = (llr_es2[0]>limit)? limit:llr_es2[0];
			llr_es2[0] = (llr_es2[0]<-limit)? -limit:llr_es2[0];

			pta = alfa2 + (Len+Mem-1)*StateNum;
			ptb = beta2 + (Len+Mem-1)*StateNum;

			dt0[0] = pta[0]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[0];
			dt0[1] = pta[1]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[4];
			dt0[2] = pta[2]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[5];
			dt0[3] = pta[3]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[1];
			dt0[4] = pta[4]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[2];
			dt0[5] = pta[5]-gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[6];
			dt0[6] = pta[6]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[7];
			dt0[7] = pta[7]-gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[3];

			dt1[0] = pta[0]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[4];
			dt1[1] = pta[1]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[0];
			dt1[2] = pta[2]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[1];
			dt1[3] = pta[3]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[5];
			dt1[4] = pta[4]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[6];
			dt1[5] = pta[5]+gs2[Len+Mem-1]-gp2[Len+Mem-1]+ptb[2];
			dt1[6] = pta[6]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[3];
			dt1[7] = pta[7]+gs2[Len+Mem-1]+gp2[Len+Mem-1]+ptb[7];

			nom = maxo_avx(maxo_avx(maxo_avx(dt1[0],dt1[1]),maxo_avx(dt1[2],dt1[3])),maxo_avx(maxo_avx(dt1[4],dt1[5]),maxo_avx(dt1[6],dt1[7])));
			denom = maxo_avx(maxo_avx(maxo_avx(dt0[0],dt0[1]),maxo_avx(dt0[2],dt0[3])),maxo_avx(maxo_avx(dt0[4],dt0[5]),maxo_avx(dt0[6],dt0[7])));
		
			llr_es2[Len+Mem-1] = nom - denom - (gs2[Len+Mem-1]<<1);
			//llr_es1[Len+Mem-1] = nom - denom;
			
			llr_es2[Len+Mem-1] = (llr_es2[Len+Mem-1]>limit)? limit:llr_es2[Len+Mem-1];
			llr_es2[Len+Mem-1] = (llr_es2[Len+Mem-1]<-limit)? -limit:llr_es2[Len+Mem-1];
		}
		/* end of decoder2 */
		/*{
			FILE *fp;

			fp = fopen("turbo_avx2.txt","at");

			for(i=0; i<Len+Mem; i++)
			{
				fprintf(fp,"%8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d\n", i, llr_es2[i], g0[i], g1[i], 
					alfa2[i*StateNum+0], alfa2[i*StateNum+1],alfa2[i*StateNum+2], alfa2[i*StateNum+3],alfa2[i*StateNum+4], alfa2[i*StateNum+5],alfa2[i*StateNum+6], alfa2[i*StateNum+7],
					beta2[i*StateNum+0], beta2[i*StateNum+1],beta2[i*StateNum+2], beta2[i*StateNum+3],beta2[i*StateNum+4], beta2[i*StateNum+5],beta2[i*StateNum+6], beta2[i*StateNum+7]);
			}

			fclose(fp);
		}*/

		deinterleave_short(llr_es2, Len, interleave_table);
	}

	for(i=0; i<Len; i++)
	{
		store_llr[i] = (double) llr_es2[i]/(double) (1<<qbits);	
	}

	deinterleave_short(llr_es1, Len, interleave_table);

	/* compute the llr_ep1 */
	for(i=0; i<Len+Mem; i++)
	{
		short nom,denom;
		short nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;	

		nom1	=	maxo_avx(alfa1[i*StateNum+0]+gs1[i]+beta1[i*StateNum+4], alfa1[i*StateNum+1]+gs1[i]+beta1[i*StateNum+0]);
		nom2	=	maxo_avx(alfa1[i*StateNum+2]-gs1[i]+beta1[i*StateNum+5], alfa1[i*StateNum+3]-gs1[i]+beta1[i*StateNum+1]);
		nom3	=	maxo_avx(alfa1[i*StateNum+4]-gs1[i]+beta1[i*StateNum+2], alfa1[i*StateNum+5]-gs1[i]+beta1[i*StateNum+6]);
		nom4	=	maxo_avx(alfa1[i*StateNum+6]+gs1[i]+beta1[i*StateNum+3], alfa1[i*StateNum+7]+gs1[i]+beta1[i*StateNum+7]);

		denom1	=	maxo_avx(alfa1[i*StateNum+0]-gs1[i]+beta1[i*StateNum+0], alfa1[i*StateNum+1]-gs1[i]+beta1[i*StateNum+4]);
		denom2	=	maxo_avx(alfa1[i*StateNum+2]+gs1[i]+beta1[i*StateNum+1], alfa1[i*StateNum+3]+gs1[i]+beta1[i*StateNum+5]);
		denom3	=	maxo_avx(alfa1[i*StateNum+4]+gs1[i]+beta1[i*StateNum+6], alfa1[i*StateNum+5]+gs1[i]+beta1[i*StateNum+2]);
		denom4	=	maxo_avx(alfa1[i*StateNum+6]-gs1[i]+beta1[i*StateNum+7], alfa1[i*StateNum+7]-gs1[i]+beta1[i*StateNum+3]);

		nom = maxo_avx(maxo_avx(nom1,nom2),maxo_avx(nom3,nom4));
		denom = maxo_avx(maxo_avx(denom1,denom2),maxo_avx(denom3,denom4));

		llr_ep1[i] = nom - denom;

		if (llr_ep1[i] > limit)
		{
			llr_ep1[i] = limit;
		}
		if (llr_ep1[i] < -limit)
		{
			llr_ep1[i] = -limit;
		}
	}

	/* compute the llr_ep2 */
	for(i=0; i<Len+Mem; i++)
	{
		short nom,denom;
		short nom1,nom2,nom3,nom4,denom1,denom2,denom3,denom4;	

		nom1	=	maxo_avx(alfa2[i*StateNum+0]+gs2[i]+beta2[i*StateNum+4], alfa2[i*StateNum+1]+gs2[i]+beta2[i*StateNum+0]);
		nom2	=	maxo_avx(alfa2[i*StateNum+2]-gs2[i]+beta2[i*StateNum+5], alfa2[i*StateNum+3]-gs2[i]+beta2[i*StateNum+1]);
		nom3	=	maxo_avx(alfa2[i*StateNum+4]-gs2[i]+beta2[i*StateNum+2], alfa2[i*StateNum+5]-gs2[i]+beta2[i*StateNum+6]);
		nom4	=	maxo_avx(alfa2[i*StateNum+6]+gs2[i]+beta2[i*StateNum+3], alfa2[i*StateNum+7]+gs2[i]+beta2[i*StateNum+7]);

		denom1	=	maxo_avx(alfa2[i*StateNum+0]-gs2[i]+beta2[i*StateNum+0], alfa2[i*StateNum+1]-gs2[i]+beta2[i*StateNum+4]);
		denom2	=	maxo_avx(alfa2[i*StateNum+2]+gs2[i]+beta2[i*StateNum+1], alfa2[i*StateNum+3]+gs2[i]+beta2[i*StateNum+5]);
		denom3	=	maxo_avx(alfa2[i*StateNum+4]+gs2[i]+beta2[i*StateNum+6], alfa2[i*StateNum+5]+gs2[i]+beta2[i*StateNum+2]);
		denom4	=	maxo_avx(alfa2[i*StateNum+6]-gs2[i]+beta2[i*StateNum+7], alfa2[i*StateNum+7]-gs2[i]+beta2[i*StateNum+3]);

		nom = maxo_avx(maxo_avx(nom1,nom2),maxo_avx(nom3,nom4));
		denom = maxo_avx(maxo_avx(denom1,denom2),maxo_avx(denom3,denom4));

		llr_ep2[i] = nom - denom;

		if (llr_ep2[i] > limit)
		{
			llr_ep2[i] = limit;
		}

		if (llr_ep2[i] < -limit)
		{
			llr_ep2[i] = -limit;
		}
	}

	for(i=0; i<Len+Mem; i++)
	{
		llr_s1[i] = llr_s1[i] + llr_es1[i] + llr_es2[i];
		llr_p1[i] = llr_p1[i] + llr_ep1[i];
		llr_p2[i] = llr_p2[i] + llr_ep2[i];
	}	

	// Tailbiting
	llr_s_tail[0]	=	llr_s1[Len];
	llr_p1_tail[0]	=	llr_p1[Len];
	llr_p2_tail[0]	=	llr_s1[Len+1];

	llr_s_tail[1]	=	llr_p1[Len+1];
	llr_p1_tail[1]	=	llr_s1[Len+2];
	llr_p2_tail[1]	=	llr_p1[Len+2];

	llr_s_tail[2]	=	llr_s1[Len];
	llr_p1_tail[2]	=	llr_p2[Len];
	llr_p2_tail[2]	=	llr_s1[Len+1];

	llr_s_tail[3]	=	llr_p2[Len+1];
	llr_p1_tail[3]	=	llr_s1[Len+2];
	llr_p2_tail[3]	=	llr_p2[Len+2];

	for(i=0; i<4; i++)
	{
		llr_s1[Len+i] = llr_s_tail[i];
		llr_p1[Len+i] = llr_p1_tail[i];
		llr_p2[Len+i] = llr_p2_tail[i];
	}

	/* decision */
	for(i=0; i<Len; i++)
	{
		dec[i] = (llr_s1[i] > 0)? 1:0;
	}

	for(i=0; i<Len+4; i++)
	{
		llr_out_rm[3*i] = llr_s1[i];
		llr_out_rm[3*i+1] = llr_p1[i];
		llr_out_rm[3*i+2] = llr_p2[i];
	}

	for(i=0; i<3*(Len+4); i++)
	{
		llr_out_rm[i] = (llr_out_rm[i]>64)? 64:llr_out_rm[i];
		llr_out_rm[i] = (llr_out_rm[i]<-64)? -64:llr_out_rm[i];			
	}

	rm_interleave_short(CodeLen, llr_out_rm, llr_out_fix, ratematch_table);

	for(i=0; i<CodeLen; i++)
	{
		llr_out[i] = (double) llr_out_fix[i]/(double) (1<<qbits);
	}

	free(llr_in_rm);
	free(llr_out_rm);
	free(llr_s1);
	free(llr_s2);
	free(llr_p1);
	free(llr_p2);
	free(llr_es1);
	free(llr_es2);
	free(llr_ep1);
	free(llr_ep2);
	free(alfa1);
	free(beta1);
	free(gs1);
	free(gp1);
	free(alfa2);
	free(beta2);
	free(gs2);
	free(gp2);
	free(llr_in_fix);
	free(llr_out_fix);
	free(g0);
	free(g1);
	free(tmp_llr_es1);
}

int Convert_LimitLLR(double *llr_in, short *llr_out, int CodeLen, const __m256i limit0_256i, const __m256i limit1_256i, short limit, int qbits)
{
	int	i;

	double	scale = (double) (1<<qbits);

	__m256d		scale_256d = _mm256_set_pd(scale,scale,scale,scale);

	for (i = 0; i+15 < CodeLen; i += 16)
	{
		__m256d		llra_in_256d,llrb_in_256d,llrc_in_256d,llrd_in_256d;

		__m128i		llra_in_128i,llrb_in_128i,llrc_in_128i,llrd_in_128i;

		__m256i		llrab_in_256i, llrcd_in_256i, llr16_in_256i;

		llra_in_256d = _mm256_load_pd(&llr_in[i]);
		llrb_in_256d = _mm256_load_pd(&llr_in[i+4]);
		llrc_in_256d = _mm256_load_pd(&llr_in[i+8]);
		llrd_in_256d = _mm256_load_pd(&llr_in[i+12]);

		llra_in_256d = _mm256_mul_pd(llra_in_256d, scale_256d);
		llrb_in_256d = _mm256_mul_pd(llrb_in_256d, scale_256d);
		llrc_in_256d = _mm256_mul_pd(llrc_in_256d, scale_256d);
		llrd_in_256d = _mm256_mul_pd(llrd_in_256d, scale_256d);

		llra_in_128i = _mm256_cvttpd_epi32(llra_in_256d);
		llrb_in_128i = _mm256_cvttpd_epi32(llrb_in_256d);
		llrc_in_128i = _mm256_cvttpd_epi32(llrc_in_256d);
		llrd_in_128i = _mm256_cvttpd_epi32(llrd_in_256d);

		llrab_in_256i = _mm256_set_m128i(llrb_in_128i,llra_in_128i);
		llrcd_in_256i = _mm256_set_m128i(llrd_in_128i,llrc_in_128i);
			
		llr16_in_256i = _mm256_packs_epi32(llrab_in_256i, llrcd_in_256i);
		llr16_in_256i = _mm256_permute4x64_epi64(llr16_in_256i, 0xd8);

		llr16_in_256i = _mm256_min_epi16(llr16_in_256i, limit0_256i);
		llr16_in_256i = _mm256_max_epi16(llr16_in_256i, limit1_256i);

		_mm256_store_si256((__m256i *)&llr_out[i], llr16_in_256i);
	}
	for (i; i < CodeLen; i++)
	{
		llr_out[i] = (short) (llr_in[i]*(1<<qbits));
		llr_out[i] = (llr_out[i] > limit)? limit : llr_out[i];
		llr_out[i] = (llr_out[i] < -limit)? -limit : llr_out[i];
	}
}

int Convert_I16_D(short *llr_in, double *llr_out, int Len, int qbits)
{
	int			i;
	__m256i		llr_I16_256i, llr0_I32_256i, llr1_I32_256i;
	__m256d		llr00_D_256d, llr01_D_256d, llr10_D_256d,llr11_D_256d;
	__m128i		llr0_I16_128i, llr1_I16_128i, llr00_I32_128i, llr01_I32_128i, llr10_I32_128i, llr11_I32_128i;
	double		scale = 1.0/(double) (1<<qbits);
	__m256d		scale_256d = _mm256_set_pd(scale,scale,scale,scale);

	//for(i=0; i+15 < Len; i += 16)
	//{
	//	llr_I16_256i = _mm256_load_si256((__m256i*)&llr_in[i]);
	//	llr0_I16_128i = _mm256_extracti128_si256(llr_I16_256i, 0);
	//	llr1_I16_128i = _mm256_extracti128_si256(llr_I16_256i, 1);

	//	llr0_I32_256i = _mm256_cvtepi16_epi32(llr0_I16_128i);
	//	llr1_I32_256i = _mm256_cvtepi16_epi32(llr1_I16_128i);

	//	llr00_I32_128i = _mm256_extracti128_si256(llr0_I32_256i, 0);
	//	llr01_I32_128i = _mm256_extracti128_si256(llr0_I32_256i, 1);
	//	llr10_I32_128i = _mm256_extracti128_si256(llr1_I32_256i, 0);
	//	llr11_I32_128i = _mm256_extracti128_si256(llr1_I32_256i, 1);

	//	llr00_D_256d = _mm256_cvtepi32_pd(llr00_I32_128i);
	//	llr01_D_256d = _mm256_cvtepi32_pd(llr01_I32_128i);
	//	llr10_D_256d = _mm256_cvtepi32_pd(llr10_I32_128i);
	//	llr11_D_256d = _mm256_cvtepi32_pd(llr11_I32_128i);

	//	llr00_D_256d = _mm256_mul_pd(llr00_D_256d, scale_256d);
	//	llr01_D_256d = _mm256_mul_pd(llr01_D_256d, scale_256d);
	//	llr10_D_256d = _mm256_mul_pd(llr10_D_256d, scale_256d);
	//	llr11_D_256d = _mm256_mul_pd(llr11_D_256d, scale_256d);
	//}
	//for(i; i<Len; i++)
	//{
	//	llr_out[i] = (double) llr_in[i]*scale;
	//}
	for(i=0; i<Len; i++)
	{
		llr_out[i] = (double) llr_in[i]*scale;
	}
}

int	GetLLRsp(short *llr_s1, short *llr_s2, short * llr_p1, short* llr_p2, short *llr_in_rm, int Len)
{
	int			i; 
	short		llr_s_tail[4], llr_p1_tail[4], llr_p2_tail[4];

	for(i=0; i<Len+4; i++)
	{
		llr_s1[i] = llr_in_rm[3*i];
		llr_s2[i] = llr_in_rm[3*i];
		llr_p1[i] = llr_in_rm[3*i+1];
		llr_p2[i] = llr_in_rm[3*i+2];
	}

	//process tail bit
	for(i=0; i<4; i++)
	{
		llr_s_tail[i] = llr_s1[Len+i];
		llr_p1_tail[i] = llr_p1[Len+i];
		llr_p2_tail[i] = llr_p2[Len+i];
	}
	llr_s1[Len]		= llr_s_tail[0];
	llr_s1[Len+1]	= llr_p2_tail[0];
	llr_s1[Len+2]	= llr_p1_tail[1];

	llr_p1[Len]		= llr_p1_tail[0];
	llr_p1[Len+1]	= llr_s_tail[1];
	llr_p1[Len+2]	= llr_p2_tail[1];

	llr_s2[Len]		= llr_s_tail[2];
	llr_s2[Len+1]	= llr_p2_tail[2];
	llr_s2[Len+2]	= llr_p1_tail[3];

	llr_p2[Len]		= llr_p1_tail[2];
	llr_p2[Len+1]	= llr_s_tail[3];
	llr_p2[Len+2]	= llr_p2_tail[3];
}

int InitGamma(short *gs1, short *gp1, short *g0, short* g1, short *llr_s1, short *llr_p1, short *llr_es2, int Len)
{
	int i;

	__m256i gs1_256i, gp1_256i, g0_256i, g1_256i, tmp0_256i, tmp1_256i;

	const __m256i neg_flag = _mm256_set1_epi16(0xff00); 

	//for (i = 0; i < Len; i += 16)
	//{
	//	tmp0_256i = _mm256_load_si256((__m256i*) &llr_s1[i]);
	//	tmp1_256i = _mm256_load_si256((__m256i*) &llr_es2[i]);
	//	tmp0_256i = _mm256_add_epi16(tmp0_256i, tmp1_256i);
	//	gs1_256i = _mm256_srai_epi16(tmp0_256i, 1);
	//	_mm256_store_si256((__m256i*) &gs1[i], gs1_256i);

	//	tmp0_256i = _mm256_load_si256((__m256i*) &llr_p1[i]);
	//	gp1_256i = _mm256_srai_epi16(tmp0_256i, 1);
	//	_mm256_store_si256((__m256i*) &gp1[i], gp1_256i);

	//	g0_256i = _mm256_add_epi16(gs1_256i, gp1_256i);
	//	g0_256i = _mm256_sign_epi16(g0_256i, neg_flag);
	//	g1_256i = _mm256_sub_epi16(gp1_256i, gs1_256i);

	//	_mm256_store_si256((__m256i*) &g0[i], g0_256i);
	//	_mm256_store_si256((__m256i*) &g1[i], g1_256i);
	//}

	for (i=0; i<Len; i++)
	{
		gs1[i] = (llr_s1[i] >> 1) + (llr_es2[i] >> 1);
		gp1[i] = (llr_p1[i] >> 1);

		g0[i] = -gs1[i] - gp1[i];
		g1[i] = -gs1[i] + gp1[i];
	}
	for (i = Len; i<Len + Mem; i++)
	{
		gs1[i] = (llr_s1[i] >> 1);
		gp1[i] = (llr_p1[i] >> 1);

		g0[i] = -gs1[i] - gp1[i];
		g1[i] = -gs1[i] + gp1[i];
	}
}

int InitGamma2(short *gs1, short *gp1, short* alfa_beta_gamma, short *llr_s1, short *llr_p1, short *llr_es2, int Len)
{
	int i;

	__m256i g0_256i, g1_256i, g00_256i, g01_256i, g10_256i, g11_256i, ga_256i, gb_256i;

	const __m256i neg_flag = _mm256_set1_epi16(0xff00); 

	for (i=0; i<(Len+Mem+1)/2; i++)
	{
		short	g0i, g1i,g0j,g1j;
		short	sgs0, sgp0, sgs1, sgp1;
		
		int		j = Len + Mem - 1 - i;
		short	*pgamma0 = alfa_beta_gamma + 32*i + 16;
		short	*pgamma1 = alfa_beta_gamma + 32*j + 16;

	/*	sgs0 = (llr_s1[i] >> 1) + (llr_es2[i] >> 1);
		sgp0 = (llr_p1[i] >> 1);*/

		sgs0 = (llr_s1[i]+llr_es2[i]) >> 1;
		sgp0 = (llr_p1[i] >> 1);

		if (j<Len)
		{
			sgs1 = (llr_s1[j]+llr_es2[j]) >> 1;
			sgp1 = (llr_p1[j] >> 1);
		}
		else
		{
			sgs1 = (llr_s1[j] >> 1);
			sgp1 = (llr_p1[j] >> 1);
		}

		gs1[i] = sgs0;
		gp1[i] = sgp0;
		gs1[j] = sgs1;
		gp1[j] = sgp1;

		g0i =  -sgs0 - sgp0;
		g1i =  -sgs0 + sgp0;

		g0j =  -sgs1 - sgp1;
		g1j =  -sgs1 + sgp1;

		//g_256i	= _mm256_set_epi16(g0j, g1j,g1j, g0j, g0j, g1j,g1j, g0j, 
		//		g0i, g0i, g1i, g1i, g1i, g1i, g0i, g0i);

		g00_256i = _mm256_set1_epi16(g0i);
		g10_256i = _mm256_set1_epi16(g0j);
		g01_256i = _mm256_set1_epi16(g1i);
		g11_256i = _mm256_set1_epi16(g1j);

		ga_256i = _mm256_blend_epi16(g00_256i, g01_256i, 0x3c);
		gb_256i = _mm256_blend_epi16(g10_256i, g11_256i, 0x66);

		g0_256i = _mm256_permute2x128_si256(ga_256i, gb_256i, 0x30);

		_mm256_store_si256 ((__m256i *) pgamma0, g0_256i);

		ga_256i = _mm256_blend_epi16(g00_256i, g01_256i, 0x66);
		gb_256i = _mm256_blend_epi16(g10_256i, g11_256i, 0x3c);

		g1_256i = _mm256_permute2x128_si256(gb_256i, ga_256i, 0x30);

		_mm256_store_si256 ((__m256i *) pgamma1, g1_256i);
	}
}

__m256i max_fun(__m256i a, __m256i b, __m256i logmap_table0, __m256i logmap_table1)
{
	__m256i		maxout, logmapdelta;

	maxout = _mm256_max_epi16(a, b);

	//// log-map approximation
	//logmapdelta = _mm256_sub_epi16(a, b);
	//logmapdelta = _mm256_abs_epi16(logmapdelta);
	//logmapdelta = _mm256_srai_epi16(logmapdelta, 0x02);
	//logmapdelta = _mm256_sub_epi16(logmap_table0, logmapdelta);
	//logmapdelta = _mm256_max_epi16(logmapdelta, logmap_table1);
	//maxout = _mm256_add_epi16(maxout, logmapdelta);

	return maxout;
}

int	Limit_Int16(short *data, int Len, short hilimit, short lolimit)
{
	int	i;
	for(i=0; i<Len; i++)
	{
		data[i]  = (data[i]>hilimit)? hilimit:data[i];
		data[i]  = (data[i]<lolimit)? lolimit:data[i];
	}
}

void turbo_decode_rm_avx(short* dec, 
						double* llr_out, 
						double* llr_in, 
						double* store_llr, 
						int Len, 
						int CodeLen,
						int iterCNum, 
						int* interleave_table,
						int* ratematch_table)
{
	int			i, s, iter_idx;
	short		llr_s_tail[4], llr_p1_tail[4], llr_p2_tail[4];
	int			inf = 4095;
	int			limit = 255;
	int			qbits = 4;
	int			bitwidth = 15;
	int			StateNum = 8;

	short		*tmp_llr;

	const __m256i limit0_256i = _mm256_set_epi64x(0x00FF00FF00FF00FF, 0x00FF00FF00FF00FF, 0x00FF00FF00FF00FF, 0x00FF00FF00FF00FF);
	const __m256i limit1_256i = _mm256_set_epi64x(0xFF01FF01FF01FF01, 0xFF01FF01FF01FF01, 0xFF01FF01FF01FF01, 0xFF01FF01FF01FF01);
	const __m256i initalfabeta_256i = _mm256_set_epi64x(0xF001F001F001F001,0xF001F001F0010000,0xF001F001F001F001,0xF001F001F0010000);

	const __m256i shf_table0 = _mm256_set_epi16(0x0706, 0x0f0e, 0x0d0c, 0x0504, 0x0302, 0x0b0a, 0x0908, 0x0100,
											   0x0d0c, 0x0b0a, 0x0504, 0x0302, 0x0f0e, 0x0908, 0x0706, 0x0100);
	const __m256i shf_table1 = _mm256_set_epi16(0x0f0e, 0x0706, 0x0504, 0x0d0c, 0x0b0a, 0x0302, 0x0100, 0x0908,
											   0x0f0e, 0x0908, 0x0706, 0x0100, 0x0d0c, 0x0b0a, 0x0504, 0x0302);
	const __m256i mshf_table0 = _mm256_set_epi16(0x0d0c, 0x0b0a, 0x0504, 0x0302, 0x0f0e, 0x0908, 0x0706, 0x0100,
											    0x0706, 0x0f0e, 0x0d0c, 0x0504, 0x0302, 0x0b0a, 0x0908, 0x0100);

	const __m256i mshf_table1 = _mm256_set_epi16(0x0f0e, 0x0908, 0x0706, 0x0100, 0x0d0c, 0x0b0a, 0x0504, 0x0302,
											   0x0f0e, 0x0706, 0x0504, 0x0d0c, 0x0b0a, 0x0302, 0x0100, 0x0908);

	const __m256i maxshf_table0 = _mm256_set_epi16(0x0706, 0x0504, 0x0302, 0x0100, 0x0f0e, 0x0d0c, 0x0b0a, 0x0908,
												0x0706, 0x0504, 0x0302, 0x0100, 0x0f0e, 0x0d0c, 0x0b0a, 0x0908);
	const __m256i maxshf_table1 = _mm256_set_epi16(0x0d0c, 0x0f0e, 0x0908, 0x0b0a, 0x0504,0x0706,0x0100,0x0302, 
												0x0d0c, 0x0f0e, 0x0908, 0x0b0a, 0x0504,0x0706,0x0100,0x0302);
	const __m256i maxshf_table2 = _mm256_set_epi16(0x0b0a, 0x0908, 0x0f0e,  0x0d0c, 0x0302, 0x0100, 0x0706, 0x0504, 
												 0x0b0a,0x0908, 0x0f0e, 0x0d0c, 0x0302, 0x0100, 0x0706, 0x0504);

	const __m256i gp_table = _mm256_set_epi16(0x00FF, 0xFF00, 0xFF00, 0x00FF, 0x00FF, 0xFF00, 0xFF00, 0x00FF,
											   0x00FF, 0x00FF, 0xFF00, 0xFF00, 0xFF00, 0xFF00, 0x00FF, 0x00FF);

	const __m256i pshf_table0 = _mm256_set_epi16(0x0d0c, 0x0908, 0x0706, 0x0302, 0x0f0e, 0x0b0a, 0x0504, 0x0100,
											   0x0706, 0x0f0e, 0x0504, 0x0d0c, 0x0b0a, 0x0302, 0x0908, 0x0100);

	const __m256i pshf_table1 = _mm256_set_epi16(0x0f0e, 0x0b0a, 0x0504, 0x0100, 0x0d0c, 0x0908, 0x0706, 0x0302,
											   0x0f0e, 0x0706, 0x0d0c, 0x0504, 0x0302, 0x0b0a, 0x0100, 0x0908);

	const __m256i logmap_table0 = _mm256_set1_epi16(0x0b); 
	const __m256i logmap_table1 = _mm256_set1_epi16(0x00); 

	// 0) llr_in_fix -> CodeLen; 
	// 1) llr_out_fix -> CodeLen;
	// 2) llr_in_rm -> 3*(Len+4);
	// 3) llr_out_rm -> 3*(Len+4);
	// 4) llr_s1_Len -> (Len+4);
	// 5) llr_s2_Len -> (Len+4);
	// 6) llr_p1_Len -> (Len+4);
	// 7) llr_p2_Len -> (Len+4);
	// 8) llr_es1_Len -> (Len+Mem)
	// 9) llr_es2_Len -> (Len+Mem)
	// 10) llr_ep1_Len -> (Len+Mem)
	// 11) llr_ep2_Len -> (Len+Mem)
	// 12) gs1 -> (Len+Mem)
	// 13) gp1 -> (Len+Mem)
	// 14) gs2 -> (Len+Mem)
	// 15) gp2 -> (Len+Mem)
	// 16) alfa_beta1 -> (Len+Mem)*sizeof(__m256i)
	// 17) alfa_beta2 -> (Len+Mem)*sizeof(__m256i)
	// 18) g0 -> (Len+Mem)
	// 19) g1 -> (Len+Mem)

	int		llr_in_fix_Len = ((CodeLen>>4)+1)<<4;
	int		llr_out_fix_Len = ((CodeLen>>4)+1)<<4;
	int		llr_in_rm_Len = (((3*(Len+4))>>4)+1)<<4;
	int		llr_out_rm_Len = (((3 * (Len + 4)) >> 4) + 1) << 4;
	int		llr_s1_Len = (((Len + 4) >> 4) + 1) << 4;
	int		llr_s2_Len = (((Len + 4) >> 4) + 1) << 4;
	int		llr_p1_Len = (((Len + 4) >> 4) + 1) << 4;
	int		llr_p2_Len = (((Len + 4) >> 4) + 1) << 4;
	int		gs1_Len = (((Len+Mem)>>4)+1)<<4;
	int		gp1_Len = (((Len+Mem)>>4)+1)<<4;
	int		gs2_Len = (((Len+Mem)>>4)+1)<<4;
	int		gp2_Len = (((Len+Mem)>>4)+1)<<4;
	int		llr_es1_Len = (((Len + Mem) >> 4) + 1) << 4;
	int		llr_es2_Len = (((Len + Mem) >> 4) + 1) << 4;
	int		alfa_beta1_Len = ((((Len+Mem)*sizeof(__m256i))>>4)+1)<<4;
	int		alfa_beta2_Len = ((((Len+Mem)*sizeof(__m256i))>>4)+1)<<4;
	int		llr_ep1_Len = (((Len + Mem) >> 4) + 1) << 4;
	int		llr_ep2_Len = (((Len + Mem) >> 4) + 1) << 4;
	int		g0_Len = (((Len+Mem)>>4)+1)<<4;
	int		g1_Len = (((Len+Mem)>>4)+1)<<4;
	int		alfa_beta_gamma_Len = (((2*(Len+Mem)*sizeof(__m256i))>>4)+1)<<4;

	int		BufferLen = llr_in_fix_Len + llr_out_fix_Len + llr_in_rm_Len + llr_out_rm_Len + llr_s1_Len + llr_s2_Len + llr_p1_Len + llr_p2_Len + gs1_Len + gp1_Len + 
		gs2_Len + gp2_Len + llr_es1_Len + llr_es2_Len + alfa_beta1_Len + alfa_beta2_Len + llr_ep1_Len + llr_ep2_Len + g0_Len + g1_Len + alfa_beta_gamma_Len;

	short		*pBuffer = (short*)_mm_malloc(BufferLen*sizeof(short), sizeof(__m256i));

	short		*llr_in_fix		= pBuffer;
	short		*llr_out_fix	= llr_in_fix + llr_in_fix_Len;
	short		*llr_in_rm		= llr_out_fix + llr_out_fix_Len;
	short		*llr_out_rm		= llr_in_rm + llr_in_rm_Len;
	short		*llr_s1 = llr_out_rm + llr_out_rm_Len;
	short		*llr_s2 = llr_s1 + llr_s1_Len;
	short		*llr_p1 = llr_s2 + llr_s2_Len;
	short		*llr_p2 = llr_p1 + llr_p1_Len;
	short		*gs1 = llr_p2 + llr_p2_Len;
	short		*gp1 = gs1 + gs1_Len;
	short		*gs2 = gp1 + gp1_Len;
	short		*gp2 = gs2 + gs2_Len;
	short		*llr_es1 = gp2 + gp2_Len;
	short		*llr_es2 = llr_es1 + llr_es1_Len;
	short		*g0 = llr_es2 + llr_es2_Len;
	short		*g1 = g0 + g0_Len;
	short		*alfa_beta1 = g1 + g1_Len;
	short		*alfa_beta2 = alfa_beta1 + alfa_beta1_Len;
	short		*llr_ep1 = alfa_beta2 + alfa_beta2_Len;
	short		*llr_ep2 = llr_ep1 + llr_ep1_Len;
	short		*alfa_beta_gamma = llr_ep2 + llr_ep2_Len;

	__m256i		*alfa_beta_gamma_256i = (__m256i*) alfa_beta_gamma;

	//printf("%d, ", BufferLen*sizeof(short));

	//tmp_llr		=	(short*) _mm_malloc(CodeLen*sizeof(short), sizeof(__m256i));

	Convert_LimitLLR(llr_in, llr_in_fix, CodeLen, limit0_256i, limit1_256i, limit, qbits);

	//rate match....
	rate_dematching_short(Len+4, CodeLen, llr_in_fix, llr_in_rm, ratematch_table);

	Limit_Int16(llr_in_rm, 3*(Len+4), limit, -limit);

	GetLLRsp(llr_s1, llr_s2, llr_p1, llr_p2, llr_in_rm, Len);

	interleave_short(llr_s2, Len, interleave_table);

	Convert_LimitLLR(store_llr, llr_es2, Len+Mem, limit0_256i, limit1_256i, limit, qbits);

	for(iter_idx=0; iter_idx<iterCNum; iter_idx++)
	{
		// decoder 1
		// initialize for alfa and beta
		_mm256_store_si256(alfa_beta_gamma_256i, initalfabeta_256i);

		/* initialize for gs and gp */
		InitGamma2(gs1, gp1, alfa_beta_gamma, llr_s1, llr_p1, llr_es2, Len);

		/* compute the alfa and beta */
		for(i=1; i<Len+Mem; i++)
		{
			int			j = Len+Mem-1-i;

			__m256i		addout0, addout1, shffout0, shffout1, maxout, nm256ab, logmapdelta;
			__m256i		g_256i, g0_256i, g1_256i, tmpab0;
			__m256i		g00_256i, g01_256i, g10_256i, g11_256i, ga_256i, gb_256i;
			__m128i		g00_128i, g01_128i, g10_128i, g11_128i;

			__m256i		*pab = alfa_beta_gamma_256i + 2*(i - 1);	

			short		nom,denom;

			tmpab0 = _mm256_load_si256(pab);
			g_256i	=	_mm256_load_si256(pab+1);

			// add
			addout0 = _mm256_add_epi16(tmpab0, g_256i);
			addout1 = _mm256_sub_epi16(tmpab0, g_256i);

			//// shuffle
			shffout0 = _mm256_shuffle_epi8(addout0, shf_table0);
			shffout1 = _mm256_shuffle_epi8(addout1, shf_table1);

			//// max function
			maxout = max_fun(shffout0, shffout1, logmap_table0, logmap_table1);

			// normalization
			if (i % 8 == 0)
			{
				__m128i		nma_128i, nmb_128i;
				__m256i		nmab_256i;

				nma_128i = _mm256_extracti128_si256(maxout, 0);
				nmb_128i = _mm256_extracti128_si256(maxout, 1);

				nma_128i = _mm_broadcastw_epi16(nma_128i);
				nmb_128i = _mm_broadcastw_epi16(nmb_128i);

				nmab_256i = _mm256_set_m128i(nmb_128i, nma_128i);
			
				maxout = _mm256_sub_epi16(maxout, nmab_256i);
			}

			//update alfa_beta1;
			_mm256_store_si256(pab+2, maxout);

			if (i >= ((Len+Mem+1)/2))
			{
				__m256i		tmpab1 = _mm256_load_si256(alfa_beta_gamma_256i + 2*(Len + Mem - i));

				__m256i		tmpshff0, tmpshff1;
				__m256i		nom_256i, denom_256i;

				tmpab1 = _mm256_permute4x64_epi64(tmpab1, 0x4E);
				
				shffout0 = _mm256_add_epi16(shffout0, tmpab1);
				shffout1 = _mm256_add_epi16(shffout1, tmpab1);

				// llres1
				tmpshff0 = _mm256_shuffle_epi8(shffout0, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shffout0);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				tmpshff0 = _mm256_shuffle_epi8(shffout1, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shffout1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				nom = *((short*) &nom_256i);
				denom = *((short*) &denom_256i);

				llr_es1[i-1] = nom - denom - (gs1[i-1]<<1);

				nom = *(((short*) &nom_256i) + 8);
				denom = *(((short*) &denom_256i) + 8);

				llr_es1[j+1] = nom - denom - (gs1[j+1]<<1);

				if (iter_idx == iterCNum - 1)
				{
					__m256i		gp_256i, shff0, shff1;

					gp_256i	=	_mm256_sign_epi16(g_256i, gp_table);

					// printf("\n");
					// shuffle and add
					shff0 = _mm256_shuffle_epi8(tmpab1, pshf_table0);
					shff1 = _mm256_shuffle_epi8(tmpab1, pshf_table1);

					shff0 = _mm256_add_epi16(shff0, tmpab0);
					shff0 = _mm256_add_epi16(shff0, gp_256i);

					shff1 = _mm256_add_epi16(shff1, tmpab0);
					shff1 = _mm256_sub_epi16(shff1, gp_256i);

					// llrep1
					tmpshff0 = _mm256_shuffle_epi8(shff0, maxshf_table0);
					tmpshff1 = _mm256_max_epi16(tmpshff0, shff0);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
					tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
					denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

					tmpshff0 = _mm256_shuffle_epi8(shff1, maxshf_table0);
					tmpshff1 = _mm256_max_epi16(tmpshff0, shff1);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
					tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
					nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

					nom = *((short*) &nom_256i);
					denom = *((short*) &denom_256i);

					llr_ep1[i-1] = nom - denom - (gp1[i-1]<<1);

					nom = *(((short*) &nom_256i) + 8);
					denom = *(((short*) &denom_256i) + 8);

					llr_ep1[j+1] = nom - denom - (gp1[j+1]<<1);
				}
			}
		}

		{
			__m256i		g_256i	=	_mm256_load_si256(alfa_beta_gamma_256i+1);

			__m256i		tmpab0 = _mm256_load_si256(alfa_beta_gamma_256i);
			__m256i		tmpab1 = _mm256_load_si256(alfa_beta_gamma_256i + 2*(Len + Mem - 1));

			__m256i		shff0, shff1, tmpshff0, tmpshff1;
			__m256i		nom_256i, denom_256i;

			short		nom,denom;

			tmpab1 = _mm256_permute4x64_epi64(tmpab1, 0x4E);

			// shuffle and add
			shff0 = _mm256_shuffle_epi8(tmpab1, mshf_table0);
			shff1 = _mm256_shuffle_epi8(tmpab1, mshf_table1);

			shff0 = _mm256_add_epi16(shff0, tmpab0);
			shff0 = _mm256_add_epi16(shff0, g_256i);

			shff1 = _mm256_add_epi16(shff1, tmpab0);
			shff1 = _mm256_sub_epi16(shff1, g_256i);

			// llres1
			tmpshff0 = _mm256_shuffle_epi8(shff0, maxshf_table0);
			tmpshff1 = _mm256_max_epi16(tmpshff0, shff0);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
			tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
			denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

			tmpshff0 = _mm256_shuffle_epi8(shff1, maxshf_table0);
			tmpshff1 = _mm256_max_epi16(tmpshff0, shff1);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
			tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
			nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

			nom = *((short*) &nom_256i);
			denom = *((short*) &denom_256i);

			llr_es1[0] = nom - denom - (gs1[0]<<1);

			nom = *(((short*) &nom_256i) + 8);
			denom = *(((short*) &denom_256i) + 8);

			llr_es1[Len+Mem-1] = nom - denom - (gs1[Len+Mem-1]<<1);

			if (iter_idx == iterCNum - 1)
			{
				__m256i		gp_256i, shff0, shff1;

				gp_256i	=	_mm256_sign_epi16(g_256i, gp_table);

				// shuffle and add
				shff0 = _mm256_shuffle_epi8(tmpab1, pshf_table0);
				shff1 = _mm256_shuffle_epi8(tmpab1, pshf_table1);

				shff0 = _mm256_add_epi16(shff0, tmpab0);
				shff0 = _mm256_add_epi16(shff0, gp_256i);

				shff1 = _mm256_add_epi16(shff1, tmpab0);
				shff1 = _mm256_sub_epi16(shff1, gp_256i);

				// llrep1
				tmpshff0 = _mm256_shuffle_epi8(shff0, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shff0);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				tmpshff0 = _mm256_shuffle_epi8(shff1, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				nom = *((short*) &nom_256i);
				denom = *((short*) &denom_256i);

				llr_ep1[0] = nom - denom - (gp1[0]<<1);

				nom = *(((short*) &nom_256i) + 8);
				denom = *(((short*) &denom_256i) + 8);

				llr_ep1[Len+Mem-1] = nom - denom - (gp1[Len+Mem-1]<<1);
			}

			//printf("%d, %d\n", llr_es1[0], llr_es1[Len+Mem-1]);
		}

		Limit_Int16(llr_es1, Len+Mem, limit, -limit);

		if (iter_idx==iterCNum-1)
		{
			Limit_Int16(llr_ep1, Len+Mem, limit, -limit);
		}

		interleave_short(llr_es1, Len, interleave_table);

		// decoder2
		// initialize for alfa2 and beta2
		_mm256_store_si256(alfa_beta_gamma_256i, initalfabeta_256i);

		/* initialize for gs2 and gp2 */
		InitGamma2(gs2, gp2, alfa_beta_gamma, llr_s2, llr_p2, llr_es1, Len);

		/* compute the alfa2 and beta2 */
		for(i=1; i<Len+Mem; i++)
		{
			int			j = Len+Mem-1-i;

			__m256i		addout0, addout1, shffout0, shffout1, maxout, nm256ab, logmapdelta;
			__m256i		g_256i, g0_256i, g1_256i, tmpab0;
			__m256i		ba_256i;

			__m256i		*pab = alfa_beta_gamma_256i + 2*(i - 1);	

			short		nom,denom;

			g_256i	=	_mm256_load_si256(pab+1);

			tmpab0 = _mm256_load_si256(pab);

			// add
			addout0 = _mm256_add_epi16(tmpab0, g_256i);
			addout1 = _mm256_sub_epi16(tmpab0, g_256i);

			// shuffle
			shffout0 = _mm256_shuffle_epi8(addout0, shf_table0);
			shffout1 = _mm256_shuffle_epi8(addout1, shf_table1);

			// max function
			maxout = max_fun(shffout0, shffout1, logmap_table0, logmap_table1);

			// normalization
			if (i % 8 == 0)
			{
				__m128i		nma_128i, nmb_128i;
				__m256i		nmab_256i;

				nma_128i = _mm256_extracti128_si256(maxout, 0);
				nmb_128i = _mm256_extracti128_si256(maxout, 1);

				nma_128i = _mm_broadcastw_epi16(nma_128i);
				nmb_128i = _mm_broadcastw_epi16(nmb_128i);

				nmab_256i = _mm256_set_m128i(nmb_128i, nma_128i);
			
				maxout = _mm256_sub_epi16(maxout, nmab_256i);
			}

			//update alfa_beta2;
			_mm256_store_si256(pab+2, maxout);

			if (i >= ((Len+Mem+1)/2))
			{
				__m256i		tmpab1 = _mm256_load_si256(alfa_beta_gamma_256i + 2*(Len + Mem - i));

				__m256i		tmpshff0, tmpshff1;
				__m256i		nom_256i, denom_256i;

				tmpab1 = _mm256_permute4x64_epi64(tmpab1, 0x4E);
				
				shffout0 = _mm256_add_epi16(shffout0, tmpab1);
				shffout1 = _mm256_add_epi16(shffout1, tmpab1);

				// llres1
				tmpshff0 = _mm256_shuffle_epi8(shffout0, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shffout0);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				tmpshff0 = _mm256_shuffle_epi8(shffout1, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shffout1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				nom = *((short*) &nom_256i);
				denom = *((short*) &denom_256i);

				llr_es2[i-1] = nom - denom - (gs2[i-1]<<1);

				nom = *(((short*) &nom_256i) + 8);
				denom = *(((short*) &denom_256i) + 8);

				llr_es2[j+1] = nom - denom - (gs2[j+1]<<1);

				if (iter_idx == iterCNum - 1)
				{
					__m256i		gp_256i, shff0, shff1;

					gp_256i	=	_mm256_sign_epi16(g_256i, gp_table);

					// printf("\n");
					// shuffle and add
					shff0 = _mm256_shuffle_epi8(tmpab1, pshf_table0);
					shff1 = _mm256_shuffle_epi8(tmpab1, pshf_table1);

					shff0 = _mm256_add_epi16(shff0, tmpab0);
					shff0 = _mm256_add_epi16(shff0, gp_256i);

					shff1 = _mm256_add_epi16(shff1, tmpab0);
					shff1 = _mm256_sub_epi16(shff1, gp_256i);

					// llrep2
					tmpshff0 = _mm256_shuffle_epi8(shff0, maxshf_table0);
					tmpshff1 = _mm256_max_epi16(tmpshff0, shff0);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
					tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
					denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

					tmpshff0 = _mm256_shuffle_epi8(shff1, maxshf_table0);
					tmpshff1 = _mm256_max_epi16(tmpshff0, shff1);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
					tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
					tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
					nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

					nom = *((short*) &nom_256i);
					denom = *((short*) &denom_256i);

					llr_ep2[i-1] = nom - denom - (gp2[i-1]<<1);

					nom = *(((short*) &nom_256i) + 8);
					denom = *(((short*) &denom_256i) + 8);

					llr_ep2[j+1] = nom - denom - (gp2[j+1]<<1);
				}
				//printf("%d, %d, %d, %d\n", i-1, j+1, llr_es1[i-1], llr_es1[j+1]);
			}
		}

		{
			__m256i		g_256i	=	_mm256_load_si256(alfa_beta_gamma_256i+1);

			__m256i		tmpab0 = _mm256_load_si256(alfa_beta_gamma_256i);
			__m256i		tmpab1 = _mm256_load_si256(alfa_beta_gamma_256i + 2*(Len + Mem - 1));

			__m256i		shff0, shff1, tmpshff0, tmpshff1;
			__m256i		nom_256i, denom_256i;

			short		nom,denom;

			tmpab1 = _mm256_permute4x64_epi64(tmpab1, 0x4E);

			// shuffle and add
			shff0 = _mm256_shuffle_epi8(tmpab1, mshf_table0);
			shff1 = _mm256_shuffle_epi8(tmpab1, mshf_table1);

			shff0 = _mm256_add_epi16(shff0, tmpab0);
			shff0 = _mm256_add_epi16(shff0, g_256i);

			shff1 = _mm256_add_epi16(shff1, tmpab0);
			shff1 = _mm256_sub_epi16(shff1, g_256i);

			// llres1
			tmpshff0 = _mm256_shuffle_epi8(shff0, maxshf_table0);
			tmpshff1 = _mm256_max_epi16(tmpshff0, shff0);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
			tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
			denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

			tmpshff0 = _mm256_shuffle_epi8(shff1, maxshf_table0);
			tmpshff1 = _mm256_max_epi16(tmpshff0, shff1);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
			tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
			tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
			nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

			nom = *((short*) &nom_256i);
			denom = *((short*) &denom_256i);

			llr_es2[0] = nom - denom - (gs2[0]<<1);

			nom = *(((short*) &nom_256i) + 8);
			denom = *(((short*) &denom_256i) + 8);

			llr_es2[Len+Mem-1] = nom - denom - (gs2[Len+Mem-1]<<1);

			if (iter_idx == iterCNum - 1)
			{
				__m256i		gp_256i, shff0, shff1;

				gp_256i	=	_mm256_sign_epi16(g_256i, gp_table);

				// shuffle and add
				shff0 = _mm256_shuffle_epi8(tmpab1, pshf_table0);
				shff1 = _mm256_shuffle_epi8(tmpab1, pshf_table1);

				shff0 = _mm256_add_epi16(shff0, tmpab0);
				shff0 = _mm256_add_epi16(shff0, gp_256i);

				shff1 = _mm256_add_epi16(shff1, tmpab0);
				shff1 = _mm256_sub_epi16(shff1, gp_256i);

				// llrep1
				tmpshff0 = _mm256_shuffle_epi8(shff0, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shff0);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				denom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				tmpshff0 = _mm256_shuffle_epi8(shff1, maxshf_table0);
				tmpshff1 = _mm256_max_epi16(tmpshff0, shff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table1);
				tmpshff1 = _mm256_max_epi16(tmpshff0, tmpshff1);
				tmpshff0 = _mm256_shuffle_epi8(tmpshff1, maxshf_table2);
				nom_256i = _mm256_max_epi16(tmpshff0, tmpshff1);

				nom = *((short*) &nom_256i);
				denom = *((short*) &denom_256i);

				llr_ep2[0] = nom - denom - (gp2[0]<<1);

				nom = *(((short*) &nom_256i) + 8);
				denom = *(((short*) &denom_256i) + 8);

				llr_ep2[Len+Mem-1] = nom - denom - (gp2[Len+Mem-1]<<1);
			}

			//printf("%d, %d\n", llr_es1[0], llr_es1[Len+Mem-1]);
		}

		Limit_Int16(llr_es2, Len+Mem, limit, -limit);

		if (iter_idx==iterCNum-1)
		{
			Limit_Int16(llr_ep2, Len+Mem, limit, -limit);
		}
		// end of decoder2 

		deinterleave_short(llr_es2, Len, interleave_table);
	}

	Convert_I16_D(llr_es2, store_llr, Len, qbits);

	deinterleave_short(llr_es1, Len, interleave_table);

	for(i=0; i<Len+Mem; i++)
	{
		llr_s1[i] = llr_s1[i] + llr_es1[i] + llr_es2[i];
		llr_p1[i] = llr_p1[i] + llr_ep1[i];
		llr_p2[i] = llr_p2[i] + llr_ep2[i];
	}	

	// Tailbiting
	llr_s_tail[0]	=	llr_s1[Len];
	llr_p1_tail[0]	=	llr_p1[Len];
	llr_p2_tail[0]	=	llr_s1[Len+1];

	llr_s_tail[1]	=	llr_p1[Len+1];
	llr_p1_tail[1]	=	llr_s1[Len+2];
	llr_p2_tail[1]	=	llr_p1[Len+2];

	llr_s_tail[2]	=	llr_s1[Len];
	llr_p1_tail[2]	=	llr_p2[Len];
	llr_p2_tail[2]	=	llr_s1[Len+1];

	llr_s_tail[3]	=	llr_p2[Len+1];
	llr_p1_tail[3]	=	llr_s1[Len+2];
	llr_p2_tail[3]	=	llr_p2[Len+2];

	for(i=0; i<4; i++)
	{
		llr_s1[Len+i] = llr_s_tail[i];
		llr_p1[Len+i] = llr_p1_tail[i];
		llr_p2[Len+i] = llr_p2_tail[i];
	}

	/* decision */
	for(i=0; i<Len; i++)
	{
		dec[i] = (llr_s1[i] > 0)? 1:0;
	}

	for(i=0; i<Len+4; i++)
	{
		llr_out_rm[3*i] = llr_s1[i];
		llr_out_rm[3*i+1] = llr_p1[i];
		llr_out_rm[3*i+2] = llr_p2[i];
	}

	Limit_Int16(llr_out_rm, 3*(Len+4), limit, -limit);

	rm_interleave_short(CodeLen, llr_out_rm, llr_out_fix, ratematch_table);

	Convert_I16_D(llr_out_fix, llr_out,CodeLen, qbits);

	_mm_free(pBuffer);
}
