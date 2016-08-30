#include<malloc.h>
#include<time.h>
#define Mem 3
#define StateNum 8
#define Inf 128
#define Limit 8
double max(double a, double b)
{
	return a > b ? a : b;
}
void MaxLogSubDecoder(double*Ls,double*Lp,double*Le,int Len,double*Output) //Len长度不含4位尾比特，不记寄存器个数(Mem)
{
	//定义alpha，beta
	int OffsetTmp_p = -Inf;
	int OffsetTmp_n = -Inf;
	double *alpha = (double*)malloc(StateNum*(Len + Mem)*sizeof(double));
	double *beta = (double*)malloc(StateNum*(Len + Mem)*sizeof(double));
	//定义四种取值不同的分支度量
	double *r0 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r_0 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r1 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r_1 = (double*)malloc((Len + Mem)*sizeof(double));
	//初始化alpha，beta
	alpha[0] = 0;
	alpha[1] = -Inf;
	alpha[2] = -Inf;
	alpha[3] = -Inf;
	alpha[4] = -Inf;
	alpha[5] = -Inf;
	alpha[6] = -Inf;
	alpha[7] = -Inf;
	beta[(Len + Mem - 1)*StateNum + 0] = 0;
	beta[(Len + Mem - 1)*StateNum + 1] = -Inf;
	beta[(Len + Mem - 1)*StateNum + 2] = -Inf;
	beta[(Len + Mem - 1)*StateNum + 3] = -Inf;
	beta[(Len + Mem - 1)*StateNum + 4] = -Inf;
	beta[(Len + Mem - 1)*StateNum + 5] = -Inf;
	beta[(Len + Mem - 1)*StateNum + 6] = -Inf;
	beta[(Len + Mem - 1)*StateNum + 7] = -Inf;
	//初始化四个分支度量
	for (int i = 0; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Le[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Le[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	for (int i = Len; i < Len + Mem; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	//计算alpha,
	OffsetTmp_n = StateNum;
	OffsetTmp_p = 0;
	for (int i = 0; i < Len + Mem-1; ++i)
	{
		alpha[OffsetTmp_n] = max(alpha[OffsetTmp_p] + r_0[i], alpha[OffsetTmp_p + 1] + r0[i]);
		alpha[OffsetTmp_n + 1] = max(alpha[OffsetTmp_p + 2] + r1[i], alpha[OffsetTmp_p + 3] + r_1[i]);
		alpha[OffsetTmp_n + 2] = max(alpha[OffsetTmp_p + 4] + r_1[i], alpha[OffsetTmp_p + 5] + r1[i]);
		alpha[OffsetTmp_n + 3] = max(alpha[OffsetTmp_p + 6] + r0[i], alpha[OffsetTmp_p + 7] + r_0[i]);
		alpha[OffsetTmp_n + 4] = max(alpha[OffsetTmp_p + 0] + r0[i], alpha[OffsetTmp_p + 1] + r_0[i]);
		alpha[OffsetTmp_n + 5] = max(alpha[OffsetTmp_p + 3] + r1[i], alpha[OffsetTmp_p + 2] + r_1[i]);
		alpha[OffsetTmp_n + 6] = max(alpha[OffsetTmp_p + 4] + r1[i], alpha[OffsetTmp_p + 5] + r_1[i]);
		alpha[OffsetTmp_n + 7] = max(alpha[OffsetTmp_p + 7] + r0[i], alpha[OffsetTmp_p + 6] + r_0[i]);
		OffsetTmp_n += StateNum;
		OffsetTmp_p += StateNum;
	}
	//计算beta
	OffsetTmp_n = (Len + Mem - 1 )* StateNum;
	OffsetTmp_p = (Len + Mem - 2)*StateNum;
	for (int i = Len + Mem-1; i > 0; --i)
	{
		beta[OffsetTmp_p] = max(beta[OffsetTmp_n] + r_0[i], beta[OffsetTmp_n+4] + r0[i]);
		beta[OffsetTmp_p+1] = max(beta[OffsetTmp_n] + r0[i], beta[OffsetTmp_n+4] + r_0[i]);
		beta[OffsetTmp_p+2] = max(beta[OffsetTmp_n+1] + r1[i], beta[OffsetTmp_n+5] + r_1[i]);
		beta[OffsetTmp_p+3] = max(beta[OffsetTmp_n+5] + r1[i], beta[OffsetTmp_n+1] + r_1[i]);
		beta[OffsetTmp_p+4] = max(beta[OffsetTmp_n+6] + r1[i], beta[OffsetTmp_n+2] + r_1[i]);
		beta[OffsetTmp_p+5] = max(beta[OffsetTmp_n+2] + r1[i], beta[OffsetTmp_n+6] + r_1[i]);
		beta[OffsetTmp_p+6] = max(beta[OffsetTmp_n+3] + r0[i], beta[OffsetTmp_n+7] + r_0[i]);
		beta[OffsetTmp_p+7] = max(beta[OffsetTmp_n+3] + r_0[i], beta[OffsetTmp_n+7] + r0[i]);
		OffsetTmp_n -= StateNum;
		OffsetTmp_p -= StateNum;
	}
	//计算外信息
	/*接收情况：原状态     输出信号    校验位     末状态    对应r
	              0           0          0          0        r_0
				  1           0          0          4        r_0
				  2           0          1          5        r_1          
				  3           0          1          1        r_1
				  4           0          1          2        r_1
				  5           0          1          6        r_1
				  6           0          0          7        r_0
				  7           0          0          3        r_0
				  0           1          1          4        r0
				  1           1          1          0        r0
				  2           1          0          1        r1          
				  3           1          0          5        r1
				  4           1          0          6        r1 
				  5           1          0          2        r1
				  6           1          1          3        r0
				  7           1          1          7        r0
    */
	for (int i = 0; i < Len + Mem; ++i)
	{
		
		Output[i] =
			max(
			max(                                     //之前为i+1错误
			max(alpha[i*StateNum + 0] + r0[i] + beta[i*StateNum + 4], alpha[i*StateNum + 1] + r0[i] + beta[i*StateNum + 0]),
			max(alpha[i*StateNum + 2] + r1[i] + beta[i*StateNum + 1], alpha[i*StateNum + 3] + r1[i] + beta[i*StateNum + 5])
			),
			max(
			max(alpha[i*StateNum + 4] + r1[i] + beta[i*StateNum + 6], alpha[i*StateNum + 5] + r1[i] + beta[i*StateNum + 2]),
			max(alpha[i*StateNum + 6] + r0[i] + beta[i*StateNum + 3], alpha[i*StateNum + 7] + r0[i] + beta[i*StateNum + 7])
			)
			)
			-
			max(
			max(
			max(alpha[i*StateNum + 0] + r_0[i] + beta[i*StateNum + 0], alpha[i*StateNum + 1] + r_0[i] + beta[i*StateNum + 4]),
			max(alpha[i*StateNum + 2] + r_1[i] + beta[i*StateNum + 5], alpha[i*StateNum + 3] + r_1[i] + beta[i*StateNum + 1])
			),
			max(
			max(alpha[i*StateNum + 4] + r_1[i] + beta[i*StateNum + 2], alpha[i*StateNum + 5] + r_1[i] + beta[i*StateNum + 6]),
			max(alpha[i*StateNum + 6] + r_0[i] + beta[i*StateNum + 7], alpha[i*StateNum + 7] + r_0[i] + beta[i*StateNum + 3])
			)
			)
			- Ls[i] - Le[i];
		if (Output[i] > Limit) Output[i]= Limit;
		else if (Output[i] < -Limit) Output[i] = -Limit;
	}
	free(alpha);
	free(beta);
	free(r0);
	free(r1);
	free(r_1);
	free(r_0);
}
void Decoder(short* dec, double*LLR_out, double*LLR_in, int Len, int CodeLen, int iterCNum, int* interleave_table, int* ratematch_table)//省去了StoreLLR
{
	double*LLR_ein = (double*)malloc((Len + 4)*sizeof(double));           //为第一个译码器输入的外信息
	double*LLR_eout = (double*)malloc((Len + Mem)*sizeof(double));        //为第二个译码器输入的外信息       
	double*LLR_in_rm = (double*)malloc(3 * (Len + 4)*sizeof(double));     //速率匹配后的信息
	double*LLR_s1 = (double*)malloc((Len + 4)*sizeof(double));            //第一个译码器的输入信息
	double*LLR_s2 = (double*)malloc((Len + 4)*sizeof(double));            //第二个译码器的输入信息
	double*LLR_p1 = (double*)malloc((Len + 4)*sizeof(double));            //第一个译码器的校验信息
	double*LLR_p2 = (double*)malloc((Len + 4)*sizeof(double));            //第二个译码器的校验信息

	rate_dematching(Len + 4, CodeLen, LLR_in, LLR_in_rm, ratematch_table);
	//Ratedematching(LLR_in, ratematch_table, LLR_in_rm, Len + 4);
	for (int i = 0; i < Len + Mem; ++i)
	{
		LLR_ein[i] = 0;
	}
	/*********************************以下处理尾比特******************************/
	for (int i = 0; i<Len + 4; ++i)                        //已经加上尾比特
	{
		LLR_s1[i] = LLR_in_rm[3 * i];               //llr_s1给SubDecoder1的系统位
		LLR_s2[i] = LLR_s1[i];               //llr_s2给SubDecoder2的系统位
		LLR_p1[i] = LLR_in_rm[3 * i + 1];             //llr_p1给SubDecoder1的校验位
		LLR_p2[i] = LLR_in_rm[3 * i + 2];             //llr_p2给SubDecoder2的校验位
	}
	//接收的尾比特与输入的尾比特要重排，有一定的映射关系

	double LLR_copys1 = LLR_s1[Len + 1];           //临时存储
	//LLR_s1[Len] = LLR_s1[Len];//llr_s_tail[0];
	LLR_s1[Len + 1] = LLR_p2[Len];//llr_p2_tail[0];   无影响
	LLR_s1[Len + 2] = LLR_p1[Len + 1];//llr_p1_tail[1]; 无影响
	
	double LLR_copyp1 = LLR_p1[Len + 2];           //临时存储
	//LLR_p1[Len] = LLR_p1[Len];//llr_p1_tail[0];              
	LLR_p1[Len + 1] = LLR_copys1;//LLR_s1[Len + 1];//llr_s_tail[1]    被覆盖 ;
	LLR_p1[Len + 2] = LLR_p2[Len + 1];//llr_p2_tail[1];       

	LLR_s2[Len] = LLR_s2[Len+2];//llr_s_tail[2];
	LLR_s2[Len + 1] = LLR_p2[Len + 2];//llr_p2_tail[2];        无影响
	LLR_s2[Len + 2] = LLR_p1[Len + 3];//llr_p1_tail[3];         无影响

	LLR_p2[Len] = LLR_copyp1;//LLR_p1[Len + 2];//llr_p1_tail[2];    被覆盖         
	LLR_p2[Len + 1] = LLR_s1[Len + 3]; //llr_s_tail[3];          无影响
	LLR_p2[Len + 2] = LLR_p2[Len + 3]; //llr_p2_tail[3];         无影响
	
	interleave_double(LLR_s2, Len, interleave_table);      //将传送给SubDecoder2的系统位交织//之后llr_s1,llr_s2就一直不变了
	/*********************************一下开始迭代计算，输出三条流的应有结果******************************/
	for (int i = 0; i < iterCNum; ++i)
	{
		//SubDecoder1
		MaxLogSubDecoder(LLR_s1, LLR_p1,LLR_ein,Len,LLR_eout);
		if (i == iterCNum - 1)
			for (int j = 0; j < Len + Mem; j++)
				LLR_p1[j] += LLR_eout[j];             //交织前的信息
		interleave_double(LLR_eout, Len, interleave_table);
		//SubDecoder2
		MaxLogSubDecoder(LLR_s2, LLR_p2, LLR_eout, Len, LLR_ein);
		if (i == iterCNum - 1)
			for (int j = 0; j < Len + Mem; j++)
			{
				LLR_p2[j] += LLR_ein[j];  //交织前的信息
				LLR_ein[j]+= LLR_s2[j] + LLR_eout[j];
			}
		else                                           
			deinterleave_double(LLR_ein, Len, interleave_table);
	}
	deinterleave_double(LLR_ein, Len, interleave_table);//最终的数据结果，此数据用于硬判决
	/*************************以下开始处理尾比特***************************/
	
	double LLR_copye = LLR_ein[Len + 2];//临时存储
	LLR_copys1 = LLR_ein[Len + 1];
	LLR_ein[Len] = LLR_ein[Len];
	LLR_ein[Len + 1] = LLR_p1[Len + 1];
	LLR_ein[Len + 2] = LLR_ein[Len];
	LLR_ein[Len + 3] = LLR_p2[Len + 1];
	
	LLR_copyp1 = LLR_p1[Len + 2];//临时存储
	//LLR_p1[Len] = LLR_p1[Len];
	LLR_p1[Len + 1] = LLR_copye;//LLR_ein[Len + 2];
	LLR_p1[Len + 2] = LLR_p2[Len];
	LLR_p1[Len + 3] = LLR_copye;//LLR_ein[Len + 2];
	
	LLR_p2[Len] = LLR_copys1;//LLR_ein[Len + 1];
	LLR_p2[Len + 1] = LLR_p2[Len];//LLR_p1[Len + 2];
	LLR_p2[Len + 3] = LLR_p2[Len + 2]; //换位
	LLR_p2[Len + 2] = LLR_copyp1;//LLR_ein[Len + 1];
	
	/********************************* 判决**********************************/
	for (int i = 0; i<Len; i++)
	{ 
		dec[i] = (LLR_ein[i] > 0) ? 1 : 0;//仅仅根据llr_s[i]进行判决，其余都不要
	}

	for (int i = 0; i<Len + 4; i++)
	{
		LLR_out[3 * i] = LLR_ein[i];
		LLR_out[3 * i + 1] = LLR_p1[i];//为什么连同校验位一同输出？
		LLR_out[3 * i + 2] = LLR_p2[i];
	}
	for (int i = 0; i<3 * (Len + 4); i++)
	{
		LLR_out[i] = (LLR_out[i]>64) ? 64 : LLR_out[i];
		LLR_out[i] = (LLR_out[i]<-64) ? -64 : LLR_out[i];
	}
	//rm_interleave_doub(Len + 4, CodeLen, LLR_out, LLR_out, ratematch_table);
	free(LLR_ein);
	free(LLR_eout);
	free(LLR_in_rm);
	free(LLR_s1);
	free(LLR_s2);
	free(LLR_p1);
	free(LLR_p2);
}
void RandomInterleavetableGenerator(int Len, int* table)
{
	int A = 0, B = 0,temp = 0;
	for (int i = 0; i < Len; ++i)
		table[i] = i;
	for (int i = 0; i < 40 * Len; ++i)
	{
		srand((int)time(NULL));
		A = (int)(rand()&Len);
		temp = table[A];
		srand((int)time(NULL));
		B = (int)(rand()&Len);
		table[A] = table[B];
		table[B] = temp;
	}
}
void SubBlockInterleavetableGenerator(int Len,int*table,int type)//真正输入的长度为K+4=Len，K为输入编码器的信息比特数
{                                                                //table长度为K+4=Len
	int ICP[32] = { 0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30,
		1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31 };
	int C = 32;  
	int count = 0;
	int R = Len % 32 == 0 ? Len / 32 : (Len / 32) + 1;
	int*matrix = (int*)malloc(C*R*sizeof(int));
	/********给矩阵初始值********/
	for (int i = 0; i < C*R - Len; ++i)
		matrix[i] = -1;
	for (int i = C*R - Len; i < C*R; ++i)
	{
		matrix[i] = count;
		++count;
	}
	if (type == 0)//type = 1对应p2
	{
		count = 0;
		for (int i = 0; i < R*C; ++i)
		{
			if (matrix[(ICP[i / R] + C*(i%R)) % (R*C)] != -1)
			{
				table[count] = matrix[(ICP[i / R] + C*(i%R)) % (R*C)];
				count++;
			}
		}
	}
	else if (type == 1)
	{
		count = 0;
		for (int i = 0; i < R*C; ++i)
		{
			if (matrix[(ICP[i / R] + C*(i%R)) % (R*C)+1] != -1)
			{
				table[count] = matrix[(ICP[i / R] + C*(i%R)+1) % (R*C)];
				count++;
			}
		}
	}
	free(matrix);
}
void BitCollectiontableGenerator(int*table_s, int*table_p1, int*table_p2,int*OutTable,int Len)//Len为每个table的Len=K+4，总体输出长度为3*（K+4）
{
	for (int i = 0; i < Len; ++i)
	{
		OutTable[i] = table_s[i];
		OutTable[Len + 2 * i] = table_p1[i];
		OutTable[Len + 2 * i + 1] = table_p2[i];
	}
}
void RatematchtableGenerator(int Len, int*table)//table长为3*len+12；Len长度为K+4
{
	int*tables = (int*)malloc(Len*sizeof(int));
	int*tablep1 = (int*)malloc(Len*sizeof(int));
	int*tablep2 = (int*)malloc(Len*sizeof(int));
	SubBlockInterleavetableGenerator(Len, tables, 0);
	SubBlockInterleavetableGenerator(Len, tablep1, 0);
	SubBlockInterleavetableGenerator(Len, tablep2, 0);
	BitCollectiontableGenerator(tables, tablep1, tablep2, table, Len);
	free(tables);
	free(tablep1);
	free(tablep2);
}
void Ratedematching(int*stream,int*table,int*rm_stream,int Len)//Len为s，p1，p2长度，为Len=K+4;out为s，p1，p2交叉传输
{
	for (int i = 0; i < Len; ++i)
	{
		rm_stream[3*table[i]] = stream[i];
		rm_stream[3 * table[Len + 2 * i] + 1] = stream[Len + 2 * i];
		rm_stream[3 * table[Len + 2 * i] + 2] = stream[Len + 2 * i+1];
	}
}
void Ratematching(int*s, int*p1, int*p2,int*table,int*out,int Len)
{
	for (int i = 0; i < Len; ++i)
	{
		out[i] = s[table[i]];
		out[Len + 2 * i] = p1[table[i]];
		out[Len + 2 * i + 1] = p2[table[i]];
	}
}