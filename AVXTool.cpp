#include<stdio.h>
#include<string.h>
#include "ippcore.h"
#include<stdlib.h>
#include<iostream>
#include <stdlib.h>
#include<time.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vml.h> 
#include <windows.h>
#include <xmmintrin.h>
#include<stack>
#include <immintrin.h> //AVX
#define INF 127
#define Mem 3
#define QuantizeType float
//New Function
__declspec(naked) unsigned __int64 GetCpuCycle(void) //一个裸函数编译器不添加任何汇编代码
{
	_asm
	{
		rdtsc  //将时间标签寄存器读入EDX:EAX寄存器
		ret
	}
}
float Mmax(float A, float B)
{
	return A > B ? A : B;
}
bool BinaryStreamGenerator(int *Result, int Len)
{
	const __int32 SEED = GetCpuCycle();
	memset(Result, 0, Len);
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MT19937, SEED);
	viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, Len,Result,0.5);
	vslDeleteStream(&stream);
	return 1;
}
bool SubBlockInterleaveTableGenerator(int*table,int Len,int type)//真正输入的长度为K+4=Len，K为输入编码器的信息比特数
{                                                                  //table长度为K+4=Len
	int ICP[32] = { 0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30,
		1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31 };
	int C = 32;
	int count = 0;
	int R = Len % 32 == 0 ? Len / 32 : (Len / 32) + 1;
	int*matrix = (int*)malloc(C*R*sizeof(int));
	if(matrix == NULL) return 0;
	/********给矩阵初始值********/
	for (int i = 0; i < C*R - Len; ++i)
		matrix[i] = -1;
	for (int i = C*R - Len; i < C*R; ++i)
	{
		matrix[i] = count;
		++count;
	}
	if (type == 0)//d0 d1
	{
		count = 0;
		for (int i = 0; i < R*C; ++i)
		{
			if (matrix[(ICP[i / R] + C*(i%R)) % (R*C)] != -1)
			{
				int tmp = (ICP[i / R] + C*(i%R)) % (R*C);
				table[count] = matrix[(ICP[i / R] + C*(i%R)) % (R*C)];
				count++;
			}
		}
	}
	else if (type == 1) //d2
	{
		count = 0;
		for (int i = 0; i < R*C; ++i)
		{
			//if (matrix[(ICP[i / R] + C*(i%R)) % (R*C) + 1] != -1)
			if (matrix[(ICP[i / R] + C*(i%R) + 1) % (R*C)] != -1)
			{
				int tmp = (ICP[i / R] + C*(i%R) + 1) % (R*C);
				//tmp = matrix[(ICP[i / R] + C*(i%R) + 1) % (R*C)];
				table[count] = matrix[(ICP[i / R] + C*(i%R) + 1) % (R*C)];
				count++;
			}
		}
	}
	free(matrix);
	return 1;
}
bool TurboCodeInternalInterleaveTableGenerator(int K)
{
	return 1;
}
bool RandomTableGenerator(int *Table, int Len)
{
	if (Len == 0 || Len == 1) return 0;
	VSLStreamStatePtr stream;
	int*Reflect = (int *)malloc(Len*sizeof(int));
	for (int i = 0; i < Len; ++i) //生成初始映射
		Table[i] = i;
	for (int i = 0; i < 10; ++i)  //如此处理10次
	{
		__int32 SEED = GetCpuCycle();
		vslNewStream(&stream, VSL_BRNG_MT19937, SEED);
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, Len, Reflect, 0, Len - 1);
		vslDeleteStream(&stream);
		for (int i = 0; i < Len; ++i)
		{
			int Swap = Reflect[i];  //与i交换的元素index
			int Tmp = Table[i];
			Table[i] = Table[Swap];
			Table[Swap] = Tmp;
		}
	}
	free(Reflect);
	return 1;
}
bool ReadTable(int*Table, int Len, std::string Type) //Type==1时为随机表，Type==0为QPP
{
	FILE*Tmp;
	std::stack<int> S;
	std::string FileName;
	//根据类别取名
	if (Type == "Random")
		FileName = "RandomTable\\RandomTable_Ele";
	else if (Type == "QPP")
		FileName = "QPPTable\\QPP_Ele";
	else if (Type == "SubBlock_M_P1")
		FileName = "SubBlockTable\\SubBlock_M_P1\\Ele";
	else if (Type == "SubBlock_P2")
		FileName = "SubBlockTable\\SubBlock_P2\\Ele";
	else
	{
		printf("Read File Error,No Such Type!\nType Name As:\nRandom\tQPP\tSubBlock_M_P1\tSubBlock_P2\n");
		return 0;
	}
	int Num = Len;
	while (Num != 0)
	{
		S.push(Num % 10);
		Num /= 10;
	}
	while (!S.empty())
	{
		char tmp = S.top() + '0';
		S.pop();
		FileName = FileName + tmp;
	}
	FileName = FileName + ".dat";
	//取名完毕，接下来打开
	int State = fopen_s(&Tmp, FileName.c_str(), "r");
	if (State != 0)//如果没有此文件，则创建一个
	{
		State = fopen_s(&Tmp, FileName.c_str(), "w");
		if (State != 0) printf("Read File Error ,Can not Write!\n");
		else if (Type == "Random") RandomTableGenerator(Table, Len);
		else if (Type == "SubBlock_M_P1") SubBlockInterleaveTableGenerator(Table, Len, 0);
		else if (Type == "SubBlock_P2")SubBlockInterleaveTableGenerator(Table, Len, 1);
		//else if (Type == "QPP")
		for (int i = 0; i < Len; ++i)
			fprintf(Tmp, "%d ", Table[i]);
		fclose(Tmp);
		State = fopen_s(&Tmp, FileName.c_str(), "r"); //再次打开
	}
	for (int i = 0; i < Len; ++i)
		fscanf_s(Tmp, "%d ", &Table[i]);
	fclose(Tmp);
	return 1;
}
bool TurboRscEncoder(int *InputStream, int Len,int *OutputEncodedMessage,int*UnCodedMessage) //Type==0为计算Parity1，此时也计算Message，Type=0反之
{
	//int ForwardPath[3] = { 1, 0, 1 };//1,1,0,1
	//int BakwardPath[3] = { 0, 1, 1 };//1,0,1,1
	bool PresentState[3] = { 0 };
	for (int i = 0; i < Len; ++i)
	{
		if (UnCodedMessage != InputStream) UnCodedMessage[i] = InputStream[i];
		bool InputBit = PresentState[1] ^ PresentState[2]^InputStream[i];//按照BakwardPath
		bool OutputBit = InputBit^PresentState[0] ^ PresentState[2];
		OutputEncodedMessage[i] = OutputBit;
		PresentState[2] = PresentState[1];
		PresentState[1] = PresentState[0];
		PresentState[0] = InputBit;
	}
	//尾比特处理，感觉不太对啊。
	
	for (int i = Len; i < Len + Mem; ++i)
	{
		bool InputBit = PresentState[1] ^ PresentState[2];//此处填充1则与原先契合
		UnCodedMessage[i] = InputBit;//InputBit变为这个
		InputBit = InputBit ^ PresentState[1] ^ PresentState[2];
		bool OutputBit = InputBit^PresentState[0] ^ PresentState[2];
		OutputEncodedMessage[i] = OutputBit;
		PresentState[2] = PresentState[1];
		PresentState[1] = PresentState[0];
		PresentState[0] = InputBit;
	}
	
	//备份一个输入为0的版本
	/*
	for (int i = Len; i < Len + Mem; ++i)
	{
		int InputBit = 0;//此处填充1则与原先契合
		UnCodedMessage[i] = InputBit;
		int OutputBit = InputBit^PresentState[0] ^ PresentState[2];
		OutputEncodedMessage[i] = OutputBit;
		PresentState[2] = PresentState[1];
		PresentState[1] = PresentState[0];
		PresentState[0] = InputBit;
	}
	*/
	return 1;
}
bool TurboEncoder(int*InputStream, int Len, int *InternalInterleaverTable, int *ParityStream, int*InterleavedParityStream,int*OutputMessageStream)
{
	TurboRscEncoder(InputStream, Len, ParityStream,OutputMessageStream); //Parity1，Type为0时计算输出x
	int*InterleavedInputMessage = (int*)malloc(Len*sizeof(int)); 
	int*InterleavedOutputMessage = (int*)malloc((Len+Mem)*sizeof(int));
	if (InterleavedOutputMessage == NULL || InterleavedInputMessage == NULL) return 0;
	if (InterleavedInputMessage == NULL)return 0;
	for (int i = 0; i < Len; ++i)//交织映射
		InterleavedInputMessage[i] = InputStream[InternalInterleaverTable[i]];
	TurboRscEncoder(InterleavedInputMessage, Len, InterleavedParityStream, InterleavedOutputMessage);//Parity2
	/*
	printf("Rsc:Message：");
	for (int i = 0; i < Len + Mem; ++i)
		printf("%d ", OutputMessageStream[i]);
	printf("\n");
	printf("R:InMessage：");
	for (int i = 0; i < Len + Mem; ++i)
		printf("%d ", InterleavedOutputMessage[i]);
	printf("\n");
	printf("Rsc:Parity1：");
	for (int i = 0; i < Len + Mem; ++i)
		printf("%d ", ParityStream[i]);
	printf("\n");
	printf("Rsc:Parity2：");
	for (int i = 0; i < Len + Mem; ++i)
		printf("%d ", InterleavedParityStream[i]);
	printf("\n");
	*/
	//处理尾比特
	int InterleavedParityStreamTail[4];
	int OutputMessageStreamTail[4];
	int ParityStreamTail[4];
	OutputMessageStreamTail[0] = OutputMessageStream[Len];
	OutputMessageStreamTail[1] = ParityStream[Len + 1];
	OutputMessageStreamTail[2] = InterleavedOutputMessage[Len];
	OutputMessageStreamTail[3] = InterleavedParityStream[Len + 1];
	ParityStreamTail[0] = ParityStream[Len];
	ParityStreamTail[1] = OutputMessageStream[Len+2];
	ParityStreamTail[2] = InterleavedParityStream[Len];
	ParityStreamTail[3] = InterleavedOutputMessage[Len+2];
	InterleavedParityStreamTail[0] = OutputMessageStream[Len + 1];
	InterleavedParityStreamTail[1] = ParityStream[Len+2];
	InterleavedParityStreamTail[2] = InterleavedOutputMessage[Len + 1];
	InterleavedParityStreamTail[3] = InterleavedParityStream[Len + 2];
	for (int i = 0; i < 4; ++i)
	{
		OutputMessageStream[Len + i] = OutputMessageStreamTail[i];
		ParityStream[Len + i] = ParityStreamTail[i];
		InterleavedParityStream[Len + i] = InterleavedParityStreamTail[i];
	}
	/*
	printf("收前Message：");
	for (int i = 0; i < Len + 4; ++i)
		printf("%d ", OutputMessageStream[i]);
	printf("\n");
	printf("收前Parity1：");
	for (int i = 0; i < Len + 4; ++i)
		printf("%d ", ParityStream[i]);
	printf("\n");
	printf("收前Parity2：");
	for (int i = 0; i < Len + 4; ++i)
		printf("%d ", InterleavedParityStream[i]);
	printf("\n");
	*/
	free(InterleavedInputMessage);
	free(InterleavedOutputMessage);
	return 1;
}
bool MessageTransmit(int*MessageInputStream, int Len, int*Output)
{
	int*ParityStream = (int*)malloc((Len + 4)*sizeof(int));
	int*InterleavedParityStream = (int*)malloc((Len + 4)*sizeof(int));
	int*MessageStream = (int*)malloc((Len + 4)*sizeof(int));
	int*InternalInterleaverTable = (int*)malloc(Len*sizeof(int));
	int*SubBlockInterleaveTable = (int*)malloc((Len + 4)*sizeof(int));
	//编码 
	ReadTable(InternalInterleaverTable, Len, "Random");
	TurboEncoder(MessageInputStream, Len, InternalInterleaverTable, ParityStream, InterleavedParityStream, MessageStream);
	//临时输出以检测结果
	/*
	for (int i = 0; i < Len + 4; ++i)
		printf("%d ", MessageStream[i]);
	printf("\n");
	for (int i = 0; i < Len + 4; ++i)
		printf("%d ", ParityStream[i]);
	printf("\n");
	for (int i = 0; i < Len + 4; ++i)
		printf("%d ", InterleavedParityStream[i]);
	printf("\n");
	*/
	//译码后三条流分别交织并收集
	//先处理d0 d1

	//SubBlockInterleaveTableGenerator(SubBlockInterleaveTable, Len + 4, 0);
	ReadTable(SubBlockInterleaveTable, Len + 4, "SubBlock_M_P1");
	for (int i = 0; i < Len + 4; ++i)
	{
		Output[i] = MessageStream[SubBlockInterleaveTable[i]];
		Output[Len + 4 + 2 * i] = ParityStream[SubBlockInterleaveTable[i]];
	}
	//再处理d2
	//SubBlockInterleaveTableGenerator(SubBlockInterleaveTable, Len + 4, 1); //此处有错,3月24日已更改
	ReadTable(SubBlockInterleaveTable, Len + 4, "SubBlock_P2");
	for (int i = 0; i < Len + 4; ++i)
		Output[Len + 4 + 2 * i + 1] = InterleavedParityStream[SubBlockInterleaveTable[i]];
		
	//RateMatching()
	//穿孔
	free(ParityStream);
	free(InterleavedParityStream);
	free(MessageStream);
	free(InternalInterleaverTable);
	free(SubBlockInterleaveTable);
	return 1;
}
bool AddGaussianNoise(int*Message, int Len, float SNR,float*ReceivedMessage) //Messsage为二进制0 1序列，输出为-1 1加上噪声
{
	float N0 = pow(10, -SNR / 10);
	float*GaussianNoise = (float*)malloc(Len*sizeof(float));
	VSLStreamStatePtr stream;
	__int32 SEED = GetCpuCycle();
	vslNewStream(&stream, VSL_BRNG_MT19937, SEED);
	vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, Len, GaussianNoise, 0, sqrt(N0));
	for (int i = 0; i < Len; ++i)
	{
		//printf("%.2f ", GaussianNoise[i]);
		ReceivedMessage[i] = 2 * ((2 * Message[i] - 1) + GaussianNoise[i]) / N0;
	}
	//printf("\n");
	vslDeleteStream(&stream);
	free(GaussianNoise);
	return 1;
}
bool NoNoise(int*Message, int Len, float SNR, float*ReceivedMessage)
{
	for (int i = 0; i < Len; ++i)
	{
		//printf("%.2f ", GaussianNoise[i]);
		ReceivedMessage[i] = (__int64)(GetCpuCycle()%8)*(2 * Message[i] - 1);
	}
	return 1;
}
bool QuantizeFloatToShort(float*input, int Len, short*output)
{
	for (int i = 0; i < Len; ++i)
	{
		
		output[i] = input[i] * 16;
		output[i] = output[i]>255 ? 255 : ((output[i] < -255) ? -255 : output[i]);
	}
		return 1;
}

bool InputDivide(float*InputMessage,int Len, float*OutputMessage,float*OutputInterleaveredMessage,float*OutputParity1,float*OutputParity2)//未完待续，输出长度为Len+Mem
{
	//RateMatching()
	//提取出三条比特流
	if (Len % 3 != 0)
	{
		printf("Error In Received Length!(Function:InputDivide)\n");
		system("pause");
	}
	float*TmpMessage = (float*)_mm_malloc((Len / 3)*sizeof(float), sizeof(__m256));
	float*TmpParity1 = (float*)_mm_malloc((Len / 3)*sizeof(float), sizeof(__m256));
	float*TmpParity2 = (float*)_mm_malloc((Len / 3)*sizeof(float), sizeof(__m256));
	float*InterleavedMessage = (float*)_mm_malloc((Len / 3)*sizeof(float), sizeof(__m256));
	float*InterleavedParity1 = (float*)_mm_malloc((Len / 3)*sizeof(float), sizeof(__m256));
	float*InterleavedParity2 = (float*)_mm_malloc((Len / 3)*sizeof(float), sizeof(__m256));
	int*Table = (int*)malloc((Len / 3)*sizeof(int));
	for (int i = 0; i < Len / 3; ++i)
	{
		TmpMessage[i] = InputMessage[i];
		TmpParity1[i] = InputMessage[Len/3 + 2 * i];
		TmpParity2[i] = InputMessage[Len/3 + 2 * i + 1];
	}
	//临时输出以检测结果
	ReadTable(Table, Len / 3, "SubBlock_M_P1");
	for (int i = 0; i < Len / 3; ++i)
	{
		InterleavedMessage[Table[i]] = TmpMessage[i];
		InterleavedParity1[Table[i]] = TmpParity1[i];
	}
	ReadTable(Table, Len / 3, "SubBlock_P2");
	for (int i = 0; i < Len / 3; ++i)
		InterleavedParity2[Table[i]] = TmpParity2[i];
	/*
	printf("尾处前Mesag：");
	for (int i = 0; i < Len / 3; ++i)
		printf("%d ", InterleavedMessage[i]>0 ? 1 : 0);
	printf("\n");
	printf("尾处前Pari1：");
	for (int i = 0; i < Len / 3; ++i)
		printf("%d ", InterleavedParity1[i]>0 ? 1 : 0);
	printf("\n");
	printf("尾处前Pari2：");
	for (int i = 0; i < Len / 3; ++i)
		printf("%d ", InterleavedParity2[i]>0 ? 1 : 0);
	printf("\n");
	*/
	for (int i = 0; i < Len/3-4; ++i)
	{
		OutputMessage[i] = InterleavedMessage[i];
		OutputParity1[i] = InterleavedParity1[i];
		OutputParity2[i] = InterleavedParity2[i];
	}
	//处理InterleavedMessage
	ReadTable(Table, Len / 3-4, "Random");
	for (int i = 0; i < Len / 3 - 4; ++i)
		OutputInterleaveredMessage[i] = OutputMessage[Table[i]];
	//接下来还原尾巴比特，总长Len+Mem=
	OutputMessage[Len/3-4] = InterleavedMessage[Len/3-4];
	OutputMessage[Len/3-3] = InterleavedParity2[Len/3-4];
	OutputMessage[Len/3-2] = InterleavedParity1[Len/3-3];
	OutputParity1[Len/3-4] = InterleavedParity1[Len/3-4];
	OutputParity1[Len/3-3] = InterleavedMessage[Len/3-3];
	OutputParity1[Len/3-2] = InterleavedParity2[Len/3-3];
	OutputParity2[Len/3-4] = InterleavedParity1[Len/3-2];
	OutputParity2[Len/3-3] = InterleavedMessage[Len/3-1];
	OutputParity2[Len/3-2] = InterleavedParity2[Len/3-1];
	OutputInterleaveredMessage[Len / 3 - 4] = InterleavedMessage[Len / 3 - 2];
	OutputInterleaveredMessage[Len / 3 - 3] = InterleavedParity2[Len / 3 - 2];
	OutputInterleaveredMessage[Len / 3 - 2] = InterleavedParity1[Len / 3 - 1];
	_mm_free(TmpMessage);
	_mm_free(TmpParity1);
	_mm_free(TmpParity2);
	_mm_free(InterleavedMessage);
	_mm_free(InterleavedParity1);
	_mm_free(InterleavedParity2);
	free(Table);
	return 1;
}
void MaxLogSubDecoderSInt(short*Ls, short*Lp, short*Le, int Len, short*Output) //Len长度不含4位尾比特，不记寄存器个数(Mem)
{
	//定义alpha，beta
	int StateNum = 8;
	int OffsetTmp_p = -INF;
	int OffsetTmp_n = -INF;
	short *alpha = (short *)malloc(StateNum*(Len)*sizeof(short));
	short *beta = (short *)malloc(StateNum*(Len)*sizeof(short));
	//定义四种取值不同的分支度量
	short *r0 = (short*)malloc((Len)*sizeof(short));
	short *r_0 = (short*)malloc((Len)*sizeof(short));
	short *r1 = (short*)malloc((Len)*sizeof(short));
	short *r_1 = (short*)malloc((Len)*sizeof(short));
	//初始化alpha，beta
	alpha[0] = 0;
	alpha[1] = -INF;
	alpha[2] = -INF;
	alpha[3] = -INF;
	alpha[4] = -INF;
	alpha[5] = -INF;
	alpha[6] = -INF;
	alpha[7] = -INF;
	beta[(Len - 1)*StateNum + 0] = 0;
	beta[(Len - 1)*StateNum + 1] = -INF;
	beta[(Len - 1)*StateNum + 2] = -INF;
	beta[(Len - 1)*StateNum + 3] = -INF;
	beta[(Len - 1)*StateNum + 4] = -INF;
	beta[(Len - 1)*StateNum + 5] = -INF;
	beta[(Len - 1)*StateNum + 6] = -INF;
	beta[(Len - 1)*StateNum + 7] = -INF;
	//初始化四个分支度量
	for (int i = 0; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Le[i] + Ls[i]);
		r_0[i] = -0.5*(Lp[i] + Le[i] + Ls[i]);
		r1[i] = 0.5*(-Lp[i] + Le[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	for (int i = Len - Mem; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Ls[i]);
		r_0[i] = -0.5*(Lp[i] + Ls[i]);
		r1[i] = 0.5*(-Lp[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	//计算alpha,
	OffsetTmp_n = StateNum;
	OffsetTmp_p = 0;
	int k = 0;
	for (int i = 0; i < Len - 1; ++i)
	{
		alpha[OffsetTmp_n] = max(alpha[OffsetTmp_p] + r_0[i], alpha[OffsetTmp_p + 1] + r0[i]);
		alpha[OffsetTmp_n + 1] = max(alpha[OffsetTmp_p + 2] + r1[i], alpha[OffsetTmp_p + 3] + r_1[i]);
		alpha[OffsetTmp_n + 2] = max(alpha[OffsetTmp_p + 4] + r_1[i], alpha[OffsetTmp_p + 5] + r1[i]);
		alpha[OffsetTmp_n + 3] = max(alpha[OffsetTmp_p + 6] + r0[i], alpha[OffsetTmp_p + 7] + r_0[i]);
		alpha[OffsetTmp_n + 4] = max(alpha[OffsetTmp_p + 0] + r0[i], alpha[OffsetTmp_p + 1] + r_0[i]);
		alpha[OffsetTmp_n + 5] = max(alpha[OffsetTmp_p + 3] + r1[i], alpha[OffsetTmp_p + 2] + r_1[i]);
		alpha[OffsetTmp_n + 6] = max(alpha[OffsetTmp_p + 4] + r1[i], alpha[OffsetTmp_p + 5] + r_1[i]);
		alpha[OffsetTmp_n + 7] = max(alpha[OffsetTmp_p + 7] + r0[i], alpha[OffsetTmp_p + 6] + r_0[i]);
		//if (i % 8 == 7)
		/*
		if (i % 8 == 0)
		{
		alpha[OffsetTmp_n] -= alpha[OffsetTmp_n];
		alpha[OffsetTmp_n + 1] -= alpha[OffsetTmp_n];
		alpha[OffsetTmp_n + 2] -= alpha[OffsetTmp_n];
		alpha[OffsetTmp_n + 3] -= alpha[OffsetTmp_n];
		alpha[OffsetTmp_n + 4] -= alpha[OffsetTmp_n];
		alpha[OffsetTmp_n + 5] -= alpha[OffsetTmp_n];
		alpha[OffsetTmp_n + 6] -= alpha[OffsetTmp_n];
		alpha[OffsetTmp_n + 7] -= alpha[OffsetTmp_n];
		}
		*/
		OffsetTmp_n += StateNum;
		OffsetTmp_p += StateNum;
		k++;
	}
	//计算beta
	k = 0;
	OffsetTmp_n = (Len - 1)* StateNum;
	OffsetTmp_p = (Len - 2)*StateNum;
	for (int i = Len - 1; i > 0; --i)
	{
		beta[OffsetTmp_p] = max(beta[OffsetTmp_n] + r_0[i], beta[OffsetTmp_n + 4] + r0[i]);
		beta[OffsetTmp_p + 1] = max(beta[OffsetTmp_n] + r0[i], beta[OffsetTmp_n + 4] + r_0[i]);
		beta[OffsetTmp_p + 2] = max(beta[OffsetTmp_n + 1] + r1[i], beta[OffsetTmp_n + 5] + r_1[i]);
		beta[OffsetTmp_p + 3] = max(beta[OffsetTmp_n + 5] + r1[i], beta[OffsetTmp_n + 1] + r_1[i]);
		beta[OffsetTmp_p + 4] = max(beta[OffsetTmp_n + 6] + r1[i], beta[OffsetTmp_n + 2] + r_1[i]);
		beta[OffsetTmp_p + 5] = max(beta[OffsetTmp_n + 2] + r1[i], beta[OffsetTmp_n + 6] + r_1[i]);
		beta[OffsetTmp_p + 6] = max(beta[OffsetTmp_n + 3] + r0[i], beta[OffsetTmp_n + 7] + r_0[i]);
		beta[OffsetTmp_p + 7] = max(beta[OffsetTmp_n + 3] + r_0[i], beta[OffsetTmp_n + 7] + r0[i]);
		k++;
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
	for (int i = 0; i < Len; ++i)
	{

		Output[i] =
			Mmax(
			Mmax(                                     //之前为i+1错误
			Mmax(alpha[i*StateNum + 0] + r0[i] + beta[i*StateNum + 4], alpha[i*StateNum + 1] + r0[i] + beta[i*StateNum + 0]),
			Mmax(alpha[i*StateNum + 2] + r1[i] + beta[i*StateNum + 1], alpha[i*StateNum + 3] + r1[i] + beta[i*StateNum + 5])
			),
			Mmax(
			Mmax(alpha[i*StateNum + 4] + r1[i] + beta[i*StateNum + 6], alpha[i*StateNum + 5] + r1[i] + beta[i*StateNum + 2]),
			Mmax(alpha[i*StateNum + 6] + r0[i] + beta[i*StateNum + 3], alpha[i*StateNum + 7] + r0[i] + beta[i*StateNum + 7])
			)
			)
			-
			Mmax(
			Mmax(
			Mmax(alpha[i*StateNum + 0] + r_0[i] + beta[i*StateNum + 0], alpha[i*StateNum + 1] + r_0[i] + beta[i*StateNum + 4]),
			Mmax(alpha[i*StateNum + 2] + r_1[i] + beta[i*StateNum + 5], alpha[i*StateNum + 3] + r_1[i] + beta[i*StateNum + 1])
			),
			Mmax(
			Mmax(alpha[i*StateNum + 4] + r_1[i] + beta[i*StateNum + 2], alpha[i*StateNum + 5] + r_1[i] + beta[i*StateNum + 6]),
			Mmax(alpha[i*StateNum + 6] + r_0[i] + beta[i*StateNum + 7], alpha[i*StateNum + 7] + r_0[i] + beta[i*StateNum + 3])
			)
			)
			- Ls[i] - Le[i];
	}
	free(alpha);
	free(beta);
	free(r0);
	free(r1);
	free(r_1);
	free(r_0);
}
bool MaxLogRscDecoderAVXFloat(float*Ls, float*Lp, float*La, int Len, float*Output,float*AlphaM,float*BetaPlusM,float*BetaMinusM)//OutPut以8对齐
{
	//AlphaM为8*Len大小
	//声明
	__m256 AlphaBetaMain, AlphaBetaSub;//用于计算前向度量，后向度量
	__m256 RMask, RLp, Rr;//用于计算分支度量
	__m256 RTmp;//用于重排
	__m256i Offset = _mm256_set_epi32(56, 48, 40, 32, 24, 16, 8, 0);
	/*****一、先正向计算Alpha*****/
	AlphaBetaMain = _mm256_set_ps(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);//AlphaBetaMain为上一次状态
	float*AlphaIndex = AlphaM;
	_mm256_store_ps(AlphaIndex, AlphaBetaMain);//将起始状态放入内存;
	AlphaIndex += 8;
	RMask = _mm256_set_ps(1, 1, -1, -1, -1, -1, 1, 1);
	for (int i = 0; i < Len-1; ++i)//注意起始和结束位置是这样的么？
	{                              //这样算完后内存为0到Len-1共Len个元素,Alpha指针指向的为i+1
		//1.先计算分支度量
		if (i < Len - Mem)
			Rr = _mm256_set1_ps(0.5*(Ls[i] + La[i]));
		else Rr = _mm256_set1_ps(0.5*Ls[i]);
		RLp = _mm256_set1_ps(0.5*Lp[i]);
		RLp = _mm256_mul_ps(RLp, RMask);
		Rr = _mm256_add_ps(Rr, RLp);
		//2.处理前向度量
		//AlphaBetaMain = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm256_add_ps(AlphaBetaMain, Rr);
		AlphaBetaMain = _mm256_permute_ps(AlphaBetaMain, 0x0000009C);
		RTmp = _mm256_permute2f128_ps(AlphaBetaMain, AlphaBetaMain, 0x00000001);
		AlphaBetaMain = _mm256_permute_ps(AlphaBetaMain, 0x0000004E);
		AlphaBetaMain = _mm256_blend_ps(AlphaBetaMain, RTmp, 0x0000003C);
		//对于相减的部分
		AlphaBetaSub = _mm256_sub_ps(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permute_ps(AlphaBetaSub, 0x000000C9);
		RTmp = _mm256_permute2f128_ps(AlphaBetaSub, AlphaBetaSub, 0x00000001);
		AlphaBetaSub = _mm256_permute_ps(AlphaBetaSub, 0x0000004E);
		AlphaBetaSub = _mm256_blend_ps(AlphaBetaSub, RTmp, 0x0000003C);
		//3.选取最大值,并更新到K+1并存入内存
		AlphaBetaMain = _mm256_max_ps(AlphaBetaMain, AlphaBetaSub);
		_mm256_store_ps(AlphaIndex, AlphaBetaMain);
		AlphaIndex += 8;
	}
	/*****二、再反向计算Beta*****/
	//注，此处不用存入内存，直接连续计算到Alpha加Beta，再存入内存,因此没有存储Beta的内存区块。
	float*BetaPlusIndex = BetaPlusM + (Len << 3);//指向最近存储的位置
	float*BetaMinusIndex = BetaMinusM + (Len  << 3);//指向最近存储的位置
	AlphaIndex = AlphaM + (Len << 3);//指向最近存储的位置
	float*OutputIndex = Output; //OutputIndex指向最近存储的位置，不同于以上(注：OutPut大小不必为8的倍数)
	AlphaBetaMain = _mm256_set_ps(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	RMask = _mm256_set_ps(1, -1, -1, 1, 1, -1, -1, 1);
	//Reference
	float *Beta = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float),sizeof(__m256));//可删去
	float*BetaIndex = Beta + (Len << 3);
	for (int i = Len - 1; i >= 0; --i)
	{
		//1.还是先计算分支度量
		if (i < Len - Mem)
			Rr = _mm256_set1_ps(0.5*(Ls[i] + La[i]));
		else Rr = _mm256_set1_ps(0.5*Ls[i]);
		RLp = _mm256_set1_ps(0.5*Lp[i]);
		RLp = _mm256_mul_ps(RLp, RMask);
		Rr = _mm256_add_ps(Rr, RLp);
		//2.计算后向度量，到计算BetaPlus和BetaAdd
		//AlphaBetaMain = _mm256_set_ps(7,6,5,4,3,2,1,0);
		AlphaBetaSub = AlphaBetaMain;
		//先加
		AlphaBetaMain = _mm256_add_ps(AlphaBetaMain, Rr);
		AlphaBetaMain = _mm256_permute_ps(AlphaBetaMain, 0x00000078);
		RTmp = _mm256_permute2f128_ps(AlphaBetaMain, AlphaBetaMain, 0x00000001);
		AlphaBetaMain = _mm256_permute_ps(AlphaBetaMain, 0x000000B1);
		AlphaBetaMain = _mm256_blend_ps(AlphaBetaMain, RTmp, 0x00000069);
		//再计算减
		AlphaBetaSub = _mm256_sub_ps(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permute_ps(AlphaBetaSub, 0x000000D2);
		RTmp = _mm256_permute2f128_ps(AlphaBetaSub, AlphaBetaSub, 0x00000001);
		AlphaBetaSub = _mm256_permute_ps(AlphaBetaSub, 0x000000B1);
		AlphaBetaSub = _mm256_blend_ps(AlphaBetaSub, RTmp, 0x00000096);
		//3.Alpha与Beta相加的步骤
		AlphaIndex -= 8;
		BetaPlusIndex -= 8;
		BetaMinusIndex -= 8;
		Rr = _mm256_load_ps(AlphaIndex); //Rr之后用不到了，当做Alpha的存储器
		RTmp = _mm256_add_ps(AlphaBetaSub, Rr);//BetaSub加上Alpha,RTmp之后也用不到了
		_mm256_store_ps(BetaMinusIndex, RTmp);
		RTmp = _mm256_add_ps(AlphaBetaMain, Rr);
		_mm256_store_ps(BetaPlusIndex, RTmp);
		/*4.计算Hmax，并得出软输出结果
		Count++;
		int E = Count & 0x7;
		if ((Count & 0x7) == 0)//满8，为1000形式时
		{
			int k = 0;
			Rr = _mm256_i32gather_ps(BetaPlusIndex, Offset, 4);
			RLp = _mm256_i32gather_ps(BetaMinusIndex, Offset, 4);
			for (k = 1; k < 8; ++k)//求出最大值（Hmax）
			{
				RTmp = _mm256_i32gather_ps(BetaPlusIndex + k, Offset, 4);//先处理BetaPlus（Rr）
				Rr = _mm256_max_ps(Rr,RTmp);
				RTmp = _mm256_i32gather_ps(BetaMinusIndex + k, Offset, 4);//再处理BetaMinus（RLp）
				RLp = _mm256_max_ps(RLp, RTmp);
			}
			//注意，OutPut为从大到小存储,且起始位置对齐
			Rr = _mm256_sub_ps(RLp, Rr);//此处为Hmax相减，完成后RLp作废
			for (int i = 0; i < Len; ++i)
				printf("%0.1f", Ls[i]);
			printf("\n");
			RLp = _mm256_load_ps(Ls+i);//此时i指向最后8个
			RTmp = _mm256_load_ps(La+i);
			RLp = _mm256_add_ps(RLp,RTmp);
			Rr = _mm256_sub_ps(Rr, RLp);
			OutputIndex -= 8;
			_mm256_store_ps(OutputIndex, Rr);
		}//if (Count & 0x111 == 0)
		*/
		//5.计算下一步的Beta
		AlphaBetaMain = _mm256_max_ps(AlphaBetaMain, AlphaBetaSub);
		BetaIndex -= 8;//可删去
		_mm256_store_ps(BetaIndex, AlphaBetaMain);//可删去
	}
	/*****三、处理最后剩下的软输出结果
	int E = Count & 0x7;
	if ((Count & 0x7) != 0)
	{
		int k = 0;
		Rr = _mm256_i32gather_ps(BetaPlusM, Offset, 4);//此时BetaPlusM==BetaPlusIndex
		RLp = _mm256_i32gather_ps(BetaMinusM, Offset, 4);
		for (k = 1; k < 8; ++k)//求出最大值（Hmax）
		{
			RTmp = _mm256_i32gather_ps(BetaPlusM + k, Offset, 4);//先处理BetaPlus（Rr）
			Rr = _mm256_max_ps(Rr, RTmp);
			RTmp = _mm256_i32gather_ps(BetaMinusM + k, Offset, 4);//再处理BetaMinus（RLp）
			RLp = _mm256_max_ps(RLp, RTmp);
		}
		Rr = _mm256_sub_ps(RLp, Rr);//此处为Hmax相减，完成后RLp作废
		RLp = _mm256_load_ps(Ls);//此时i指向最后8个
		RTmp = _mm256_load_ps(La);
		RLp = _mm256_add_ps(RLp, RTmp);
		Rr = _mm256_sub_ps(Rr, RLp);
		_mm256_store_ps(Output, Rr);
	}
	*/
	/*****三、处理软输出结果*****/
	BetaPlusIndex = BetaPlusM;
	BetaMinusIndex = BetaMinusM;
	OutputIndex = Output;
	for (int Count= 0; Count < Len; Count+=8)
	{
		Rr = _mm256_i32gather_ps(BetaPlusIndex, Offset, 4);
		RLp = _mm256_i32gather_ps(BetaMinusIndex, Offset, 4);
		for (int k = 1; k < 8; ++k)//求出最大值（Hmax）
		{
			RTmp = _mm256_i32gather_ps(BetaPlusIndex + k, Offset, 4);//先处理BetaPlus（Rr）
			Rr = _mm256_max_ps(Rr, RTmp);
			RTmp = _mm256_i32gather_ps(BetaMinusIndex + k, Offset, 4);//再处理BetaMinus（RLp）
			RLp = _mm256_max_ps(RLp, RTmp);
		}
		//注意，OutPut为从大到小存储,且起始位置对齐
		//Rr = _mm256_sub_ps(RLp, Rr);//此处为Hmax相减，完成后RLp作废
		Rr = _mm256_sub_ps(Rr, RLp); 
		RLp = _mm256_load_ps(Ls+Count);//此时i指向最后8个
		RTmp = _mm256_load_ps(La+Count);
		RLp = _mm256_add_ps(RLp, RTmp);
		Rr = _mm256_sub_ps(Rr, RLp);
		_mm256_store_ps(OutputIndex, Rr);
		OutputIndex += 8;
		BetaPlusIndex += 64;
		BetaMinusIndex += 64;
	}
	return 1;
}
inline bool MaxLogRscDecoderSSEIntD0(short*Ls, short*Lp, short*La, int Len, short*Output, __m128i*BetaM)
{//注Common大小为8*__m256i
	__m128i AlphaBetaMain, AlphaBetaSub, Offset, RL1, RL0;
	__m256i R1, R2, R3, R4, R5, R6, R7, Rr, RTmp;
	__m128i ShuffleAlphaBetaPlus = _mm_setr_epi8(8, 9, 0, 1, 2, 3, 10, 11, 12, 13, 4, 5, 6, 7, 14, 15);//低位为Plus，高位为Sub
	__m128i ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 8, 9, 10, 11, 2, 3, 4, 5, 12, 13, 14, 15, 6, 7);
	__m256i ShuffleLshift64Invert = _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
	__m256i ShuffleLshift32 = _mm256_setr_epi8(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3);
	__m256i ShuffleLshift16 = _mm256_setr_epi8(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1);
	short L1;
	short L0;
	__m128i*BetaIndex = BetaM + Len;//指向最后一个元素的下一个位置
	__m128i*RrIndex = (__m128i*)(Ls + Len - 8);
	__m128i*RaIndex = (__m128i*)(La + Len - 8);
	__m128i*RpIndex = (__m128i*)(Lp + Len - 8);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	_mm_store_si128(--BetaIndex, AlphaBetaMain);
	//初始化Rr
	RL1 = _mm_loadu_si128(RrIndex);
	Offset = _mm_loadu_si128(RaIndex);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_loadu_si128(RpIndex);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
	//重新定义Rr值
	RrIndex = (__m128i*)(Ls + ((Len >> 3) << 3) - 8);
	RaIndex = (__m128i*)(La + ((Len >> 3) << 3) - 8);
	RpIndex = (__m128i*)(Lp + ((Len >> 3) << 3) - 8);
	/*先计算后向度量*/
	for (int i = Len - 1; i>0; --i)
	{
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		/*先得出分支度量*/
		//L1 = (Ls[i] + La[i] + Lp[i]) >> 1;
		//L0 = L1 - Lp[i];
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//归一化
		if ((i & 0x7) == 0)//相当于每8次运算中有一次归一化
		{
			//更新
			RL1 = _mm_load_si128(RrIndex--);
			Offset = _mm_load_si128(RaIndex--);
			RL1 = _mm_add_epi16(RL1, Offset);
			Offset = _mm_load_si128(RpIndex--);
			RL0 = _mm_sub_epi16(RL1, Offset);
			RL1 = _mm_add_epi16(RL1, Offset);
			RL0 = _mm_srai_epi16(RL0, 1);
			RL1 = _mm_srai_epi16(RL1, 1);
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
			//用AVX2指令实现
			Offset = _mm_broadcastw_epi16(AlphaBetaMain);
			AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		}
		_mm_store_si128(--BetaIndex, AlphaBetaMain);
	}

	/*再计算前向度量并计算似然比*/
	__m128i*LsIndex = (__m128i*)Ls;
	__m128i*LaIndex = (__m128i*)La;
	__m128i*OutputIndex = (__m128i*)Output;
	BetaIndex = BetaM;
	RrIndex = (__m128i*)Ls;
	RaIndex = (__m128i*)La;
	RpIndex = (__m128i*)Lp;
	ShuffleLshift64Invert = _mm256_setr_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
	ShuffleAlphaBetaPlus = _mm_setr_epi8(2, 3, 4, 5, 10, 11, 12, 13, 0, 1, 6, 7, 8, 9, 14, 15);
	ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 6, 7, 8, 9, 14, 15, 2, 3, 4, 5, 10, 11, 12, 13);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	//初始化Rr
	RL1 = _mm_load_si128(RrIndex++);
	Offset = _mm_load_si128(RaIndex++);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_load_si128(RpIndex++);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	for (int i = 0; i < Len; ++i)
	{
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//先得出分支度量
		//L1 = (Ls[i] + La[i] + Lp[i]) >> 1;
		//L0 = L1 - Lp[i];
		//Offset = _mm_set_epi16(L1, L1, L0, L0, L0, L0, L1, L1);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(BetaIndex++);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		/*******处理Hmax*******/
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		Rr = _mm256_max_epi16(Rr, RTmp);
		switch (i & 0x7)
		{
		case 0:R1 = Rr; break;//低位为Plus
		case 1:R2 = Rr; break;
		case 2:R3 = Rr; break;
		case 3:R4 = Rr; break;
		case 4:R5 = Rr; break;
		case 5:R6 = Rr; break;
		case 6:R7 = Rr; break;
		case 7:
		{
				  //更新Rr
				  RL1 = _mm_load_si128(RrIndex++);
				  Offset = _mm_load_si128(RaIndex++);
				  RL1 = _mm_add_epi16(RL1, Offset);
				  Offset = _mm_load_si128(RpIndex++);
				  RL0 = _mm_sub_epi16(RL1, Offset);
				  RL1 = _mm_add_epi16(RL1, Offset);
				  RL0 = _mm_srai_epi16(RL0, 1);
				  RL1 = _mm_srai_epi16(RL1, 1);
				  //归一化
				  Offset = _mm_broadcastw_epi16(AlphaBetaMain);
				  AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
				  //此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
				  R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
				  R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
				  R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
				  R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
				  R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
				  R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
				  R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
				  //Plus与Minus相减,到此Rr的低八位数值全一样
				  Offset = _mm256_extractf128_si256(R1, 1);
				  Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
				  //减去La
				  AlphaBetaSub = _mm_load_si128(LaIndex++);
				  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
				  //存入内存
				  _mm_store_si128(OutputIndex++, Offset);
				  break;
		}
		}
	}
	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去Ls+La
		AlphaBetaSub = _mm_load_si128(LsIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_load_si128(LaIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		_mm_store_si128(OutputIndex++, Offset);
	}
	return 1;
}
inline bool MaxLogRscDecoderSSEIntD1(short*Ls, short*Lp, int Len, short*Output, __m128i*BetaM)
{//注Common大小为8*__m256i
	__m128i AlphaBetaMain, AlphaBetaSub, Offset, RL1, RL0;
	__m256i R1, R2, R3, R4, R5, R6, R7, Rr, RTmp;
	__m128i ShuffleAlphaBetaPlus = _mm_setr_epi8(8, 9, 0, 1, 2, 3, 10, 11, 12, 13, 4, 5, 6, 7, 14, 15);//低位为Plus，高位为Sub
	__m128i ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 8, 9, 10, 11, 2, 3, 4, 5, 12, 13, 14, 15, 6, 7);
	__m256i ShuffleLshift64Invert = _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
	__m256i ShuffleLshift32 = _mm256_setr_epi8(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3);
	__m256i ShuffleLshift16 = _mm256_setr_epi8(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1);
	short L1;
	short L0;
	__m128i*BetaIndex = BetaM + Len;//指向最后一个元素的下一个位置
	__m128i*RrIndex = (__m128i*)(Ls + Len - 8);
	__m128i*RpIndex = (__m128i*)(Lp + Len - 8);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	_mm_store_si128(--BetaIndex, AlphaBetaMain);
	//初始化Rr
	RL1 = _mm_loadu_si128(RrIndex);
	Offset = _mm_loadu_si128(RpIndex);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
	//重新定义Rr值
	RrIndex = (__m128i*)(Ls + ((Len >> 3) << 3) - 8);
	RpIndex = (__m128i*)(Lp + ((Len >> 3) << 3) - 8);
	/*先计算后向度量*/
	for (int i = Len - 1; i>0; --i)
	{
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		/*先得出分支度量*/
		//L1 = (Ls[i] + La[i] + Lp[i]) >> 1;
		//L0 = L1 - Lp[i];
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//归一化
		if ((i & 0x7) == 0)//相当于每8次运算中有一次归一化
		{
			//更新
			RL1 = _mm_load_si128(RrIndex--);
			Offset = _mm_load_si128(RpIndex--);
			RL0 = _mm_sub_epi16(RL1, Offset);
			RL1 = _mm_add_epi16(RL1, Offset);
			RL0 = _mm_srai_epi16(RL0, 1);
			RL1 = _mm_srai_epi16(RL1, 1);
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
			//用AVX2指令实现
			Offset = _mm_broadcastw_epi16(AlphaBetaMain);
			AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		}
		_mm_store_si128(--BetaIndex, AlphaBetaMain);
	}

	/*再计算前向度量并计算似然比*/
	__m128i*LsIndex = (__m128i*)Ls;
	__m128i*OutputIndex = (__m128i*)Output;
	BetaIndex = BetaM;
	RrIndex = (__m128i*)Ls;
	RpIndex = (__m128i*)Lp;
	ShuffleLshift64Invert = _mm256_setr_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
	ShuffleAlphaBetaPlus = _mm_setr_epi8(2, 3, 4, 5, 10, 11, 12, 13, 0, 1, 6, 7, 8, 9, 14, 15);
	ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 6, 7, 8, 9, 14, 15, 2, 3, 4, 5, 10, 11, 12, 13);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	//初始化Rr
	RL1 = _mm_load_si128(RrIndex++);
	Offset = _mm_load_si128(RpIndex++);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	for (int i = 0; i < Len; ++i)
	{
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//先得出分支度量
		//L1 = (Ls[i] + La[i] + Lp[i]) >> 1;
		//L0 = L1 - Lp[i];
		//Offset = _mm_set_epi16(L1, L1, L0, L0, L0, L0, L1, L1);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(BetaIndex++);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		/*******处理Hmax*******/
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		Rr = _mm256_max_epi16(Rr, RTmp);
		switch (i & 0x7)
		{
		case 0:R1 = Rr; break;//低位为Plus
		case 1:R2 = Rr; break;
		case 2:R3 = Rr; break;
		case 3:R4 = Rr; break;
		case 4:R5 = Rr; break;
		case 5:R6 = Rr; break;
		case 6:R7 = Rr; break;
		case 7:
		{
				  //更新Rr
				  RL1 = _mm_load_si128(RrIndex++);
				  Offset = _mm_load_si128(RpIndex++);
				  RL0 = _mm_sub_epi16(RL1, Offset);
				  RL1 = _mm_add_epi16(RL1, Offset);
				  RL0 = _mm_srai_epi16(RL0, 1);
				  RL1 = _mm_srai_epi16(RL1, 1);
				  //归一化
				  Offset = _mm_broadcastw_epi16(AlphaBetaMain);
				  AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
				  //此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
				  R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
				  R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
				  R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
				  R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
				  R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
				  R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
				  R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
				  //Plus与Minus相减,到此Rr的低八位数值全一样
				  Offset = _mm256_extractf128_si256(R1, 1);
				  Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
				  //减去Ls+La
				  AlphaBetaSub = _mm_load_si128(LsIndex++);
				  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
				  //存入内存
				  _mm_store_si128(OutputIndex++, Offset);
				  break;
		}
		}
	}
	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去Ls+La
		AlphaBetaSub = _mm_load_si128(LsIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		_mm_store_si128(OutputIndex++, Offset);
	}
	return 1;
}
inline bool MaxLogRscDecoderAVXIntDulD1(short*Ls1, short*Ls2, short*Lp1, short*Lp2,int Len, short*Output1, short*Output2, __m256i*BetaM, __m256i*CommenPlusM, __m256i*CommenMinusM)
{//注Common大小为8*__m256i
	__m256i AlphaBetaMain, AlphaBetaSub, Rr, Offset, RMask;
	__m256i OffsetTable = _mm256_set_epi32(112, 96, 80, 64, 48, 32, 16, 0);
	__m256i ShuffleAlphaPlus = _mm256_setr_epi32(1, 2, 5, 6, 0, 3, 4, 7);
	__m256i ShuffleAlphaMinus = _mm256_setr_epi32(0, 3, 4, 7, 1, 2, 5, 6);
	__m256i ShuffleBetaPlus = _mm256_setr_epi32(4, 0, 1, 5, 6, 2, 3, 7);
	__m256i ShuffleBetaMinus = _mm256_setr_epi32(0, 4, 5, 1, 2, 6, 7, 3);
	__m256i ShuffleTableBeforeStore = _mm256_setr_epi8(0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15);
	short L1_1;
	short L1_0;
	short L2_1;
	short L2_0;
	__m256i*BetaIndex = BetaM + Len;//指向最后一个元素的下一个位置
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);
	_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	/*先计算后向度量*/
	for (int i = Len - 1; i>0; --i)
	{
		AlphaBetaSub = AlphaBetaMain;
		/*先得出分支度量*/
		L1_1 = (Ls1[i] + Lp1[i]) >> 1;
		L1_0 = L1_1 - Lp1[i];
		L2_1 = (Ls2[i] + Lp2[i]) >> 1;
		L2_0 = L2_1 - Lp2[i];
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1);
		//对于相加的部分
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//归一化

		if ((i & 0x7) == 0x7)//相当于每8次运算中有一次归一化
		{
			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);
		}

		_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	}
	/*再计算前向度量并计算似然比*/
	__m128i*Ls1Index = (__m128i*)Ls1;
	__m128i*Ls2Index = (__m128i*)Ls2;
	__m256i*CommonPlusIndex = CommenPlusM;
	__m256i*CommonMinusIndex = CommenMinusM;
	__m128i*Output1Index = (__m128i*)Output1;
	__m128i*Output2Index = (__m128i*)Output2;
	int*CountPlusIndex;
	int*CountMinusIndex;
	BetaIndex = BetaM;
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);

	for (int i = 0; i < Len; ++i)
	{
		AlphaBetaSub = AlphaBetaMain;
		//先得出分支度量
		L1_1 = (Ls1[i] + Lp1[i]) >> 1;
		L1_0 = L1_1 - Lp1[i];
		L2_1 = (Ls2[i] + Lp2[i]) >> 1;
		L2_0 = L2_1 - Lp2[i];
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1);
		//先计算加
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleAlphaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleAlphaMinus);
		//再加上Alpha
		Offset = _mm256_load_si256(BetaIndex++);//加一次为256bit
		Rr = _mm256_add_epi16(AlphaBetaMain, Offset); //此处用Rr代替Main
		Offset = _mm256_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//存入内存
		_mm256_store_si256(CommonPlusIndex++, Rr);//加一次为256bit
		_mm256_store_si256(CommonMinusIndex++, Offset);
		//再存入公共存储区,以下Rr，AlphaBetaSub和Offset可使用
		//直接更新下一次
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		if ((i & 0x7) == 0x7)  //i%8=7时,则此时要计算得出一次OutPut并归一化
		{
			//归一化
			
			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);

			//生成输出信息，公共存储空间归0
			CommonPlusIndex = CommenPlusM;
			CommonMinusIndex = CommenMinusM;
			CountPlusIndex = (int*)CommenPlusM;
			CountMinusIndex = (int*)CommenMinusM;
			Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			
			for (int k = 0; k < 7; ++k)
			{
				Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
				Rr = _mm256_max_epi16(Offset, Rr);//Plus
				Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
				AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
			}
			//Plus减去Minus并重排
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
			Rr = _mm256_permute4x64_epi64(Rr, 0x8D);//11011000
			//减去Ls+La
			AlphaBetaSub = _mm256_loadu2_m128i(Ls2Index++, Ls1Index++);
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
			
		}//if ((i & 0x111) == 0x111)
	}
	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		CountPlusIndex = (int*)CommenPlusM;
		CountMinusIndex = (int*)CommenMinusM;
		Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
		AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
		for (int k = 0; k < 7; ++k)
		{
			Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			Rr = _mm256_max_epi16(Offset, Rr);//Plus
			Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
		}
		//Plus减去Minus并重排
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
		Rr = _mm256_permute4x64_epi64(Rr, 0x8D);
		//减去Ls+La
		AlphaBetaSub = _mm256_loadu2_m128i(Ls2Index++, Ls1Index++);
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
	}
	return 1;
}
inline bool MaxLogRscDecoderAVXIntDulD0(short*Ls1, short*Ls2, short*Lp1, short*Lp2, short*La1, short*La2, int Len, short*Output1, short*Output2, __m256i*BetaM, __m256i*CommenPlusM, __m256i*CommenMinusM)
{//注Common大小为8*__m256i
	__m256i AlphaBetaMain, AlphaBetaSub, Rr, Offset, RMask;
	__m256i OffsetTable = _mm256_set_epi32(112, 96, 80, 64, 48, 32, 16, 0);
	__m256i ShuffleAlphaPlus = _mm256_setr_epi32(1, 2, 5, 6, 0, 3, 4, 7);
	__m256i ShuffleAlphaMinus = _mm256_setr_epi32(0, 3, 4, 7, 1, 2, 5, 6);
	__m256i ShuffleBetaPlus = _mm256_setr_epi32(4, 0, 1, 5, 6, 2, 3, 7);
	__m256i ShuffleBetaMinus = _mm256_setr_epi32(0, 4, 5, 1, 2, 6, 7, 3);
	__m256i ShuffleTableBeforeStore = _mm256_setr_epi8(0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15);
	short L1_1;
	short L1_0;
	short L2_1;
	short L2_0;
	__m256i*BetaIndex = BetaM + Len;//指向最后一个元素的下一个位置
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);
	_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	/*先计算后向度量*/
	for (int i = Len - 1; i>0; --i)
	{
		AlphaBetaSub = AlphaBetaMain;
		/*先得出分支度量*/
		if (i < Len - Mem)
		{
			L1_1 = (Ls1[i] + La1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + La2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		else
		{
			L1_1 = (Ls1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1);
		//对于相加的部分
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//归一化

		if ((i & 0x7) == 0x7)//相当于每8次运算中有一次归一化
		{
			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);
		}

		_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	}
	/*再计算前向度量并计算似然比*/
	__m128i*Ls1Index = (__m128i*)Ls1;
	__m128i*Ls2Index = (__m128i*)Ls2;
	__m128i*La1Index = (__m128i*)La1;
	__m128i*La2Index = (__m128i*)La2;
	__m256i*CommonPlusIndex = CommenPlusM;
	__m256i*CommonMinusIndex = CommenMinusM;
	__m128i*Output1Index = (__m128i*)Output1;
	__m128i*Output2Index = (__m128i*)Output2;
	int*CountPlusIndex;
	int*CountMinusIndex;
	BetaIndex = BetaM;
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);

	for (int i = 0; i < Len; ++i)
	{
		AlphaBetaSub = AlphaBetaMain;
		//先得出分支度量
		if (i < Len - Mem)
		{
			L1_1 = (Ls1[i] + La1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + La2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		else
		{
			L1_1 = (Ls1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1);
		//先计算加
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleAlphaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleAlphaMinus);
		//再加上Alpha
		Offset = _mm256_load_si256(BetaIndex++);//加一次为256bit
		Rr = _mm256_add_epi16(AlphaBetaMain, Offset); //此处用Rr代替Main
		Offset = _mm256_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//存入内存
		_mm256_store_si256(CommonPlusIndex++, Rr);//加一次为256bit
		_mm256_store_si256(CommonMinusIndex++, Offset);
		//再存入公共存储区,以下Rr，AlphaBetaSub和Offset可使用
		//直接更新下一次
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		if ((i & 0x7) == 0x7)  //i%8=7时,则此时要计算得出一次OutPut并归一化
		{
			//归一化

			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);

			//生成输出信息，公共存储空间归0
			CommonPlusIndex = CommenPlusM;
			CommonMinusIndex = CommenMinusM;
			CountPlusIndex = (int*)CommenPlusM;
			CountMinusIndex = (int*)CommenMinusM;
			Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			for (int k = 0; k < 7; ++k)
			{
				Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
				Rr = _mm256_max_epi16(Offset, Rr);//Plus
				Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
				AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
			}
			//Plus减去Minus并重排
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
			Rr = _mm256_permute4x64_epi64(Rr, 0x8D);//11011000
			//减去La

			AlphaBetaSub = _mm256_loadu2_m128i(La2Index++, La1Index++);//此处每加一次为128bit
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
		}//if ((i & 0x111) == 0x111)
	}
	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		CountPlusIndex = (int*)CommenPlusM;
		CountMinusIndex = (int*)CommenMinusM;
		Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
		AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
		for (int k = 0; k < 7; ++k)
		{
			Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			Rr = _mm256_max_epi16(Offset, Rr);//Plus
			Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
		}
		//Plus减去Minus并重排
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
		Rr = _mm256_permute4x64_epi64(Rr, 0x8D);
		//减去La
		AlphaBetaSub = _mm256_loadu2_m128i(La2Index++, La1Index++);//此处每加一次为128bit
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
	}
	return 1;
}
inline bool MaxLogRscDecoderAVXIntDul(short*Ls1, short*Ls2, short*Lp1, short*Lp2, short*La1, short*La2, int Len, short*Output1, short*Output2, __m256i*BetaM, __m256i*CommenPlusM,__m256i*CommenMinusM)
{//注Common大小为8*__m256i
	__m256i AlphaBetaMain, AlphaBetaSub, Rr, Offset,RMask;
	__m256i OffsetTable = _mm256_set_epi32(112, 96, 80, 64, 48, 32, 16, 0);
	__m256i ShuffleAlphaPlus = _mm256_setr_epi32(1, 2, 5, 6, 0, 3, 4, 7);
	__m256i ShuffleAlphaMinus = _mm256_setr_epi32(0, 3, 4, 7, 1, 2, 5, 6);
	__m256i ShuffleBetaPlus = _mm256_setr_epi32(4, 0, 1, 5, 6, 2, 3, 7);
	__m256i ShuffleBetaMinus = _mm256_setr_epi32(0, 4, 5, 1, 2, 6, 7, 3);
	__m256i ShuffleTableBeforeStore = _mm256_setr_epi8(0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15);
	short L1_1;
	short L1_0;
	short L2_1;
	short L2_0;
	__m256i*BetaIndex = BetaM+Len;//指向最后一个元素的下一个位置
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);
	_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	/*先计算后向度量*/
	for (int i = Len-1; i>0; --i)
	{
		AlphaBetaSub = AlphaBetaMain;
		/*先得出分支度量*/
		if (i < Len - Mem)
		{
			L1_1 =(Ls1[i] + La1[i] + Lp1[i])>>1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + La2[i] + Lp2[i])>>1;
			L2_0 = L2_1 - Lp2[i];
		}
		else
		{
			L1_1 = (Ls1[i] + Lp1[i])>>1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 =( Ls2[i] + Lp2[i])>>1;
			L2_0 = L2_1 - Lp2[i];
		}
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1);
		//对于相加的部分
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//归一化
		
		if ((i & 0x7) == 0x7)//相当于每8次运算中有一次归一化
		{
			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);
		}
		
		_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	}
	/*再计算前向度量并计算似然比*/
	__m128i*Ls1Index = (__m128i*)Ls1;
	__m128i*Ls2Index = (__m128i*)Ls2;
	__m128i*La1Index = (__m128i*)La1;
	__m128i*La2Index = (__m128i*)La2;
	__m256i*CommonPlusIndex = CommenPlusM;
	__m256i*CommonMinusIndex = CommenMinusM;
	__m128i*Output1Index = (__m128i*)Output1;
	__m128i*Output2Index = (__m128i*)Output2;
	int*CountPlusIndex;
	int*CountMinusIndex;
	BetaIndex = BetaM;
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);

	for (int i = 0; i < Len; ++i)
	{
		AlphaBetaSub = AlphaBetaMain;
		//先得出分支度量
		if (i < Len - Mem)
		{
			L1_1 = (Ls1[i] + La1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + La2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		else
		{
			L1_1 = (Ls1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1);
		//先计算加
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleAlphaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleAlphaMinus);
		//再加上Alpha
		Offset = _mm256_load_si256(BetaIndex++);//加一次为256bit
		Rr = _mm256_add_epi16(AlphaBetaMain, Offset); //此处用Rr代替Main
		Offset = _mm256_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//存入内存
		_mm256_store_si256(CommonPlusIndex++, Rr);//加一次为256bit
		_mm256_store_si256(CommonMinusIndex++, Offset);
		//再存入公共存储区,以下Rr，AlphaBetaSub和Offset可使用
        //直接更新下一次
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		if ((i & 0x7) == 0x7)  //i%8=7时,则此时要计算得出一次OutPut并归一化
		{
			//归一化
			
			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);
			
			//生成输出信息，公共存储空间归0
			CommonPlusIndex = CommenPlusM; 
			CommonMinusIndex = CommenMinusM;
			CountPlusIndex = (int*)CommenPlusM;
			CountMinusIndex = (int*)CommenMinusM;
			Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			for (int k = 0; k < 7; ++k)
			{
				Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
				Rr= _mm256_max_epi16(Offset, Rr);//Plus
				Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
				AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
			}
			//Plus减去Minus并重排
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
			Rr = _mm256_permute4x64_epi64(Rr, 0x8D);//11011000
			//减去Ls+La

			Offset = _mm256_loadu2_m128i(La2Index++, La1Index++);//此处每加一次为128bit
			AlphaBetaSub = _mm256_loadu2_m128i(Ls2Index++, Ls1Index++);
			AlphaBetaSub = _mm256_add_epi16(AlphaBetaSub, Offset);
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
		}//if ((i & 0x111) == 0x111)
	}
	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		CountPlusIndex = (int*)CommenPlusM;
		CountMinusIndex = (int*)CommenMinusM;
		Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
		AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
		for (int k = 0; k < 7; ++k)
		{
			Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			Rr = _mm256_max_epi16(Offset, Rr);//Plus
			Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
		}
		//Plus减去Minus并重排
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
		Rr = _mm256_permute4x64_epi64(Rr, 0x8D);
		//减去Ls+La
		Offset = _mm256_loadu2_m128i(La2Index++, La1Index++);//此处每加一次为128bit
		AlphaBetaSub = _mm256_loadu2_m128i(Ls2Index++, Ls1Index++);
		AlphaBetaSub = _mm256_add_epi16(AlphaBetaSub, Offset);
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
	}
	return 1;
}
inline bool MaxLogRscDecoderAVXIntDulV2(short*Ls1, short*Ls2, short*Lp1, short*Lp2, short*La1, short*La2, int Len, short*Output1, short*Output2, __m256i*BetaM, __m256i*CommenPlusM, __m256i*CommenMinusM)
{//注Common大小为8*__m256i
	__m256i AlphaBetaMain, AlphaBetaSub, Rr, Offset, R1, R2, R3, R4, R5, R6, R7, R8;
	__m256i OffsetTable = _mm256_set_epi32(112, 96, 80, 64, 48, 32, 16, 0);
	__m256i ShuffleAlphaBetaPlus = _mm256_setr_epi32(4, 0, 1, 5, 6, 2, 3, 7);
	__m256i ShuffleAlphaBetaMinus = _mm256_setr_epi32(0, 4, 5, 1, 2, 6, 7, 3);
	__m256i ShuffleTableBeforeStore = _mm256_setr_epi8(0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15);
	short L1_1;
	short L1_0;
	short L2_1;
	short L2_0;
	__m256i*BetaIndex = BetaM + Len;//指向最后一个元素的下一个位置
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);
	_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	/*先计算后向度量*/
	for (int i = Len - 1; i>0; --i)
	{
		AlphaBetaSub = AlphaBetaMain;
		/*先得出分支度量*/
		if (i < Len - Mem)
		{
			L1_1 = (Ls1[i] + La1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + La2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		else
		{
			L1_1 = (Ls1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1);
		//对于相加的部分
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//归一化
		if ((i & 0x7) == 0x7)//相当于每8次运算中有一次归一化
		{
			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);
		}

		_mm256_store_si256(--BetaIndex, AlphaBetaMain);
	}
	/*再计算前向度量并计算似然比*/
	__m128i*Ls1Index = (__m128i*)Ls1;
	__m128i*Ls2Index = (__m128i*)Ls2;
	__m128i*La1Index = (__m128i*)La1;
	__m128i*La2Index = (__m128i*)La2;
	__m256i*CommonPlusIndex = CommenPlusM;
	__m256i*CommonMinusIndex = CommenMinusM;
	__m128i*Output1Index = (__m128i*)Output1;
	__m128i*Output2Index = (__m128i*)Output2;
	int*CountPlusIndex;
	int*CountMinusIndex;
	BetaIndex = BetaM;
	AlphaBetaMain = _mm256_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, 0, 0);
	ShuffleAlphaBetaPlus = _mm256_setr_epi32(1, 2, 5, 6, 0, 3, 4, 7);
	ShuffleAlphaBetaMinus = _mm256_setr_epi32(0, 3, 4, 7, 1, 2, 5, 6);
	for (int i = 0; i < Len; ++i)
	{
		AlphaBetaSub = AlphaBetaMain;
		//先得出分支度量
		if (i < Len - Mem)
		{
			L1_1 = (Ls1[i] + La1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + La2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		else
		{
			L1_1 = (Ls1[i] + Lp1[i]) >> 1;
			L1_0 = L1_1 - Lp1[i];
			L2_1 = (Ls2[i] + Lp2[i]) >> 1;
			L2_0 = L2_1 - Lp2[i];
		}
		Rr = _mm256_set_epi16(L1_1, L2_1, L1_1, L2_1, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_0, L2_0, L1_1, L2_1, L1_1, L2_1);
		//先计算加
		AlphaBetaMain = _mm256_add_epi16(AlphaBetaMain, Rr);//相加
		AlphaBetaMain = _mm256_permutevar8x32_epi32(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm256_sub_epi16(AlphaBetaSub, Rr);
		AlphaBetaSub = _mm256_permutevar8x32_epi32(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Alpha
		Offset = _mm256_load_si256(BetaIndex++);//加一次为256bit
		Rr = _mm256_add_epi16(AlphaBetaMain, Offset); //此处用Rr代替Main
		Offset = _mm256_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//存入寄存器
		_mm256_store_si256(CommonPlusIndex++, Rr);//加一次为256bit
		_mm256_store_si256(CommonMinusIndex++, Offset);
		//再存入公共存储区,以下Rr，AlphaBetaSub和Offset可使用
		//直接更新下一次
		AlphaBetaMain = _mm256_max_epi16(AlphaBetaMain, AlphaBetaSub);
		if ((i & 0x7) == 0x7)  //i%8=7时,则此时要计算得出一次OutPut并归一化
		{
			//归一化

			Offset = _mm256_broadcastd_epi32(_mm256_castsi256_si128(AlphaBetaMain));
			AlphaBetaMain = _mm256_sub_epi16(AlphaBetaMain, Offset);

			//生成输出信息，公共存储空间归0
			CommonPlusIndex = CommenPlusM;
			CommonMinusIndex = CommenMinusM;
			CountPlusIndex = (int*)CommenPlusM;
			CountMinusIndex = (int*)CommenMinusM;
			Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			for (int k = 0; k < 7; ++k)
			{
				Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
				Rr = _mm256_max_epi16(Offset, Rr);//Plus
				Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
				AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
			}
			//Plus减去Minus并重排
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
			Rr = _mm256_permute4x64_epi64(Rr, 0x8D);//11011000
			//减去Ls+La
			Offset = _mm256_loadu2_m128i(La2Index++, La1Index++);//此处每加一次为128bit
			AlphaBetaSub = _mm256_loadu2_m128i(Ls2Index++, Ls1Index++);
			AlphaBetaSub = _mm256_add_epi16(AlphaBetaSub, Offset);
			Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
			_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
		}//if ((i & 0x111) == 0x111)
	}
	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		CountPlusIndex = (int*)CommenPlusM;
		CountMinusIndex = (int*)CommenMinusM;
		Rr = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
		AlphaBetaSub = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
		for (int k = 0; k < 7; ++k)
		{
			Offset = _mm256_i32gather_epi32(CountPlusIndex++, OffsetTable, 2);//加一次为32bit
			Rr = _mm256_max_epi16(Offset, Rr);//Plus
			Offset = _mm256_i32gather_epi32(CountMinusIndex++, OffsetTable, 2);
			AlphaBetaSub = _mm256_max_epi16(Offset, AlphaBetaSub);//Minus
		}
		//Plus减去Minus并重排
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		Rr = _mm256_shuffle_epi8(Rr, ShuffleTableBeforeStore);
		Rr = _mm256_permute4x64_epi64(Rr, 0x8D);
		//减去Ls+La
		Offset = _mm256_loadu2_m128i(La2Index++, La1Index++);//此处每加一次为128bit
		AlphaBetaSub = _mm256_loadu2_m128i(Ls2Index++, Ls1Index++);
		AlphaBetaSub = _mm256_add_epi16(AlphaBetaSub, Offset);
		Rr = _mm256_sub_epi16(Rr, AlphaBetaSub);
		_mm256_storeu2_m128i(Output2Index++, Output1Index++, Rr);
	}
	return 1;
}
double TurboDecoderAVX(float*Ls, float*Lp0, float*Lp1, float*LIns, int *Table, float*Output, int Len, int InterNum)//Len为加上Mem的
{
	float*AlphaM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*BetaPlusM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*BetaMinusM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*Le1 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	float*Le0 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	float*La1 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	memset(Output, 0, Len*sizeof(float));
	//All System Go
	LARGE_INTEGER Start, End, Freq;
	QueryPerformanceCounter(&Start);
	QueryPerformanceFrequency(&Freq);
	for (int i = 1; i <= InterNum; ++i)
	{
		MaxLogRscDecoderAVXFloat(Ls, Lp0, Output, Len, Le0, AlphaM, BetaPlusM, BetaMinusM);
		//MaxLogSubDecoderAno(Ls, Lp0, Output, Len, Le0);
		for (int k = 0; k < Len - Mem; ++k)
			La1[k] = Le0[Table[k]];
		for (int k = Len - Mem; k < Len; ++k)
			La1[k] = Le0[k];
		MaxLogRscDecoderAVXFloat(LIns, Lp1, La1, Len, Le1, AlphaM, BetaPlusM, BetaMinusM);
		//MaxLogSubDecoderAno(LIns, Lp1, La1, Len, Le1);
		for (int k = 0; k < Len - Mem; ++k)
			Output[Table[k]] = Le1[k];
		for (int k = Len - Mem; k < Len; ++k)
			Output[k] = Le1[k];
	}
	QueryPerformanceCounter(&End);
	//
	_mm_free(AlphaM);
	_mm_free(BetaPlusM);
	_mm_free(BetaMinusM);
	_mm_free(Le1);
	_mm_free(La1);
	_mm_free(Le0);
	return (double)(End.QuadPart - Start.QuadPart) / (double)Freq.QuadPart;
}
double TurboDecoderAVXIntDul(short*Ls0, short*Ls1, short*Lp0_1, short*Lp1_1, short*Lp0_2, short*Lp1_2, short*LIns0, short*LIns1, int*Table, short*Output0, short*Output1, int Len, int InterNum)
{//译码器号_译码器内支路号
	__m256i*AlphaM = (__m256i*)_mm_malloc(Len*sizeof(__m256), sizeof(__m256));
	__m256i*BetaPlusM = (__m256i*)_mm_malloc(8*sizeof(__m256), sizeof(__m256));
	__m256i*BetaMinusM = (__m256i*)_mm_malloc(8*sizeof(__m256), sizeof(__m256));
	short*Le1_0 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	short*Le0_0 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	short*La0 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	short*Le1_1 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	short*Le0_1 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	short*La1 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	memset(Output0, 0, Len*sizeof(short));
	memset(Output1, 0, Len*sizeof(short));
	//All System Go
	LARGE_INTEGER Start, End, Freq;
	QueryPerformanceCounter(&Start);
	QueryPerformanceFrequency(&Freq);
	for (int i = 1; i <= InterNum; ++i)
	{
		MaxLogRscDecoderAVXIntDulV2(Ls0,Ls1,Lp0_1,Lp1_1,Output0,Output1, Len, Le0_0,Le1_0,AlphaM, BetaPlusM, BetaMinusM);
		//MaxLogSubDecoderSInt(Ls0, Lp0_1, Output0, Len, Le0_0);
		//MaxLogSubDecoderSInt(Ls1, Lp1_1, Output1, Len, Le1_0);
		for (int k = 0; k < Len - Mem; ++k)
		{
			La0[k] = Le0_0[Table[k]];
			La1[k] = Le1_0[Table[k]];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			La0[k] = Le0_0[k];
			La1[k] = Le1_0[k];
		}
		MaxLogRscDecoderAVXIntDulV2(LIns0, LIns1, Lp0_2, Lp1_2, La0, La1,Len, Le0_1,Le1_1, AlphaM, BetaPlusM, BetaMinusM);
		//MaxLogSubDecoderSInt(LIns0, Lp0_2, La0, Len, Le0_1);
		//MaxLogSubDecoderSInt(LIns1, Lp1_2, La1, Len, Le1_1);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output0[Table[k]] = Le0_1[k];
			Output1[Table[k]] = Le1_1[k];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output0[k] = Le0_1[k];
			Output1[k] = Le1_1[k];
		}
	}
	QueryPerformanceCounter(&End);
	//
	_mm_free(AlphaM);
	_mm_free(BetaPlusM);
	_mm_free(BetaMinusM);
	_mm_free(La1);
	_mm_free(La0);
	_mm_free(Le1_0);
	_mm_free(Le0_0);
	_mm_free(Le1_1);
	_mm_free(Le0_1);
	return (double)(End.QuadPart - Start.QuadPart) / (double)Freq.QuadPart;
}
double TurboDecoderSSEIntV2(short*Ls, short*Lp0, short*Lp1, short*LIns, int*Table, short*Output, int Len, int InterNum)
{
	//译码器号_译码器内支路号
	__m128i*AlphaM = (__m128i*)_mm_malloc(Len*sizeof(__m128i), sizeof(__m128i));
	short*La = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m128i));
	memset(Output, 0, Len*sizeof(short));
	//All System Go
	LARGE_INTEGER Start, End, Freq;
	QueryPerformanceCounter(&Start);
	QueryPerformanceFrequency(&Freq);
	for (int i = 1; i <= InterNum; ++i)
	{
		MaxLogRscDecoderSSEIntD0(Ls, Lp0, Output, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[k] = La[Table[k]];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = LIns[k];
		}
		MaxLogRscDecoderSSEIntD1(Output, Lp1, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[Table[k]] = La[k];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = La[k];
		}
	}
	QueryPerformanceCounter(&End);
	//
	_mm_free(AlphaM);
	_mm_free(La);
	return (double)(End.QuadPart - Start.QuadPart) / (double)Freq.QuadPart;
}
double TurboDecoderAVXIntDulV2(short*Ls0, short*Ls1, short*Lp0_1, short*Lp1_1, short*Lp0_2, short*Lp1_2, short*LIns0, short*LIns1, int*Table, short*Output0, short*Output1, int Len, int InterNum)
{
	//译码器号_译码器内支路号
	__m256i*AlphaM = (__m256i*)_mm_malloc(Len*sizeof(__m256), sizeof(__m256));
	__m256i*BetaPlusM = (__m256i*)_mm_malloc(8 * sizeof(__m256), sizeof(__m256));
	__m256i*BetaMinusM = (__m256i*)_mm_malloc(8 * sizeof(__m256), sizeof(__m256));
	short*La0 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	short*La1 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
	memset(Output0, 0, Len*sizeof(short));
	memset(Output1, 0, Len*sizeof(short));
	//All System Go
	LARGE_INTEGER Start, End, Freq;
	QueryPerformanceCounter(&Start);
	QueryPerformanceFrequency(&Freq);
	for (int i = 1; i <= InterNum; ++i)
	{
		MaxLogRscDecoderAVXIntDulD0(Ls0, Ls1, Lp0_1, Lp1_1, Output0, Output1, Len, La0, La1, AlphaM, BetaPlusM, BetaMinusM);
		
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output0[k] = La0[Table[k]];
			Output1[k] = La1[Table[k]];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output0[k] = LIns0[k];
			Output1[k] = LIns1[k];
		}
		
		MaxLogRscDecoderAVXIntDulD1(Output0, Output1, Lp0_2, Lp1_2, Len, La0, La1, AlphaM, BetaPlusM, BetaMinusM);
		
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output0[Table[k]] = La0[k];
			Output1[Table[k]] = La1[k];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output0[k] = La0[k];
			Output1[k] = La1[k];
		}
		
	}
	QueryPerformanceCounter(&End);
	_mm_free(AlphaM);
	_mm_free(BetaPlusM);
	_mm_free(BetaMinusM);
	_mm_free(La1);
	_mm_free(La0);
	return (double)(End.QuadPart - Start.QuadPart) / (double)Freq.QuadPart;
}
//Reference Function
void MaxLogSubDecoderOri(float*Ls, float*Lp, float*Le, int Len, float*Output) //Len长度不含4位尾比特，不记寄存器个数(Mem)
{
	//定义alpha，beta
	int StateNum = 8;
	int OffsetTmp_p = -INF;
	int OffsetTmp_n = -INF;
	float *alpha = (float *)malloc(StateNum*(Len)*sizeof(float));
	float *beta = (float *)malloc(StateNum*(Len)*sizeof(float));
	//定义四种取值不同的分支度量
	double *r0 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r_0 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r1 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r_1 = (double*)malloc((Len + Mem)*sizeof(double));
	//初始化alpha，beta
	alpha[0] = 0;
	alpha[1] = -INF;
	alpha[2] = -INF;
	alpha[3] = -INF;
	alpha[4] = -INF;
	alpha[5] = -INF;
	alpha[6] = -INF;
	alpha[7] = -INF;
	beta[(Len - 1)*StateNum + 0] = 0;
	beta[(Len - 1)*StateNum + 1] = -INF;
	beta[(Len - 1)*StateNum + 2] = -INF;
	beta[(Len - 1)*StateNum + 3] = -INF;
	beta[(Len - 1)*StateNum + 4] = -INF;
	beta[(Len - 1)*StateNum + 5] = -INF;
	beta[(Len - 1)*StateNum + 6] = -INF;
	beta[(Len - 1)*StateNum + 7] = -INF;
	//初始化四个分支度量
	for (int i = 0; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Le[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Le[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	for (int i = Len - Mem; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	//计算alpha,
	OffsetTmp_n = StateNum;
	OffsetTmp_p = 0;
	for (int i = 0; i < Len - 1; ++i)
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
	OffsetTmp_n = (Len - 1)* StateNum;
	OffsetTmp_p = (Len - 2)*StateNum;
	for (int i = Len - 1; i > 0; --i)
	{
		beta[OffsetTmp_p] = max(beta[OffsetTmp_n] + r_0[i], beta[OffsetTmp_n + 4] + r0[i]);
		beta[OffsetTmp_p + 1] = max(beta[OffsetTmp_n] + r0[i], beta[OffsetTmp_n + 4] + r_0[i]);
		beta[OffsetTmp_p + 2] = max(beta[OffsetTmp_n + 1] + r1[i], beta[OffsetTmp_n + 5] + r_1[i]);
		beta[OffsetTmp_p + 3] = max(beta[OffsetTmp_n + 5] + r1[i], beta[OffsetTmp_n + 1] + r_1[i]);
		beta[OffsetTmp_p + 4] = max(beta[OffsetTmp_n + 6] + r1[i], beta[OffsetTmp_n + 2] + r_1[i]);
		beta[OffsetTmp_p + 5] = max(beta[OffsetTmp_n + 2] + r1[i], beta[OffsetTmp_n + 6] + r_1[i]);
		beta[OffsetTmp_p + 6] = max(beta[OffsetTmp_n + 3] + r0[i], beta[OffsetTmp_n + 7] + r_0[i]);
		beta[OffsetTmp_p + 7] = max(beta[OffsetTmp_n + 3] + r_0[i], beta[OffsetTmp_n + 7] + r0[i]);
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
	for (int i = 0; i < Len; ++i)
	{

		Output[i] =
			Mmax(
			Mmax(                                     //之前为i+1错误
			Mmax(alpha[i*StateNum + 0] + r0[i] + beta[i*StateNum + 4], alpha[i*StateNum + 1] + r0[i] + beta[i*StateNum + 0]),
			Mmax(alpha[i*StateNum + 2] + r1[i] + beta[i*StateNum + 1], alpha[i*StateNum + 3] + r1[i] + beta[i*StateNum + 5])
			),
			Mmax(
			Mmax(alpha[i*StateNum + 4] + r1[i] + beta[i*StateNum + 6], alpha[i*StateNum + 5] + r1[i] + beta[i*StateNum + 2]),
			Mmax(alpha[i*StateNum + 6] + r0[i] + beta[i*StateNum + 3], alpha[i*StateNum + 7] + r0[i] + beta[i*StateNum + 7])
			)
			)
			-
			Mmax(
			Mmax(
			Mmax(alpha[i*StateNum + 0] + r_0[i] + beta[i*StateNum + 0], alpha[i*StateNum + 1] + r_0[i] + beta[i*StateNum + 4]),
			Mmax(alpha[i*StateNum + 2] + r_1[i] + beta[i*StateNum + 5], alpha[i*StateNum + 3] + r_1[i] + beta[i*StateNum + 1])
			),
			Mmax(
			Mmax(alpha[i*StateNum + 4] + r_1[i] + beta[i*StateNum + 2], alpha[i*StateNum + 5] + r_1[i] + beta[i*StateNum + 6]),
			Mmax(alpha[i*StateNum + 6] + r_0[i] + beta[i*StateNum + 7], alpha[i*StateNum + 7] + r_0[i] + beta[i*StateNum + 3])
			)
			)
			- Ls[i] - Le[i];
	}
	free(alpha);
	free(beta);
	free(r0);
	free(r1);
	free(r_1);
	free(r_0);
}
void MaxLogSubDecoderSIntAno(short*Ls, short*Lp, short*Le, int Len, short*Output) //Len长度不含4位尾比特，不记寄存器个数(Mem)
{
	//定义alpha，beta
	int StateNum = 8;
	int OffsetTmp_p = -INF;
	int OffsetTmp_n = -INF;
	short *alphaplus = (short *)malloc(StateNum*(Len)*sizeof(short));
	short *alpha = (short *)malloc(StateNum*(Len)*sizeof(short));
	short *alphaminus = (short *)malloc(StateNum*(Len)*sizeof(short));
	short *beta = (short *)malloc(StateNum*(Len)*sizeof(short));
	//定义四种取值不同的分支度量
	short *r0 = (short*)malloc((Len)*sizeof(short));
	short *r_0 = (short*)malloc((Len)*sizeof(short));
	short *r1 = (short*)malloc((Len)*sizeof(short));
	short *r_1 = (short*)malloc((Len)*sizeof(short));
	//初始化alpha，beta
	alpha[0] = 0;
	alpha[1] = -INF;
	alpha[2] = -INF;
	alpha[3] = -INF;
	alpha[4] = -INF;
	alpha[5] = -INF;
	alpha[6] = -INF;
	alpha[7] = -INF;
	beta[(Len - 1)*StateNum + 0] = 0;
	beta[(Len - 1)*StateNum + 1] = -INF;
	beta[(Len - 1)*StateNum + 2] = -INF;
	beta[(Len - 1)*StateNum + 3] = -INF;
	beta[(Len - 1)*StateNum + 4] = -INF;
	beta[(Len - 1)*StateNum + 5] = -INF;
	beta[(Len - 1)*StateNum + 6] = -INF;
	beta[(Len - 1)*StateNum + 7] = -INF;
	//初始化四个分支度量
	for (int i = 0; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Le[i] + Ls[i]);
		r_0[i] = -0.5*(Lp[i] + Le[i] + Ls[i]);
		r1[i] = 0.5*(-Lp[i] + Le[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	for (int i = Len - Mem; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Ls[i]);
		r_0[i] = -0.5*(Lp[i] + Ls[i]);
		r1[i] = 0.5*(-Lp[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	//计算beta
	OffsetTmp_n = (Len - 1)* StateNum;
	OffsetTmp_p = (Len - 2)*StateNum;
	for (int i = Len - 1; i > 0; --i)
	{
		beta[OffsetTmp_p] = max(beta[OffsetTmp_n] + r_0[i], beta[OffsetTmp_n + 4] + r0[i]);
		beta[OffsetTmp_p + 1] = max(beta[OffsetTmp_n] + r0[i], beta[OffsetTmp_n + 4] + r_0[i]);
		beta[OffsetTmp_p + 2] = max(beta[OffsetTmp_n + 1] + r1[i], beta[OffsetTmp_n + 5] + r_1[i]);
		beta[OffsetTmp_p + 3] = max(beta[OffsetTmp_n + 5] + r1[i], beta[OffsetTmp_n + 1] + r_1[i]);
		beta[OffsetTmp_p + 4] = max(beta[OffsetTmp_n + 6] + r1[i], beta[OffsetTmp_n + 2] + r_1[i]);
		beta[OffsetTmp_p + 5] = max(beta[OffsetTmp_n + 2] + r1[i], beta[OffsetTmp_n + 6] + r_1[i]);
		beta[OffsetTmp_p + 6] = max(beta[OffsetTmp_n + 3] + r0[i], beta[OffsetTmp_n + 7] + r_0[i]);
		beta[OffsetTmp_p + 7] = max(beta[OffsetTmp_n + 3] + r_0[i], beta[OffsetTmp_n + 7] + r0[i]);
		OffsetTmp_n -= StateNum;
		OffsetTmp_p -= StateNum;
	}
	//计算alpha,
	OffsetTmp_n = StateNum;
	OffsetTmp_p = 0;
	for (int i = 0; i < Len - 1; ++i)
	{
		alphaplus[OffsetTmp_p] = alpha[OffsetTmp_p + 1] + r0[i]+beta[OffsetTmp_p];
		alphaplus[OffsetTmp_p + 1] = alpha[OffsetTmp_p + 2] + r1[i] + beta[OffsetTmp_p+1];
		alphaplus[OffsetTmp_p + 2] = alpha[OffsetTmp_p + 5] + r1[i] + beta[OffsetTmp_p+2];
		alphaplus[OffsetTmp_p + 3] = alpha[OffsetTmp_p + 6] + r0[i] + beta[OffsetTmp_p+3];
		alphaplus[OffsetTmp_p + 4] = alpha[OffsetTmp_p + 0] + r0[i] + beta[OffsetTmp_p+4];
		alphaplus[OffsetTmp_p + 5] = alpha[OffsetTmp_p + 3] + r1[i] + beta[OffsetTmp_p+5];
		alphaplus[OffsetTmp_p + 6] = alpha[OffsetTmp_p + 4] + r1[i] + beta[OffsetTmp_p+6];
		alphaplus[OffsetTmp_p + 7] = alpha[OffsetTmp_p + 7] + r0[i] + beta[OffsetTmp_p+7];

		alphaminus[OffsetTmp_p] = alpha[OffsetTmp_p] + r_0[i] + beta[OffsetTmp_p];
		alphaminus[OffsetTmp_p + 1] = alpha[OffsetTmp_p + 3] + r_1[i] + beta[OffsetTmp_p+1];
		alphaminus[OffsetTmp_p + 2] = alpha[OffsetTmp_p + 4] + r_1[i] + beta[OffsetTmp_p+2];
		alphaminus[OffsetTmp_p + 3] = alpha[OffsetTmp_p + 7] + r_0[i] + beta[OffsetTmp_p+3];
		alphaminus[OffsetTmp_p + 4] = alpha[OffsetTmp_p + 1] + r_0[i] + beta[OffsetTmp_p+4];
		alphaminus[OffsetTmp_p + 5] = alpha[OffsetTmp_p + 2] + r_1[i] + beta[OffsetTmp_p+5];
		alphaminus[OffsetTmp_p + 6] = alpha[OffsetTmp_p + 5] + r_1[i] + beta[OffsetTmp_p+6];
		alphaminus[OffsetTmp_p + 7] = alpha[OffsetTmp_p + 6] + r_0[i] + beta[OffsetTmp_p+7];

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
	alphaplus[OffsetTmp_p] = alpha[OffsetTmp_p + 1] + r0[Len-1] + beta[OffsetTmp_p];
	alphaplus[OffsetTmp_p + 1] = alpha[OffsetTmp_p + 2] + r1[Len - 1] + beta[OffsetTmp_p + 1];
	alphaplus[OffsetTmp_p + 2] = alpha[OffsetTmp_p + 5] + r1[Len - 1] + beta[OffsetTmp_p + 2];
	alphaplus[OffsetTmp_p + 3] = alpha[OffsetTmp_p + 6] + r0[Len - 1] + beta[OffsetTmp_p + 3];
	alphaplus[OffsetTmp_p + 4] = alpha[OffsetTmp_p + 0] + r0[Len - 1] + beta[OffsetTmp_p + 4];
	alphaplus[OffsetTmp_p + 5] = alpha[OffsetTmp_p + 3] + r1[Len - 1] + beta[OffsetTmp_p + 5];
	alphaplus[OffsetTmp_p + 6] = alpha[OffsetTmp_p + 4] + r1[Len - 1] + beta[OffsetTmp_p + 6];
	alphaplus[OffsetTmp_p + 7] = alpha[OffsetTmp_p + 7] + r0[Len - 1] + beta[OffsetTmp_p + 7];

	alphaminus[OffsetTmp_p] = alpha[OffsetTmp_p] + r_0[Len - 1] + beta[OffsetTmp_p];
	alphaminus[OffsetTmp_p + 1] = alpha[OffsetTmp_p + 3] + r_1[Len - 1] + beta[OffsetTmp_p + 1];
	alphaminus[OffsetTmp_p + 2] = alpha[OffsetTmp_p + 4] + r_1[Len - 1] + beta[OffsetTmp_p + 2];
	alphaminus[OffsetTmp_p + 3] = alpha[OffsetTmp_p + 7] + r_0[Len - 1] + beta[OffsetTmp_p + 3];
	alphaminus[OffsetTmp_p + 4] = alpha[OffsetTmp_p + 1] + r_0[Len - 1] + beta[OffsetTmp_p + 4];
	alphaminus[OffsetTmp_p + 5] = alpha[OffsetTmp_p + 2] + r_1[Len - 1] + beta[OffsetTmp_p + 5];
	alphaminus[OffsetTmp_p + 6] = alpha[OffsetTmp_p + 5] + r_1[Len - 1] + beta[OffsetTmp_p + 6];
	alphaminus[OffsetTmp_p + 7] = alpha[OffsetTmp_p + 6] + r_0[Len - 1] + beta[OffsetTmp_p + 7];
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
	for (int i = 0; i < Len; ++i)
	{

		Output[i] =
			Mmax(
			Mmax(
			Mmax(alphaplus[i*StateNum + 0], alphaplus[i*StateNum + 1]),
			Mmax(alphaplus[i*StateNum + 2], alphaplus[i*StateNum + 3])
			),
			Mmax(
			Mmax(alphaplus[i*StateNum + 4], alphaplus[i*StateNum + 5]),
			Mmax(alphaplus[i*StateNum + 6], alphaplus[i*StateNum + 7])
			)
			)
			-
			Mmax(
			Mmax(                                     //之前为i+1错误
			Mmax(alphaminus[i*StateNum + 0], alphaminus[i*StateNum + 1]),
			Mmax(alphaminus[i*StateNum + 2], alphaminus[i*StateNum + 3])
			),
			Mmax(
			Mmax(alphaminus[i*StateNum + 4], alphaminus[i*StateNum + 5]),
			Mmax(alphaminus[i*StateNum + 6], alphaminus[i*StateNum + 7])
			)
			)
			- Ls[i] - Le[i];
	}
	free(alpha);
	free(beta);
	free(r0);
	free(r1);
	free(r_1);
	free(r_0);
}
void MaxLogSubDecoderSIntNor(short*Ls, short*Lp, short*Le, int Len, short*Output) //Len长度不含4位尾比特，不记寄存器个数(Mem)
{
	//定义alpha，beta
	int StateNum = 8;
	int OffsetTmp_p = -INF;
	int OffsetTmp_n = -INF;
	short *alpha = (short *)malloc(StateNum*(Len)*sizeof(short));
	short *beta = (short *)malloc(StateNum*(Len)*sizeof(short));
	//定义四种取值不同的分支度量
	short *r0 = (short*)malloc((Len + Mem+1)*sizeof(short));
	short *r_0 = (short*)malloc((Len + Mem+1)*sizeof(short));
	short *r1 = (short*)malloc((Len + Mem+1)*sizeof(short));
	short *r_1 = (short*)malloc((Len + Mem+1)*sizeof(short));
	//初始化alpha，beta
	alpha[0] = 0;
	alpha[1] = -INF;
	alpha[2] = -INF;
	alpha[3] = -INF;
	alpha[4] = -INF;
	alpha[5] = -INF;
	alpha[6] = -INF;
	alpha[7] = -INF;
	beta[(Len - 1)*StateNum + 0] = 0;
	beta[(Len - 1)*StateNum + 1] = -INF;
	beta[(Len - 1)*StateNum + 2] = -INF;
	beta[(Len - 1)*StateNum + 3] = -INF;
	beta[(Len - 1)*StateNum + 4] = -INF;
	beta[(Len - 1)*StateNum + 5] = -INF;
	beta[(Len - 1)*StateNum + 6] = -INF;
	beta[(Len - 1)*StateNum + 7] = -INF;
	//初始化四个分支度量
	for (int i = 0; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Le[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Le[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	for (int i = Len - Mem; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Ls[i]);
		r_1[i] = -r1[i];
	}
	//计算alpha,
	OffsetTmp_n = StateNum;
	OffsetTmp_p = 0;
	int k = 0;
	for (int i = 0; i < Len - 1; ++i)
	{
		alpha[OffsetTmp_n] = max(alpha[OffsetTmp_p] + r_0[i], alpha[OffsetTmp_p + 1] + r0[i]);
		alpha[OffsetTmp_n + 1] = max(alpha[OffsetTmp_p + 2] + r1[i], alpha[OffsetTmp_p + 3] + r_1[i]);
		alpha[OffsetTmp_n + 2] = max(alpha[OffsetTmp_p + 4] + r_1[i], alpha[OffsetTmp_p + 5] + r1[i]);
		alpha[OffsetTmp_n + 3] = max(alpha[OffsetTmp_p + 6] + r0[i], alpha[OffsetTmp_p + 7] + r_0[i]);
		alpha[OffsetTmp_n + 4] = max(alpha[OffsetTmp_p + 0] + r0[i], alpha[OffsetTmp_p + 1] + r_0[i]);
		alpha[OffsetTmp_n + 5] = max(alpha[OffsetTmp_p + 3] + r1[i], alpha[OffsetTmp_p + 2] + r_1[i]);
		alpha[OffsetTmp_n + 6] = max(alpha[OffsetTmp_p + 4] + r1[i], alpha[OffsetTmp_p + 5] + r_1[i]);
		alpha[OffsetTmp_n + 7] = max(alpha[OffsetTmp_p + 7] + r0[i], alpha[OffsetTmp_p + 6] + r_0[i]);
		if ((OffsetTmp_n / StateNum) % 8 == 0)
		{
			alpha[OffsetTmp_n] -= alpha[OffsetTmp_n];
			alpha[OffsetTmp_n + 1] -= alpha[OffsetTmp_n];
			alpha[OffsetTmp_n + 2] -= alpha[OffsetTmp_n];
			alpha[OffsetTmp_n + 3] -= alpha[OffsetTmp_n];
			alpha[OffsetTmp_n + 4] -= alpha[OffsetTmp_n];
			alpha[OffsetTmp_n + 5] -= alpha[OffsetTmp_n];
			alpha[OffsetTmp_n + 6] -= alpha[OffsetTmp_n];
			alpha[OffsetTmp_n + 7] -= alpha[OffsetTmp_n];
		}
		k++;
		OffsetTmp_n += StateNum;
		OffsetTmp_p += StateNum;
	}
	//计算beta
	
	OffsetTmp_n = (Len - 1)* StateNum;
	OffsetTmp_p = (Len - 2)*StateNum;
	k = 0;;
	for (int i = Len - 1; i > 0; --i)
	{
		beta[OffsetTmp_p] = max(beta[OffsetTmp_n] + r_0[i], beta[OffsetTmp_n + 4] + r0[i]);
		beta[OffsetTmp_p + 1] = max(beta[OffsetTmp_n] + r0[i], beta[OffsetTmp_n + 4] + r_0[i]);
		beta[OffsetTmp_p + 2] = max(beta[OffsetTmp_n + 1] + r1[i], beta[OffsetTmp_n + 5] + r_1[i]);
		beta[OffsetTmp_p + 3] = max(beta[OffsetTmp_n + 5] + r1[i], beta[OffsetTmp_n + 1] + r_1[i]);
		beta[OffsetTmp_p + 4] = max(beta[OffsetTmp_n + 6] + r1[i], beta[OffsetTmp_n + 2] + r_1[i]);
		beta[OffsetTmp_p + 5] = max(beta[OffsetTmp_n + 2] + r1[i], beta[OffsetTmp_n + 6] + r_1[i]);
		beta[OffsetTmp_p + 6] = max(beta[OffsetTmp_n + 3] + r0[i], beta[OffsetTmp_n + 7] + r_0[i]);
		beta[OffsetTmp_p + 7] = max(beta[OffsetTmp_n + 3] + r_0[i], beta[OffsetTmp_n + 7] + r0[i]);
		if ((OffsetTmp_p / StateNum) % 8 == 0)
		{
			beta[OffsetTmp_p] -= beta[OffsetTmp_p];
			beta[OffsetTmp_p + 1] -= beta[OffsetTmp_p];
			beta[OffsetTmp_p + 2] -= beta[OffsetTmp_p];
			beta[OffsetTmp_p + 3] -= beta[OffsetTmp_p];
			beta[OffsetTmp_p + 4] -= beta[OffsetTmp_p];
			beta[OffsetTmp_p + 5] -= beta[OffsetTmp_p];
			beta[OffsetTmp_p + 6] -= beta[OffsetTmp_p];
			beta[OffsetTmp_p + 7] -= beta[OffsetTmp_p];
		}
		k++;
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
	for (int i = 0; i < Len; ++i)
	{

		Output[i] =
			Mmax(
			Mmax(                                     //之前为i+1错误
			Mmax(alpha[i*StateNum + 0] + r0[i] + beta[i*StateNum + 4], alpha[i*StateNum + 1] + r0[i] + beta[i*StateNum + 0]),
			Mmax(alpha[i*StateNum + 2] + r1[i] + beta[i*StateNum + 1], alpha[i*StateNum + 3] + r1[i] + beta[i*StateNum + 5])
			),
			Mmax(
			Mmax(alpha[i*StateNum + 4] + r1[i] + beta[i*StateNum + 6], alpha[i*StateNum + 5] + r1[i] + beta[i*StateNum + 2]),
			Mmax(alpha[i*StateNum + 6] + r0[i] + beta[i*StateNum + 3], alpha[i*StateNum + 7] + r0[i] + beta[i*StateNum + 7])
			)
			)
			-
			Mmax(
			Mmax(
			Mmax(alpha[i*StateNum + 0] + r_0[i] + beta[i*StateNum + 0], alpha[i*StateNum + 1] + r_0[i] + beta[i*StateNum + 4]),
			Mmax(alpha[i*StateNum + 2] + r_1[i] + beta[i*StateNum + 5], alpha[i*StateNum + 3] + r_1[i] + beta[i*StateNum + 1])
			),
			Mmax(
			Mmax(alpha[i*StateNum + 4] + r_1[i] + beta[i*StateNum + 2], alpha[i*StateNum + 5] + r_1[i] + beta[i*StateNum + 6]),
			Mmax(alpha[i*StateNum + 6] + r_0[i] + beta[i*StateNum + 7], alpha[i*StateNum + 7] + r_0[i] + beta[i*StateNum + 3])
			)
			)
			- Ls[i] - Le[i];
	}
	free(alpha);
	free(beta);
	free(r0);
	free(r1);
	free(r_1);
	free(r_0);
}
void MaxLogSubDecoderAno(float*Ls, float*Lp, float*Le, int Len, float*Output) //Len长度不含4位尾比特，不记寄存器个数(Mem)
{
	//定义alpha，beta
	int StateNum = 8;
	int OffsetTmp_p = -INF;
	int OffsetTmp_n = -INF;
	float *alpha = (float *)malloc(StateNum*(Len)*sizeof(float));
	float *beta = (float *)malloc(StateNum*(Len)*sizeof(float));
	float *betaplus = (float *)malloc(StateNum*(Len)*sizeof(float));
	float *betaminus = (float *)malloc(StateNum*(Len)*sizeof(float));
	//定义四种取值不同的分支度量
	double *r0 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r_0 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r1 = (double*)malloc((Len + Mem)*sizeof(double));
	double *r_1 = (double*)malloc((Len + Mem)*sizeof(double));
	//初始化alpha，beta
	alpha[0] = 0;
	alpha[1] = -INF;
	alpha[2] = -INF;
	alpha[3] = -INF;
	alpha[4] = -INF;
	alpha[5] = -INF;
	alpha[6] = -INF;
	alpha[7] = -INF;
	beta[(Len - 1)*StateNum + 0] = 0;
	beta[(Len - 1)*StateNum + 1] = -INF;
	beta[(Len - 1)*StateNum + 2] = -INF;
	beta[(Len - 1)*StateNum + 3] = -INF;
	beta[(Len - 1)*StateNum + 4] = -INF;
	beta[(Len - 1)*StateNum + 5] = -INF;
	beta[(Len - 1)*StateNum + 6] = -INF;
	beta[(Len - 1)*StateNum + 7] = -INF;
	//初始化四个分支度量
	for (int i = 0; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Le[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Le[i] + Ls[i]);
		r_1[i] = -r1[i];
	}

	for (int i = Len - Mem; i < Len; ++i)
	{
		r0[i] = 0.5*(Lp[i] + Ls[i]);
		r_0[i] = -r0[i];
		r1[i] = 0.5*(-Lp[i] + Ls[i]);
		r_1[i] = -r1[i];
	}

	//计算alpha,
	OffsetTmp_n = StateNum;
	OffsetTmp_p = 0;
	for (int i = 0; i < Len - 1; ++i)
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
	OffsetTmp_n = (Len - 1)* StateNum;
	OffsetTmp_p = (Len - 2)*StateNum;
	for (int i = Len - 1; i > 0; --i)
	{

		betaplus[OffsetTmp_n] = beta[OffsetTmp_n + 4] + r0[i] + alpha[OffsetTmp_n];
		betaplus[OffsetTmp_n + 1] = beta[OffsetTmp_n] + r0[i] + alpha[OffsetTmp_n + 1];
		betaplus[OffsetTmp_n + 2] = beta[OffsetTmp_n + 1] + r1[i] + alpha[OffsetTmp_n + 2];
		betaplus[OffsetTmp_n + 3] = beta[OffsetTmp_n + 5] + r1[i] + alpha[OffsetTmp_n + 3];
		betaplus[OffsetTmp_n + 4] = beta[OffsetTmp_n + 6] + r1[i] + alpha[OffsetTmp_n + 4];
		betaplus[OffsetTmp_n + 5] = beta[OffsetTmp_n + 2] + r1[i] + alpha[OffsetTmp_n + 5];
		betaplus[OffsetTmp_n + 6] = beta[OffsetTmp_n + 3] + r0[i] + alpha[OffsetTmp_n + 6];
		betaplus[OffsetTmp_n + 7] = beta[OffsetTmp_n + 7] + r0[i] + alpha[OffsetTmp_n + 7];

		betaminus[OffsetTmp_n] = beta[OffsetTmp_n] + r_0[i] + alpha[OffsetTmp_n];
		betaminus[OffsetTmp_n + 1] = beta[OffsetTmp_n + 4] + r_0[i] + alpha[OffsetTmp_n + 1];
		betaminus[OffsetTmp_n + 2] = beta[OffsetTmp_n + 5] + r_1[i] + alpha[OffsetTmp_n + 2];
		betaminus[OffsetTmp_n + 3] = beta[OffsetTmp_n + 1] + r_1[i] + alpha[OffsetTmp_n + 3];
		betaminus[OffsetTmp_n + 4] = beta[OffsetTmp_n + 2] + r_1[i] + alpha[OffsetTmp_n + 4];
		betaminus[OffsetTmp_n + 5] = beta[OffsetTmp_n + 6] + r_1[i] + alpha[OffsetTmp_n + 5];
		betaminus[OffsetTmp_n + 6] = beta[OffsetTmp_n + 7] + r_0[i] + alpha[OffsetTmp_n + 6];
		betaminus[OffsetTmp_n + 7] = beta[OffsetTmp_n + 3] + r_0[i] + alpha[OffsetTmp_n + 7];


		beta[OffsetTmp_p] = max(beta[OffsetTmp_n] + r_0[i], beta[OffsetTmp_n + 4] + r0[i]);
		beta[OffsetTmp_p + 1] = max(beta[OffsetTmp_n] + r0[i], beta[OffsetTmp_n + 4] + r_0[i]);
		beta[OffsetTmp_p + 2] = max(beta[OffsetTmp_n + 1] + r1[i], beta[OffsetTmp_n + 5] + r_1[i]);
		beta[OffsetTmp_p + 3] = max(beta[OffsetTmp_n + 5] + r1[i], beta[OffsetTmp_n + 1] + r_1[i]);
		beta[OffsetTmp_p + 4] = max(beta[OffsetTmp_n + 6] + r1[i], beta[OffsetTmp_n + 2] + r_1[i]);
		beta[OffsetTmp_p + 5] = max(beta[OffsetTmp_n + 2] + r1[i], beta[OffsetTmp_n + 6] + r_1[i]);
		beta[OffsetTmp_p + 6] = max(beta[OffsetTmp_n + 3] + r0[i], beta[OffsetTmp_n + 7] + r_0[i]);
		beta[OffsetTmp_p + 7] = max(beta[OffsetTmp_n + 3] + r_0[i], beta[OffsetTmp_n + 7] + r0[i]);

		OffsetTmp_n -= StateNum;
		OffsetTmp_p -= StateNum;
	}

	betaplus[OffsetTmp_n] = beta[OffsetTmp_n + 4] + r0[0] + alpha[OffsetTmp_n];
	betaplus[OffsetTmp_n + 1] = beta[OffsetTmp_n] + r0[0] + alpha[OffsetTmp_n + 1];
	betaplus[OffsetTmp_n + 2] = beta[OffsetTmp_n + 1] + r1[0] + alpha[OffsetTmp_n + 2];
	betaplus[OffsetTmp_n + 3] = beta[OffsetTmp_n + 5] + r1[0] + alpha[OffsetTmp_n + 3];
	betaplus[OffsetTmp_n + 4] = beta[OffsetTmp_n + 6] + r1[0] + alpha[OffsetTmp_n + 4];
	betaplus[OffsetTmp_n + 5] = beta[OffsetTmp_n + 2] + r1[0] + alpha[OffsetTmp_n + 5];
	betaplus[OffsetTmp_n + 6] = beta[OffsetTmp_n + 3] + r0[0] + alpha[OffsetTmp_n + 6];
	betaplus[OffsetTmp_n + 7] = beta[OffsetTmp_n + 7] + r0[0] + alpha[OffsetTmp_n + 7];

	betaminus[OffsetTmp_n] = beta[OffsetTmp_n] + r_0[0] + alpha[OffsetTmp_n];
	betaminus[OffsetTmp_n + 1] = beta[OffsetTmp_n + 4] + r_0[0] + alpha[OffsetTmp_n + 1];
	betaminus[OffsetTmp_n + 2] = beta[OffsetTmp_n + 5] + r_1[0] + alpha[OffsetTmp_n + 2];
	betaminus[OffsetTmp_n + 3] = beta[OffsetTmp_n + 1] + r_1[0] + alpha[OffsetTmp_n + 3];
	betaminus[OffsetTmp_n + 4] = beta[OffsetTmp_n + 2] + r_1[0] + alpha[OffsetTmp_n + 4];
	betaminus[OffsetTmp_n + 5] = beta[OffsetTmp_n + 6] + r_1[0] + alpha[OffsetTmp_n + 5];
	betaminus[OffsetTmp_n + 6] = beta[OffsetTmp_n + 7] + r_0[0] + alpha[OffsetTmp_n + 6];
	betaminus[OffsetTmp_n + 7] = beta[OffsetTmp_n + 3] + r_0[0] + alpha[OffsetTmp_n + 7];

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
	for (int i = 0; i < Len; ++i)
	{

		Output[i] =
			Mmax(
			Mmax(
			Mmax(betaplus[i*StateNum + 0], betaplus[i*StateNum + 1]),
			Mmax(betaplus[i*StateNum + 2], betaplus[i*StateNum + 3])
			),
			Mmax(
			Mmax(betaplus[i*StateNum + 4], betaplus[i*StateNum + 5]),
			Mmax(betaplus[i*StateNum + 6], betaplus[i*StateNum + 7])
			)
			)
			-
			Mmax(
			Mmax(                                     //之前为i+1错误
			Mmax(betaminus[i*StateNum + 0], betaminus[i*StateNum + 1]),
			Mmax(betaminus[i*StateNum + 2], betaminus[i*StateNum + 3])
			),
			Mmax(
			Mmax(betaminus[i*StateNum + 4], betaminus[i*StateNum + 5]),
			Mmax(betaminus[i*StateNum + 6], betaminus[i*StateNum + 7])
			)
			)
			- Ls[i] - Le[i];
	}
	free(alpha);
	free(betaplus);
	free(betaminus);
	free(beta);
	free(r0);
	free(r1);
	free(r_1);
	free(r_0);
}
bool TurboDecoderNorOri(float*Ls, float*Lp0, float*Lp1, float*LIns, int *Table, float*Output, int Len, int InterNum)//Len为加上Mem的
{
	float*AlphaM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*BetaPlusM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*BetaMinusM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*Le1 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	float*Le0 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	float*La1 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	memset(Output, 0, Len*sizeof(float));
	//All System Go
	for (int i = 1; i <= InterNum; ++i)
	{
		//MaxLogRscDecoder(Ls, Lp0, Output, Len, Le0,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderOri(Ls, Lp0, Output, Len, Le0);
		for (int k = 0; k < Len - Mem; ++k)
			La1[k] = Le0[Table[k]];
		for (int k = Len - Mem; k < Len; ++k)
			La1[k] = Le0[k];
		//MaxLogRscDecoder(Ls, Lp1, La1, Len,Le1,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderOri(LIns, Lp1, La1, Len, Le1);
		for (int k = 0; k < Len - Mem; ++k)
			Output[Table[k]] = Le1[k];
		for (int k = Len - Mem; k < Len; ++k)
			Output[k] = Le1[k];
	}
	//
	_mm_free(AlphaM);
	_mm_free(BetaPlusM);
	_mm_free(BetaMinusM);
	_mm_free(Le1);
	_mm_free(La1);
	_mm_free(Le0);
	return 1;
}
bool TurboDecoderNorAno(float*Ls, float*Lp0, float*Lp1, float*LIns, int *Table, float*Output, int Len, int InterNum)//Len为加上Mem的
{
	float*AlphaM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*BetaPlusM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*BetaMinusM = (float*)_mm_malloc(Len*sizeof(__m256)*sizeof(float), sizeof(__m256));
	float*Le1 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	float*Le0 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	float*La1 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));
	memset(Output, 0, Len*sizeof(float));
	//All System Go
	for (int i = 1; i <= InterNum; ++i)
	{
		//MaxLogRscDecoder(Ls, Lp0, Output, Len, Le0,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderAno(Ls, Lp0, Output, Len, Le0);
		for (int k = 0; k < Len - Mem; ++k)
			La1[k] = Le0[Table[k]];
		for (int k = Len - Mem; k < Len; ++k)
			La1[k] = Le0[k];
		//MaxLogRscDecoder(Ls, Lp1, La1, Len,Le1,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderAno(LIns, Lp1, La1, Len, Le1);
		for (int k = 0; k < Len - Mem; ++k)
			Output[Table[k]] = Le1[k];
		for (int k = Len - Mem; k < Len; ++k)
			Output[k] = Le1[k];
	}
	//
	_mm_free(AlphaM);
	_mm_free(BetaPlusM);
	_mm_free(BetaMinusM);
	_mm_free(Le1);
	_mm_free(La1);
	_mm_free(Le0);
	return 1;
}
bool TurboDecoderNorSInt(short*Ls, short*Lp0, short*Lp1, short*LIns, int *Table, short*Output, int Len, int InterNum)//Len为加上Mem的
{
	short*Le1 = (short*)malloc(Len*sizeof(short));
	short*Le0 = (short*)malloc(Len*sizeof(short));
	short*La1 = (short*)malloc(Len *sizeof(short));
	memset(Output, 0, Len*sizeof(float));
	//_mm_free(Le0);
	//All System Go
	for (int i = 1; i <= InterNum; ++i)
	{
		//MaxLogRscDecoder(Ls, Lp0, Output, Len, Le0,AlphaM,BetaPlusM,BetaMinusM);
		//Le0 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
		MaxLogSubDecoderSInt(Ls, Lp0, Output, Len, Le0);
		for (int k = 0; k < Len - Mem; ++k)
			La1[k] = Le0[Table[k]];
		for (int k = Len - Mem; k < Len; ++k)
			La1[k] = Le0[k];
		//_mm_free(Le0);
		//MaxLogRscDecoder(Ls, Lp1, La1, Len,Le1,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderSInt(LIns, Lp1, La1, Len, Le1);
		for (int k = 0; k < Len - Mem; ++k)
			Output[Table[k]] = Le1[k];
		for (int k = Len - Mem; k < Len; ++k)
			Output[k] = Le1[k];
	}
	//
	free(Le1);
	free(La1);
	free(Le0);
	return 1;
}
bool TurboDecoderNorSIntAno(short*Ls, short*Lp0, short*Lp1, short*LIns, int *Table, short*Output, int Len, int InterNum)//Len为加上Mem的
{
	short*Le1 = (short*)malloc(Len*sizeof(short));
	short*Le0 = (short*)malloc(Len*sizeof(short));
	short*La1 = (short*)malloc(Len *sizeof(short));
	memset(Output, 0, Len*sizeof(float));
	//_mm_free(Le0);
	//All System Go
	for (int i = 1; i <= InterNum; ++i)
	{
		//MaxLogRscDecoder(Ls, Lp0, Output, Len, Le0,AlphaM,BetaPlusM,BetaMinusM);
		//Le0 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256));
		MaxLogSubDecoderSIntAno(Ls, Lp0, Output, Len, Le0);
		for (int k = 0; k < Len - Mem; ++k)
			La1[k] = Le0[Table[k]];
		for (int k = Len - Mem; k < Len; ++k)
			La1[k] = Le0[k];
		//_mm_free(Le0);
		//MaxLogRscDecoder(Ls, Lp1, La1, Len,Le1,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderSIntAno(LIns, Lp1, La1, Len, Le1);
		for (int k = 0; k < Len - Mem; ++k)
			Output[Table[k]] = Le1[k];
		for (int k = Len - Mem; k < Len; ++k)
			Output[k] = Le1[k];
	}
	//
	free(Le1);
	free(La1);
	free(Le0);
	return 1;
}
bool TurboDecoderNorSIntNor(short*Ls, short*Lp0, short*Lp1, short*LIns, int *Table, short*Output, int Len, int InterNum)//Len为加上Mem的
{
	short*Le1 = (short*)malloc(Len*sizeof(short));
	short*Le0 = (short*)malloc(Len*sizeof(short));
	short*La1 = (short*)malloc(Len *sizeof(short));
	memset(Output, 0, Len*sizeof(float));
	//All System Go
	for (int i = 1; i <= InterNum; ++i)
	{
		//MaxLogRscDecoder(Ls, Lp0, Output, Len, Le0,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderSIntNor(Ls, Lp0, Output, Len, Le0);
		for (int k = 0; k < Len - Mem; ++k)
			La1[k] = Le0[Table[k]];
		for (int k = Len - Mem; k < Len; ++k)
			La1[k] = Le0[k];
		//MaxLogRscDecoder(Ls, Lp1, La1, Len,Le1,AlphaM,BetaPlusM,BetaMinusM);
		MaxLogSubDecoderSIntNor(LIns, Lp1, La1, Len, Le1);
		for (int k = 0; k < Len - Mem; ++k)
			Output[Table[k]] = Le1[k];
		for (int k = Len - Mem; k < Len; ++k)
			Output[k] = Le1[k];
	}
	//
	free(Le1);
	free(La1);
	free(Le0);
	return 1;
}
//FinalVersion
inline bool MaxLogRscDecoderSSEInt(short*Ls, short*Lp, short*La, int Len, short*Output, __m128i*BetaM)
{//注Common大小为8*__m256i
	__m128i AlphaBetaMain, AlphaBetaSub, Offset, RL1, RL0;
	__m256i R1, R2, R3, R4, R5, R6, R7, Rr, RTmp;
	__m128i ShuffleAlphaBetaPlus = _mm_setr_epi8(8, 9, 0, 1, 2, 3, 10, 11, 12, 13, 4, 5, 6, 7, 14, 15);//低位为Plus，高位为Sub
	__m128i ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 8, 9, 10, 11, 2, 3, 4, 5, 12, 13, 14, 15, 6, 7);
	__m256i ShuffleLshift64Invert = _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
	__m256i ShuffleLshift32 = _mm256_setr_epi8(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3);
	__m256i ShuffleLshift16 = _mm256_setr_epi8(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1);
	short L1;
	short L0;
	__m128i*BetaIndex = BetaM;//指向最后一个元素的下一个位置
	__m256i*BetaIndex256 = (__m256i*)BetaM;
	__m128i*RrIndex = (__m128i*)(Ls + Len - 8);
	__m128i*RaIndex = (__m128i*)(La + Len - 8);
	__m128i*RpIndex = (__m128i*)(Lp + Len - 8);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	_mm_store_si128(BetaIndex++, AlphaBetaMain);
	//初始化Rr
	RL1 = _mm_loadu_si128(RrIndex);
	Offset = _mm_loadu_si128(RaIndex);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_loadu_si128(RpIndex);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
	//重新定义Rr值
	RrIndex = (__m128i*)(Ls + ((Len >> 3) << 3) - 8);
	RaIndex = (__m128i*)(La + ((Len >> 3) << 3) - 8);
	RpIndex = (__m128i*)(Lp + ((Len >> 3) << 3) - 8);
	/*先计算后向度量*/
	for (int i = Len - 1; i>0; --i)
	{
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//RL1=_mm_bsrli_si128(RL1, 2);//00111001
		//RL0 = _mm_bsrli_si128(RL0, 2);
		/*先得出分支度量*/
		//L1 = (Ls[i] + La[i] + Lp[i]) >> 1;
		//L0 = L1 - Lp[i];
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//归一化
		if ((i & 0x7) == 0)//相当于每8次运算中有一次归一化
		{
			//更新
			RL1 = _mm_load_si128(RrIndex--);
			Offset = _mm_load_si128(RaIndex--);
			RL1 = _mm_add_epi16(RL1, Offset);
			Offset = _mm_load_si128(RpIndex--);
			RL0 = _mm_sub_epi16(RL1, Offset);
			RL1 = _mm_add_epi16(RL1, Offset);
			RL0 = _mm_srai_epi16(RL0, 1);
			RL1 = _mm_srai_epi16(RL1, 1);
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
			//用AVX2指令实现归一化
			Offset = _mm_broadcastw_epi16(AlphaBetaMain);
			AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		}
		_mm_store_si128(BetaIndex++, AlphaBetaMain);
	}
	/*再计算前向度量并计算似然比*/
	__m128i*LsIndex = (__m128i*)Ls;
	__m128i*LaIndex = (__m128i*)La;
	__m128i*OutputIndex = (__m128i*)Output;
	BetaIndex = BetaM + Len;
	RrIndex = (__m128i*)Ls;
	RaIndex = (__m128i*)La;
	RpIndex = (__m128i*)Lp;
	ShuffleLshift64Invert = _mm256_setr_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
	ShuffleAlphaBetaPlus = _mm_setr_epi8(2, 3, 4, 5, 10, 11, 12, 13, 0, 1, 6, 7, 8, 9, 14, 15);
	ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 6, 7, 8, 9, 14, 15, 2, 3, 4, 5, 10, 11, 12, 13);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	//初始化Rr
	RL1 = _mm_load_si128(RrIndex++);
	Offset = _mm_load_si128(RaIndex++);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_load_si128(RpIndex++);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	for (int i = 0; i < Len; ++i)
	{
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
		//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		RL1 = _mm_bsrli_si128(RL1, 2);//00111001
		RL0 = _mm_bsrli_si128(RL0, 2);
		//先得出分支度量
		//L1 = (Ls[i] + La[i] + Lp[i]) >> 1;
		//L0 = L1 - Lp[i];
		//Offset = _mm_set_epi16(L1, L1, L0, L0, L0, L0, L1, L1);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		/*******处理Hmax*******/
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		Rr = _mm256_max_epi16(Rr, RTmp);
		switch (i & 0x7)
		{
		case 0:R1 = Rr; break;//低位为Plus
		case 1:R2 = Rr; break;
		case 2:R3 = Rr; break;
		case 3:R4 = Rr; break;
		case 4:R5 = Rr; break;
		case 5:R6 = Rr; break;
		case 6:R7 = Rr; break;
		case 7:
		{
				  //更新Rr
				  RL1 = _mm_load_si128(RrIndex++);
				  Offset = _mm_load_si128(RaIndex++);
				  RL1 = _mm_add_epi16(RL1, Offset);
				  Offset = _mm_load_si128(RpIndex++);
				  RL0 = _mm_sub_epi16(RL1, Offset);
				  RL1 = _mm_add_epi16(RL1, Offset);
				  RL0 = _mm_srai_epi16(RL0, 1);
				  RL1 = _mm_srai_epi16(RL1, 1);
				  //归一化
				  Offset = _mm_broadcastw_epi16(AlphaBetaMain);
				  AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
				  //此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
				  R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
				  R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
				  R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
				  R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
				  R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
				  R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
				  R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
				  //Plus与Minus相减,到此Rr的低八位数值全一样
				  Offset = _mm256_extractf128_si256(R1, 1);
				  Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
				  //减去Ls+La
				  AlphaBetaSub = _mm_load_si128(LsIndex++);
				  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
				  AlphaBetaSub = _mm_load_si128(LaIndex++);
				  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
				  //存入内存
				  _mm_store_si128(OutputIndex++, Offset);
				  break;
		}
		}
	}
	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去Ls+La
		AlphaBetaSub = _mm_load_si128(LsIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_load_si128(LaIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		_mm_store_si128(OutputIndex++, Offset);
	}
	return 1;
}
double TurboDecoderSSEInt(short*Ls, short*Lp0, short*Lp1, short*LIns, int*Table, short*Output, int Len, int InterNum)
{//译码器号_译码器内支路号
	//__m128i*AlphaM = (__m128i*)_mm_malloc(Len*sizeof(__m128i), sizeof(__m128i));
	__m128i*AlphaM = (__m128i*)_mm_malloc(Len*sizeof(__m128i), sizeof(__m256i));
	short*La = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m128i));
	memset(Output, 0, Len*sizeof(short));
	//All System Go
	LARGE_INTEGER Start, End, Freq;
	QueryPerformanceCounter(&Start);
	QueryPerformanceFrequency(&Freq);
	for (int i = 1; i <= InterNum; ++i)
	{
		MaxLogRscDecoderSSEInt(Ls, Lp0, Output, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[k] = La[Table[k]];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = 0;
		}
		MaxLogRscDecoderSSEInt(LIns, Lp1, Output, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[Table[k]] = La[k];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = 0;
		}
	}
	QueryPerformanceCounter(&End);
	//
	_mm_free(AlphaM);
	_mm_free(La);
	return (double)(End.QuadPart - Start.QuadPart) / (double)Freq.QuadPart;
}
inline bool MaxLogRscDecoderSSEIntSpecify(short*Ls, short*Lp, short*La, int Len, short*Output, __m128i*BetaM)
{//注Common大小为8*__m256i
	__m128i AlphaBetaMain, AlphaBetaSub, Offset, RL1, RL0, RTable;
	__m256i R1, R2, R3, R4, R5, R6, R7, Rr, RTmp;
	__m128i ShuffleAlphaBetaPlus = _mm_setr_epi8(8, 9, 0, 1, 2, 3, 10, 11, 12, 13, 4, 5, 6, 7, 14, 15);//低位为Plus，高位为Sub
	__m128i ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 8, 9, 10, 11, 2, 3, 4, 5, 12, 13, 14, 15, 6, 7);
	__m256i ShuffleLshift64Invert = _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
	__m256i ShuffleLshift32 = _mm256_setr_epi8(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3);
	__m256i ShuffleLshift16 = _mm256_setr_epi8(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1);
	short L1;
	short L0;
	__m128i*BetaIndex = BetaM;//指向最后一个元素的下一个位置
	__m128i*RrIndex;
	__m128i*RaIndex;
	__m128i*RpIndex;
	/*计算后向度量*/
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	_mm_store_si128(BetaIndex++, AlphaBetaMain);
	//当末尾不规整时
	int i = Len - 1;
	if ((Len & 0x7) != 0)
	{
		RrIndex = (__m128i*)(Ls + Len - 8);
		RaIndex = (__m128i*)(La + Len - 8);
		RpIndex = (__m128i*)(Lp + Len - 8);
		//初始化Rr
		RL1 = _mm_loadu_si128(RrIndex);
		Offset = _mm_loadu_si128(RaIndex);
		RL1 = _mm_add_epi16(RL1, Offset);
		Offset = _mm_loadu_si128(RpIndex);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
		for (; i > 0; --i)
		{
			//归一化
			if ((i & 0x7) == 7)//相当于每8次运算中有一次归一化
				break;
			Offset = _mm_broadcastw_epi16(RL1);
			AlphaBetaSub = _mm_broadcastw_epi16(RL0);
			Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
			//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
			AlphaBetaSub = AlphaBetaMain;
			//对于相加的部分
			AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
			AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
			//对于相减的部分
			AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
			AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
			//取最大值并存AlphaM
			AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
			_mm_store_si128(BetaIndex++, AlphaBetaMain);
		}
	}
	//再处理8位规整的
	//重新定义Rr值
	RrIndex = (__m128i*)(Ls + ((Len >> 3) << 3));
	RaIndex = (__m128i*)(La + ((Len >> 3) << 3));
	RpIndex = (__m128i*)(Lp + ((Len >> 3) << 3));
	//更新
	RL1 = _mm_load_si128(--RrIndex);
	Offset = _mm_load_si128(--RaIndex);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_load_si128(--RpIndex);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));//与下面合并
	//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
	//用AVX2指令实现归一化
	Offset = _mm_broadcastw_epi16(AlphaBetaMain);
	AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
	for (; i>7; i -= 8)
	{
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理7
		//Offset = _mm_broadcastw_epi16(RL1);
		//AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		Offset = _mm_blend_epi16(RL1, RL0, 0xFE);//11111110
		RTable = _mm_setr_epi8(0, 1, 14, 15, 14, 15, 0, 1, 0, 1, 14, 15, 14, 15, 0, 1);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理6
		Offset = _mm_blend_epi16(RL1, RL0, 0x01);//00000001
		RTable = _mm_setr_epi8(2, 3, 0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理5
		Offset = _mm_blend_epi16(RL1, RL0, 0x02);//00000010
		RTable = _mm_setr_epi8(4, 5, 2, 3, 2, 3, 4, 5, 4, 5, 2, 3, 2, 3, 4, 5);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理4
		Offset = _mm_blend_epi16(RL1, RL0, 0x04);//00000100
		RTable = _mm_setr_epi8(6, 7, 4, 5, 4, 5, 6, 7, 6, 7, 4, 5, 4, 5, 6, 7);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理3
		Offset = _mm_blend_epi16(RL1, RL0, 0x08);//00001000
		RTable = _mm_setr_epi8(8, 9, 6, 7, 6, 7, 8, 9, 8, 9, 6, 7, 6, 7, 8, 9);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理2
		Offset = _mm_blend_epi16(RL1, RL0, 0x10);//00010000
		RTable = _mm_setr_epi8(10, 11, 8, 9, 8, 9, 10, 11, 10, 11, 8, 9, 8, 9, 10, 11);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理1
		Offset = _mm_blend_epi16(RL1, RL0, 0x20);//00100000
		RTable = _mm_setr_epi8(12, 13, 10, 11, 10, 11, 12, 13, 12, 13, 10, 11, 10, 11, 12, 13);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理0
		Offset = _mm_blend_epi16(RL1, RL0, 0x40);//01000000
		RTable = _mm_setr_epi8(14, 15, 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);

		//归一化后存入内存，并更新
		//if ((i & 0x7) == 0)//相当于每8次运算中有一次归一化
		RL1 = _mm_load_si128(--RrIndex);
		Offset = _mm_load_si128(--RaIndex);
		RL1 = _mm_add_epi16(RL1, Offset);
		Offset = _mm_load_si128(--RpIndex);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));//与下面合并
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//用AVX2指令实现归一化
		Offset = _mm_broadcastw_epi16(AlphaBetaMain);
		AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);
	}

	//处理最后7位
	//初始化Rr
	//RL1 = _mm_loadu_si128((__m128i*)Ls);
	//Offset = _mm_loadu_si128((__m128i*)La);
	//RL1 = _mm_add_epi16(RL1, Offset);
	//Offset = _mm_loadu_si128((__m128i*)Lp);
	//RL0 = _mm_sub_epi16(RL1, Offset);
	//RL1 = _mm_add_epi16(RL1, Offset);
	//RL0 = _mm_srai_epi16(RL0, 1);
	//RL1 = _mm_srai_epi16(RL1, 1);
	//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
	for (; i > 0; --i)
	{
		//归一化
		//if ((i & 0x7) == 0x7)//相当于每8次运算中有一次归一化
		//	break;
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);
	}

	/*再计算前向度量并计算似然比*/
	__m128i*LsIndex = (__m128i*)Ls;
	__m128i*LaIndex = (__m128i*)La;
	__m128i*OutputIndex = (__m128i*)Output;
	BetaIndex = BetaM + Len;
	RrIndex = (__m128i*)Ls;
	RaIndex = (__m128i*)La;
	RpIndex = (__m128i*)Lp;
	ShuffleLshift64Invert = _mm256_setr_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
	ShuffleAlphaBetaPlus = _mm_setr_epi8(2, 3, 4, 5, 10, 11, 12, 13, 0, 1, 6, 7, 8, 9, 14, 15);
	ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 6, 7, 8, 9, 14, 15, 2, 3, 4, 5, 10, 11, 12, 13);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	//初始化Rr
	RL1 = _mm_load_si128(RrIndex++);
	Offset = _mm_load_si128(RaIndex++);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_load_si128(RpIndex++);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	i = 8;
	//开始时处理规整的部分
	for (; i <= Len; i += 8)
	{
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理0
		//Offset = _mm_broadcastw_epi16(RL1);
		//AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		//Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
		//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_set_epi16(L1, L1, L0, L0, L0, L0, L1, L1);
		Offset = _mm_blend_epi16(RL1, RL0, 0xFE);//11111110
		RTable = _mm_setr_epi8(0, 1, 0, 1, 14, 15, 14, 15, 14, 15, 14, 15, 0, 1, 0, 1);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R1 = _mm256_max_epi16(Rr, RTmp);

		//处理1
		Offset = _mm_blend_epi16(RL1, RL0, 0x01);//00000001
		RTable = _mm_setr_epi8(2, 3, 2, 3, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R2 = _mm256_max_epi16(Rr, RTmp);

		//处理2
		Offset = _mm_blend_epi16(RL1, RL0, 0x02);//00000010
		RTable = _mm_setr_epi8(4, 5, 4, 5, 2, 3, 2, 3, 2, 3, 2, 3, 4, 5, 4, 5);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R3 = _mm256_max_epi16(Rr, RTmp);

		//处理3
		Offset = _mm_blend_epi16(RL1, RL0, 0x04);//00000100
		RTable = _mm_setr_epi8(6, 7, 6, 7, 4, 5, 4, 5, 4, 5, 4, 5, 6, 7, 6, 7);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R4 = _mm256_max_epi16(Rr, RTmp);

		//处理4
		Offset = _mm_blend_epi16(RL1, RL0, 0x08);//00001000
		RTable = _mm_setr_epi8(8, 9, 8, 9, 6, 7, 6, 7, 6, 7, 6, 7, 8, 9, 8, 9);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R5 = _mm256_max_epi16(Rr, RTmp);

		//处理5
		Offset = _mm_blend_epi16(RL1, RL0, 0x10);//00010000
		RTable = _mm_setr_epi8(10, 11, 10, 11, 8, 9, 8, 9, 8, 9, 8, 9, 10, 11, 10, 11);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R6 = _mm256_max_epi16(Rr, RTmp);

		//处理6
		Offset = _mm_blend_epi16(RL1, RL0, 0x20);//00100000
		RTable = _mm_setr_epi8(12, 13, 12, 13, 10, 11, 10, 11, 10, 11, 10, 11, 12, 13, 12, 13);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R7 = _mm256_max_epi16(Rr, RTmp);

		//处理7
		Offset = _mm_blend_epi16(RL1, RL0, 0x40);//01000000
		RTable = _mm_setr_epi8(14, 15, 14, 15, 12, 13, 12, 13, 12, 13, 12, 13, 14, 15, 14, 15);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		Rr = _mm256_max_epi16(Rr, RTmp);

		//更新Rr
		RL1 = _mm_load_si128(RrIndex++);
		Offset = _mm_load_si128(RaIndex++);
		RL1 = _mm_add_epi16(RL1, Offset);
		Offset = _mm_load_si128(RpIndex++);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//归一化
		Offset = _mm_broadcastw_epi16(AlphaBetaMain);
		AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		//此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去Ls+La
		AlphaBetaSub = _mm_load_si128(LsIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_load_si128(LaIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		if ((i & 0x3F) == 0)//纯粹为了减小运行次数
		{
			AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
			Offset = _mm_min_epi16(Offset, AlphaBetaSub);
			AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
			Offset = _mm_max_epi16(Offset, AlphaBetaSub);
		}
		_mm_store_si128(OutputIndex++, Offset);
	}


	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		//更新Rr
		//RL1 = _mm_load_si128(--RrIndex);
		//Offset = _mm_load_si128(--RaIndex);
		//RL1 = _mm_add_epi16(RL1, Offset);
		//Offset = _mm_load_si128(--RpIndex);
		//RL0 = _mm_sub_epi16(RL1, Offset);
		//RL1 = _mm_add_epi16(RL1, Offset);
		//RL0 = _mm_srai_epi16(RL0, 1);
		//RL1 = _mm_srai_epi16(RL1, 1);
		//ShuffleLshift16 = _mm256_setr_epi8(14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理最后的剩余部分
		i -= 8;
		for (; i < Len; ++i)
		{
			Offset = _mm_broadcastw_epi16(RL1);
			AlphaBetaSub = _mm_broadcastw_epi16(RL0);
			Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
			AlphaBetaSub = AlphaBetaMain;
			//先计算加
			AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
			AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
			//再计算相减
			AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
			AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
			//再加上Beta
			Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
			Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
			Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
			//直接更新下一次
			AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
			//处理Hmax
			Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
			//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
			Rr = _mm256_max_epi16(Rr, RTmp);
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
			Rr = _mm256_max_epi16(Rr, RTmp);
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
			Rr = _mm256_max_epi16(Rr, RTmp);

			switch (i & 0x7)
			{
			case 0:R1 = Rr; break;//低位为Plus
			case 1:R2 = Rr; break;
			case 2:R3 = Rr; break;
			case 3:R4 = Rr; break;
			case 4:R5 = Rr; break;
			case 5:R6 = Rr; break;
			case 6:R7 = Rr; break;
			case 7:
			{
					  //更新Rr
					  RL1 = _mm_load_si128(RrIndex++);
					  Offset = _mm_load_si128(RaIndex++);
					  RL1 = _mm_add_epi16(RL1, Offset);
					  Offset = _mm_load_si128(RpIndex++);
					  RL0 = _mm_sub_epi16(RL1, Offset);
					  RL1 = _mm_add_epi16(RL1, Offset);
					  RL0 = _mm_srai_epi16(RL0, 1);
					  RL1 = _mm_srai_epi16(RL1, 1);
					  //归一化
					  Offset = _mm_broadcastw_epi16(AlphaBetaMain);
					  AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
					  //此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
					  R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
					  R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
					  R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
					  R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
					  R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
					  R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
					  R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
					  //Plus与Minus相减,到此Rr的低八位数值全一样
					  Offset = _mm256_extractf128_si256(R1, 1);
					  Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
					  //减去Ls+La
					  AlphaBetaSub = _mm_load_si128(LsIndex++);
					  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
					  AlphaBetaSub = _mm_load_si128(LaIndex++);
					  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
					  //存入内存
					  // if ((i & 0xFF) == 0xFF)
					  // {
					  //  AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
					  //  Offset = _mm_min_epi16(Offset, AlphaBetaSub);
					  // AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
					  // Offset = _mm_max_epi16(Offset, AlphaBetaSub);
					  //}
					  _mm_store_si128(OutputIndex++, Offset);
			}
			}
		}
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去Ls+La
		AlphaBetaSub = _mm_load_si128(LsIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_load_si128(LaIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
		Offset = _mm_min_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
		Offset = _mm_max_epi16(Offset, AlphaBetaSub);
		_mm_store_si128(OutputIndex++, Offset);
	}
	return 1;
}
double TurboDecoderSSEIntSpecify(short*Ls, short*Lp0, short*Lp1, short*LIns, int*Table, short*Output, int Len, int InterNum)
{//译码器号_译码器内支路号
	//__m128i*AlphaM = (__m128i*)_mm_malloc(Len*sizeof(__m128i), sizeof(__m128i));
	__m128i*AlphaM = (__m128i*)_mm_malloc(Len*sizeof(__m128i), sizeof(__m128i));
	short*La = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m128i));
	memset(Output, 0, Len*sizeof(short));
	//All System Go
	LARGE_INTEGER Start, End, Freq;
	QueryPerformanceCounter(&Start);
	QueryPerformanceFrequency(&Freq);
	for (int i = 1; i <= InterNum; ++i)
	{
		MaxLogRscDecoderSSEIntSpecify(Ls, Lp0, Output, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[k] = La[Table[k]];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = 0;
		}
		MaxLogRscDecoderSSEIntSpecify(LIns, Lp1, Output, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[Table[k]] = La[k];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = 0;
		}
	}
	QueryPerformanceCounter(&End);
	//
	_mm_free(AlphaM);
	_mm_free(La);
	return (double)(End.QuadPart - Start.QuadPart) / (double)Freq.QuadPart;
}
inline bool MaxLogRscDecoderSSEIntSpecifyD0(short*Ls, short*Lp, short*La, int Len, short*Output, __m128i*BetaM)
{//注Common大小为8*__m256i
	__m128i AlphaBetaMain, AlphaBetaSub, Offset, RL1, RL0, RTable;
	__m256i R1, R2, R3, R4, R5, R6, R7, Rr, RTmp;
	__m128i ShuffleAlphaBetaPlus = _mm_setr_epi8(8, 9, 0, 1, 2, 3, 10, 11, 12, 13, 4, 5, 6, 7, 14, 15);//低位为Plus，高位为Sub
	__m128i ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 8, 9, 10, 11, 2, 3, 4, 5, 12, 13, 14, 15, 6, 7);
	__m256i ShuffleLshift64Invert = _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
	__m256i ShuffleLshift32 = _mm256_setr_epi8(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3);
	__m256i ShuffleLshift16 = _mm256_setr_epi8(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1);
	short L1;
	short L0;
	__m128i*BetaIndex = BetaM;//指向最后一个元素的下一个位置
	__m128i*RrIndex;
	__m128i*RaIndex;
	__m128i*RpIndex;
	/*计算后向度量*/
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	_mm_store_si128(BetaIndex++, AlphaBetaMain);
	//当末尾不规整时
	int i = Len - 1;
	if ((Len & 0x7) != 0)
	{
		RrIndex = (__m128i*)(Ls + Len - 8);
		RaIndex = (__m128i*)(La + Len - 8);
		RpIndex = (__m128i*)(Lp + Len - 8);
		//初始化Rr
		RL1 = _mm_loadu_si128(RrIndex);
		Offset = _mm_loadu_si128(RaIndex);
		RL1 = _mm_add_epi16(RL1, Offset);
		Offset = _mm_loadu_si128(RpIndex);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
		for (; i > 0; --i)
		{
			//归一化
			if ((i & 0x7) == 7)//相当于每8次运算中有一次归一化
				break;
			Offset = _mm_broadcastw_epi16(RL1);
			AlphaBetaSub = _mm_broadcastw_epi16(RL0);
			Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
			//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
			AlphaBetaSub = AlphaBetaMain;
			//对于相加的部分
			AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
			AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
			//对于相减的部分
			AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
			AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
			//取最大值并存AlphaM
			AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
			_mm_store_si128(BetaIndex++, AlphaBetaMain);
		}
	}
	//再处理8位规整的
	//重新定义Rr值
	RrIndex = (__m128i*)(Ls + ((Len >> 3) << 3));
	RaIndex = (__m128i*)(La + ((Len >> 3) << 3));
	RpIndex = (__m128i*)(Lp + ((Len >> 3) << 3));
	//更新
	RL1 = _mm_load_si128(--RrIndex);
	Offset = _mm_load_si128(--RaIndex);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_load_si128(--RpIndex);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));//与下面合并
	//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
	//用AVX2指令实现归一化
	Offset = _mm_broadcastw_epi16(AlphaBetaMain);
	AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
	for (; i>7; i -= 8)
	{
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理7
		//Offset = _mm_broadcastw_epi16(RL1);
		//AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		Offset = _mm_blend_epi16(RL1, RL0, 0xFE);//11111110
		RTable = _mm_setr_epi8(0, 1, 14, 15, 14, 15, 0, 1, 0, 1, 14, 15, 14, 15, 0, 1);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理6
		Offset = _mm_blend_epi16(RL1, RL0, 0x01);//00000001
		RTable = _mm_setr_epi8(2, 3, 0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理5
		Offset = _mm_blend_epi16(RL1, RL0, 0x02);//00000010
		RTable = _mm_setr_epi8(4, 5, 2, 3, 2, 3, 4, 5, 4, 5, 2, 3, 2, 3, 4, 5);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理4
		Offset = _mm_blend_epi16(RL1, RL0, 0x04);//00000100
		RTable = _mm_setr_epi8(6, 7, 4, 5, 4, 5, 6, 7, 6, 7, 4, 5, 4, 5, 6, 7);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理3
		Offset = _mm_blend_epi16(RL1, RL0, 0x08);//00001000
		RTable = _mm_setr_epi8(8, 9, 6, 7, 6, 7, 8, 9, 8, 9, 6, 7, 6, 7, 8, 9);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理2
		Offset = _mm_blend_epi16(RL1, RL0, 0x10);//00010000
		RTable = _mm_setr_epi8(10, 11, 8, 9, 8, 9, 10, 11, 10, 11, 8, 9, 8, 9, 10, 11);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理1
		Offset = _mm_blend_epi16(RL1, RL0, 0x20);//00100000
		RTable = _mm_setr_epi8(12, 13, 10, 11, 10, 11, 12, 13, 12, 13, 10, 11, 10, 11, 12, 13);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理0
		Offset = _mm_blend_epi16(RL1, RL0, 0x40);//01000000
		RTable = _mm_setr_epi8(14, 15, 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);

		//归一化后存入内存，并更新
		//if ((i & 0x7) == 0)//相当于每8次运算中有一次归一化
		RL1 = _mm_load_si128(--RrIndex);
		Offset = _mm_load_si128(--RaIndex);
		RL1 = _mm_add_epi16(RL1, Offset);
		Offset = _mm_load_si128(--RpIndex);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));//与下面合并
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//用AVX2指令实现归一化
		Offset = _mm_broadcastw_epi16(AlphaBetaMain);
		AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);
	}

	//处理最后7位
	//初始化Rr
	//RL1 = _mm_loadu_si128((__m128i*)Ls);
	//Offset = _mm_loadu_si128((__m128i*)La);
	//RL1 = _mm_add_epi16(RL1, Offset);
	//Offset = _mm_loadu_si128((__m128i*)Lp);
	//RL0 = _mm_sub_epi16(RL1, Offset);
	//RL1 = _mm_add_epi16(RL1, Offset);
	//RL0 = _mm_srai_epi16(RL0, 1);
	//RL1 = _mm_srai_epi16(RL1, 1);
	//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
	for (; i > 0; --i)
	{
		//归一化
		//if ((i & 0x7) == 0x7)//相当于每8次运算中有一次归一化
		//	break;
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);
	}

	/*再计算前向度量并计算似然比*/
	__m128i*LsIndex = (__m128i*)Ls;
	__m128i*LaIndex = (__m128i*)La;
	__m128i*OutputIndex = (__m128i*)Output;
	BetaIndex = BetaM + Len;
	RrIndex = (__m128i*)Ls;
	RaIndex = (__m128i*)La;
	RpIndex = (__m128i*)Lp;
	ShuffleLshift64Invert = _mm256_setr_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
	ShuffleAlphaBetaPlus = _mm_setr_epi8(2, 3, 4, 5, 10, 11, 12, 13, 0, 1, 6, 7, 8, 9, 14, 15);
	ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 6, 7, 8, 9, 14, 15, 2, 3, 4, 5, 10, 11, 12, 13);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	//初始化Rr
	RL1 = _mm_load_si128(RrIndex++);
	Offset = _mm_load_si128(RaIndex++);
	RL1 = _mm_add_epi16(RL1, Offset);
	Offset = _mm_load_si128(RpIndex++);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	i = 8;
	//开始时处理规整的部分
	for (; i <= Len; i += 8)
	{
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理0
		//Offset = _mm_broadcastw_epi16(RL1);
		//AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		//Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
		//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_set_epi16(L1, L1, L0, L0, L0, L0, L1, L1);
		Offset = _mm_blend_epi16(RL1, RL0, 0xFE);//11111110
		RTable = _mm_setr_epi8(0, 1, 0, 1, 14, 15, 14, 15, 14, 15, 14, 15, 0, 1, 0, 1);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R1 = _mm256_max_epi16(Rr, RTmp);

		//处理1
		Offset = _mm_blend_epi16(RL1, RL0, 0x01);//00000001
		RTable = _mm_setr_epi8(2, 3, 2, 3, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R2 = _mm256_max_epi16(Rr, RTmp);

		//处理2
		Offset = _mm_blend_epi16(RL1, RL0, 0x02);//00000010
		RTable = _mm_setr_epi8(4, 5, 4, 5, 2, 3, 2, 3, 2, 3, 2, 3, 4, 5, 4, 5);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R3 = _mm256_max_epi16(Rr, RTmp);

		//处理3
		Offset = _mm_blend_epi16(RL1, RL0, 0x04);//00000100
		RTable = _mm_setr_epi8(6, 7, 6, 7, 4, 5, 4, 5, 4, 5, 4, 5, 6, 7, 6, 7);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R4 = _mm256_max_epi16(Rr, RTmp);

		//处理4
		Offset = _mm_blend_epi16(RL1, RL0, 0x08);//00001000
		RTable = _mm_setr_epi8(8, 9, 8, 9, 6, 7, 6, 7, 6, 7, 6, 7, 8, 9, 8, 9);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R5 = _mm256_max_epi16(Rr, RTmp);

		//处理5
		Offset = _mm_blend_epi16(RL1, RL0, 0x10);//00010000
		RTable = _mm_setr_epi8(10, 11, 10, 11, 8, 9, 8, 9, 8, 9, 8, 9, 10, 11, 10, 11);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R6 = _mm256_max_epi16(Rr, RTmp);

		//处理6
		Offset = _mm_blend_epi16(RL1, RL0, 0x20);//00100000
		RTable = _mm_setr_epi8(12, 13, 12, 13, 10, 11, 10, 11, 10, 11, 10, 11, 12, 13, 12, 13);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R7 = _mm256_max_epi16(Rr, RTmp);

		//处理7
		Offset = _mm_blend_epi16(RL1, RL0, 0x40);//01000000
		RTable = _mm_setr_epi8(14, 15, 14, 15, 12, 13, 12, 13, 12, 13, 12, 13, 14, 15, 14, 15);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		Rr = _mm256_max_epi16(Rr, RTmp);

		//更新Rr
		RL1 = _mm_load_si128(RrIndex++);
		Offset = _mm_load_si128(RaIndex++);
		RL1 = _mm_add_epi16(RL1, Offset);
		Offset = _mm_load_si128(RpIndex++);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//归一化
		Offset = _mm_broadcastw_epi16(AlphaBetaMain);
		AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		//此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去La
		//AlphaBetaSub = _mm_load_si128(LsIndex++);
		//Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_load_si128(LaIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		if ((i & 0x3F) == 0)//纯粹为了减小运行次数
		{
			AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
			Offset = _mm_min_epi16(Offset, AlphaBetaSub);
			AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
			Offset = _mm_max_epi16(Offset, AlphaBetaSub);
		}
		_mm_store_si128(OutputIndex++, Offset);
	}


	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		//更新Rr
		//RL1 = _mm_load_si128(--RrIndex);
		//Offset = _mm_load_si128(--RaIndex);
		//RL1 = _mm_add_epi16(RL1, Offset);
		//Offset = _mm_load_si128(--RpIndex);
		//RL0 = _mm_sub_epi16(RL1, Offset);
		//RL1 = _mm_add_epi16(RL1, Offset);
		//RL0 = _mm_srai_epi16(RL0, 1);
		//RL1 = _mm_srai_epi16(RL1, 1);
		//ShuffleLshift16 = _mm256_setr_epi8(14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理最后的剩余部分
		i -= 8;
		for (; i < Len; ++i)
		{
			Offset = _mm_broadcastw_epi16(RL1);
			AlphaBetaSub = _mm_broadcastw_epi16(RL0);
			Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
			AlphaBetaSub = AlphaBetaMain;
			//先计算加
			AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
			AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
			//再计算相减
			AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
			AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
			//再加上Beta
			Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
			Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
			Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
			//直接更新下一次
			AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
			//处理Hmax
			Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
			//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
			Rr = _mm256_max_epi16(Rr, RTmp);
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
			Rr = _mm256_max_epi16(Rr, RTmp);
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
			Rr = _mm256_max_epi16(Rr, RTmp);

			switch (i & 0x7)
			{
			case 0:R1 = Rr; break;//低位为Plus
			case 1:R2 = Rr; break;
			case 2:R3 = Rr; break;
			case 3:R4 = Rr; break;
			case 4:R5 = Rr; break;
			case 5:R6 = Rr; break;
			case 6:R7 = Rr; break;
			case 7:
			{
					  //更新Rr
					  RL1 = _mm_load_si128(RrIndex++);
					  Offset = _mm_load_si128(RaIndex++);
					  RL1 = _mm_add_epi16(RL1, Offset);
					  Offset = _mm_load_si128(RpIndex++);
					  RL0 = _mm_sub_epi16(RL1, Offset);
					  RL1 = _mm_add_epi16(RL1, Offset);
					  RL0 = _mm_srai_epi16(RL0, 1);
					  RL1 = _mm_srai_epi16(RL1, 1);
					  //归一化
					  Offset = _mm_broadcastw_epi16(AlphaBetaMain);
					  AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
					  //此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
					  R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
					  R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
					  R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
					  R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
					  R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
					  R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
					  R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
					  //Plus与Minus相减,到此Rr的低八位数值全一样
					  Offset = _mm256_extractf128_si256(R1, 1);
					  Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
					  //减去La
					  AlphaBetaSub = _mm_load_si128(LaIndex++);
					  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
					  //存入内存
					  // if ((i & 0xFF) == 0xFF)
					  // {
					  //  AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
					  //  Offset = _mm_min_epi16(Offset, AlphaBetaSub);
					  // AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
					  // Offset = _mm_max_epi16(Offset, AlphaBetaSub);
					  //}
					  _mm_store_si128(OutputIndex++, Offset);
			}
			}
		}
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去La
		AlphaBetaSub = _mm_load_si128(LaIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
		Offset = _mm_min_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
		Offset = _mm_max_epi16(Offset, AlphaBetaSub);
		_mm_store_si128(OutputIndex++, Offset);
	}
	return 1;
}
inline bool MaxLogRscDecoderSSEIntSpecifyD1(short*Ls, short*Lp, int Len, short*Output, __m128i*BetaM)
{//注Common大小为8*__m256i
	__m128i AlphaBetaMain, AlphaBetaSub, Offset, RL1, RL0, RTable;
	__m256i R1, R2, R3, R4, R5, R6, R7, Rr, RTmp;
	__m128i ShuffleAlphaBetaPlus = _mm_setr_epi8(8, 9, 0, 1, 2, 3, 10, 11, 12, 13, 4, 5, 6, 7, 14, 15);//低位为Plus，高位为Sub
	__m128i ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 8, 9, 10, 11, 2, 3, 4, 5, 12, 13, 14, 15, 6, 7);
	__m256i ShuffleLshift64Invert = _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
	__m256i ShuffleLshift32 = _mm256_setr_epi8(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3);
	__m256i ShuffleLshift16 = _mm256_setr_epi8(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1);
	short L1;
	short L0;
	__m128i*BetaIndex = BetaM;//指向最后一个元素的下一个位置
	__m128i*RrIndex;
	__m128i*RpIndex;
	/*计算后向度量*/
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	_mm_store_si128(BetaIndex++, AlphaBetaMain);
	//当末尾不规整时
	int i = Len - 1;
	if ((Len & 0x7) != 0)
	{
		RrIndex = (__m128i*)(Ls + Len - 8);
		RpIndex = (__m128i*)(Lp + Len - 8);
		//初始化Rr
		RL1 = _mm_loadu_si128(RrIndex);
		//RL1 = _mm_add_epi16(RL1, Offset);
		Offset = _mm_loadu_si128(RpIndex);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
		for (; i > 0; --i)
		{
			//归一化
			if ((i & 0x7) == 7)//相当于每8次运算中有一次归一化
				break;
			Offset = _mm_broadcastw_epi16(RL1);
			AlphaBetaSub = _mm_broadcastw_epi16(RL0);
			Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
			//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
			AlphaBetaSub = AlphaBetaMain;
			//对于相加的部分
			AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
			AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
			//对于相减的部分
			AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
			AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
			//取最大值并存AlphaM
			AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
			_mm_store_si128(BetaIndex++, AlphaBetaMain);
		}
	}
	//再处理8位规整的
	//重新定义Rr值
	RrIndex = (__m128i*)(Ls + ((Len >> 3) << 3));
	RpIndex = (__m128i*)(Lp + ((Len >> 3) << 3));
	//更新
	RL1 = _mm_load_si128(--RrIndex);
	Offset = _mm_load_si128(--RpIndex);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));//与下面合并
	//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
	//用AVX2指令实现归一化
	Offset = _mm_broadcastw_epi16(AlphaBetaMain);
	AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
	for (; i>7; i -= 8)
	{
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理7
		//Offset = _mm_broadcastw_epi16(RL1);
		//AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		Offset = _mm_blend_epi16(RL1, RL0, 0xFE);//11111110
		RTable = _mm_setr_epi8(0, 1, 14, 15, 14, 15, 0, 1, 0, 1, 14, 15, 14, 15, 0, 1);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理6
		Offset = _mm_blend_epi16(RL1, RL0, 0x01);//00000001
		RTable = _mm_setr_epi8(2, 3, 0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理5
		Offset = _mm_blend_epi16(RL1, RL0, 0x02);//00000010
		RTable = _mm_setr_epi8(4, 5, 2, 3, 2, 3, 4, 5, 4, 5, 2, 3, 2, 3, 4, 5);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理4
		Offset = _mm_blend_epi16(RL1, RL0, 0x04);//00000100
		RTable = _mm_setr_epi8(6, 7, 4, 5, 4, 5, 6, 7, 6, 7, 4, 5, 4, 5, 6, 7);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理3
		Offset = _mm_blend_epi16(RL1, RL0, 0x08);//00001000
		RTable = _mm_setr_epi8(8, 9, 6, 7, 6, 7, 8, 9, 8, 9, 6, 7, 6, 7, 8, 9);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理2
		Offset = _mm_blend_epi16(RL1, RL0, 0x10);//00010000
		RTable = _mm_setr_epi8(10, 11, 8, 9, 8, 9, 10, 11, 10, 11, 8, 9, 8, 9, 10, 11);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理1
		Offset = _mm_blend_epi16(RL1, RL0, 0x20);//00100000
		RTable = _mm_setr_epi8(12, 13, 10, 11, 10, 11, 12, 13, 12, 13, 10, 11, 10, 11, 12, 13);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);

		//处理0
		Offset = _mm_blend_epi16(RL1, RL0, 0x40);//01000000
		RTable = _mm_setr_epi8(14, 15, 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);

		//归一化后存入内存，并更新
		//if ((i & 0x7) == 0)//相当于每8次运算中有一次归一化
		RL1 = _mm_load_si128(--RrIndex);
		Offset = _mm_load_si128(--RpIndex);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));//与下面合并
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//用AVX2指令实现归一化
		Offset = _mm_broadcastw_epi16(AlphaBetaMain);
		AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);
	}

	//处理最后7位
	//初始化Rr
	//RL1 = _mm_loadu_si128((__m128i*)Ls);
	//Offset = _mm_loadu_si128((__m128i*)La);
	//RL1 = _mm_add_epi16(RL1, Offset);
	//Offset = _mm_loadu_si128((__m128i*)Lp);
	//RL0 = _mm_sub_epi16(RL1, Offset);
	//RL1 = _mm_add_epi16(RL1, Offset);
	//RL0 = _mm_srai_epi16(RL0, 1);
	//RL1 = _mm_srai_epi16(RL1, 1);
	//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift64Invert));
	//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift64Invert));
	for (; i > 0; --i)
	{
		//归一化
		//if ((i & 0x7) == 0x7)//相当于每8次运算中有一次归一化
		//	break;
		Offset = _mm_broadcastw_epi16(RL1);
		AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x66);//01100110
		RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_set_epi16(L1, L0, L0, L1, L1, L0, L0, L1);
		AlphaBetaSub = AlphaBetaMain;
		//对于相加的部分
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//对于相减的部分
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//取最大值并存AlphaM
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		_mm_store_si128(BetaIndex++, AlphaBetaMain);
	}

	/*再计算前向度量并计算似然比*/
	__m128i*LsIndex = (__m128i*)Ls;
	__m128i*OutputIndex = (__m128i*)Output;
	BetaIndex = BetaM + Len;
	RrIndex = (__m128i*)Ls;
	RpIndex = (__m128i*)Lp;
	ShuffleLshift64Invert = _mm256_setr_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);
	ShuffleAlphaBetaPlus = _mm_setr_epi8(2, 3, 4, 5, 10, 11, 12, 13, 0, 1, 6, 7, 8, 9, 14, 15);
	ShuffleAlphaBetaMinus = _mm_setr_epi8(0, 1, 6, 7, 8, 9, 14, 15, 2, 3, 4, 5, 10, 11, 12, 13);
	AlphaBetaMain = _mm_set_epi16(-INF, -INF, -INF, -INF, -INF, -INF, -INF, 0);
	//初始化Rr
	RL1 = _mm_load_si128(RrIndex++);
	Offset = _mm_load_si128(RpIndex++);
	RL0 = _mm_sub_epi16(RL1, Offset);
	RL1 = _mm_add_epi16(RL1, Offset);
	RL0 = _mm_srai_epi16(RL0, 1);
	RL1 = _mm_srai_epi16(RL1, 1);
	i = 8;
	//开始时处理规整的部分
	for (; i <= Len; i += 8)
	{
		RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理0
		//Offset = _mm_broadcastw_epi16(RL1);
		//AlphaBetaSub = _mm_broadcastw_epi16(RL0);
		//Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
		//RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//Offset = _mm_set_epi16(L1, L1, L0, L0, L0, L0, L1, L1);
		Offset = _mm_blend_epi16(RL1, RL0, 0xFE);//11111110
		RTable = _mm_setr_epi8(0, 1, 0, 1, 14, 15, 14, 15, 14, 15, 14, 15, 0, 1, 0, 1);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R1 = _mm256_max_epi16(Rr, RTmp);

		//处理1
		Offset = _mm_blend_epi16(RL1, RL0, 0x01);//00000001
		RTable = _mm_setr_epi8(2, 3, 2, 3, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R2 = _mm256_max_epi16(Rr, RTmp);

		//处理2
		Offset = _mm_blend_epi16(RL1, RL0, 0x02);//00000010
		RTable = _mm_setr_epi8(4, 5, 4, 5, 2, 3, 2, 3, 2, 3, 2, 3, 4, 5, 4, 5);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R3 = _mm256_max_epi16(Rr, RTmp);

		//处理3
		Offset = _mm_blend_epi16(RL1, RL0, 0x04);//00000100
		RTable = _mm_setr_epi8(6, 7, 6, 7, 4, 5, 4, 5, 4, 5, 4, 5, 6, 7, 6, 7);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R4 = _mm256_max_epi16(Rr, RTmp);

		//处理4
		Offset = _mm_blend_epi16(RL1, RL0, 0x08);//00001000
		RTable = _mm_setr_epi8(8, 9, 8, 9, 6, 7, 6, 7, 6, 7, 6, 7, 8, 9, 8, 9);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R5 = _mm256_max_epi16(Rr, RTmp);

		//处理5
		Offset = _mm_blend_epi16(RL1, RL0, 0x10);//00010000
		RTable = _mm_setr_epi8(10, 11, 10, 11, 8, 9, 8, 9, 8, 9, 8, 9, 10, 11, 10, 11);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R6 = _mm256_max_epi16(Rr, RTmp);

		//处理6
		Offset = _mm_blend_epi16(RL1, RL0, 0x20);//00100000
		RTable = _mm_setr_epi8(12, 13, 12, 13, 10, 11, 10, 11, 10, 11, 10, 11, 12, 13, 12, 13);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		R7 = _mm256_max_epi16(Rr, RTmp);

		//处理7
		Offset = _mm_blend_epi16(RL1, RL0, 0x40);//01000000
		RTable = _mm_setr_epi8(14, 15, 14, 15, 12, 13, 12, 13, 12, 13, 12, 13, 14, 15, 14, 15);
		Offset = _mm_shuffle_epi8(Offset, RTable);
		AlphaBetaSub = AlphaBetaMain;
		//先计算加
		AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
		AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
		//再计算相减
		AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
		AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
		//再加上Beta
		Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
		Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
		Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
		//直接更新下一次
		AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
		//处理Hmax
		Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
		//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
		Rr = _mm256_max_epi16(Rr, RTmp);
		RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
		Rr = _mm256_max_epi16(Rr, RTmp);

		//更新Rr
		RL1 = _mm_load_si128(RrIndex++);
		Offset = _mm_load_si128(RpIndex++);
		RL0 = _mm_sub_epi16(RL1, Offset);
		RL1 = _mm_add_epi16(RL1, Offset);
		RL0 = _mm_srai_epi16(RL0, 1);
		RL1 = _mm_srai_epi16(RL1, 1);
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//归一化
		Offset = _mm_broadcastw_epi16(AlphaBetaMain);
		AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
		//此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去Ls+La
		AlphaBetaSub = _mm_load_si128(LsIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		if ((i & 0x3F) == 0)//纯粹为了减小运行次数
		{
			AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
			Offset = _mm_min_epi16(Offset, AlphaBetaSub);
			AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
			Offset = _mm_max_epi16(Offset, AlphaBetaSub);
		}
		_mm_store_si128(OutputIndex++, Offset);
	}


	if ((Len & 0x7) != 0)//如果不是8的倍数则要多处理一次软输出
	{
		//更新Rr
		//RL1 = _mm_load_si128(--RrIndex);
		//Offset = _mm_load_si128(--RaIndex);
		//RL1 = _mm_add_epi16(RL1, Offset);
		//Offset = _mm_load_si128(--RpIndex);
		//RL0 = _mm_sub_epi16(RL1, Offset);
		//RL1 = _mm_add_epi16(RL1, Offset);
		//RL0 = _mm_srai_epi16(RL0, 1);
		//RL1 = _mm_srai_epi16(RL1, 1);
		//ShuffleLshift16 = _mm256_setr_epi8(14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);
		//RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
		//处理最后的剩余部分
		i -= 8;
		for (; i < Len; ++i)
		{
			Offset = _mm_broadcastw_epi16(RL1);
			AlphaBetaSub = _mm_broadcastw_epi16(RL0);
			Offset = _mm_blend_epi16(Offset, AlphaBetaSub, 0x3C);//00111100
			RL1 = _mm_shuffle_epi8(RL1, _mm256_castsi256_si128(ShuffleLshift16));
			RL0 = _mm_shuffle_epi8(RL0, _mm256_castsi256_si128(ShuffleLshift16));
			AlphaBetaSub = AlphaBetaMain;
			//先计算加
			AlphaBetaMain = _mm_add_epi16(AlphaBetaMain, Offset);//相加
			AlphaBetaMain = _mm_shuffle_epi8(AlphaBetaMain, ShuffleAlphaBetaPlus);//重排
			//再计算相减
			AlphaBetaSub = _mm_sub_epi16(AlphaBetaSub, Offset);
			AlphaBetaSub = _mm_shuffle_epi8(AlphaBetaSub, ShuffleAlphaBetaMinus);
			//再加上Beta
			Offset = _mm_load_si128(--BetaIndex);//加一次为128bit
			Rr = _mm256_add_epi16(_mm256_castsi128_si256(AlphaBetaMain), _mm256_castsi128_si256(Offset)); //!!!!此处用Rr代替Main
			Offset = _mm_add_epi16(AlphaBetaSub, Offset); //此处用Offset代替Sub
			//直接更新下一次
			AlphaBetaMain = _mm_max_epi16(AlphaBetaMain, AlphaBetaSub);
			//处理Hmax
			Rr = _mm256_set_m128i(Offset, _mm256_castsi256_si128(Rr));
			//对于每个寄存器进行三次max和三次shuffle得出Hmax，外加一次相减得出
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift64Invert);
			Rr = _mm256_max_epi16(Rr, RTmp);
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift32);
			Rr = _mm256_max_epi16(Rr, RTmp);
			RTmp = _mm256_shuffle_epi8(Rr, ShuffleLshift16);
			Rr = _mm256_max_epi16(Rr, RTmp);

			switch (i & 0x7)
			{
			case 0:R1 = Rr; break;//低位为Plus
			case 1:R2 = Rr; break;
			case 2:R3 = Rr; break;
			case 3:R4 = Rr; break;
			case 4:R5 = Rr; break;
			case 5:R6 = Rr; break;
			case 6:R7 = Rr; break;
			case 7:
			{
					  //更新Rr
					  RL1 = _mm_load_si128(RrIndex++);
					  Offset = _mm_load_si128(RpIndex++);
					  RL0 = _mm_sub_epi16(RL1, Offset);
					  RL1 = _mm_add_epi16(RL1, Offset);
					  RL0 = _mm_srai_epi16(RL0, 1);
					  RL1 = _mm_srai_epi16(RL1, 1);
					  //归一化
					  Offset = _mm_broadcastw_epi16(AlphaBetaMain);
					  AlphaBetaMain = _mm_sub_epi16(AlphaBetaMain, Offset);
					  //此时进行8个寄存器的blend和存储R1、R2......R7、Rr，
					  R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
					  R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
					  R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
					  R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
					  R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
					  R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
					  R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
					  //Plus与Minus相减,到此Rr的低八位数值全一样
					  Offset = _mm256_extractf128_si256(R1, 1);
					  Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
					  //减去Ls+La
					  AlphaBetaSub = _mm_load_si128(LsIndex++);
					  Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
					  //存入内存
					  // if ((i & 0xFF) == 0xFF)
					  // {
					  //  AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
					  //  Offset = _mm_min_epi16(Offset, AlphaBetaSub);
					  // AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
					  // Offset = _mm_max_epi16(Offset, AlphaBetaSub);
					  //}
					  _mm_store_si128(OutputIndex++, Offset);
			}
			}
		}
		R1 = _mm256_blend_epi16(R1, R2, 0x0202);//00000010
		R1 = _mm256_blend_epi16(R1, R3, 0x0404);//00000100
		R1 = _mm256_blend_epi16(R1, R4, 0x0808);//00001000
		R1 = _mm256_blend_epi16(R1, R5, 0x1010);//00010000
		R1 = _mm256_blend_epi16(R1, R6, 0x2020);//00100000
		R1 = _mm256_blend_epi16(R1, R7, 0x4040);//01000000
		R1 = _mm256_blend_epi16(R1, Rr, 0x8080);//10000000
		//Plus与Minus相减,到此Rr的低八位数值全一样
		Offset = _mm256_extractf128_si256(R1, 1);
		Offset = _mm_sub_epi16(_mm256_castsi256_si128(R1), Offset);
		//减去Ls+La
		AlphaBetaSub = _mm_load_si128(LsIndex++);
		Offset = _mm_sub_epi16(Offset, AlphaBetaSub);
		//存入内存
		AlphaBetaSub = _mm_set_epi16(1024, 1024, 1024, 1024, 1024, 1024, 1024, 1024);
		Offset = _mm_min_epi16(Offset, AlphaBetaSub);
		AlphaBetaSub = _mm_set_epi16(-1024, -1024, -1024, -1024, -1024, -1024, -1024, -1024);
		Offset = _mm_max_epi16(Offset, AlphaBetaSub);
		_mm_store_si128(OutputIndex++, Offset);
	}
	return 1;
}
double TurboDecoderSSEIntSpecifyV2(short*Ls, short*Lp0, short*Lp1, int*Table, short*Output, int Len, int InterNum)
{
	//译码器号_译码器内支路号
	__m128i*AlphaM = (__m128i*)_mm_malloc(Len*sizeof(__m128i), sizeof(__m128i));
	short*La = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m128i));
	memset(Output, 0, Len*sizeof(short));
	//All System Go
	LARGE_INTEGER Start, End, Freq;
	QueryPerformanceCounter(&Start);
	QueryPerformanceFrequency(&Freq);
	for (int i = 1; i <= InterNum; ++i)
	{
		MaxLogRscDecoderSSEIntSpecifyD0(Ls, Lp0, Output, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[k] = La[Table[k]];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = 0;
		}
		MaxLogRscDecoderSSEIntSpecifyD1(Output, Lp1, Len, La, AlphaM);
		for (int k = 0; k < Len - Mem; ++k)
		{
			Output[Table[k]] = La[k];
		}
		for (int k = Len - Mem; k < Len; ++k)
		{
			Output[k] = La[k];
		}
	}
	QueryPerformanceCounter(&End);
	//
	_mm_free(AlphaM);
	_mm_free(La);
	return (double)(End.QuadPart - Start.QuadPart) / (double)Freq.QuadPart;
}
//TestZone
bool HardDecision(float*Input, int Len, short*Output)
{
	for (int i = 0; i < Len; ++i)
		Output[i] = Input[i] >= 0 ? 1 : 0;
	return 1;
}
bool HardDecision(short*Input, int Len, int*Output)
{
	for (int i = 0; i < Len; ++i)
		Output[i] = Input[i] >= 0 ? 1 : 0;
	return 1;
}
bool HardDecision(short*Input, int Len, short*Output)
{
	for (int i = 0; i < Len; ++i)
		Output[i] = Input[i] >= 0 ? 1 : 0;
	return 1;
}
int ErrorNum(int*Input, int *OriginMessage, int Len)
{
	int Sum = 0;
	for (int i = 0; i < Len; ++i)
	if (Input[i] != OriginMessage[i]) ++Sum;
	return Sum;
}
int ErrorNum(int*Input, short *OriginMessage, int Len)
{
	int Sum = 0;
	for (int i = 0; i < Len; ++i)
	if (Input[i] != OriginMessage[i]) ++Sum;
	return Sum;
}
int ErrorNum(short*Input, int *OriginMessage, int Len)
{
	int Sum = 0;
	for (int i = 0; i < Len; ++i)
	if (Input[i] != OriginMessage[i]) ++Sum;
	return Sum;
}
int ErrorNum(short*Input, short *OriginMessage, int Len)
{
	int Sum = 0;
	for (int i = 0; i < Len; ++i)
	if (Input[i] != OriginMessage[i]) ++Sum;
	return Sum;
}
int main()
{
	//多线程测速度
	
	
	int Len = 6144;      //交织长度
	int InterNum = 4;    //迭代次数
	int PNum = 640;      //总Turbo码数
	int TLen = PNum*Len;//总解码比特数
	int**InputMessageArray = (int**)malloc(PNum*sizeof(int*));
	int**OutputMessageArray = (int**)malloc(PNum*sizeof(int*));
	float**ReceivedMessageArray = (float**)malloc(PNum*sizeof(float*));
	float**DivideMessageArray = (float**)malloc(PNum*sizeof(float*));
	float**DivideParity2Array = (float**)malloc(PNum*sizeof(float*));
	float**DivideParity1Array = (float**)malloc(PNum*sizeof(float*));
	float**DivideIMessageArray = (float**)malloc(PNum*sizeof(float*));
	short**DivideMessageSArray = (short**)malloc(PNum*sizeof(short*));
	short**DivideParityS2Array = (short**)malloc(PNum*sizeof(short*));
	short**DivideParityS1Array = (short**)malloc(PNum*sizeof(short*));
	short**DivideIMessageSArray = (short**)malloc(PNum*sizeof(short*));
	short**SoftOutputSArray = (short**)malloc(PNum*sizeof(short*));

	double*BER = (double*)malloc(PNum*sizeof(double));
	int*Table = (int*)malloc((Len + Mem)*sizeof(int));
	ReadTable(Table, Len, "Random");
	Sleep(10);
	//生成信息
    #pragma omp parallel for
	for (int j = 0; j < PNum; ++j)
	{
		InputMessageArray[j] = (int*)malloc(Len*sizeof(int));
		OutputMessageArray[j] = (int*)malloc(3 * (Len + 4)*sizeof(int));
		ReceivedMessageArray[j] = (float*)malloc(3 * (Len + 4)*sizeof(float));
		DivideMessageArray[j] = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m128i));
		DivideParity2Array[j] = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m128i));
		DivideParity1Array[j] = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m128i));
		DivideIMessageArray[j] = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m128i));
		DivideMessageSArray[j] = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
		DivideParityS2Array[j] = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
		DivideParityS1Array[j] = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
		DivideIMessageSArray[j] = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
		SoftOutputSArray[j] = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256i));

		BinaryStreamGenerator(InputMessageArray[j], Len);//产生序列
		MessageTransmit(InputMessageArray[j], Len, OutputMessageArray[j]);//发送Turbo码
		AddGaussianNoise(OutputMessageArray[j], 3 * (Len + 4), 3, ReceivedMessageArray[j]);//加入高斯白噪声
		InputDivide(ReceivedMessageArray[j], 3 * (Len + 4), DivideMessageArray[j], DivideIMessageArray[j], DivideParity1Array[j], DivideParity2Array[j]);//解码率匹配

		//量化
		QuantizeFloatToShort(DivideMessageArray[j], Len + Mem, DivideMessageSArray[j]);
		QuantizeFloatToShort(DivideIMessageArray[j], Len + Mem, DivideIMessageSArray[j]);
		QuantizeFloatToShort(DivideParity2Array[j], Len + Mem, DivideParityS2Array[j]);
		QuantizeFloatToShort(DivideParity1Array[j], Len + Mem, DivideParityS1Array[j]);
	}
	Sleep(10);
	//大家一起来解码
	DWORD Begin = GetTickCount();
	#pragma omp parallel for
	for (int j = 0; j < PNum; ++j)
		TurboDecoderSSEIntSpecify(DivideMessageSArray[j], DivideParityS1Array[j], DivideParityS2Array[j], DivideIMessageSArray[j], Table, SoftOutputSArray[j], Len + Mem, InterNum);
		//TurboDecoderSSEIntSpecifyV2(DivideMessageSArray[j], DivideParityS1Array[j], DivideParityS2Array[j], Table, SoftOutputSArray[j], Len + Mem, InterNum);
	DWORD End = GetTickCount();
	double TimeInMs = End - Begin;
	printf("时间(Ms)：%lf\n", TimeInMs);
	printf("速度(Mbps)：%lf\n", TLen / (TimeInMs / 1000)/(1024*1024));
	for (int j = 0; j < PNum; ++j)
	{
		free(InputMessageArray[j]);
		free(OutputMessageArray[j]);
		free(ReceivedMessageArray[j]);
		_mm_free(DivideMessageArray[j]);
		_mm_free(DivideParity2Array[j]);
		_mm_free(DivideParity1Array[j]);
		_mm_free(DivideIMessageArray[j]);
		_mm_free(DivideMessageSArray[j]);
		_mm_free(DivideParityS2Array[j]);
		_mm_free(DivideParityS1Array[j]);
		_mm_free(DivideIMessageSArray[j]);
		_mm_free(SoftOutputSArray[j]);
	}
	free(InputMessageArray);
	free(OutputMessageArray);
	free(ReceivedMessageArray);
	free(DivideMessageArray);
	free(DivideParity2Array);
	free(DivideParity1Array);
	free(DivideIMessageArray);
	free(DivideMessageSArray);
	free(DivideParityS2Array);
	free(DivideParityS1Array);
	free(DivideIMessageSArray);
	free(SoftOutputSArray);
	free(BER);
	free(Table);
	return 1;
	

	//测误码率

	/*
	int Len = 6144; //交织长度
	int Rate = 3;  //码率 3为1/3,2为1/2
	int InterNum = 4; //迭代次数
	int Cycle = 200;  //帧个数
	int*InputMessage1 = (int*)malloc(Len*sizeof(int));
	int*OutputMessage1 = (int*)malloc(3 * (Len + 4)*sizeof(int));
	float*ReceivedMessage1 = (float*)malloc(3 * (Len + 4)*sizeof(float));

	float*DivideMessage1 = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m256));
	float*DivideParity12 = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m256));
	float*DivideParity11 = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m256));
	float*DivideIMessage1 = (float*)_mm_malloc((Len + Mem)*sizeof(float), sizeof(__m256));
	float*SoftOutput1 = (float*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(float), sizeof(__m256));

	short*DivideMessageS1 = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
	short*DivideParityS12 = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
	short*DivideParityS11 = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
	short*DivideIMessageS1 = (short*)_mm_malloc((Len + Mem)*sizeof(short), sizeof(__m128i));
	short*SoftOutputS1 = (short*)_mm_malloc(8 * (Len / 8 + 1)*sizeof(short), sizeof(__m256i));
	int*Table = (int*)malloc((Len + Mem)*sizeof(int));
	int TLen = Len * Cycle;
	int RateLen = (Len + 4)*Rate;
	double TimeAVX;
	int ErrorNumberAVX1;
	int FrameErrorAVX;
	int ErrorTmpAVX;
	double L_DB = -3.1; //低SNR
	double H_DB = 2;    //最高SNR
	for (double db = H_DB; db >= L_DB; db = db - 0.2)
	{
		TimeAVX = 0;
		FrameErrorAVX = 0;
		ErrorNumberAVX1 = 0;
		for (int i = 0; i < Cycle; ++i)
		{
			BinaryStreamGenerator(InputMessage1, Len);

			MessageTransmit(InputMessage1, Len, OutputMessage1);//输出穿孔
			AddGaussianNoise(OutputMessage1, 3 * (Len + 4), db, ReceivedMessage1);
			//穿孔的元素为0
			for (int i = RateLen; i < 3 * (Len + 4); ++i)
				ReceivedMessage1[i] = 0;
			InputDivide(ReceivedMessage1, 3 * (Len + 4), DivideMessage1, DivideIMessage1, DivideParity11, DivideParity12);
			ReadTable(Table, Len, "Random");
			QuantizeFloatToShort(DivideMessage1, Len + Mem, DivideMessageS1);
			QuantizeFloatToShort(DivideIMessage1, Len + Mem, DivideIMessageS1);
			QuantizeFloatToShort(DivideParity12, Len + Mem, DivideParityS12);
			QuantizeFloatToShort(DivideParity11, Len + Mem, DivideParityS11);
			//TimeAVX += TurboDecoderSSEIntSpecify(DivideMessageS1, DivideParityS11, DivideParityS12, DivideIMessageS1, Table, SoftOutputS1, Len + Mem, InterNum);
			TimeAVX += TurboDecoderSSEIntSpecifyV2(DivideMessageS1, DivideParityS11, DivideParityS12, Table, SoftOutputS1, Len + Mem, InterNum);
			HardDecision(SoftOutputS1, Len, SoftOutputS1);
			ErrorTmpAVX = ErrorNum(SoftOutputS1, InputMessage1, Len);
			if (ErrorTmpAVX > 0) FrameErrorAVX++;
			ErrorNumberAVX1 += ErrorTmpAVX;
		}
		printf("dB：%lfdB\n", db);
		printf("耗时：%lfMbps\n", ((double)TLen / TimeAVX) / (1024 * 1024));
		printf("误比特率：%lf\n", (double)(ErrorNumberAVX1) / (TLen));
		printf("误帧率：%lf\n", (double)(FrameErrorAVX) / (Cycle));
		printf("=========================================\n");
	}
	free(OutputMessage1);
	free(ReceivedMessage1);
	free(InputMessage1);
	_mm_free(DivideMessage1);
	_mm_free(DivideMessageS1);
	_mm_free(SoftOutput1);
	_mm_free(DivideParity11);
	_mm_free(DivideParity12);
	_mm_free(DivideIMessage1);
	_mm_free(DivideParityS11);
	_mm_free(DivideParityS12);
	_mm_free(DivideIMessageS1);
	return 1;
	*/
}
