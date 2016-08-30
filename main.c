#include <stdio.h>
#include <stdlib.h>
#include<malloc.h>
#include <math.h>
#include <mkl.h>
#include <mkl_vml.h> 
#include <windows.h>
#include <xmmintrin.h>
#include <immintrin.h> //AVX
#include "TurboCode.h"
#include "CYHTurbo.h"

int main()
{
	/***参数定义并分配空间***/
	short int Len; //信息长度
	short int CodeLen; //块长度
	short int IterCNum;//迭代次数
	short int Cycle;//Turbo帧个数
	FILE*InterleaveTable;
	FILE*MatchTable;
	int*BinaryStream = (int*)_mm_malloc(Len*sizeof(int), sizeof(__m256i));
	int*TransStream = (int*)_mm_malloc(CodeLen*sizeof(int), sizeof(__m256i));
	float*ReceivedStream = (float*)_mm_malloc(CodeLen*sizeof(float), sizeof(__m256));
	int*OutBinaryStream = (int*)_mm_malloc(Len*sizeof(int), sizeof(__m256i));

	float SNR_L;//最低信噪比
	float SNR_H;//最高信噪比
	float SNR_Pace;//信噪比步进长度
	while (1)
	{
		/***初始化变量***/
		/***输入参数（信息长度，交织次数，信噪比范围，信噪比步进，Turbo帧个数）***/
		/***产生随机比特串***/
		/***交织映射表、速率匹配映射表生成***/

	}
	return 0;
}
