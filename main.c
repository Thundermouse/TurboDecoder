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
	/***�������岢����ռ�***/
	short int Len; //��Ϣ����
	short int CodeLen; //�鳤��
	short int IterCNum;//��������
	short int Cycle;//Turbo֡����
	FILE*InterleaveTable;
	FILE*MatchTable;
	int*BinaryStream = (int*)_mm_malloc(Len*sizeof(int), sizeof(__m256i));
	int*TransStream = (int*)_mm_malloc(CodeLen*sizeof(int), sizeof(__m256i));
	float*ReceivedStream = (float*)_mm_malloc(CodeLen*sizeof(float), sizeof(__m256));
	int*OutBinaryStream = (int*)_mm_malloc(Len*sizeof(int), sizeof(__m256i));

	float SNR_L;//��������
	float SNR_H;//��������
	float SNR_Pace;//����Ȳ�������
	while (1)
	{
		/***��ʼ������***/
		/***�����������Ϣ���ȣ���֯����������ȷ�Χ������Ȳ�����Turbo֡������***/
		/***����������ش�***/
		/***��֯ӳ�������ƥ��ӳ�������***/

	}
	return 0;
}
