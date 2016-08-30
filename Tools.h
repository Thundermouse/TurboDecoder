#include<stdio.h>
#include<malloc.h>
int IsTable(int*table, int k)
{
	int*test = (int*)malloc(k*sizeof(int));
	int flag = 1;
	for (int i = 0; i < k; ++i)
		test[i] = 0;
	for (int i = 0; i < k; ++i)
		test[table[i]] = 1;
	for (int i = 0; i < k;++i)
		if (test[i] == 0)
		{
			flag = 0;
			break;
		}
	return flag;
}
int CmpVector(int*table1, int*table2, int k)
{
	int flag = 1;
	for (int i = 0; i < k;++i)
		if (table1[i] != table2[i])
		{
		flag = 0;
		break;
		}
	return flag;
}