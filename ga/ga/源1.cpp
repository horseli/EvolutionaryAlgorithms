#include<iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
const int population = 100;//��Ⱥ��С
const int num = 60000;//����cd
const double cross_possbility = 0.7;//�������
const double muta_possbility = 0.07;//�������
const int popu_new = population*(cross_possbility + muta_possbility);//����ͱ������������Ⱥ
const int length = 28;//������볤��
struct gen//����ṹ
{
	unsigned int gens[30][28] = { 0 };//Ⱦɫ�����
	double inf[30];//ʵ��ֵ
	double suit;//��Ӧ��ֵ
};
gen group[population];//��Ⱥ
gen group_new[177] = { 0 };//
int w[10];
void f(int *w)
{
	for (int i = 0; i<10; i++)
	{
		if (i<5)
			w[i] = 0;
		else w[i] = 1;
		cout << w[i];
	}
}
int main()
{
	int n[10] = { 0 };
	f(n);
	for (int i = 0; i<10; i++)
	{
		cout << w[i] << endl;
	}
	system("pause");
}