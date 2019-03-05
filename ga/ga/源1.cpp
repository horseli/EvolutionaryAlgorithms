#include<iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
const int population = 100;//种群大小
const int num = 60000;//代数cd
const double cross_possbility = 0.7;//交叉概率
const double muta_possbility = 0.07;//变异概率
const int popu_new = population*(cross_possbility + muta_possbility);//交叉和变异产生的新种群
const int length = 28;//基因编码长度
struct gen//基因结构
{
	unsigned int gens[30][28] = { 0 };//染色体编码
	double inf[30];//实际值
	double suit;//适应度值
};
gen group[population];//种群
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