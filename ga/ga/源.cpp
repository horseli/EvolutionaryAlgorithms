#include<iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
//是否要饮用
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
void cross();//交叉
void mutation();//变异
void initialize();//种群初始化
void selection1();//选择方法轮盘赌
void selection2();//选择方法锦标赛
double evaluation(double *x);//适应度函数评价
void randgens(unsigned int w[30][28]);//产生随机的二进制编码
double convert2to10(unsigned int w[28]);//二进制转十进制
//double convert10to2(int x);//十进制转二进制
void cross_two(gen a, gen b,int setcross,int n); //选定2个染色体进行交叉
void muta_one(gen a,int setmuta,int n);//选定的染色体进行变异操作

void selection1()
{
	srand((unsigned)time(NULL));
	int min=group_new[0].suit;
	/*for (int i = 0; i < 177; i++)
	{
		if (group_new[i].suit < min)
			min = group_new[i].suit;
	}*/
	int total=0;
	double suitation[177];//轮盘赌的适应值转换
	for (int i = 0; i < 177; i++)
	{
		suitation[i] = 1 / (group_new[i].suit + 1);
	}
	for (int i = 0; i < 177; i++)
	{
		total += suitation[i];
	}
	int poss[177];
	poss[0] = suitation[0] / total;
	for (int i = 1; i < 177; i++)
	{
		poss[i] = poss[i - 1] + suitation[i];
	}
	for (int j = 0; j < population; j++)
	{
		int temp = rand();
		for (int i = 0; i < 177; i++)
		{
			if (temp < poss[0] && i == 0 || i>0 && temp > poss[i] && temp < poss[i + 1])
			{
				group[j] = group_new[i];
					break;
			}
		}
	}
}
double evaluation(double *x)
{
	double value=0;
	for (int i = 0; i < 30; i++)
	{
		value += x[i] * x[i];
	}
	return value;
}
void muta_one(gen a, int setmuta, int n)
{
	srand((unsigned)time(NULL));
	group_new[n] = a;
	for (int i = 0; i < 30; i++)
	{
		group_new[n].gens[i][setmuta] = (a.gens[i][setmuta] + 1) % 2;
		group_new[n].inf[i] = convert2to10(group_new[n].gens[i]);
	}
	group_new[n].suit = evaluation(group_new[n].inf);
}
void cross_two(gen a, gen b,int setcross,int n)
{
	/*group_new[n] = a;
	group_new[n+1] = b;*/
	for (int j = 0; j < 30; j++)
	{
		for (int i = setcross; i < 28; i++)
		{
			int temp = group_new[n].gens[j][i];
			group_new[n].gens[j][i] = group_new[n + 1].gens[j][i];
			group_new[n + 1].gens[j][i] = temp;
			/*group_new[n].inf = convert2to10(group_new[n].gens);
			group_new[n].suit = evaluation(group_new[n].inf);
			group_new[n + 1].inf = convert2to10(group_new[n].gens);
			group_new[n + 1].suit = evaluation(group_new[n].inf);*/
		}
		
			group_new[n].inf[j] = convert2to10(group_new[n].gens[j]);//交叉后的新值
			group_new[n + 1].inf[j] = convert2to10(group_new[n].gens[j]);
		}
		group_new[n].suit = evaluation(group_new[n].inf);
		group_new[n + 1].suit = evaluation(group_new[n].inf);
	
	//

}
double convert2to10(unsigned int w[28])
{
	double x=0;
	for (int i = 0; i < 28; i++)
	{
		x += w[i] * pow(2, i);
	}
	x = (x - 134217728) / pow(10, 6);
	return x;
}
void randgens(unsigned int w[30][28])
{
	srand((unsigned)time(NULL));
	
		unsigned int a[28];
		for (int j = 0; j < 30; j++)
		{
			for (int i = 0; i < 28; i++)
			{
				if (rand() > 0.5)
					a[i] = 1;
				else a[i] = 0;
			}

			for (int i = 0; i < 28; i++)
			{
				w[j][i] = a[i];//编码赋值
			}
		}
}
void initialize()
{	
	srand((unsigned)time(NULL));
	for (int i = 0; i < population; i++)
	{
	    randgens(group[i].gens);
		for (int j = 0; j < 30; j++)
		{
			group[i].inf[j] = convert2to10(group[i].gens[j]);
		}
		group[i].suit = evaluation(group[i].inf);
	}
}
void cross()
{
	srand((unsigned)time(NULL));
	int num = population*cross_possbility;
	int*p = new int[num];
	for (int i = 0; i < num; i++)
	{
		p[i] = rand() * 100 ;
	}
	int setCross = (rand() * 100)%28;//交叉部位
	for (int i = 0; i < num; i = i + 2)
	{
		cross_two(group[p[i]], group[p[i + 1]],setCross, i+100);		
	}
}
void mutation()
{
	srand((unsigned)time(NULL));
	int num = population*muta_possbility;
	int*p = new int[num];
	for (int i = 0; i < num; i++)
	{
		p[i] = rand() * 100;
	}
	int setMuta = (rand() * 100) % 29;//确定变异位置 28位二进制数
	for (int i = 0; i < num; i++)
	{
		muta_one(group[p[i]], setMuta,i+170);
	}
}

