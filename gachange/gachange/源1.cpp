#include<iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
//是否要饮用
const int population = 100;//种群大小
const int num = 60000;//代数cd
const int Size = 30;//维度 
const double cross_possbility = 0.7;//交叉概率
const double muta_possbility = 0.07;//变异概率
const int popu_new = population*(cross_possbility + muta_possbility);//交叉和变异产生的新种群
//const int length = 24;//基因编码长度
//const int total_len = 720;//总长度 
const double high = 8;
const double low = -8;
struct gen//基因结构
{
	//unsigned int gens[30][24];//染色体编码
	double inf[30];//实际值
	double suit;//适应度值
};
gen group[population];//种群
					  //gen group_new[177];//
void cross();//交叉
void mutation();//变异
void initialize();//种群初始化
void selection1();//选择方法轮盘赌
void selection2();//选择方法锦标赛
double evaluation(double *x);//适应度函数评价
void randgens(unsigned int w[30]);//产生随机的二进制编码
//double convert2to10(unsigned int w[24]);//二进制转十进制
										//double convert10to2(int x);//十进制转二进制
void cross_two(int setcross, int n, int m); //选定2个染色体进行交叉
void muta_one(int i, int j, int k);//选定的染色体进行变异操作

void selection2()
{

}
void selection1()
{

	///	int min=group_new[0].suit;
	/*for (int i = 0; i < 177; i++)
	{
	if (group_new[i].suit < min)
	min = group_new[i].suit;
	}*/
	double total = 0;
	double suitation[population];//轮盘赌的适应值转换
	for (int i = 0; i <population; i++)
	{

		suitation[i] = 1 / (group[i].suit + 1);
		//cout<<group_new[i].suit<<"g"<<suitation[i]<<"d \n";
	}
	for (int i = 0; i < population; i++)
	{
		total += suitation[i];
		//cout<<total<<"s";
	}
	double poss[population];
	poss[0] = suitation[0] / total;
	for (int i = 1; i < population; i++)//累积概率 
	{
		poss[i] = poss[i - 1] + suitation[i] / total;
		//cout<<poss[i]<<" ";//
	}
	gen group_temp[population];
	for (int j = 0; j < population; j++)
	{
		double temp = (double)(rand()) / RAND_MAX;
		//cout<<temp<<"s\n";
		for (int i = 0; i < population; i++)
		{
			//cout<<i<<"k";
			if (temp < poss[0] && i == 0 || i<176 && i>0 && temp > poss[i - 1] && temp < poss[i] || i == 176 && poss[175]<temp)
			{
				group_temp[j] = group[i];
				//cout<<"s"<<i<<"s\n";
				break;
			}

		}
	}
	for (int i = 0; i < population; i++)
	{
		group[i] = group_temp[i];
	}
	//输出最小值 
	double min = group[0].suit;
	int k = 0;
	for (int i = 1; i<100; i++)
	{
		if (min>group[i].suit)
		{
			min = group[i].suit;
			k = i;
		}
	}
	//cout<<min<<" ";
	/*for(int i=0;i<30;i++){

	cout<<group[i].inf[i]<<" ";
	}
	cout<<"\t";*/
	printf("%10.6f", min);
	cout << endl;
}
double evaluation(double *x)
{
	double value = 0;
	for (int i = 0; i < 30; i++)
	{
		value += x[i] * x[i];
	}
	//cout<<value<<"t";//
	return value;
}
void muta_one(int i, int j)
{

	//group_new[n] = a;
	/*for (int i = 0; i < 30; i++)
	{
	group_new[n].gens[i][setmuta] = (a.gens[i][setmuta] + 1) % 2;
	//group_new[n].inf[i] = convert2to10(group_new[n].gens[i]);
	//group_new[n].gens[i][setmuta2] = (a.gens[i][setmuta2] + 1) % 2;
	group_new[n].inf[i] = convert2to10(group_new[n].gens[i]);
	}*/

	group[i].inf[j] = ((double)(rand() % 1000) / 1000.0)*(high - low) + low;
	//group[i].inf[j] = convert2to10(group[i].gens[j]);
	//group[i].suit = evaluation(group[i].inf);
}
void cross_two(int setcross, int n1, int n2)
{
	//group_new[n] = a;
	//group_new[n + 1] = b;
	//int s1 = setcross / length;
	int s2 = setcross%Size;
	for (int i = s2; i < Size; i++)//setcross后的第一列
	{
		int temp = group[n1].inf[s1];
		group[n1].inf[s1] = group[n2].inf[s1];
		group[n2].inf[s1] = temp;
	}
	

}
double convert2to10(unsigned int w[24])
{
	double x = 0;
	for (int i = 0; i < 24; i++)
	{
		x += w[i] * pow(2, 23 - i);
	}
	x = (x - 8388608) / pow(10, 6);
	return x;
}
void randgens(unsigned int w[30][24])
{


	unsigned int a[24];
	for (int j = 0; j < 30; j++)
	{
		for (int i = 0; i < 24; i++)
		{
			if ((double)rand() / RAND_MAX > 0.5)
				a[i] = 1;
			else a[i] = 0;

		}

		for (int i = 0; i < 24; i++)
		{
			w[j][i] = a[i];//编码赋值
		}
	}
}
void initialize()
{

	for (int i = 0; i < population; i++)
	{
		randgens(group[i].gens);
		for (int j = 0; j < 30; j++)
		{
			group[i].inf[j] = convert2to10(group[i].gens[j]);
		}
		group[i].suit = evaluation(group[i].inf);
	}

	/*for(int i=0;i<30;i++)
	{
	group[0].gens[i][0]=0;
	for(int j=1;j<24;j++)
	{
	group[0].gens[i][j]=1;
	}

	}
	for (int j = 0; j < 30; j++)
	{
	group[0].inf[j] = convert2to10(group[0].gens[j]);
	}
	group[0].suit = evaluation(group[0].inf);//

	for (int i = 1; i < population; i++)
	{
	randgens(group[i].gens);
	for (int j = 0; j < 30; j++)
	{
	group[i].inf[j] = convert2to10(group[i].gens[j]);
	}
	group[i].suit = evaluation(group[i].inf);
	}*/
}
void cross()
{

	int num = population*cross_possbility;
	int*p = new int[num];
	int i = 0;
	while (i < num)
	{
		int temp = 0;
		p[i] = rand() % 100;
		/*for(int j=0;j<i;j++)
		{
		if(p[i]==p[j])
		{
		temp=1;
		break;}
		}
		if(temp==1)
		continue;*/

		i++;
	}
	//for(int i=0;i<70;i++)
	//{
	//	cout<<p[i]<<" ";
	//}//
	for (int i = 0; i < num; i = i + 2)
	{
		int setCross = rand() % (length);//交叉部位
		cross_two(setCross, p[i], p[i + 1]);
	}
}
void mutation()
{

	//int num = population*muta_possbility;
	//int*p = new int[num];
	//int temp;
	for (int i = 0; i<population; i++)
	{
		for (int j = 0; j<Size; j++)
		{
			for (int k = 0; k<length; k++)
			{
				if (rand() / RAND_MAX<muta_possbility)
					muta_one(i, j, k);
			}
		}
	}
	/*for (int i = 0; i < num; i++)
	{
	p[i] = rand()%100;
	}

	//cout<<setMuta<<"p";//确定变异位置 28位二进制数
	for (int i = 0; i < num; i++)
	{
	int setMuta1 = (rand()) % (length);
	//int setMuta2 = (rand()) % (length);
	muta_one(group[p[i]], setMuta1,i+170);
	}*/
}
int main()
{

	srand((unsigned)time(NULL));
	//cout<<RAND_MAX<<endl;
	initialize();
	//cout<<"xx"<<endl;
	/*for (int i = 0; i < population; i++)
	{
	for (int j = 0; j < 30; j++)
	{
	for (int k = 0; k < 28; k++)
	{
	cout << group[i].gens[j][k];
	}
	cout << " ";
	}
	cout << endl;
	}
	for(int i=0;i<population;i++)
	{
	group_new[i]=group[i];
	}
	*/
	for (int in = 0; in<num; in++)//60000代 
	{

		/*for (int i = 0; i<population; i++)
		{
		group_new[i] = group[i];
		//cout<<group[i].suit<<" ";
		}*
		cout << group[0].suit << "ss";
		for (int i = 0; i<177; i++)
		{
		if (i == 100)
		cout << "@@@";
		cout << group[i].suit << " ";
		}*/
		cross();
		//cout<<"xx"<<endl;
		mutation();
		//cout<<"yy\n";

		//cout<<"xx"<<endl;
		//printf("%.6f",convert2to10(group[1].gens[1]));
		selection1();
		//printf("%f.6",6.6666666);
	}
	system("pause");
}