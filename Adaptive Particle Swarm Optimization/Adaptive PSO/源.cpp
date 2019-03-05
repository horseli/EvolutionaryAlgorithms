#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>
#include<time.h>
#include<fstream>
#include <iomanip>
using namespace std;


#define NVARS 30      //维度
#define GROUP 30
#define MAXGENS 3000

const double high = 100;
const double low = -100;
const double vhigh = 15;
const double vlow = -15;
const double learnRateMax = 1.0;
const double learnRateMin = 0.1;
const double FES = 200000;
double C1 = 2;
double C2 = 2;
double w = 0.9;
double f,s1f,s2f,s3f, s4f;
int period = 1;//阶段1,2,3,4
double acceleRate;
double learningRate;
int g;
double d[GROUP];
double dmin, dmax;
double dg;
int fes;

struct particle {
	double v[NVARS];
	double x[NVARS];

	double fitness;
	double xpbest[NVARS];
	double fitbest;
};

struct particle population[GROUP + 1];

void initialize();
double randval(double, double);
void evaluate();
void evaluate_one(int);
void updatevx();
double meanDistance(int n);   /*计算平均距离d*/
double evaFactor();/*evalu factor f*/
void classify();/*classify into the four moods*/
void parameterControl();/*w,c1,c2的改变*/
double gaussrand(double mu, double deta);/*高斯扰动用于精英学习*/
void elitistLearning();

void initialize()
{
	for (int i = 0; i < GROUP; i++)
	{
		population[i].fitness = 0;
		for (int j = 0; j < NVARS; j++)
		{
			/*initailize gbest*/
			population[i].x[j] = randval(high, low);
			population[i].xpbest[j] = population[i].x[j];
			population[i].v[j] = randval(vhigh, vlow);
			population[i].fitness += pow(population[i].x[j], 2);
		}
		population[i].fitbest = population[i].fitness;
	}
	/*find gbest*/
	int min = 0;
	for (int i = 1; i < GROUP; i++)
	{
		if (population[min].fitness > population[i].fitness)
			min = i;
	}
	population[GROUP].fitness = population[min].fitness;
	for (int i = 0; i < NVARS; i++)
	{
		population[GROUP].x[i] = population[min].x[i];
	}
}

double randval(double low, double high)
{

	double val;
	val = ((double)(rand() % 1000) / 1000.0)*(high - low) + low;
	val = low + (high - low)*rand()*1.0 / RAND_MAX;//这个更好
	return(val);
}
void evaluate_one(int i) {
	population[i].fitness = 0;
	for (int j = 0; j < NVARS; j++)
	{
		population[i].fitness += pow(population[i].x[j], 2);
	}
	fes++;
}
void evaluate()
{

	int min = GROUP;
	for (int i = 0; i < GROUP; i++)
	{
		/*fitness evalutate*/
		
		evaluate_one(i);
		/*update pbest*/
		if (population[i].fitness < population[i].fitbest)
		{
			population[i].fitbest = population[i].fitness;
			for (int k = 0; k < NVARS; k++)
			{
				population[i].xpbest[k] = population[i].x[k];
			}
		}
		/*find gbest*/
		if (population[i].fitness < population[GROUP].fitness)
			min = i;
	}
	/*update gbest*/
	if (min != GROUP)
	{
		population[GROUP].fitness = population[min].fitness;
		for (int k = 0; k < NVARS; k++)
		{
			population[GROUP].x[k] = population[min].x[k];
			population[GROUP].v[k] = population[min].v[k];
			population[GROUP].xpbest[k] = population[min].xpbest[k];
		}
	}
}

void updatevx()
{
	double randr1, randr2;
	for (int i = 0; i < GROUP; i++)
	{
		for (int j = 0; j < NVARS; j++)
		{
			randr1 = rand() % 1000 / 1000.0;
			randr2 = rand() % 1000 / 1000.0;
			/*update velocity*/
			population[i].v[j] = population[i].v[j] * w + C1*randr1*(population[i].xpbest[j] - population[i].x[j]) +
				C2*randr2*(population[GROUP].x[j] - population[i].x[j]);
			/*if cross the border*/
			if (population[i].v[j] > vhigh)
				population[i].v[j] = vhigh;
			if (population[i].v[j] <vlow)
				population[i].v[j] = vlow;
			/*update x*/
			population[i].x[j] += population[i].v[j];
			if (population[i].x[j] > high)
				population[i].x[j] = high;
			if (population[i].x[j] <low)
				population[i].x[j] = low;
		}
	}
}

void methodrandw()
{
	srand(time(NULL));

	initialize();
	for (int i = 0; i < MAXGENS; i++)
	{
		w = randval(0.4, 0.6);
		updatevx();
		evaluate();

		//cout << i << " " << population[GROUP].fitness << endl;
	}
	cout << population[GROUP].fitness << endl;
}

void methodLDM()
{
	srand(time(NULL));

	initialize();
	double wmax = 0.9;
	double wmin = 0.4;
	for (int i = 0; i < MAXGENS; i++)
	{
		w = wmax - (wmax - wmin) / MAXGENS*i;
		updatevx();
		evaluate();

		//cout << i << " " << population[GROUP].fitness << endl;
	}
	cout << population[GROUP].fitness << endl;
}
void methodConcave()
{
	srand(time(NULL));

	initialize();
	double wmax = 0.9;
	double wmin = 0.4;
	for (int i = 0; i < MAXGENS; i++)
	{
		w = wmax - (wmax - wmin)* pow((double(i) / MAXGENS), 2);
		updatevx();
		evaluate();

		//cout << i << " " << population[GROUP].fitness << endl;
	}
	cout << population[GROUP].fitness << endl;
}

void methodConvex()
{
	srand(time(NULL));

	initialize();
	double wmax = 0.9;
	double wmin = 0.4;
	for (int i = 0; i < MAXGENS; i++)
	{
		w = wmin + (wmax - wmin)* pow((double(i) / MAXGENS) - 1, 2);
		updatevx();
		evaluate();
	}
	cout << population[GROUP].fitness << endl;
}

double meanDistance(int n) {
	double d=0;
	for (int i = 0; i < GROUP; i++) {
		double xd = 0;
		if (i == n)
			continue;
		for (int j = 0; j < NVARS; j++) {
			
			xd += pow(population[n].x[j] - population[i].x[j], 2);

		}
		d += sqrt(xd);
	}
	d = d / (GROUP - 1);
	return d;
}
double evaFactor() {
	
	
	for (int i = 0; i < GROUP; i++)
	{
		d[i] = meanDistance(i);
		if (i == 0)
			dmax = dmin = d[i];/*初始化最大和最小值*/
		if (d[i] > dmax)
			dmax = d[i];
		if (d[i] < dmin)
			dmin = d[i];
	}
	dg = meanDistance(GROUP);
	if (dmin > dg)
		dmin = dg;
	f = (dg - dmin) / (dmax - dmin);
	return f;
}
/*singleton and rulebase to fuzzclassfication*/
/*确定的分类 直接分  不确定的分类考虑前后*/
void classify() {
	if (f < 0.2)
		period = 3;
	else if (f < 0.3) {    //2,3
		if (period == 1)
			period = 2;//前一阶段为1时改变，否则不变
		/*else if (f < 0.5 / 3)
			period = 3;
		else period = 2;*/
		else if (period == 4)
			period = 2;
	}
	else if (f < 0.4) {
		period = 2;
	}
	else if (f < 0.6) {    //2,1
		if (period == 4)
			period = 1;
		/*else if (f < 0.5)
			period = 2;
		else period = 1;*/
		else if (period == 3)
			period = 1;
	}
	else if (f < 0.7) {
		period = 1;
	}
	else if (f < 0.8) {    //1,4
		if (period == 3)
			period = 4;
		/*else if (f < 11.5 / 15)
			period = 1;
		else period = 4;*/
		else if (period == 2)
			period = 4;
	}
	else period = 4;
		
}
/*w and c1 c2*/
void parameterControl() {
	/* w */
	w = 1 / (1 + 1.5*exp(-2.6*f));
	/* c1 and c2 */
	acceleRate = randval(0.05, 0.1);
	double increMent;
	switch (period)
	{
	case 1:
		increMent = randval(0, acceleRate);
		C1 += increMent;
		//increMent = randval(0, acceleRate);
		C2 -= increMent;
		break;
	case 2:
		increMent = randval(0, acceleRate/2);
		C1 += increMent;
		//increMent = randval(0, acceleRate / 2);
		C2 -=increMent;
		break;
	case 3:
		increMent = randval(0, acceleRate / 2);
		C1 += increMent;
		//increMent = randval(0, acceleRate / 2);
		C2 += increMent;
		break;
	case 4:
		increMent = randval(0, acceleRate);
		C1 -= increMent;
		//increMent = randval(0, acceleRate);
		C2 += increMent;
		break;
	default:
		break;
	}
	/*单个越界操作*/
	if (C1 > 2.5)
		C1 = 2.5;
	if (C1 < 1.5)
		C1 = 1.5;
	if (C2 > 2.5)
		C2 = 2.5;
	if (C2 < 1.5)
		C2 = 1.5;
	/*和越界操作*/
	if (C1 + C2 > 4 || C1 + C2 < 3) {
		C1 = 4 * C1 / (C1 + C2);
		C2 = 4 * C2 / (C1 + C2);
	}
}
double gaussrand(double mu,double deta)  //期望为E，方差为V
{
	double f;
	double x = randval(low, high);
	//cout << x << endl;
	f = 1 / (sqrt(2 * 3.1415926)*deta)*exp(-pow((x - mu), 2) / (2 * pow(deta, 2)));
	return f;
}

void elitistLearning() {
	int d = abs(rand()) % NVARS;
	learningRate = learnRateMax - (learnRateMax - learnRateMin)*double(fes) / FES;
	double Pd = population[GROUP].x[d];
	double gbest = population[GROUP].fitness;
	/*最大和最小值取问题的限制*/
	double xmax = high;
	double xmin = low;
	/*double xmax = population[0].x[d];
	double xmin = xmax;
	for (int i = 1; i < GROUP; i++) {
		if (xmax < population[i].x[d])
			xmax = population[i].x[d];
		if (xmin > population[i].x[d])
			xmin = population[i].x[d];
	}*/
	population[GROUP].x[d] = population[GROUP].x[d] + (xmax - xmin)*gaussrand(0, learningRate);
	evaluate_one(GROUP);
	/*not better*/
	int m=0;
	for (int i = 1; i < GROUP; i++) {
		if (population[m].fitness < population[i].fitness)
			m = i;
	}
	if (population[GROUP].fitness > gbest)
	{
		
		population[m] = population[GROUP];
		

		population[GROUP].fitness = gbest;		
		population[GROUP].x[d] = Pd;
		population[m].fitbest = population[GROUP].fitness;
	}
}
void apso() {
	for (fes=0; fes < FES; fes++) {
		evaFactor();
		classify();
		parameterControl();
		updatevx();
		elitistLearning();
		evaluate();
	}
	cout << population[GROUP].fitness << endl;
}
void pso() {
	for (fes = 0; fes < FES; fes++) {
		/*evaFactor();
		classify();
		parameterControl();*/
		//w = randval(0.4, 0.6);
		w = 0.5;
		updatevx();
		//elitistLearning();
		evaluate();
	}
	cout << population[GROUP].fitness << endl;
}
int main()
{


	srand(time(NULL));
	
	for (int i = 0; i < 30; i++) {
		initialize();
		int fes = 0;
		//acceleRate = randval(0.05, 0.1);
		apso();
		//for (; fes < FES; fes++) {
		//	//methodrandw();
		//	evaFactor();
		//	classify();
		//	parameterControl();
		//	updatevx();
		//	//elitistLearning();
		//	evaluate();
		//	//cout << population[GROUP].fitness << endl;
		//	//cout << setw(10) << population[GROUP].fitness <<setw(10)<<w << setw(10) << C1  <<setw(10)<< C2 <<setw(15)<<dg<<setw(15)<<f<< setw(10) <<period<< endl;
		//}
		//cout << population[GROUP].fitness << endl;
	}
	//methodrandw();
	/*for (int i = 0; i < 100; i++)
	{
		methodLDM();
	}
	for (int i = 0; i < 100; i++)
	{
		methodrandw();
	}
	for (int i = 0; i < 100; i++)
	{
		methodConcave();
	}
	for (int i = 0; i < 100; i++)
	{
		methodConvex();
	}*/
	/*for (int i = 0; i < 100; i++)
	{
		cout << gaussrand(0, 1) << endl;
	}*/
	system("pause");
}