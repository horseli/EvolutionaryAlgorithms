#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>
#include<time.h>
#include<fstream>
using namespace std;

#define C1 2
#define C2 2
#define NVARS 10      //维度
#define GROUP 30
#define MAXGENS 3000

const double high=100;
const double low=-100;
const double vhigh = 15;
const double vlow = -15;
double w=0.9;

struct particle {
	double v[NVARS];
	double x[NVARS];

	double fitness;
	double xpbest[NVARS];
	double fitbest;
};

struct particle population[GROUP+1];

void initialize();
double randval(double, double);
void evaluate();
void evaluate_one(particle);
void updatevx();

void initialize()
{
	for (int i = 0; i < GROUP; i++)
	{
		population[i].fitness = 0;
		for (int j = 0; j < NVARS; j++)
		{
			/*initailize gbest*/
			population [i].x[j] = randval(high, low);
			population[i].xpbest[j] = population[i].x[j];
			population[i].v[j] = randval(vhigh,vlow);
			population[i].fitness += pow(population[i].x[j],2);
		}
		population[i].fitbest=population[i].fitness;
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

void evaluate()
{
	
	int min=GROUP;
	for (int i = 0; i < GROUP; i++)
	{
		/*fitness evalutate*/
		population[i].fitness = 0;
		for (int j = 0; j < NVARS; j++)
		{
			population[i].fitness += pow(population[i].x[j], 2);
		}
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
			population[GROUP].x[k] = population[min] .x[k];
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
			randr1 = rand()%1000 / 1000.0;
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

		cout << i << " " << population[GROUP].fitness << endl;
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

		cout << i << " " << population[GROUP].fitness << endl;
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
		w = wmax - (wmax - wmin)* pow((double(i)/ MAXGENS),2);
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
		w = wmin +(wmax - wmin)* pow((double(i) / MAXGENS)-1, 2);
		updatevx();
		evaluate();
	}
	cout << population[GROUP].fitness << endl;
}

int main()
{
	//methodLDM();
	methodrandw();
	/*for (int i = 0; i < 100; i++)
	{
		methodLDM();*/
	/*}
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
	system("pause");
}