#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>
#include<time.h>
using namespace std;

#define POPSIZE 100               /* population size */
#define MAXGENS 60000             /* max. number of generations */
#define NVARS 30                 /* no. of problem variables */
#define PXOVER 0.7               /* probability of crossover */
#define PMUTATION 0.07           /* probability of mutation */
#define TRUE 1
#define FALSE 0

int generation;                  /* current generation no. */
int cur_best;

struct genotype /* genotype (GT), a member of the population */
{
	double gene[NVARS];        /* a string of variables */
	double fitness;            /* GT's fitness */
	double upper[NVARS];       /* GT's variables upper bound */
	double lower[NVARS];       /* GT's variables lower bound */
	double rfitness;           /* relative fitness */
	double cfitness;           /* cumulative fitness */
};

struct genotype population[POPSIZE + 1];    /* population */
struct genotype newpopulation[POPSIZE + 1]; /* new population; */
void initialize(void);
double randval(double, double);
void evaluate(void);
void keep_the_best(void);
void elitist(void);
void select(void);
void crossover(void);
void Xover(int, int);
void swap(double *, double *);
void mutate(void);
void report(void);


void initialize(void)
{
	double lbound = -8;//上下界
	double ubound = 8;
	for (int j = 0; j < POPSIZE; j++)
	{

		population[j].fitness = 0;
		population[j].rfitness = 0;
		population[j].cfitness = 0;
		for (int i = 0; i < NVARS; i++)
		{
			population[j].lower[i] = lbound;
			population[j].upper[i] = ubound;
			population[j].gene[i] = randval(population[j].lower[i],
				population[j].upper[i]);
			//cout<<population[j].gene[i]<<"##";///
		}
		//cout<<population[j].fitness<<"%";
	}
}

double randval(double low, double high)
{
	double val;
	val = ((double)(rand() % 1000) / 1000.0)*(high - low) + low;
	return(val);
}

void evaluate(void)
{
	for (int i = 0; i < POPSIZE; i++)
	{
		for (int j = 0; j < NVARS; j++)
		{
			population[i].fitness += population[i].gene[j] * population[i].gene[j];
		}
		//cout<<population[i].fitness<<"%";
	}
}

void keep_the_best()//找到最优解 
{
	int mem;
	int i;
	cur_best = 0; /* stores the index of the best individual */
	population[POPSIZE].fitness = population[cur_best].fitness;
	for (mem = 0; mem < POPSIZE; mem++)
	{
		if (population[mem].fitness < population[POPSIZE].fitness)//找最小值
		{
			cur_best = mem;
			population[POPSIZE].fitness = population[mem].fitness;
		}
	}
	/* once the best member in the population is found, copy the genes */
	for (i = 0; i < NVARS; i++)
		population[POPSIZE].gene[i] = population[cur_best].gene[i];
}

void elitist()
{
	int i;
	double best, worst;             // best(max) and worst(min) fitness values 
	int best_mem, worst_mem;		// indexes of the best and worst member 

	best = population[0].fitness;
	worst = population[0].fitness;
	best_mem = worst_mem = 0;
	for (i = 1; i<POPSIZE; i++) {
		if (population[i].fitness>best) {
			best = population[i].fitness;
			best_mem = i;
			//cout<<"***"<<best;
		}
		if (population[i].fitness<worst) {
			worst = population[i].fitness;
			worst_mem = i;
			//cout<<"^^^"<<worst;
		}
		//cout<<"^^^"<<population[i].fitness;
	}
	if (worst < population[POPSIZE].fitness)
	{
		//cout<<"^^^^^^^";
		population[POPSIZE] = population[worst_mem];//替换最小值(min)
		for (int j = 0; j < NVARS; j++)
		{
			population[POPSIZE].gene[j] = population[worst_mem].gene[j];
		}
	}
	else {
		population[best_mem] = population[POPSIZE];
		for (int j = 0; j < NVARS; j++)
		{
			population[best_mem].gene[j] = population[POPSIZE].gene[j];
		}
	}
}

void select(void)
{
	int mem, i, j, k;
	double sum = 0;
	double p;

	/* find total fitness of the population */
	for (mem = 0; mem < POPSIZE; mem++)
	{
		sum += 1 / (population[mem].fitness + 1);
	}

	/* calculate relative fitness */
	for (mem = 0; mem < POPSIZE; mem++)
	{
		population[mem].rfitness = (1 / (population[mem].fitness + 1)) / sum;

	}
	population[0].cfitness = population[0].rfitness;

	/* calculate cumulative fitness */
	for (mem = 1; mem < POPSIZE; mem++)
	{
		population[mem].cfitness = population[mem - 1].cfitness +
			population[mem].rfitness;
		//cout<<population[mem].cfitness<<"&&";
	}

	/* finally select survivors using cumulative fitness. */

	for (i = 0; i < POPSIZE; i++)
	{
		p = rand() % 1000 / 1000.0;
		//cout<<p<<"$";
		if (p < population[0].cfitness)
			newpopulation[i] = population[0];
		else
		{
			for (j = 0; j < POPSIZE; j++)
				if (p >= population[j].cfitness &&
					p<population[j + 1].cfitness)
				{
					newpopulation[i] = population[j + 1];
					//cout<<j+1<<"@@@";
				}
		}
	}
	/* once a new population is created, copy it back */

	for (i = 0; i < POPSIZE; i++)
		population[i] = newpopulation[i];
}

void crossover(void)
{
	int i, mem, one;
	int first = 0; /* count of the number of members chosen */
	double x;

	for (mem = 0; mem < POPSIZE; ++mem)
	{
		x = rand() % 1000 / 1000.0;
		if (x < PXOVER)
		{
			++first;
			if (first % 2 == 0)//找到两个交配一次 
				Xover(one, mem);
			else
				one = mem;
		}
	}
}

void Xover(int one, int two)
{
	int i;
	int point; /* crossover point */

			   /* select crossover point */
	if (NVARS > 1)//数据维度大小分类讨论 
	{
		if (NVARS == 2)
			point = 1;
		else
			point = (rand() % (NVARS - 1)) + 1;//至少有一个可以交换

		for (i = 0; i < point; i++)
			swap(&population[one].gene[i], &population[two].gene[i]);

	}
}

void swap(double*x, double*y)
{
	double temp = *x;
	*x = *y;
	*y = temp;
}

void mutate(void)
{
	double temp;
	for (int i = 0; i < POPSIZE; i++)
	{
		for (int j = 0; j < NVARS; j++)
		{
			temp = rand() % 1000 / 1000.0;
			//cout<<temp<<"&";
			if (temp < PMUTATION)
			{
				population[i].gene[j] = randval(population[i].lower[j], population[i].upper[j]);
			}
		}
	}
}

void report(void)
{
	int i;
	double best_val;            /* best population fitness */
	double avg;                 /* avg population fitness */
	double stddev;              /* std. deviation of population fitness */
	double sum_square;          /* sum of square for std. calc */
	double square_sum;          /* square of sum for std. calc */
	double sum;                 /* total population fitness */

	sum = 0.0;
	sum_square = 0.0;

	for (i = 0; i < POPSIZE; i++)
	{
		sum += population[i].fitness;
		sum_square += population[i].fitness * population[i].fitness;
	}

	avg = sum / (double)POPSIZE;
	square_sum = avg * avg * POPSIZE;
	stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
	best_val = population[POPSIZE].fitness;
}

int main()
{
	srand(time(NULL));
	generation = 0;
	initialize();
	evaluate();
	keep_the_best();
	while (generation<MAXGENS)
	{
		generation++;
		select();             //选择操作
		crossover();       //交配操作
		mutate();            //变异操作
							 //report();
		evaluate();
		elitist();
		cout << population[POPSIZE].fitness << endl;
	}

}