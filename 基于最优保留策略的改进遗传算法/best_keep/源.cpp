//第二次修改 求最大值
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>
#include<time.h>
#include<fstream>
using namespace std;

#define POPSIZE 50               /* population size */
#define MAXGENS 30000             /* max. number of generations */
#define NVARS 30                 /* no. of problem variables */
#define PXOVER 0.9               /* probability of crossover */
#define PMUTATION 0.09           /* probability of mutation */
#define TRUE 1
#define FALSE 0

int generation;                  /* current generation no. */
int cur_best;
int Lflag = 0;
int f = 0;
int besttime = 0;
double mmin[NVARS];
double mmax[NVARS];

ofstream file;


struct genotype /* genotype (GT), a member of the population */
{
	double gene[NVARS];        /* a string of variables */
	double fitness;            /* GT's fitness */
	double upper[NVARS];       /* GT's variables upper bound */
	double lower[NVARS];       /* GT's variables lower bound */
	double rfitness;           /* relative fitness */
	double cfitness;           /* cumulative fitness */
	double fit;//求最小值时 fitness的转换
};
struct genotype population[POPSIZE + 1];    /* population */
struct genotype newpopulation[POPSIZE + 1]; /* new population; */
void initialize();
double randval(double, double);
void evaluate();
void keep_the_best();
void elitist();
void select();
void crossover();
void Xover(int, int);
void swap(double *, double *);
void mutate();
void report();
double find_min(int);
double find_max(int);
double evaluate_best(double *x);

void initialize()
{
	int i, j;
	double lbound, ubound;
	remove("output.dat");
	/* initialize variables within the bounds */
	for (i = 0; i < NVARS; i++) {
		lbound = -2;
		ubound = 2;
		for (j = 0; j <= POPSIZE; j++) {
			population[j].fit = 0;
			population[j].rfitness = 0;
			population[j].cfitness = 0;
			population[j].lower[i] = lbound;
			population[j].upper[i] = ubound;
			population[j].gene[i] = randval(population[j].lower[i],
				population[j].upper[i]);
		}
	}
}
double randval(double low, double high)
{
	double val;
	val = ((double)(rand() % 1000) / 1000.0)*(high - low) + low;
	val = low + (high - low)*rand()*1.0 / RAND_MAX;//这个更好
	return(val);
}

void keep_the_best()
{
	int mem; int i;
	cur_best = 0; /* stores the index of the best individual */
	population[POPSIZE].fitness = population[cur_best].fitness;
	for (mem = 0; mem < POPSIZE; mem++)
	{
		if (population[mem].fitness > population[POPSIZE].fitness)
		{
			cur_best = mem;
			population[POPSIZE].fitness = population[mem].fitness;
			population[POPSIZE].fit = population[mem].fit;
		}
	}
	/* once the best member in the population is found, copy the genes */
	cout << cur_best << " ";//
	for (i = 0; i < NVARS; i++)
		population[POPSIZE].gene[i] = population[cur_best].gene[i];
}

void evaluate()
{
	int mem;
	int i;
	double x[NVARS + 1];
	for (mem = 0; mem < POPSIZE; mem++)
	{
		population[mem].fit = 0;
		for (i = 0; i < NVARS; i++)
		{
			//x[i + 1] = population[mem].gene[i];
			population[mem].fit += population[mem].gene[i] * population[mem].gene[i];
		}  //根据

		/*for (int i = 0; i < NVARS; i++)
		{
			int m = population[mem].gene[i] + 0.5;
			population[mem].fit += m*m;
		}*/

		/*double y1 = 0;
		double y2 = 1;
		for (int i = 0; i < NVARS; i++)
		{
			y1 += abs(population[mem].gene[i]);
			y2 *= abs(population[mem].gene[i]);
		}
		population[mem].fit = y1 + y2;*/
		population[mem].fitness = 1 / (population[mem].fit + 1);
		//cout << population[mem].fitness << "&&";
	}

}

void elitist()
{
	int i;  double best, worst;             /* best and worst fitness values */
	int best_mem, worst_mem; /* indexes of the best and worst member */
	best = population[0].fitness;  worst = population[0].fitness;
	if (population[POPSIZE].fitness == 0)
		population[POPSIZE].fitness = best;

	for (i = 0; i < POPSIZE - 1; ++i)
	{
		if (population[i].fitness > population[i + 1].fitness)
		{
			if (population[i].fitness >= best)   //越大越好
			{
				best = population[i].fitness;
				best_mem = i;
			}
			if (population[i + 1].fitness <= worst)
			{
				worst = population[i + 1].fitness;
				worst_mem = i + 1;
			}
		}
		else
		{
			if (population[i].fitness <= worst)
			{
				worst = population[i].fitness;
				worst_mem = i;
			}
			if (population[i + 1].fitness >= best)
			{
				best = population[i + 1].fitness;
				best_mem = i + 1;
			}
		}
	}
	if (best >= population[POPSIZE].fitness)
	{
		for (i = 0; i < NVARS; i++)
			population[POPSIZE].gene[i] = population[best_mem].gene[i];
		population[POPSIZE].fitness = population[best_mem].fitness;
		population[POPSIZE].fit = population[best_mem].fit;
	}
	else
	{
		for (i = 0; i < NVARS; i++)
			population[worst_mem].gene[i] = population[POPSIZE].gene[i];
		population[worst_mem].fitness = population[POPSIZE].fitness;
		population[worst_mem].fit = population[POPSIZE].fit;
	}
}

int findworst()
{
	int t=2;
	for (int i = 3; i < POPSIZE; i++)
	{
		if (population[i].fitness < population[t].fitness)
			t = i;
	}
	return t;
}

void select()//最优保留改进的选择
{
	int mem, i, j, k; double sum = 0;  double p;
	/* find total fitness of the population */
	for (mem = 2; mem < POPSIZE - 1; mem++)
		sum += population[mem].fitness;

	/* calculate relative fitness */
	for (mem = 2; mem < POPSIZE - 1; mem++)
		population[mem].rfitness = population[mem].fitness / sum;

	population[2].cfitness = population[2].rfitness;
	/* calculate cumulative fitness */
	for (mem = 3; mem < POPSIZE - 1; mem++)
	{
		population[mem].cfitness = population[mem - 1].cfitness +
			population[mem].rfitness;
	}
	//找到最小值
	int min = findworst();

	/* finally select survivors using cumulative fitness. */
	//for (i = 2; i < POPSIZE - 1; i++)
	//{
	//	p = rand() % 1000 / 1000.0; // p=rand()*1.0/RAND_MAX;
	//	if (p < population[2].cfitness)
	//		newpopulation[i] = population[2];
	//	else
	//	{
	//		for (j = 2; j < POPSIZE - 1; j++)
	//			if (p >= population[j].cfitness &&
	//				p < population[j + 1].cfitness)
	//			{
	//				
	//				newpopulation[i] = population[j + 1];
	//			}
	//		if ((j + 1) == min)//放弃重选
	//			continue;
	//	}

	//}
	i = 2;
	while (i < POPSIZE - 1)
	{
		p = rand() % 1000 / 1000.0; // p=rand()*1.0/RAND_MAX;
		if (p < population[2].cfitness)
			newpopulation[i] = population[2];
		else
		{
			for (j = 2; j < POPSIZE - 1; j++)
				if (p >= population[j].cfitness &&
					p < population[j + 1].cfitness)
				{

					newpopulation[i] = population[j + 1];
				}
			if ((j + 1) == min)//放弃重选
				continue;
		}
		i++;
	}
	/* once a new population is created, copy it back */

	for (i = 2; i < POPSIZE - 1; i++)
		population[i] = newpopulation[i];
	/*keep the best in 1 2 and n*/
	population[0] = population[POPSIZE];
	population[1] = population[POPSIZE];
	population[POPSIZE - 1] = population[POPSIZE];
	//记录之前的最大最小值
	for (int a = 0; a < NVARS; a++)
	{
		mmin[a] = find_min(a);
		//cout << mmin[a]<<endl;
	}
	for (int a = 0; a < NVARS; a++)
		mmax[a] = find_max(a);
}

void crossover(void)//改进后的变异
{
	int i, mem, one;
	int first = 0; /* count of the number of members chosen */
	double x;

	for (mem = 2; mem < POPSIZE; ++mem)
	{
		x = rand() % 1000 / 1000.0;
		if (x < PXOVER)
		{
			++first;
			if (first % 2 == 0)
				Xover(one, mem);
			else
				one = mem;
		}
	}
	/*the second one must take part in the crossover*/
	one = rand() % 99 + 1;
	mem = 1;
	Xover(one, mem);
}

void Xover(int one, int two)
{
	int i;
	int point; /* crossover point */

			   /* select crossover point */
	if (NVARS > 1)
	{
		if (NVARS == 2)
			point = 1;
		else
			point = (rand() % (NVARS - 1)) + 1;

		for (i = 0; i < point; i++)
			swap(&population[one].gene[i], &population[two].gene[i]);
	}
}

void swap(double *x, double *y)
{
	double temp;

	temp = *x;
	*x = *y;
	*y = temp;

}

void mutate()
{
	int i, j;
	double lbound, hbound, min, max;
	double x;
	genotype temp = population[0];
	/*muta for the best*/

	for (i = Lflag; i != Lflag + NVARS; i++)
	{
		/*min = find_min(i % (NVARS));
		max = find_max(i % (NVARS));*/
		min = mmin[i % (NVARS)];
		max = mmax[i % (NVARS)];
		temp.gene[i % (NVARS)] = randval(min, max);
		temp.fitness = evaluate_best(temp.gene);
		file.open("muta.txt", ios::app);
		file << generation << '\t' << besttime << min << '  ' << max << '  ' << temp.gene[i % (NVARS)] << "\n";
		//cout << "@\n";
		if (temp.fitness > population[0].fitness)
		{
			population[0].gene[i % (NVARS)] = temp.gene[i % (NVARS)];
			//population[0].fitness=
			Lflag = i % (NVARS);
			cout << i % (NVARS) << "%%%%%%%%%%%";//
												 //fprintf(f, "%.6f   ,%.6f", i % (NVARS), generation);
			besttime++;

			file << "success" << generation << '\t' << besttime << min << '  ' << max << '  ' << temp.gene[i % (NVARS)] << "\n";

			break;
		}
		
		//cout<<Lflag<<"#";
		//cout << i << "$";//
	}
	file.close();
	/*mutua from the second one*/
	for (i = 1; i < POPSIZE; i++)
		for (j = 0; j < NVARS; j++)
		{
			x = rand() % 1000 / 1000.0;
			if (x < PMUTATION)
			{
				/* find the bounds on the variable to be mutated */
				lbound = population[i].lower[j];
				hbound = population[i].upper[j];
				population[i].gene[j] = randval(lbound, hbound);//popsize中的upper无定义
																/*if (population[i].gene[j] == 0)
																{
																cout << "[" << i << "," << j << "]" ;
																cout << lbound << "()" << hbound << endl;
																}*/
			}
		}
	
	//cout << i << "&";//
}

double evaluate_best(double *gene)
{
	double fit = 0;
	for (int i = 0; i < NVARS; i++)
	{
		//x[i + 1] = population[mem].gene[i];
		fit += gene[i] * gene[i];
	}  //根据
	fit = 1 / (fit + 1);
	return fit;
}

double find_min(int n)
{
	double min = population[0].gene[n];
	for (int i = 1; i < POPSIZE; i++)
	{
		if (min > population[i].gene[n])
		{
			min = population[i].gene[n];
		}
	}
	if (min < 0 && 1.2*min>-2)
		min *= 1.2;
	if (min > 0 && 0.8*min<2)
		min *= 0.8;
	return min;
}

double find_max(int n)
{
	double max = population[0].gene[n];
	for (int i = 1; i < POPSIZE; i++)
	{
		if (max < population[i].gene[n])
		{
			max = population[i].gene[n];
		}
	}
	double min = find_min(n);
	if (min > 0 && 1.2*max<2)
		max *= 1.2;
	if (min < 0 && 0.8*max>-2)
		max *= 0.8;
	return max;
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

int main(void)
{
	srand(time(NULL)); //增加这么一句，以保证每次产生随机数的时候都不同
	int i;
	generation = 0;
	initialize();
	evaluate();
	keep_the_best();
	cout << population[POPSIZE].fit << endl;
	while (generation<10000)
	{
		generation++;
		select();             //选择操作
		crossover();       //交配操作
		mutate();            //变异操作
							 //report();
		evaluate();
		elitist();
		cout << generation << " ";
		/*for (int i = 0; i < 10; i++)
		cout << population[0].gene[i] << "^";*/
		printf("%.26f\n", population[0].fit);
	}

	system("pause");
}



