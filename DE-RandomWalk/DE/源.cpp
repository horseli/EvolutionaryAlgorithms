#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>
#include<time.h>
#include<fstream>
using namespace std;

#define POPSIZE 20               /* population size */
#define FES 300000             /* max. number of generations */
#define NVARS 30                 /* no. of problem variables */
#define CR 0.1               /* probability of crossover */
#define F 0.1				/* amplified */
#define TRUE 1
#define FALSE 0
#define LOW -10
#define HIGH 10
double RW = 0.1;
int generation;
struct De
{
	double demen[NVARS];
	double fitness;

};
De population[POPSIZE];
De newpop[POPSIZE];
int cur_best;		/*best index*/
int r1, r2;
int fes;
void initialize();
double randval(double, double);
void evaluate_one(De&);
void evaluate(bool);         /*use bool to judge whether evaluate population or newpop*/
int randint(int low, int high);
void muta(De&);
void crossover(int);
void cycle();               /*the cycle for muat and crossover in a run*/
void select();

void initialize() {
	/*intial population*/
	for (int i = 0; i < POPSIZE; i++) {
		for (int j = 0; j < NVARS; j++) {
			population[i].demen[j] = randval(LOW, HIGH);
		}
		population[i].fitness = 0;
	}
	evaluate(TRUE);
	/*find cur_best*/	
	for (int i = 0; i < POPSIZE; i++) {
		if (population[cur_best].fitness > population[i].fitness) {
			cur_best = i;
		}
	}
	/*initial newpop fitness*/
	for (int i = 0; i < POPSIZE; i++) {
		newpop[i].fitness = 0;
	}
}
double randval(double low, double high) {
	return low + (high - low)*rand()*1.0 / RAND_MAX;
}
int randint(int low, int high) {
	return low + abs(rand()) % NVARS;     /*0-29+0*/
}
void evaluate_one(De&p) {
	p.fitness = 0;

	/*for (int i = 0; i < NVARS; i++) {
		p.fitness += pow(p.demen[i],2);
	}*/
	for (int i = 0; i < NVARS-1; i++)
	{
		p.fitness += 100*pow((p.demen[i+1]-pow(p.demen[i],2)),2)+pow(p.demen[i]-1, 2);
	}
	fes++;
}
void evaluate(bool ol) {
	if (ol == TRUE) {
		for (int i = 0; i < POPSIZE; i++) {
			evaluate_one(population[i]);
		}
	}
	if (ol == FALSE) {
		for (int i = 0; i < POPSIZE; i++) {
			evaluate_one(newpop[i]);
		}
	}
}
void muta(De&p) {
	r1 = randint(0, POPSIZE);
	r2 = randint(0, POPSIZE);
	for (int i = 0; i < NVARS; i++) {
		p.demen[i] = population[cur_best].demen[i] + F*(population[r1].demen[i] - population[r2].demen[i]);
		/*if overflow*/
		if (p.demen[i] > HIGH)
			p.demen[i] = HIGH;
		if (p.demen[i] < LOW)
			p.demen[i] = LOW;
	}
}
int rn = randint(0, NVARS);
void crossover(int n) {
	double rd;
	RW = 0.1 - 0.099*fes /FES;
	for (int i = 0; i < NVARS; i++) {
		rd = randval(0, 1);
		//cout << rn << endl;
		if (rd < CR || i == rn)
			continue;
		else if (randval(0, 1) < RW)
			newpop[n].demen[i] = randval(LOW, HIGH);
		else newpop[n].demen[i] = population[n].demen[i];
		/*else if (rd > CR && i != rn) {
			newpop[n].demen[i] = population[n].demen[i];
}*/
	}
}
void cycle() {
	for (int i = 0; i < POPSIZE; i++) {
		muta(newpop[i]);
		crossover(i);
	}
}
void select() {
	/*evaluate newpop*/
	evaluate(FALSE);
	/*select*/
	for (int i = 0; i < POPSIZE; i++) {
		if (newpop[i].fitness < population[i].fitness) {
			population[i].fitness = newpop[i].fitness;
			for (int j = 0; j < NVARS; j++) {
				population[i].demen[j] = newpop[i].demen[j];
			}
		}
	}
	//evaluate(TRUE);
	/*update cur_best*/
	for (int i = 0; i < POPSIZE; i++) {
		if (population[cur_best].fitness > population[i].fitness) {
			cur_best = i;
		}
	}
	/*times count*/
	//fes++;
}
int main() {
	srand(time(NULL));
	initialize();
	cout << population[cur_best].fitness << endl;
	while (fes < 300000) {
		cycle();
		select();
		/*for (int i = 0; i < 10; i++)
			cout << population[2].demen[i] << "^^";*/
		cout << population[2].fitness << endl;
		generation++;
	}
	/*initialize();
	for (int i = 0; i < POPSIZE; i++) {
		cout << population[i].fitness << endl;
	}*/
	system("pause");
}
