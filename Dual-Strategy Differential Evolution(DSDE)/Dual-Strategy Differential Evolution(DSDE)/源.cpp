#include<iostream>
#include<cmath>
#include<time.h>
using namespace std;
#define epslo0 pow(10,-2)//accuracy level
#define epslo1 pow(10,-3)
#define epslo2 pow(10,-4)
#define Max_FEs 5*pow(10,4)
#define N 80//population size
#define T 40//archive technique
#define F 0.5//amplification factor
#define CR 0.9//crossover rate
#define RUNtotal 51//run times
#define NVARS 2//dimention
#define LOW -1.1
#define HIGH 1.1
#define LOWM 4//interval of M in specie clustering
#define HIGHM 20
#define max_its 200
#define con_its 30
#define lamida 0.9
#define shPreference -5
#define Theta pow(10,-4)
struct particle {
	double x[NVARS];
	double fintness;
};
struct popPart{
	particle*pop;
	int len;
};
struct popClust {
	popPart*pt;
	int len;
};

void evaluateSinglef1(particle&);
void evalute();
double randval(double, double);
int randint(int low, int high);
void initialization();
void sort(particle*,int);
void sort(double*p, int n);
void findNearest(particle*, int);//找到最近的
void specieClustering();
void dualMutation();
void mutation(const popPart,int);//排序的 和在pop2中的起始位置
void crossover(int);
void APCClustering();
bool findMax(double, double);
void proSelection();
void archiveTechnique();
particle pop0[N];
popClust pop1;
particle pop2[2 * N]; //for selection
particle pop3[N];//for crossover
int fes = 0;
int M;//划分的种群大小
int popBest;//pop0中最优的下标
int record[N];//记录pop0[i]是否在
double r0[2 * N][2 * N];//responsibility
double r1[2 * N][2 * N];
double a0[2 * N][2 * N];//avalibility
double a1[2 * N][2 * N];
double s[2 * N][2 * N];//distance
int sc[N];//记录是否改变
particle archive[10000];
int archivelen = 0;
int pop2len = 0;
void initialization() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < NVARS; j++) {
			pop0[i].x[j] = randval(LOW, HIGH);
		}
		evaluateSinglef1(pop0[i]);
		sc[i] = 0;//intialize sci
	}
}
double randval(double low, double high) {
	return low + (high - low)*rand()*1.0 / RAND_MAX;
}
int randint(int low, int high) {
	return rand()%(high-low+1)+low;
}
void evaluateSinglef1(particle&p) {//f1为一元函数
	double x = p.x[0];
	/*for (int i = 0; i < N; i++) {
		double x = p.x[0];
		if (x >= 0 && x < 2.5)
			p.fintness = 80 * (2.5 - x);
		else if (x < 5)
			p.fintness = 64 * (x - 2.5);
		else if (x < 7.5)
			p.fintness = 64 * (7.5 - x);
		else if (x < 12.5)
			p.fintness = 28 * (x - 7.5);
		else if (x < 17.5)
			p.fintness = 28 * (17.5 - x);
		else if (x < 22.5)
			p.fintness = 32 * (x - 17.5);
		else if (x < 27.5)
			p.fintness = 32 * (27.5 - x);
		else p.fintness = 80 * (x - 27.5);
	}*/

	//p.fintness = pow(sin(5 * 3.1415926*p.x[0]), 6);
	//p.fintness = exp(-2 * log(2)*((x - 0.08) / 0.854)*((x - 0.08) / 0.854))*pow(sin(5 * 3.1415926*(pow(x, 0.75) - 0.05)), 6);
	double y = p.x[1];
	//p.fintness = 200 - pow(x*x + y - 11, 2) - pow(x + y*y - 7, 2);
	p.fintness = -4*((4 - 2.1*x*x + pow(x, 4)*1.0 / 3)*x*x + x*y + (4 * y*y - 4)*y*y);
	fes++;
}
void sort(particle*p,int n) {//从大到小
	for (int i = 1; i < n; i++) {
		int r = 0;
		for (int j = 0; j < n-i; j++) {
			if (p[j].fintness < p[j + 1].fintness) {
				particle temp = p[j];
				p[j] = p[j + 1];
				p[j + 1] = temp;
				r = 1;
			}
		}
		if (r == 0)
			break;
	}
}
void sort(double*p, int n) {//从小到大排列
	int *index = new int[n];//下标
	for (int i = 0; i < n; i++) {
		index[i] = i;
	}
	for (int i = 1; i < n; i++) {
		int r = 0;
		for (int j = 0; j < n - i; j++) {
			if (p[j] > p[j + 1]) {
				double temp = p[j];
				p[j] = p[j + 1];
				p[j + 1] = temp;
				r = 1;
				int t = index[j];
				index[j] = index[j + 1];
				index[j + 1] = t;
			}
		}
		if (r == 0)
			break;
	}

	for (int i = 0; i < n; i++)
		p[i] = index[i];//返回下标
}
void findNearest(particle*p, int len) {
		double distance[N];
		for (int j = 0; j < N; j++) {
			distance[j] = 0;
			for(int k=0;k<NVARS;k++)
			distance[j] += pow((pop0[j].x[k] - pop0[popBest].x[k]),2);//到本身的距离是最小的 包括在内
		}
		sort(distance, N);//返回的distance为排序后的下标
		int dex = -1;
		for (int i = 0; i < N; i++) {
			if (record[int(distance[i])])continue;
			dex++;
			record[(int)distance[i]] = 1;//已被移除
			p[dex] = pop0[(int)distance[i]];//移动
			if (dex == len - 1)
				break;
		}
}
void specieClustering() {
	M = randint( LOWM ,HIGHM );
	int groN = N / M;
	pop1.len =groN ;
	pop1.pt = new popPart[groN];
	//sort(pop0, N);
	for (int i = 0; i < N; i++)//记录重置为0
		record[i] = 0;
	for (int i = 0; i < groN; i++) {
		if (i == groN - 1){
			pop1.pt[i].pop = new particle[M + N%M];//最后一个种群数量多
			pop1.pt[i].len = M + N%M;
		}
		else {
			pop1.pt[i].pop = new particle[M]; 
			pop1.pt[i].len = M;
		}
		popBest = -1;
		for (int j = 0; j < N; j++) {//find popbest
			if (record[j])
				continue;
			if (popBest == -1)
				popBest = j;
			if (pop0[popBest].fintness < pop0[j].fintness)
				popBest = j;
		}
		findNearest(pop1.pt[i].pop, pop1.pt[i].len);
		/*tiaosi*/
		//for(int i=0;i<pop1.len;i++)
		cout << "&&&&&&&&" << endl;
		for (int x = 0; x < pop1.pt[i].len; x++)
			cout << pop1.pt[i].pop[x].x[0]<<" ";
		cout << endl;
	}
}
void mutation(const popPart p, int n) {
	double*dis = new double[p.len];
	for (int i = 0; i < p.len; i++) {
		dis[i] = p.pop[i].fintness;
	}
	sort(dis, p.len);//dis变成了下标 从小到大排列
	/*tiaosi*/
	cout << endl;
	for (int i = 0; i < p.len; i++)
		cout << p.pop[(int)dis[i]].fintness << "^";
	for (int i = 0; i < p.len; i++) {
		if (i >= p.len / 2)//lbest
		{
			int bes = int(dis[p.len - 1]);
			int r1 = randint(0, p.len - 1);
			int r2 = randint(0, p.len - 1);
			for (int j = 0; j < NVARS; j++) {
				int ind = int(dis[i]);
				pop2[i+n].x[j] = p.pop[bes].x[j]+ F*(p.pop[r1].x[j] - p.pop[r2].x[j]);
				if (pop2[i + n].x[j] > HIGH)//越界处理
					pop2[i + n].x[j] = HIGH;
				if (pop2[i + n].x[j] < LOW)
					pop2[i + n].x[j] = LOW;
			}
			crossover(i + n);
		}
		else {
			int r1 = randint(0, p.len - 1);
			int r2 = randint(0, p.len - 1);
			int r3 = randint(0, p.len - 1);
			for (int j = 0; j < NVARS; j++) {
				int ind = int(dis[i]);
				pop2[i + n].x[j] = p.pop[ind].x[j] + rand()/double(RAND_MAX)*(p.pop[r1].x[j] - p.pop[ind].x[j]) + F*(p.pop[r2].x[j] - p.pop[r3].x[j]);//rand()不是0到一容易越界？
				if (pop2[i + n].x[j] > HIGH)
					pop2[i + n].x[j] = HIGH;
				if (pop2[i + n].x[j] < LOW)
					pop2[i + n].x[j] = LOW;
			}
		}

		evaluateSinglef1(pop2[i+n]);//适应值计算
	}
}
void crossover(int n) {
	double rd;
	int rn = randint(0, NVARS);
	for (int i = 0; i < NVARS; i++) {
		rd = randval(0, 1);
		if (rd < CR || i == rn) {}//不变
		else pop2[n].x[i] = pop0[n].x[i];
	}
}
void dualMutation() {
	int n = 0;
	for (int i = 0; i < pop1.len; i++) {
		mutation(pop1.pt[i], n);
		n += pop1.pt[i].len;
	}
	for (int i = N; i < 2*N; i++) {///////////80
		pop2[i] = pop0[i-N];
	}
}
bool findMax(double a, double b) {
	if (a > b)
		return true;
	else return false;
}
void APCClustering() {
	double dislen[2 * N*(2 * N - 1)];
	double median;
	for (int i = 0; i < 2 * N; i++)//initialization
		for (int j = i; j < 2 * N; j++) {
			a0[i][j] = 0;
			a0[j][i] = 0;
			if (i != j) {
				s[i][j] = 0;   //r[k][k]?
				for (int k = 0; k < NVARS; k++)
					s[i][j] -= pow(pop2[i].x[k] - pop2[j].x[k], 2);
				s[j][i] = s[i][j];
			}
			else s[i][j] = shPreference;
			//r0[i][j] = s[i][j];
		}
	//s[i][i]赋予中位数
	int le = 0;
	for (int i = 0; i < 2 * N; i++) {
		for (int j = 0; j < 2 * N; j++) {
			if (i == j)
				continue;
			dislen[le] = s[i][j];
			le++;
		}
	}
	for (int i = 1; i < le; i++) {
		int re = 0;
		for (int j = 0; j < le - i; j++) {
			if (dislen[j] < dislen[j + 1])
			{
				double temp = dislen[j];
				dislen[j] = dislen[j + 1];
				dislen[j + 1] = temp;
				re = 1;
			}
		}
		if (re == 0)
			break;
	}
	median = dislen[2 * N*(2 * N - 1) / 2];
	for (int i = 0; i < 2 * N; i++)
		s[i][i] = median;

	int mits = 0;
	int recordChange = 0;
	int  cits = 0;
	int *exemplar = new int[2*N];
	int len = 0;
	int *exeTemp = new int[2*N];
	int lenTemp = 0;
	int lenlast = 0;
	while (mits < max_its&&cits < con_its) {
		recordChange = 0;
		exemplar = new int[2*N];
		len = 0;
		/*find max to r*/
		for (int i = 0; i < 2 * N; i++) {
			for (int k = 0; k < 2 * N; k++) {
				double max;
				int max_index;
				if (i == k)
					max = s[i][0];
				else max = a0[i][0] + s[i][0];//初始化最大值
				max_index = 0;
				for (int j = 0; j < 2 * N; j++) {
					if (i == k) {
						if (j == k) {
							if (j == 0) {
								max = s[i][1];
								max_index = 1;
							}
							continue;
						}
						if (max < s[i][j]) {//find max
							max = s[i][j];
							max_index = j;
						}
					}
					else {
						if (j == k) {
							if (j == 0) {
								max = a0[i][1] + s[i][1];
								max_index = 1;
							}
							continue;
						}
						if (max < a0[i][j] + s[i][j]) {//find max
							max = a0[i][j] + s[i][j];
							max_index = j;
						}
					}
				}
				r0[i][k] = s[i][k] - max;//r更新
			}
		}
		//cout << r0[0][0] << endl;
		/*update r*/
		if (mits >= 1) {//除第一次外都要更新
			for (int i = 0; i < 2 * N; i++)
				for (int j = 0; j < 2 * N; j++) {
					r0[i][j] = (1 - lamida)*r0[i][j] + lamida*r1[i][j];
				}
		}
		/*find max to a*/
		for (int i = 0; i < 2 * N; i++) {
			for (int k = 0; k < 2 * N; k++) {
				double min = 0;
				for (int j = 0; j < 2 * N; j++) {
					if (j == k || j == i)
						continue;
					min += 0 > r0[j][k] ? 0 : r0[j][k];
				}
				if (i != k)
					a0[i][k] = 0 < r0[k][k] + min ? 0 : (r0[i][k] + min);//a更新
																		 /*to a[k,k]*/
				else {
					min += 0 > r0[i][k] ? 0 : r0[i][k];
					a0[k][k] = min;
				}//a更新
			}
		}

		/*update a with lamida*/
		if (mits >= 1) {//除第一次外都要更新
			for (int i = 0; i < 2 * N; i++)
				for (int j = 0; j < 2 * N; j++) {
					//r0[i][j] = (1 - lamida)*r0[i][j] + lamida*r1[i][j];
					a0[i][j] = (1 - lamida)*a0[i][j] + lamida*a1[i][j];
				}
		}
		for (int i = 0; i < 2 * N; i++)//保存上一次的值
			for (int j = 0; j < 2 * N; j++) {
				r1[i][j] = r0[i][j];
				a1[i][j] = a0[i][j];
			}

		/*find the exemplar*/
		for (int i = 0; i < 2 * N; i++) {
			int maxk = 0;
			for (int k = 0; k < 2 * N; k++) {
				if (r0[i][maxk] + a0[i][maxk] < r0[i][k] + a0[i][k])
					maxk = k;
			}
			int kk = 0;
			for (int s = 0; s < len; s++) {
				if (exemplar[s] == maxk)
					kk = 1;
			}
			if (kk == 0) {
				exemplar[len] = maxk;
				len++;
				//cout << len << endl;
				//初始化exeTemp
				if (mits == 0) {
					exeTemp[lenTemp] = maxk;
					lenTemp++;
				}
			}

		}
		//记录是否改变
		if (len != lenTemp&&mits)
			recordChange = 1;
		else if (mits)for (int j = 0; j < len; j++) {//除第一次外都要比较
			if (exemplar[j] != exeTemp[j]) {
				recordChange = 1;
				break;
			}
		}
		//exeTemp 记录上一次的exemplar
		if (mits) {
			lenTemp = 0;
			for (int i = 0; i < len; i++) {//除第一次和相同都要赋值传递给temp
				exeTemp[i] = exemplar[i];
				lenTemp++;
			}
		}
		if (recordChange == 0)
			cits++;
		else cits = 0;
		mits++;
		/*taiosi*/
		//cout << len << endl;
		/*for (int i = 0; i < 10; i++)
		cout << r0[0][i] << " ";
		cout << endl<<endl;*/
		//cout << len <<mits<< endl;
		//cout << len << endl;
		/*	if (mits>50)
		if (len > lenlast)
		break;
		lenlast = len;*/
	}
	//根据距离进行分类各个point属于哪个exemplar
	int recordd[2 * N];
	for (int i = 0; i < 2 * N; i++)
		recordd[i] = 0;//记录是否被取

	pop1.len = len;
	pop1.pt = new popPart[len];//subpopulation new
	for (int i = 0; i < len; i++) {
		pop1.pt[i].pop = new particle[2 * N];
		pop1.pt[i].pop[0] = pop2[exemplar[i]];//赋初值
		pop1.pt[i].len = 1;
		recordd[exemplar[i]] = 1;
	}
	//int x = 0;
	//int y = 0;
	for (int i = 0; i < 2 * N; i++) {
		int choice = 0;
		for (int j = 0; j < len; j++) {
			if (r0[i][exemplar[j]] + a0[i][exemplar[j]] > r0[i][exemplar[choice]] + a0[i][choice])//find max to cluster
				choice = j;
		}
		/*	for (int k = 0; k < pop1.pt[choice].len; k++) {
		if(pop2[i].x)
		}*/
		if (recordd[i])//自身已经添加
		{
			//y++;
			continue;
		}
		pop1.pt[choice].pop[pop1.pt[choice].len] = pop2[i];//add
		//x++;
		pop1.pt[choice].len++;
		recordd[i] = 1;
	}
}
void proSelection() {
	double fmin, fmax;
	double*p = new double[pop1.len];
	int*lent = new int[pop1.len];
	for (int i = 0; i < pop1.len; i++) {
		sort(pop1.pt[i].pop, pop1.pt[i].len);
		if (i == 0) {
			fmin = pop1.pt[i].pop[0].fintness;
			fmax = pop1.pt[i].pop[0].fintness;
		}
		else {
			if (fmin > pop1.pt[i].pop[0].fintness)
				fmin = pop1.pt[i].pop[0].fintness;
			if (fmax < pop1.pt[i].pop[0].fintness)
				fmax = pop1.pt[i].pop[0].fintness;
		}
	}
	/*tiaosi*/
	//cout << "*"<<pop1.pt[2].len<<endl;
	//for (int i = 0; i < pop1.pt[2].len; i++) {
	//	cout << pop1.pt[2].pop[i].fintness << " ";
	//}
	for (int i = 0; i < pop1.len; i++)//定义选择概率
	{
		p[i] = (pop1.pt[i].pop[0].fintness - fmin + Theta) / (fmax - fmin + Theta);
		lent[i] = 0;//定义0
	}
	int N0 = 0, i = 0;
	while (N0 < N) {
		if (pop1.pt[i].len>lent[i]) {
			if (double(rand()) / RAND_MAX < p[i]) {
				//particle temp = pop0[N0];
				pop0[N0] = pop1.pt[i].pop[lent[i]];//insert
				lent[i]++;
				N0++;
			}
		}
		i++;
		if (i >= pop1.len)
			i = 0;
	}
}
void archiveTechnique() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (pop0[i].fintness == pop2[j+N].fintness) {
				sc[i]++;
				double temp = sc[i];
				sc[i] = sc[j];
				sc[j] = temp;
				break;
			}
		}
	}
	for (int i = 0; i < N; i++) {
		if (sc[i] >= T) {
		/*	popBest = i;
			for (int i = 0; i < N; i++)
				record[i] = 0;
			findNearest()*/
			double distance[N];
			for (int j = 0; j < N; j++) {
				for (int k = 0; k<NVARS; k++)
					distance[j] += pow((pop0[j].x[k] - pop0[i].x[k]), 2);
			}
			sort(distance, N);//返回的distance为排序后的下标
			int lenM = 0;
			for (int j = 0; j < N; j++) {//store and reintialize
				if (pop0[(int)distance[j]].fintness > pop0[j].fintness)
					continue;
				if (lenM == M+1)
					break;
				archive[archivelen] = pop0[j];
				archivelen++;
				lenM++;
				for (int k = 0; j < NVARS; j++) {//reintialize
					pop0[j].x[k] = randval(LOW, HIGH);
				}
				evaluateSinglef1(pop0[j]);
				sc[j] = 0;//intialize sci
			}
		}
	}
}
void findOptim() {
	for (int i = 0; i < archivelen; i++) {
		if (archive[i].x[0] < epslo1)
			cout << 0 << " ";
		if(archive[i].x[0]>30-epslo2) {
			cout <<5<< " ";
		}
	}
	//for (int i = 0; i < N; i++) {
	//	if (pop0[i].x[0] < 0.1+epslo0&&pop0[i].x[0]>0.1-epslo2)
	//		cout << -1 << " ";
	//	if (pop0[i].x[0]>(30 - epslo0)&&pop0[i].x[0]<30) {
	//		cout << -5 << " ";
	//	}
	//}
	cout << "\n$$$$$$" << endl;
	/*for (int i = 0; i < N; i++) {
		if (pop0[i].x[0] < 0.1 + epslo2&&pop0[i].x[0]>0.1 - epslo2)
			cout << -1 << " ";
		if (pop0[i].x[0]>(0.3 - epslo2) && pop0[i].x[0]<0.3+epslo2) 
			cout << -2 << " ";
		if (pop0[i].x[0]>(0.5 - epslo2) && pop0[i].x[0]<0.5 + epslo2) 
			cout << -3 << " ";
		if (pop0[i].x[0]>(0.7 - epslo2) && pop0[i].x[0]<0.7 + epslo2) 
			cout << -4 << " ";
		if (pop0[i].x[0]>(0.9 - epslo2) && pop0[i].x[0]<0.9 + epslo2) 
			cout << -5 << " ";
		}*/
	/*for (int i = 0; i < N; i++) {
		if ();
	}*/
	
}
void FullCycle() {
	initialization();
	while (fes < Max_FEs) {
		specieClustering();
		dualMutation();
		APCClustering();
		proSelection();
		//archiveTechnique();
		findOptim();
		cout << fes<<" "<<pop1.len;
		cout << endl;
	}
	//system("pasue");
}
int main() {
	srand(time(NULL));
	FullCycle();
	system("pause");
}