#include "problem.h"
#include "utils.h"
#include "random.h"
#include "timer.h"
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cstdio>
using namespace std;


#define USAGE_STRING "Usage: ./ls SEED (vns4|ins) INSTANCE P[0] P[1] ... P[n-1]"


//variables used in depLocalSearch.cpp
static unsigned int seed;			//seed
static int* gbest;					//global best so far
static FitnessType fgbest;			//fitness of the global best so far
static int nfes;					//number of evaluations performed
static int nfesFoundAt;				//nfes when the global best has been found
static unsigned long timeFoundAt;	//milliseconds when the global best has been found
static int* tmpint;					//temporary memory
static int* lsHistory;				//local search history (only one perm)
static bool lsmode;					//flag indicating if local search is running or not
static int permByteSize;			//permutation size in bytes (used for memcpy)
static int nls;						//number of local searches performed
static int nImprovingls;			//number of improving local searches
static FitnessType totImprovingls;	//total delta of improving in local searches
static int nfesls;					//evaluations spent in local searches
static int ngen;					//number of generations performed
static int nrestarts;				//number of restart performed
static bool gbestls;				//true if the gbest has been obtained in a local search
static int* x;						//my variable for x genotype
static FitnessType fx;				//my variable for x fitness
static unsigned long totalTime;		//total time at the end
static char lstype[100];			//local search type
static int* inx;					//initial genotype
static FitnessType infx;			//initial fitness
static int improvingSteps;			//number of improving steps

//some functions
void initializeVariables();
void finalizeVariables();
inline bool termination();
inline void updateGbest(int* x, FitnessType fx);
inline void ins(int* x, int i, int j);
void outputBefore();
void outputAfter();

//include code
#include "onlineprint.h"
#include "depLocalSearch.cpp"

//main function
int main(int argc, char** argv) {
	//check command line parameters
	if (argc<3) {
		cerr << USAGE_STRING << endl;
		return EXIT_FAILURE;
	}
	//read seed
	sscanf(argv[1],"%u",&seed);
	//read ls type
	strcpy(lstype,argv[2]);
	//read instance
	readInstance(argv[3]);
	//read seed solution
	x = new int[n];
	for (int i=0; i<n; i++)
		x[i] = atoi(argv[i+4]);
	//set local search function pointer
	void (*localSearch)(int*,FitnessType&);
	if (strcmp(lstype,"vns4")==0)
		localSearch = localSearch_vns4;
	else if (strcmp(lstype,"ins")==0)
		localSearch = localSearch_ins;
	else {
		cerr << "LS TYPE IS UNKNOWN!!!" << endl;
		return EXIT_FAILURE;
	}
	//compute fitness of the seed solution
	fx = eval(x);
	//initialize static variables for local search
	initializeVariables();
	//copy x,fx in inx,infx
	memcpy(inx,x,sizeof(int)*n);
	infx = fx;
	//initialize rng
	if (!seed)				//if zero
		seed = randSeed();	//from /dev/random
	initRand(seed);
	//print output before running
	outputBefore();
	//init timer
	setTimer();
	//run the local search
	localSearch(x,fx);
	//get total time
	totalTime = getTimer();
	//print output after running
	outputAfter();
	//finalize memory
	finalizeVariables();
	//return success
	return EXIT_SUCCESS;
	//done
}


//initialize variables
void initializeVariables() {
	gbest = new int[n];
	fgbest = fx;
	nfes = 0;
	nfesFoundAt = 0;
	timeFoundAt = 0;
	tmpint = new int[2*n];	//the implementation of local search requires this amount of memory
	lsHistory = new int[n];
	for (int i=0; i<n; i++)
		lsHistory[i] = -1;	//something not possible
	lsmode = true;
	permByteSize = sizeof(int)*n;
	nls = 0;
	nImprovingls = 0;
	totImprovingls = 0;
	nfesls = 0;
	ngen = 0;
	nrestarts = 0;
	gbestls = false;
	inx = new int[n];
	improvingSteps = 0;
}


//finalize variables
void finalizeVariables() {
	delete[] x;
	delete[] gbest;
	delete[] tmpint;
	delete[] lsHistory;
	delete[] inx;
}


//always false (needed by the code in depLocalSearch.cpp)
bool termination() {
	return false;
}


//update global best
inline void updateGbest(int* x, FitnessType fx) {
#ifdef MINIMIZATION
	if (fx<fgbest) {
#else
	if (fx>fgbest) {
#endif
		fgbest = fx;
		memcpy(gbest,x,permByteSize);
		nfesFoundAt = nfes;
		timeFoundAt = getTimer();
		improvingSteps++;
#ifdef ONLINEPRINT
		bestPrint();
#endif
#ifdef MYDEBUG
		if (!permValid(gbest,n)) {
			cout<<"updateGbest gbest"<<endl;
			exit(1);
		}
#endif
	}
}


//insertion move
void ins(int* x, int i, int j) {
	//perform insertion (i,j) on x assuming nothing on i and j
	//(i,j) means "move job at pos i to pos j"
	int t = x[i];
	if (i<j) //forward
		memmove(x+i,x+i+1,sizeof(int)*(j-i));
	else //backward
		memmove(x+j+1,x+j,sizeof(int)*(i-j));
	x[j] = t;
}

//write the output before the run
void outputBefore() {
	//only cerr
	cerr << "LS type =\t\t" << lstype << endl;
	cerr << "Starting fitness =\t" << fx << endl;
}

//write the ouput after the run
void outputAfter() {
	//cerr
	char stime[5000];
	cerr << "Final fitness =\t\t" << fx << endl;
	millis2str(totalTime,stime);
	cerr << "Total time =\t\t" << stime << endl;
	cerr << "Total nfes =\t\t" << nfes << endl;
	millis2str(timeFoundAt,stime);
	cerr << "Time found at =\t\t" << stime << endl;
	cerr << "Nfes found at =\t\t" << nfesFoundAt << endl;
	cerr << "Improving steps =\t" << improvingSteps << endl;
	cerr << "Fitness improvement =\t" << totImprovingls << endl;
	//cout
	char sperm1[5000],sperm2[5000];
	cout << "lstype,seedrng,inst,n,startsol,startfit,finalsol,finalfit,totalTime,totalNfes,timeFoundAt,nfesFoundAt,improvingSteps,fitImpr" << endl;
	perm2str(inx,n,sperm1);
	perm2str(x,n,sperm2);
#ifdef FIT_REAL
	fprintf(stdout,"%s,%u,%s,%d,%s,%.10lf,%s,%.10lf,%u,%d,%u,%d,%d,%.10lf\n",lstype,seed,instance,n,sperm1,infx,sperm2,fx,totalTime,nfes,timeFoundAt,nfesFoundAt,improvingSteps,totImprovingls);
#else
	fprintf(stdout,"%s,%u,%s,%d,%s,%d,%s,%d,%u,%d,%u,%d,%d,%d\n",lstype,seed,instance,n,sperm1,infx,sperm2,fx,totalTime,nfes,timeFoundAt,nfesFoundAt,improvingSteps,totImprovingls);
#endif
}

