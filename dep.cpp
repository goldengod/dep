#include "global.h"
#include "dep.h"
#include "problem.h"
#include "random.h"
#include "timer.h"
#include <cstring>   //memcpy, memset, memmove
#include <climits>   //INT_MAX
#include <unistd.h>  //alarm
#include <cstdio>    //FILE, fprintf
#include <iostream>  //cerr, endl
#include <cstdlib>   //exit, EXIT_FAILURE
using namespace std;

#ifdef MYDEBUG
#include "utils.h"
#endif


//population size definition (declared in dep.h)
int np;

//initial scale factor definition (declared in dep.h)
double finit;

//alpha value definition (declared in dep.h)
double alpha;

//heu and ls value definition (declared in dep.h)
int heu;
int ls;

//frfactor definition (declared in dep.h)
double frfactor;

//min/max value definitions for scale factor f (used by jde rule and declared in dep.h)
double fmin;
double fmax;

//global best definition (declared in dep.h)
int* gbest;
int fgbest;
int nfesFoundAt;
int stageFoundAt;

//number of evaluations and generation definition (declared in dep.h)
int nfes;
int ngen;

//execution time definition (declared in dep.h)
unsigned long execTime; //for time save/resume scheme see save_resume.cpp

//number of restarts performed (declared in dep.h)
int nrestarts;
int nforcedrestarts; //forced

//restarts statistics definition (declared in dep.h)
int minStageLength;
int maxStageLength;
double avgStageLength;
int improvingStages;

//local search statistics definition (declared in dep.h)
int nfesls;
int nls;
int nImprovingls;
int totImprovingls;
bool gbestls;

//flag for saving execution state/memory (declared in dep.h)
bool haveToSave;

//number of improving steps performed definitino (declared in dep.h)
int improvingSteps;
int lsImprovingSteps;

//inner population definition (Iliffe/display style)
static int** x;
static int* fx;
static int** ix; //inverse population

//inner temporary population definition (Iliffe/display style)
static int** y1;
static int* fy1;
static int** y2;
static int* fy2;

//inner storage for main and temporary population (Iliffe/display style)
static int* ps; //x-y1-y2-ix interleaved //population storage
static int* fs; //fx-fy1-fy2          //fitness storage

//inner scale factor main and temporary memory
static double* sfx;
static double* sfy;
static double* sfs; //sfx-sfy

//inner temporary memory for diffMutation/crossover
static int* tmpint;  //length = n^2 + 3*n (1 matrix for cycles + 1 sorting seq. + 1 perm.) //EXCHANGES!!!

//inner bool flag indicating if fitnesses are all equal
static bool sameFitness;

//inner generation when last restart happened (to count restart statistics)
static int lastRestart;

//inner fitness at the begin of the current stage (to count the number of improving stages)
static int fgbestAtStageStart;

//inner permutation size in bytes (used for memcpy)
static int permByteSize;

//inner local search history (only one perm)
static int* lsHistory;

//inner flag indicating local search is running or not
static bool lsmode;

//inner variable used to manage execTime in case of save/resume
static unsigned long savedTime = 0; //for time save/resume scheme see save_resume.cpp

//inner variables to manage forced restart
static int nfesWhenToForceRestart;
static int forcedRestartPeriod;

//inner variable for diameter
static int diameter;


#ifdef GFC
//Branchless min between x and y, bits has to be sizeof(int)*8-1 (see http://aggregate.org/MAGIC)
#define FAST_MIN(x,y,bits) ((x)+((((y)-(x))>>(bits))&((y)-(x))))
#endif

//relative deviation of y wrt x
#define REL_DEV(y,x) (((y)-(x))/(double)(x))


//inner functions definitions
void depSave();
inline void updateGbest(int* x, int fx);
inline bool termination();
void popInit();
bool popEvolve();
inline void randmergess(int* s, int& l, int* x, double f=1.0);//EXCHANGES!!! (LIMITS) F>1
inline void randss(int* s, int& l, int* x, double f=1.0);     //EXCHANGES!!! (LIMITS) F<1
void diffMutation(int i);
inline void tpii(int* ch, int* fa, int* mo, int* ifa, int c1, int c2);
void crossover(int i);
void selection(int i);
inline void restartStatistics();
bool popRestart();
void vns4(int* x, int& fx);
#ifdef GFC
bool intStep(int* x, int& fx, bool first);
#else
bool intStep(int* x, int& fx);
#endif
void insStep(int* x, int& fx);
#ifdef GFC
inline void fwdins(int* x, int i, int j);
inline void bwdins(int* x, int i, int j);
inline void bwdins2(int* x, int i, int j);
#endif
inline void ins(int* x, int i, int j);
bool popForcedRestart();


//online printing functions
#ifdef ONLINEPRINT
#include "onlineprint.h"
#endif


#include "save_resume.cpp"


//set dep default parameters
void depDefaultParameters() {
//the manifest constants TFT and MAKESPAN are mutually exclusive
#ifdef TFT
	np = 100;
	finit = 0.5;
	alpha = 0.01;
	heu = 0; //sarebbe 1 ma metto 0 perche' devo passare il file euristica al main
	ls = B_LS;
	frfactor = 0.25;
	fmin = 0.1; //as standard jde rule
	fmax = 1.0; //as standard jde rule
#endif
#ifdef MAKESPAN
	np = 20;
	finit = 0.5;
	alpha = 0.01;
	heu = 0; //sarebbe 1 ma metto 0 perche' devo passare il file euristica al main
	ls = L_LS;
	frfactor = 0.25;
	fmin = 0.1; //as standard jde rule
	fmax = 1.0; //as standard jde rule
#endif
}


//main dep function
void dep() {
	//init timer
	setTimer();
	//init population
	popInit();
	//set a flag to distinguish termination inside popRestart (due to local search) or popEvolve
	bool lsTermination = false; //needed to correctly handle restart statistics
	//evolution loop
	do {
		//some printing
#ifdef ONLINEPRINT
		genPrint();
#endif
		//check and save dep execution state/memory
		if (haveToSave)
			depSave();
		//check and restart population
		if (popRestart() || popForcedRestart()) { //considering the short circuit
			//if true, maxnfes was exhausted during local search, so break evolution loop
			lsTermination = true;
			break;
		}
		//evolve population
	} while (popEvolve());
	//compute final execution time
	execTime = getTimer() + savedTime; //for time save/resume scheme see save_resume.cpp
	//adjust ngen and restart statistics (but only if the termination doesnt happen in a local search)
	if (lsTermination) {
		ngen++;
		nrestarts++; //let restartStatistics() believes that there was a restart
		restartStatistics();
		nrestarts--; //undo the dummy restart
	}
	//some final printing
#ifdef ONLINEPRINT
	genPrint();
#endif
	//done
}


//allocate dep memory
void depAlloc() {
	//gbest memory
	gbest = new int[n];
	//main and temporary population memory (Iliffe/display style with x[i]-y1[i]-y2[i]-ix[i])
	x = new int*[np];
	y1 = new int*[np];
	y2 = new int*[np];
	ix = new int*[np];
	ps = new int[4*np*n];
	for (int i=0; i<np; i++) {
		x[i] = ps+(i*4*n);
		y1[i] = x[i]+n;
		y2[i] = y1[i]+n;
		ix[i] = y2[i]+n;
	}
	fs = new int[3*np];
	fx = fs;
	fy1 = fx+np;
	fy2 = fy1+np;
	//main and temporary scale factor memory
	sfs = new double[2*np];
	sfx = sfs;
	sfy = sfx+np;
	//temporary diffMutation/crossover memory
	tmpint = new int[n*n+3*n]; //1 matrix cycles + 1 exchanges seq. + 1 perms.  //EXCHANGES!!!
	//set the permutation byte size
	permByteSize = sizeof(int)*n;
	//lsHistory
	lsHistory = new int[n];
	memset(lsHistory,0,permByteSize); //all zeros
	//init to false the save flag
	haveToSave = false;
	//init the diameter
	diameter = n-1; //VALID FOR EXCHANGES!!!
	//done
}


//deallocate dep memory
void depFree() {
	delete[] gbest;
	delete[] x;
	delete[] y1;
	delete[] y2;
	delete[] ps;
	delete[] fs;
	delete[] ix;
	delete[] sfs;
	delete[] tmpint;
	delete[] lsHistory;
}


//update global best
inline void updateGbest(int* x, int fx) {
	if (fx<fgbest) {
		fgbest = fx;
		memcpy(gbest,x,permByteSize);
		nfesFoundAt = nfes;
		stageFoundAt = nrestarts;
		gbestls = lsmode;
		nfesWhenToForceRestart = nfes + forcedRestartPeriod;
		improvingSteps++;
		if (lsmode)
			lsImprovingSteps++;
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


//check if max nfes has been exceeded
inline bool termination() {
	return nfes>=maxnfes || (maxTime>0 && getTimer()>=maxTime); //also time in milliseconds!!!
}


//initialize population
void popInit() {
	//init nfes and fgbest
	nfes = 0;
	fgbest = INT_MAX; //+inf
	//random population initialization apart one individual from heuristic
	for (int i=0; i<np; i++) {
		if (i==0 && heu==1)
			memcpy(x[0],heup,permByteSize); //heuristic permutation
		else
			prand(n,x[i]); //random
		fx[i] = eval(x[i]);
		nfes++;
		updateGbest(x[i],fx[i]);
		sfx[i] = finit;
		for (int k=0; k<n; k++) //ix[i] = inverse of x[i]
			ix[i][x[i][k]] = k;
#ifdef MYDEBUG
		if (!permValid(x[i],n)) {
			cout<<"popInit x["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(ix[i],n)) {
			cout<<"popInit ix["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//init ngen, nrestarts, samefitness
	ngen = 1; //initialization was the first generation
	nrestarts = 0;
	nforcedrestarts = 0;
	sameFitness = false;
	//init restart statistics variables
	lastRestart = 0;
	minStageLength = maxnfes; //it's impossible to be more
	maxStageLength = 0; //it's impossible to be less
	avgStageLength = 0.;
	fgbestAtStageStart = INT_MAX; //for sure at stage end it will be better
	improvingStages = 0;
	//init local search statistics variables
	nfesls = 0;
	nls = 0;
	nImprovingls = 0;
	totImprovingls = 0;
	//init lsmode to false
	lsmode = false;
	//init variables for forced restart
	nfesWhenToForceRestart = forcedRestartPeriod = maxnfes*frfactor;
	//init improving steps
	improvingSteps = lsImprovingSteps = 0;
	//done
}


//evolve population
bool popEvolve() {
	//differential mutation, crossover, evaluation, update best, check termination
	int i;
	for (i=0; i<np; i++) {
		//generate one mutant
		diffMutation(i);
		//generate two offspring
		crossover(i);
		//evaluate first offspring
		fy1[i] = eval(y1[i]);
		nfes++;
		updateGbest(y1[i],fy1[i]);
		//check termination
		if (termination())
			return false;
		//evaluate second offspring
		fy2[i] = eval(y2[i]);
		nfes++;
		updateGbest(y2[i],fy2[i]);
		//check termination
		if (termination())
			return false;
		//population loop done
	}
	//selection and compute sameFitness
	sameFitness = true;
	for (i=0; i<np; i++) {
		selection(i);
		sameFitness &= !(fx[i]-fx[0]); //it's the same of "sameFitness = sameFitness && fx[i]==fx[0]"
	}
	//increment generation counter
	ngen++;
	//termination has not been triggered so return true
	return true;
	//done
}

//USED FOR F>1 (NO CHECK INSIDE)
//s/l = sequence that bring x towards a 1 cycle permutation (SWAP_k(i,j) is s[2*k]<->s[2*k+1])
//f is used to bound the insertions chain
void randmergess(int* s, int& l, int* x, double f) {
	//some variables
	int nc,nv,nexc,i,j,t,lim,r;
	double tempDouble;
	static bool* v = (bool*)(tmpint+n); //alive with x,c,lc (reuse part of s not used by x)
	static int* lc = tmpint+(2*n);      //alive with x,s,c,v (memory just after v)
	static int* c = tmpint+(3*n);       //alive with x,v,c,lc (last chunk of memory of size n^2)
	/* PHASE1 / Cycles decomposition of x (**ALL** cycles) / outputs will be:
		- nc is the TOTAL number of cycles of x
		- c is a (n*n)-matrix s.t. c[i] contains the items of i-th cycle of x
		- lc is a n-vector s.t. lc[i] is the length of c[i]
		- nexc is the (double of the) total number of exchanges for UNsorting x (it is SUM_{i=0}^{nc-1}(lc[i]*(n-lc[i])))
		- COST: O(n)
	*/
	//initialize v to all false and other counters to 0 (also the first cycle length)
	memset(v,0,sizeof(bool)*n); //all false
	nc = nexc = nv = i = lc[0] = 0;
	//scan all the items in the order given by the permutation (first open cycle is nc=0)
	for (;;) {
		//insert current index in the currently open cycle (and increase its length)
		c[nc*n+(lc[nc]++)] = i; //in matrix form = c[nc][lc[nc]++] = i;
		//set the visited flag of the i-th index and increase nv
		v[i] = true;
		nv++;
		//if all indexes have been visited, then close the cycle and exit the loop
		if (nv==n) {
			nexc += lc[nc]*(n-lc[nc]);
			nc++;
			break;
		}
		//try to move index by following the cycle
		i = x[i];
		//check if cycle has to be closed
		if (v[i]) { //same as if (i==c[nc*n]) { //matrix form i==c[nc][0]
			//close the cycle
			nexc += lc[nc]*(n-lc[nc]);
			nc++;
			//open next cycle by setting/resetting its length
			lc[nc] = 0;
			//move to next unvisited index (no need to move circularly, leftmost unv. index ensured)
			while (v[++i]);
		} //end-if
	} //end-infinite-for
#ifdef MYDEBUG
	int debugSum2 = 0;
	for (int debugi=0; debugi<nc; debugi++)
		debugSum2 += lc[debugi];
	if (debugSum2!=n) {
		cout << "sum of cycles lengths does not match n" << endl;
		exit(1);
	}
#endif
	/* PHASE2 / Sorting sequence of x towards a 1-cycle permutation / outputs will be:
		- l (==f*(n-nc)-(n-nc)), i.e. the length of a decomposition (sorting sequence) of x
		- s, i.e. a sequence of exchanges that moves x towards one of a 1 cycle permutation
		- does not touch x
		- COST: O(n^2)
	*/
	//compute lim: "n-nc" are the exchanges from E to X; thus I need "max{f*(n-nc),diameter}" exchanges; since "n-nc" are already in X the lim is "max{f*(n-nc),diameter} - (n-nc)"
	lim = tempDouble = f*(n-nc);	//main equation part
	if (lim<tempDouble)				//ceil rounding ok
		lim++;						//ceil rounding ok
	if (lim>diameter)				//max part
		lim = diameter;				//max part
	lim -= n-nc;					//discount part
#ifdef MYDEBUG //check if lim+(n-nc)<=diameter
if (lim+(n-nc)>diameter) {
	cout << "problem with lim in randmergess" << endl;
	exit(1);
}
#endif
	//initialize sequence length
	l = 0;
	//loop until lim exchanges have been done
	while (l<lim) {
#ifdef MYDEBUG //check if nexc is ok
int debugSum1 = 0;
for (int debugi=0; debugi<nc; debugi++)
	debugSum1 += lc[debugi]*(n-lc[debugi]);
if (nexc!=debugSum1) {
	cout << "problem with nexc in randmergess" << endl;
	exit(1);
}
#endif
		//1st roulette wheel to peak cycle i (each cycle "i" has probability "(lc[i]*(n-lc[i]))/nexc")
		r = irand(nexc);
		for (i=0; i<nc; i++) {
			t = lc[i]*(n-lc[i]);	//#exchanges for cycle i
			if (r<t)                //k is the chosen cycle
				break;
			r -= t;                 //update r to implement roulette wheel
		}
		//2nd roulette wheel to peak cycle j (each cycle "j" has weight "lc[j] if j!=i, 0 if j==i" and the sum of weights is "n-lc[i]")
		r = irand(n-lc[i]);
		for (j=0; j<nc; j++) {
			t = j!=i ? lc[j] : 0;	//#exchanges for cycle j
			if (r<t)				//h is the chosen cycle
				break;
			r -= t;					//update r to implement roulette wheel
		}
#ifdef MYDEBUG //cycle i and j must be different
if (i==j) {
	cout << "randmergess: cycles i and j are equal and this is not possible!!!" << endl;
	exit(1);
}
#endif
		//store the exchange indexes in s (by selecting uniformly random indexes from cycles i and j) and increase l (note that EXC(x,y)==EXC(y,x))
		s[2*l] = c[i*n+irand(lc[i])];
		s[2*l+1] = c[j*n+irand(lc[j])];
		l++;
		//update-part1 nexc for next iteration (subtract "lc[i]*(n-lc[i]) + lc[j]*(n-lc[j])")
		nexc -= lc[i]*(n-lc[i]) + lc[j]*(n-lc[j]);
		//make i the longest cycle between cycles i and j
		if (lc[i]<lc[j]) {
			t = i;
			i = j;
			j = t;
		}
		//merge cycles i and j by appending items in cycle j to cycle i, and increase the length of cycle i
		memcpy(c+(i*n+lc[i]),c+(j*n),sizeof(int)*lc[j]);
		lc[i] += lc[j];
		//copy last cycle in cycle j (so cycle j is removed) and decrease the number of cycles
		nc--;
		memcpy(c+(j*n),c+(nc*n),sizeof(int)*lc[nc]);
		lc[j] = lc[nc];
		//update-part2 nexc for next iteration (add the exchanges for the new cycle i, i.e., add "lc[i]*(n-lc[i])")
		nexc += lc[i]*(n-lc[i]);
		//goto next iteration
	} //end of exchanges loop
#ifdef MYDEBUG
	if (l!=lim) { //limited version
		cout << "in randmergess l does not match lim" << endl;
		exit(1);
	}
#endif
	//done
}







//USED FOR F<1 (NO CHECK INSIDE)
//s/l = sequence that sorts x using exchanges (SWAP_k(i,j) is s[2*k]<->s[2*k+1])
//f is used to bound the insertions chain (it will be long l=ceil(f*original_l))
void randss(int* s, int& l, int* x, double f) {
	//some variables
	int i,nc,vc,nv,nexc,r,e,p1,p2,t,lc1,lc2,lim; //lim used in limited version
	double tempDouble; //temp variable for ceil rounding
	static bool* v = (bool*)(tmpint+n); //alive with x,c,lc (reuse part of s not used by x)
	static int* lc = tmpint+(2*n);      //alive with x,s,c,v (memory just after v)
	static int* c = tmpint+(3*n);       //alive with x,v,c,lc (last chunk of memory of size n^2)
	/* PHASE1 / Cycles decomposition of x (store only cycles **GREATER THAN 1**) / outputs will be:
		- nc is the TOTAL number of cycles of x
		- vc is the number of cycles of x GREATER THAN 1 (so the length of c)
		- c is a (n*n)-matrix s.t. c[i] contains the items of i-th cycle (g.t. 1) of x
		- lc is a n-vector s.t. lc[i] is the length of c[i]
		- nexc is the total number of exchanges for sorting x (it is SUM_{i=0}^{vc-1}(lc[i] choose 2))
		- COST: O(n)
	*/
	//initialize v to all false and other counters to 0 (also the first cycle length)
	memset(v,0,sizeof(bool)*n); //all false
	nc = vc = nexc = nv = i = lc[0] = 0;
	//scan all the items in the order given by the permutation (first open cycle is nc=0)
	for (;;) {
		//insert current index in the currently open cycle (and increase its length)
		c[vc*n+(lc[vc]++)] = i; //in matrix form = c[vc][lc[vc]++] = i;
		//set the visited flag of the i-th index and increase nv
		v[i] = true;
		nv++;
		//if all indexes have been visited, then close the cycle and exit the loop
		if (nv==n) {
			if (lc[vc]>1) { //update nexc and increase vc iff its length is greater than 1
				nexc += lc[vc]*(lc[vc]-1)/2; //(lc[vc] choose 2)
				vc++;
			}
			nc++;
			break;
		}
		//try to move index by following the cycle
		i = x[i];
		//check if cycle has to be closed
		if (v[i]) { //same as if (i==c[vc*n]) { //matrix form i==c[vc][0]
			//close the cycle
			if (lc[vc]>1) {
				nexc += lc[vc]*(lc[vc]-1)/2; //(lc[vc] choose 2)
				vc++;
			}
			nc++;
			//open next cycle by setting/resetting its length
			lc[vc] = 0;
			//move to next unvisited index (no need to move circularly, leftmost unv. index ensured)
			while (v[++i]);
		} //end-if
	} //end-infinite-for
#ifdef MYDEBUG
	int debugSum2 = 0;
	for (int debugi=0; debugi<vc; debugi++)
		debugSum2 += lc[debugi];
	debugSum2 += (nc-vc);
	if (debugSum2!=n) {
		cout << "sum of cycles lengths does not match n" << endl;
		exit(1);
	}
#endif
	/* PHASE2 / Sorting sequence of x / outputs will be:
		- l (==n-nc), i.e. the length of a decomposition (sorting sequence) of x (NO MORE BECAUSE OF LIMITED VERSION)
		- s, i.e. a sorting sequence of x such that:
			x^-1 = EXC(s[0],s[1]) * ... * EXC(s[2*(l-1)],s[2*(l-1)+1]) = PROD_{i=0}^{l-1} EXC(s[2*i],s[2*i+1])
			x = EXC(s[2*(l-1)],s[2*(l-1)+1]) * ... * EXC(s[0],s[1]) = PROD_{i=l-1}^{0} EXC(s[2*i],s[2*i+1])
		- does not sort x
		- COST: O(n^2)
	*/
	//initialize decomposition length to 0 (it will be n-nc at the end, but in lim. version f*(n-nc))
	l = 0;
	//compute limit for limited version
	//lim = f*(n-nc) + .5; //ceil rounding //(THERE WAS A BUG, it is rounding and not ceil rounding)
	lim = tempDouble = f*(n-nc);  //ceil rounding ok
	if (lim<tempDouble)           //ceil rounding ok
		lim++;                    //ceil rounding ok
	//loop until there are n cycles of length 1
	//while (vc>0) { //same as "while (l<n-nc)"  //version without limit
	while (l<lim) {                              //version with limit
#ifdef MYDEBUG
		int debugSum = 0;
		for (int debugi=0; debugi<vc; debugi++)
			debugSum += lc[debugi]*(lc[debugi]-1)/2;
		if (debugSum!=nexc) {
			cout << "sum of (lc[i] choose 2) does not match nexc" << endl;
			exit(1);
		}
#endif
		//select a random cycle i using roulette wheel where lc[i]s are the slice sizes
		r = irand(nexc);
		for (i=0; i<vc; i++) {
			e = lc[i]*(lc[i]-1)/2;  //#exchanges in cycle i
			if (r<e)                //i is the chosen cycle
				break;
			r -= e;                 //update r to implement roulette wheel
		}
		//decrease nexc by the exchanges of cycle i
		nexc -= lc[i]*(lc[i]-1)/2; //(lc[i] choose 2)
		//select two random indexes in cycle i
		twoRandIndices(lc[i],p1,p2);
		//make p1<p2 (note that EXC(x,y)==EXC(y,x))
		if (p1>p2) {
			t = p1;
			p1 = p2;
			p2 = t;
		}
		//store the exchange indexes in s and increase l (note that EXC(x,y)==EXC(y,x))
		s[2*l] = c[n*i+p1];
		s[2*l+1] = c[n*i+p2];
		l++;
		//two new cycles: C'=[0,p1]+(p2,lc[i]) and C''=(p1,p2]
		//compute length of C' (lc1) and C'' (lc2)
		lc1 = p1-p2+lc[i];   //len[0,p1] + len(p2,lc[i]) = p1+1 + lc[i]-p2-1 = p1-p2+lc[i]
		lc2 = p2-p1;         //len(p1,p2]
		//4 cases:
		//(1) lc1>1 and lc2>1: update cycle i with C' and add cycle C'' at the end
		//(2) lc1>1 and lc2=1: update cycle i with C'
		//(3) lc1=1 and lc2>1: update cycle i with C''
		//(4) lc1=1 and lc2=1: remove cycle i by moving last cycle in position i
		if (lc1>1 && lc2>1) { //case 1
			//add cycle C'' at the last position and increment vc
			memcpy(c+n*vc,c+n*i+p1+1,sizeof(int)*lc2);
			lc[vc] = lc2;
			vc++;
			//update cycle i with C' by moving the slice (p2,lc[i]) (if not empty) after p1
			t = lc[i]-p2-1;   //len(p2,lc[i])
			if (t>0)
				memmove(c+n*i+p1+1,c+n*i+p2+1,sizeof(int)*t);
			lc[i] = lc1;
			//increase nexc by the exchanges of C' and C''
			nexc += (lc1*(lc1-1)+lc2*(lc2-1))/2; //(lc1 choose 2)+(lc2 choose 2)
			//case 1 done
		} else if (lc1>1 && lc2==1) { //case 2
			//update cycle i with C' by moving the slice (p2,lc[i]) (if not empty) after p1
			t = lc[i]-p2-1;   //len(p2,lc[i])
			if (t>0)
				memmove(c+n*i+p1+1,c+n*i+p2+1,sizeof(int)*t);
			lc[i] = lc1;
			//increase nexc by the exchanges in C'
			nexc += lc1*(lc1-1)/2; //(lc1 choose 2)
			//case 2 done
		} else if (lc2>1) { //case 3
			//update cycle i with C'' by moving the slice (p1,p2] at position 0
			memmove(c+n*i,c+n*i+p1+1,sizeof(int)*lc2);
			lc[i] = lc2;
			//increase nexc by the exchanges in C''
			nexc += lc2*(lc2-1)/2; //(lc2 choose 2)
			//case 3 done
		} else { //case 4
			//move last cycle (vc-1)th in position i and decrease vc (if cycle i is not the last)
			vc--;
			if (i!=vc) {
				memcpy(c+n*i,c+n*vc,sizeof(int)*lc[vc]);
				lc[i] = lc[vc];
			}
			//case 4 done
		} //end-cases
	} //end-while
#ifdef MYDEBUG
	if (l!=lim) { //limited version (prima era l!=n-nc)
		cout << "in randss l does not match lim" << endl;
		exit(1);
	}
#endif
	//done
}



//differential mutation of individual i
void diffMutation(int i) { //EXCHANGES!!!
	//DE/rand/1: y1[i] = x[r0] + F * (x[r1] - x[r2])
	//initialize variables
	int r0,r1,r2,*p,*p0,*pp,*ixx,lss,j,exci,excj,t;
	static int* ss = tmpint;	//temp memory (max length 2*n)	//EXCHANGES
	static int* z = tmpint;		//tempmem (n) alive with ss so z DESTROYED ON RANDSS or RANDMERGESS //EXCHANGES!!!
	//(1) get 3 different indexes r0,r1,r2 different also from i
	threeRandIndicesDiffFrom(np,i,r0,r1,r2);
	//(2) compute scale factor and truncation bound using jde rule (need before randss for limit version)
	sfy[i] = urand()<.1 ? fmin+(fmax-fmin)*urandi() : sfx[i];
	//(3) check and handle the easy case F=1
	if (sfy[i]==1.) {
		//begin case F=1
		//compute y1[i] = x[r0] + (x[r1] - x[r2]) = x[r0] * x[r2]^-1 * x[r1]
		p = y1[i];
		p0 = x[r0];
		ixx = ix[r2];
		pp = x[r1];
		for (j=0; j<n; j++)
			p[j] = p0[ixx[pp[j]]];
#ifdef MYDEBUG
		int myp[n];
		for (j=0; j<n; j++)
			myp[j] = x[r0][ix[r2][x[r1][j]]];
		for (j=0; j<n; j++)
			if (y1[i][j]!=myp[j]) {
				cout << "Problem with the case F=1" << endl;
				exit(1);
			}
#endif
		//nothing else to do (this is an easy case) so return
		return;
		//end case F=1
	}
	//(4) distinguish the 2 cases F<1 and F>1
	if (sfy[i]<1.) {
		//begin case F<1: build z, randss on z, compute starting y
		//(5a) build z that will be input for randss (consider x[r1]-x[r2] = x[r2]^-1*x[r1] = sort_seq.(x[r1]^-1*x[r2])) thus z=x[r1]^-1*x[r2]
		ixx = ix[r1];
		p = x[r2];
		for (j=0; j<n; j++)
			z[j] = ixx[*p++]; //same as ixx[p[j]] == ixx[x[r2][j]]
		//(6a) ss,lss = randss(z) using limited version (sfy[i])
		randss(ss,lss,z,sfy[i]);
		//(7a) compute starting y, i.e., copy x[r0] on y1[i]
		memcpy(y1[i],x[r0],permByteSize);
		//end case F<1
	} else {
		//begin case F>1: build z, randss on z, compute starting y
		//(5b) build z that will be input for randmergess (consider I have to go from x[r1]-x[r2] = x[r2]^-1*x[r1] towards a 1cycle permutation) thus z=x[r2]^-1*x[r1]
		//(7b anticipated) compute starting y, i.e., assign to y1[i] the composition x[r0]*x[r2]^-1*x[r1] == x[r0]*z
		ixx = ix[r2];
		p = x[r1];
		p0 = x[r0];
		pp = y1[i];
		for (j=0; j<n; j++) {
			z[j] = ixx[*p++]; //same as ixx[p[j]] == ixx[x[r1][j]]
			pp[j] = p0[z[j]]; // same as y1[i][j] = x[r0][z[j]];
		}
		//(6b) ss,lss = randmergess(z) using limited version (sfy[i])
		randmergess(ss,lss,z,sfy[i]);
		//end case F>1
	}
	//(8) apply all the exchanges in ss to y1[i], since we are using limited version of randss
	p = y1[i];
	pp = ss; //... since ss is a static variable
	for (j=0; j<lss; j++) {
		exci = *pp++;
		excj = *pp++;
		t = p[exci];
		p[exci] = p[excj];
		p[excj] = t;
	}
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"diffMutation y1["<<i<<"]"<<endl;
		exit(1);
	}
#endif   
	//done	
}


//crossover TPII
inline void tpii(int* ch, int* fa, int* mo, int* ifa, int c1, int c2) {
	//child takes [c1,c2] from father, the other from mother using the order in mother
	//ifa is the inverse of fa
	//variables
	int k,t;
	//father part (between c1 and c2, both included)
	memcpy(ch+c1,fa+c1,sizeof(int)*(c2-c1+1));
	//left part (before c1) from mother (with j initialization to zero)
	//j = 0; //... needed by classical impl commented
	for (k=0; k<c1; k++) {
		while (ifa[(t=*mo++)]>=c1 && ifa[t]<=c2); //... same of what is commented below
		//do {
		//   t = *mo++; //... same as t = mo[j++];
		//} while (ifa[t]>=c1 && ifa[t]<=c2);
		*ch++ = t; //... same as ch[k] = t;
	}
	//right part (after c2) from mother (j was already set to the right position)
	ch += c2-c1; //... now ch is &ch_original[c2]
	for (k=c2+1; k<n; k++) {
		while (ifa[(t=*mo++)]>=c1 && ifa[t]<=c2); //... same of what is commented below
		//do {
		//   t = *mo++; //... same as t = mo[j++];
		//} while (ifa[t]>=c1 && ifa[t]<=c2);
		*++ch = t; //... same as ch[k] = t; //note that ch before the for was set to &ch_original[c2]
	}
	//done
}


//crossover of individual
void crossover(int i) {
	//TWO POINT CROSSOVER OF AGA: it produces two sons y1[i], y2[i] from the parents x[i], y1[i]
	//a son takes [c1,c2] from father, the other from mother using the order in mother
	//variables
	static int* tch = tmpint;     //temporary child for offspring 2
	static int* imut = tmpint+n;  //mutant inverse for offspring 2
	int c1,c2,t,*p;
	//two cut points c1<=c2 (also equal cutpoints)
	c1 = irand(n);
	c2 = irand(n);
	if (c1>c2) {
		t = c1;
		c1 = c2;
		c2 = t;
	}
	//(1) x[i] is father, y1[i] is mother, y2[i] becomes son
	tpii(y2[i],x[i],y1[i],ix[i],c1,c2);
#ifdef MYDEBUG
	if (!permValid(y2[i],n)) {
		cout<<"crossover y2["<<i<<"]"<<endl;
		exit(1);
	}
#endif
	//(2) y1[i] is father, x[i] is mother, y1[i] becomes son (computing father inverse and using a temp array)
	p = y1[i];
	for (t=0; t<n; t++)
		imut[*p++] = t; //... same as imut[p[t]] = t;
#ifdef MYDEBUG
	if (!permValid(imut,n)) {
		cout<<"crossover imut"<<endl;
		exit(1);
	}
#endif
	tpii(tch,y1[i],x[i],imut,c1,c2);
	memcpy(y1[i],tch,permByteSize);
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"crossover y1["<<i<<"]"<<endl;
		exit(1);
	}
#endif
	//done
}


//selection of individual i
void selection(int i) {
	//variables
	int* t;
	//crisp pre-selection between y1[i]/y2[i] (exchange pointers without copy)
	if (fy2[i]<fy1[i]) { //tie favors y1[i]
		fy1[i] = fy2[i];
		t = y1[i];
		y1[i] = y2[i];
		y2[i] = t;
#ifdef MYDEBUG
		if (!permValid(y1[i],n)) {
			cout<<"selection1 y1["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(y2[i],n)) {
			cout<<"selection1 y2["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//alpha-selection between y1[i]/x[i] (exchange pointers without copy)
	//the manifest constants TFT and MAKESPAN are mutually exclusive
#ifdef TFT
	if ( fy1[i]<fx[i] || urand()<alpha-REL_DEV(fy1[i],fx[i]) ) { //tie favors, in some way, x[i]
#endif
#ifdef MAKESPAN
	if ( fy1[i]<=fx[i] || urand()<alpha-REL_DEV(fy1[i],fx[i]) ) { //tie favors, in some way, y1[i]
#endif
		fx[i] = fy1[i];
		t = y1[i];
		y1[i] = x[i];
		x[i] = t; //t is now the new individual
		sfx[i] = sfy[i];  //jde rule
		//since there was replacement, update ix[i] by computing the inverse
		int* inv = ix[i];
		for (int k=0; k<n; k++)
			inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
#ifdef MYDEBUG
		if (!permValid(y1[i],n)) {
			cout<<"selection2 y1["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(x[i],n)) {
			cout<<"selection2 x["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(ix[i],n)) {
			cout<<"selection2 ix["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//done
}


inline void restartStatistics() {
	//compute restart statistics (call after increasing nrestarts)
	int stageLength = ngen-lastRestart;
	if (stageLength<minStageLength)
		minStageLength = stageLength;
	if (stageLength>maxStageLength)
		maxStageLength = stageLength;
	avgStageLength += (stageLength-avgStageLength)/nrestarts; //moving average (see wikipedia)
	lastRestart = ngen;
	if (fgbest<fgbestAtStageStart)
		improvingStages++;
	fgbestAtStageStart = fgbest;
	//done
}


//restart population if restart has been triggered
bool popRestart() {
	//if all fitnesses are not the same do nothing
	if (!sameFitness)
		return false; //false since maxnfes was not exhausted
	//the first (it's a best) individual is keeped and the other are randomized
	for (int i=1; i<np; i++) {
		prand(n,x[i]);
		fx[i] = INT_MAX;  //+infinity
		sfx[i] = finit;
		//since there was replacement, update ix[i] by computing inverse
		int* inv = ix[i];
		int* t = x[i];
		for (int k=0; k<n; k++)
			inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
#ifdef MYDEBUG
		if (!permValid(x[i],n)) {
			cout<<"restart x["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(ix[i],n)) {
			cout<<"restart ix["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//anyway reset sfx[0]
	sfx[0] = finit;
	//local search (vns4) on x[0]
	if (ls==B_LS) { //baldwin
		memcpy(tmpint,x[0],permByteSize);
		int ft = fx[0];
		vns4(tmpint,ft);
	} else if (ls==L_LS) {//lamarck
		vns4(x[0],fx[0]);
		for (int k=0; k<n; k++)
			ix[0][x[0][k]] = k;
	}
	//update nfesWhenToForceRestart, thus if normal restart don't do forcedrestart
	nfesWhenToForceRestart = nfes + forcedRestartPeriod;
	//increase restarts counter
	nrestarts++;
	//restart statistics
	restartStatistics();
	//return true/false if budget of nfes was exhausted in the local search
	return termination();
	//done
}


//vns4 local search: full interchange ls (1st imp. style) + one step insert ls (best imp. style)
void vns4(int* x, int& fx) {
	//if x is in the history, return
	int i;
	for (i=0; i<n; i++)
		if (lsHistory[i]!=x[i])
			break;
	if (i==n)
		return;
	//set lsmode
	lsmode = true;
	//save the fitness of the seed for statistics
	int fseed = fx;
	//perform full interchange + insertion step
	int fprev;
	do {
		fprev = fx;
#ifdef GFC
		bool firstInt = true; //used the first time to call computeGFC
		while (intStep(x,fx,firstInt))   //interchange local search (using gfc)
			firstInt = false;
#else
		while (intStep(x,fx));           //interchange local search (without gfc)
#endif
		if (termination()) //break the do-while if termination
			break;
		insStep(x,fx);                   //insertion step
		if (termination()) //break the do-while if termination
			break;
	} while (fx<fprev);
	//save x in the history
	memcpy(lsHistory,x,permByteSize);
	//update local search statistics
	nls++;
	if (fx<fseed) {
		nImprovingls++;
		totImprovingls += fseed-fx;
	}
	//unset lsmode
	lsmode = false;
	//done
}


#ifdef GFC
//VERSIONI DI INTSTEP E INSSTEP CHE USANO GFC!!!!!!!!!!!!!!
bool intStep(int* x, int& fx, bool first) {
	//GFC VERSION
	//one step of interchange local search
	//1st improvement style using a random permutation
	//variables
	const static int bits = sizeof(int)*8-1;
	static int* r = tmpint+n; //since tmpint is used to copy x in the case of B_LS
	static int nm1 = n-1;
	int i,j,ii,jj,t,fseed;
	//compute gfc data of x BUT ONLY IF FIRST application in vns4!!!
	if (first)
		computeGFC(x);
	//set fitness of the seed and get a random permutation
	fseed = fx;
	prand(n,r);
	//scan for interchange moves
	for (ii=0; ii<nm1; ii++) {
		for (jj=ii+1; jj<n; jj++) {
			//get i,j from the random permutation
			i = r[ii];
			j = r[jj];
			//interchange jobs at positions i and j in x
			t = x[i];
			x[i] = x[j];
			x[j] = t; //now t is equal to x[j]
			//evaluate x using gfc, update gbest and check termination
			fx = evalGFC(x,FAST_MIN(i,j,bits));
#ifdef MYDEBUG
			if (fx!=eval(x)) {
				cerr<<"intStep evaluation mismatch"<<endl;
				exit(1);
			}
#endif
			nfes++;
			nfesls++;
			updateGbest(x,fx);
			if (termination())
				return false;
			//if better UPDATE gfc data and return true, x and fx are modified
			if (fx<fseed) {
				evalUpdateGFC(x,FAST_MIN(i,j,bits));
				return true;
			}
			//if here, not better, so reset the interchange above by redoing it
			x[j] = x[i];
			x[i] = t; //remember that t was x[j]
			//done
		}
	}
	//if here, no return above, so fx has not been improved, thus reset fx and return false
	fx = fseed;
	return false;
	//done
}


void insStep(int* x, int& fx) {
	//GFC VERSION
	//one step of insertion local search
	//best improvement style, so no random permutation
	//(i,j) means "move job at pos i to pos j"
	//variables
	int i,j,fbest,bi,bj,ft;
	//no need for computing gfc data since intStep did it!!!
	//computeGFC(x)
	//init fbest to fx of the seed
	fbest = fx;
	//scan for insertion by doing rightmost insertions before
	for (i=n-1; i>=0; i--) {
		for (j=i+1; j<n; j++) { //do (i,j) and (j,i)
			//(i,j), since i<j for sure this is a forward insertion
			fwdins(x,i,j);
			//evaluate x using gfc, update gbest and check termination
			ft = evalGFC(x,i); //since i<j
#ifdef MYDEBUG
			if (ft!=eval(x)) {
				cerr<<"insStep evaluation mismatch"<<endl;
				exit(1);
			}
#endif
			nfes++;
			nfesls++;
			updateGbest(x,ft);
			if (termination())
				return;
			//update fbest,bi,bj if better ft is better
			if (ft<fbest) { //with <= I simulate < with the scan from 0 to n
				fbest = ft;
				bi = i;
				bj = j;
			}
			//before doing (j,i) check if it's case that (j,i)!=(i,j)
			if (i+1!=j) { //!= means <
				//perform two times (j,i): one to undo (i,j) and one to do (j,i)
				bwdins2(x,j,i);
				//evaluate x using gfc, update gbest and check termination
				ft = evalGFC(x,i); //since i<j
				nfes++;
				nfesls++;
				updateGbest(x,ft);
				if (termination())
					return;
				//update fbest,bi,bj if better ft is better
				if (ft<fbest) { //with <= I simulate < with the scan from 0 to n
					fbest = ft;
					bi = j; //inverted since the insertion is (j,i)
					bj = i; //inverted since the insertion is (j,i)
				}
				//undo (j,i) by doing (i,j)
				fwdins(x,i,j);
			} else //in the case that (j,i) was not performed, simply undo (i,j) by doing (j,i)
				bwdins(x,j,i);
			//done both (i,j) and (j,i)
		}
	}
	//if there was an improvement wrt the seed solution, update the solution by doing (bi,bj) insertion
	if (fbest<fx) {
		ins(x,bi,bj);
		fx = fbest;
	}
	//done
}

#else

//VERSIONI DI INTSTEP E INSSTEP CHE ***NON*** USANO GFC!!!!!!!!!!!!!!
bool intStep(int* x, int& fx) {
	//NO-GFC VERSION
	//one step of interchange local search
	//1st improvement style using a random permutation
	//variables
	static int* r = tmpint+n; //since tmpint is used to copy x in the case of B_LS
	static int nm1 = n-1;
	int i,j,ii,jj,t,fseed;
	//set fitness of the seed and get a random permutation
	fseed = fx;
	prand(n,r);
	//scan for interchange moves
	for (ii=0; ii<nm1; ii++) {
		for (jj=ii+1; jj<n; jj++) {
			//get i,j from the random permutation
			i = r[ii];
			j = r[jj];
			//interchange jobs at positions i and j in x
			t = x[i];
			x[i] = x[j];
			x[j] = t; //now t is equal to x[j]
			//evaluate x, update gbest and check termination
			fx = eval(x);
			nfes++;
			nfesls++;
			updateGbest(x,fx);
			if (termination())
				return false;
			//if better return true, x and fx are modified
			if (fx<fseed)
				return true;
			//if here, not better, so reset the interchange above by redoing it
			x[j] = x[i];
			x[i] = t; //remember that t was x[j]
			//done
		}
	}
	//if here, no return above, so fx has not been improved, thus reset fx and return false
	fx = fseed;
	return false;
	//done
}


void insStep(int* x, int& fx) {
	//NO-GFC VERSION
	//one step of insertion local search
	//best improvement style, so no random permutation
	//(i,j) means "move job at pos i to pos j"
	//variables
	int i,j,fbest,bi,bj,ft;
	//init fbest to fx of the seed
	fbest = fx;
	//scan for insertion moves in the same order of GFC version
	for (i=n-1; i>=0; i--) {
		for (j=i+1; j<n; j++) {
			ins(x,i,j);
			ft = eval(x);
			nfes++;
			nfesls++;
			updateGbest(x,ft);
			if (termination())
				return;
			if (ft<fbest) {
				fbest = ft;
				bi = i;
				bj = j;
			}
			ins(x,j,i);
			if (i+1!=j) {
				ins(x,j,i);
				ft = eval(x);
				nfes++;
				nfesls++;
				updateGbest(x,ft);
				if (termination())
					return;
				if (ft<fbest) {
					fbest = ft;
					bi = j; //note that j and i are reversed here!!!
					bj = i; //note that j and i are reversed here!!!
				}
				ins(x,i,j);
			}
		}
	}
	/*
	for (i=n-1; i>=0; i--) {
		for (j=n-1; j>=0; j--) {
			if (i!=j && i!=j+1) { //(i,i+1)==(i+1,i) so avoid to redo
				//perform the insertion (i,j) on x
				ins(x,i,j);
				//evaluate x, update gbest and check termination
				ft = eval(x);
				nfes++;
				nfesls++;
				updateGbest(x,ft);
				if (termination())
					return;
				//update fbest,bi,bj if better ft is better
				if (ft<fbest) {
					fbest = ft;
					bi = i;
					bj = j;
				}
				//undo the insertion (i,j) on x by doing the insertion (j,i) on x
				ins(x,j,i);
			}
		}
	}
	*/
	//if there was an improvement wrt the seed solution, update the solution by doing (bi,bj) insertion
	if (fbest<fx) {
		ins(x,bi,bj);
		fx = fbest;
	}
	//done
}
#endif


#ifdef GFC
inline void fwdins(int* x, int i, int j) {
	//assume i<j and perform insertion (i,j) on x (this is a forward insertion)
	int t = x[i];
	memmove(x+i,x+i+1,sizeof(int)*(j-i));
	x[j] = t;
}


inline void bwdins(int* x, int i, int j) {
	//assume i>j and perform insertion (i,j) on x (this is a backward insertion)
	int t = x[i];
	memmove(x+j+1,x+j,sizeof(int)*(i-j));
	x[j] = t;
}


inline void bwdins2(int* x, int i, int j) {
	//assume i>j+1 and perform insertion (i,j) two times on x (these are two backward insertions)
	//note that, if i==j+1 then x remains unchanged but this case is not handled here
	int t1 = x[i-1];
	int t2 = x[i];
	memmove(x+j+2,x+j,sizeof(int)*(i-j-1));
	x[j] = t1;
	x[j+1] = t2;
}
#endif


inline void ins(int* x, int i, int j) {
	//perform insertion (i,j) on x assuming nothing on i and j
	//(i,j) means "move job at pos i to pos j"
	int t = x[i];
	if (i<j) //forward
		memmove(x+i,x+i+1,sizeof(int)*(j-i));
	else //backward
		memmove(x+j+1,x+j,sizeof(int)*(i-j));
	x[j] = t;
}


bool popForcedRestart() {
	//check if restart has to be forced, return if not the case
	if (nfes<nfesWhenToForceRestart)
		return false; //false since maxnfes was not exhausted
	//find best individual
	int ibest = 0;
	int i;
	for (i=1; i<np; i++)
		if (fx[i]<fx[0])
			ibest = i;
	//restart population except ibest
	for (i=0; i<np; i++) {
		if (i!=ibest) {
			prand(n,x[i]);
			fx[i] = INT_MAX; //+inf
			//update ix
			int* inv = ix[i];
			int* t = x[i];
			for (int k=0; k<n; k++)
				inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
#ifdef MYDEBUG
			if (!permValid(x[i],n)) {
				cout<<"forcedRestart x["<<i<<"]"<<endl;
				exit(1);
			}
			if (!permValid(ix[i],n)) {
				cout<<"forcedRestart ix["<<i<<"]"<<endl;
				exit(1);
			}
#endif
		}
		sfx[i] = finit; //finit also for the best
	}
	//lamarckian local search (vns4) on x[ibest]
	vns4(x[ibest],fx[ibest]);
	for (int k=0; k<n; k++)
		ix[ibest][x[ibest][k]] = k;
	//update nfesWhenToForceRestart
	nfesWhenToForceRestart = nfes + forcedRestartPeriod;
	//update statistics
	nrestarts++;
	nforcedrestarts++;
	restartStatistics();
	//return true/false if budget of nfes was exhausted in the local search
	return termination();
	//done
}

