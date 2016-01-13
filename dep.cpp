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
static int* tmpint;  //length = n*(n-1)/2 + 2*n (one sorting sequence + two permutations)

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
inline void randbs(int* s, int& l, int* x, int* ainv, int nainv);
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
	tmpint = new int[n*(n-1)/2+2*n]; //one sorting sequence + two permutations
	//set the permutation byte size
	permByteSize = sizeof(int)*n;
	//lsHistory
	lsHistory = new int[n];
	memset(lsHistory,0,permByteSize); //all zeros
	//init to false the save flag
	haveToSave = false;
	//init the diameter
	diameter = n*(n-1)/2; //VALID FOR ADJACENT SWAPS!!!
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


//s/l = sequence that sorts x, ainv/nainv = precomputed adjacent inversions of x
inline void randbs(int* s, int& l, int* x, int* ainv, int nainv) {
	int i,j,t;
	//initialize sequence length to 0
	l = 0;
	//sorting loop
	while (nainv>0) {
		//randomly select an adjacent inversion (i,j)=(i,i+1)
		t = irand(nainv);
		i = ainv[t];
		j = i+1;
		//remove (i,j), i.e. i, from the list of adjacent inversions, i.e. ainv
		ainv[t] = ainv[--nainv];
		//add previous inversion, i.e. (i-1,i)=(t,i), to ainv if the case
		t = i-1;
		if (i>0 && x[t]<x[i] && x[t]>x[j])
			ainv[nainv++] = t;
		//add next inversion, i.e. (i,i+1)=(j,t), to ainv if the case
		t = j+1;
		if (t<n && x[t]<x[i] && x[t]>x[j])
			ainv[nainv++] = j;
		//apply the swap (i,j) to x
		t = x[i];
		x[i] = x[j];
		x[j] = t;
		//save the swap in s (only i)
		s[l++] = i;
	}
	//done
}


//differential mutation of individual i
void diffMutation(int i) {
	//DE/rand/1: y1[i] = x[r0] + F * (x[r1] - x[r2])
	//(0) initialize variables
	int r0,r1,r2,t,k,j,w,*p,*ixx,*pt,*pp,*p0;
	double tempDouble;                           //temp variable for ceil rounding
	static int* ss = tmpint;                     //use the temp memory (sorting sequence)
	static int* z = tmpint+(n*(n-1)/2);          //use the temp memory (just at right of ss)
	static int* ainv = z+n;                      //use the temp memory (just at right of z)
//int mydebug;//debug
	//(1) get 3 different indexes r0,r1,r2 different also from i
	threeRandIndicesDiffFrom(np,i,r0,r1,r2);
	//(2) compute scale factor using jde rule
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
		//begin case F<1: build z with its adjacent inversions, randbs on z, compute k, compute starting y
		//(5a) build z with its adjacent inversions ainv/k (k is #ainv)
		//z = x[r1]-x[r2] = x[r2]^-1*x[r1] = sort_seq(x[r1]^-1*x[r2]) (I already know x[r1]^-1=ix[r1])
		ixx = ix[r1];
		p = x[r2];
		*z = ixx[*p];           //first outside the for ... it's the same of "z[0] = ixx[p[0]]"
		k = 0;                  //init no. of adjacent inversions
		pt = ainv;              //due because ainv is static
		for (j=1; j<n; j++) {   //from the second to the last inside the for (so the if is ok)
			z[j] = ixx[*++p];    //"ixx[*++p]" is the same of "ixx[p[j]]"
			if (z[j-1]>z[j]) {
				*pt++ = j-1;      //considering also next line, it is the same of "ainv[k++]=j-1"
				k++;
			}
		}
		//(6a) ss,w = sort_sequence(z) using randbs sorting loop
		randbs(ss,w,z,ainv,k);
		//(7a) compute k, i.e. where to truncate ss, basing on w and sfy[i] and using ceil rounding
		k = tempDouble = sfy[i] * w;  //tempDouble needed to implement ceil rounding
		if (k<tempDouble)             //this is to implement ceil rounding
			k++;
		//(8a) compute starting y: y1[i] becomes a copy of x[r0]
		memcpy(y1[i],x[r0],permByteSize);
		//end case F<1
	} else {
		//begin case F>1: build z with its adjacent inversions, compute starting y, randbs on z, compute k
		//(5b) build z with its adjacent inversions ainv/k (k is #ainv)
		//DEC(R) = SORT(E->X) + SORT(X->R) = X + SORT(X->R) = X + SORT(R^-1*X) = X + SORT(R*X)
		//      thus the z to sort has to be R*(x[r1]-x[r2]) = R * x[r2]^-1 * x[r1]
		//      note: the multiplication at left by R can be simulated, i.e. (R*X)[i] = R[X[i]] = n-1-X[i]
		//(8b - anticipated) together compute also starting y, i.e. x[r0] * x[r2]^-1 * x[r1]
		ixx = ix[r2];
		p = x[r1];
		*z = n-1-ixx[*p];       //first outside the for ... it's the same of "z[0] = n-1-ixx[p[0]]"
		p0 = x[r0];             //for starting y (8b)
		pp = y1[i];             //for starting y (8b)
		*pp = p0[ixx[*p]];      //for starting y (8b), 1st outside for, equal to y1[i][0]=x[r0][ixx[x[r1][0]]]
		k = 0;                  //init no. of adjacent inversions
		pt = ainv;              //due because ainv is static
		for (j=1; j<n; j++) {   //from the second to the last inside the for (so the if is ok)
			z[j] = n-1-ixx[*++p];//"ixx[*++p]" is the same of "ixx[p[j]]"
			pp[j] = p0[ixx[*p]]; //for starting y (8b), equal to "y1[i][j] = x[r0][ixx[x[r1][j]]]"
			if (z[j-1]>z[j]) {
				*pt++ = j-1;      //considering also next line, it is the same of "ainv[k++]=j-1"
				k++;
			}
		}
		//(6b) ss,w = sort_sequence(z) using randbs sorting loop
		randbs(ss,w,z,ainv,k);
#ifdef MYDEBUG
		//I have to check that "x[r2]^-1 * x[r1] * ss = reverse_identity"
		int myp[n];
		for (j=0; j<n; j++)
			myp[j] = ix[r2][x[r1][j]];
		for (j=0; j<w; j++) {
			t = myp[ss[j]];
			myp[ss[j]] = myp[ss[j]+1];
			myp[ss[j]+1] = t;
		}
		for (j=0; j<n; j++)
			if (myp[j]!=n-1-j) {
				cout << "Problem with case F>1" << endl;
				exit(1);
			}
#endif
		//(8b) compute k, i.e. where to truncate ss, basing on w and sfy[i] and using ceil rounding
		//in this case (F>1), L=D-W, K = min{ceil(F*L),D}-L = min{ceil(F*(D-W)),D}-D+W
		k = tempDouble = sfy[i]*(diameter-w);  //tempDouble needed to implement ceil rounding
		if (k<tempDouble)                      //this is to implement ceil rounding
			k++;
		if (k>diameter)                        //this is the minimum part of the formula above
			k = diameter;
		k += w - diameter;                     //this is the last "-D+W" part
//mydebug=w;//debug
#ifdef MYDEBUG
		if (k<0 || k>w) {
			cout << "k is not ok in the case F>1" << endl;
			exit(1);
		}
#endif
		//end case F>1
	}
	//(9) apply the first k entries of ss to y1[i]
	p = y1[i];
	pt = ss;             //due beecause ss is a static variable
	for (j=0; j<k; j++) {
		w = *pt++;        //"*pt++" is the same of "ss[j]"
		t = p[w];         //this line and the next two implement the adjacent swap
		p[w] = p[w+1];
		p[w+1] = t;
	}
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"diffMutation y1["<<i<<"]"<<endl;
		exit(1);
	}
	//to check the following please uncomment the two lines with the variable "mydebug" above
	//check the correct number of generators
	//int myz[n];    //put x[r1]^-1 * x[r2] in myz
	//for (j=0; j<n; j++)
	//   myz[j] = ix[r1][x[r2][j]];
	//int myai[n];   //put the adj. inv. of myz in myai
	//int mynai = 0;
	//for (j=1; j<n; j++)
	//   if (myz[j-1]>myz[j])
	//      myai[mynai++] = j-1;
	//int myss[n*n]; //put the sort. seq. of myz in myss
	//int mylss;
	//randbs(myss,mylss,myz,myai,mynai);
	//int mygen;     //compute the no. of generators to check
	//double mytd;
	//mygen = mytd = sfy[i]*mylss;
	//if (mygen<mytd)
	//   mygen++;
	//if (sfy[i]<1.) { //check if everything ok - case f<1
	//   if (k!=mygen) {
	//      cout << "Number of generators mismatch when F<1" << endl;
	//      exit(1);
	//   }
	//} else {       //check if everything ok - case f>1
	//   //L=D-W, K = min{ceil(F*L),D}-L = min{ceil(F*(D-W)),D}-D+W
	//   if (mygen>diameter)
	//      mygen = diameter;
	//   int myl2 = mydebug; //length of 2nd part
	//   int myl1 = diameter-myl2; //length of 1st part
	//   if (myl1!=mylss) {
	//      cout << "Problem with lengths when F>1" << endl;
	//      exit(1);
	//   }
	//   if (mylss+k!=mygen) {
	//      cout << "Number of generators mismatch when F>1" << endl;
	//      exit(1);
	//   }
	//}
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

