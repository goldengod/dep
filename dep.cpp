#include "global.h"
#include "dep.h"
#include "problem.h"
#include "random.h"
#include "timer.h"
#include "utils.h"
#include "loglookup.h"
#include <cstring>   //memcpy, memset, memmove, strcpy, strcmp
#include <climits>   //INT_MAX
#include <unistd.h>  //alarm
#include <cstdio>    //FILE, fprintf
#include <iostream>  //cerr, endl
#include <cstdlib>   //exit, EXIT_FAILURE
using namespace std;

#ifdef MYDEBUG
#include "utils.h"
#endif

//definitions of input variables declared in dep.h
int np;						//population size
double finit;				//initial scale factor
double crinit;				//initial crossover probability (needed for obxcr crossover)
double alpha;				//alpha value (needed for alpha selection)
int heu;					//0=use heuristic in initialization, 1=do not use (needed for randheu initialization)
int ls;						//0=no local search, 1=baldwinian local search, 2=lamarkian local search (needed for restarts with local search)
double frfactor;			//used for forced restart (forced restart factor (0,1+) - restart forced if no gbest update for maxnfes*frfactor nfes)
double ffmin;				//min value for F parameter on jde rule
double ffmax;				//max value for F parameter on jde rule
double crmin;				//min value for CR parameter on jde rule (needed for obxcr crossover)
double crmax;				//max value for CR parameter on jde rule (needed for obxcr crossover)
char sgenerators[10];		//generating set in differential mutation
char sinitialization[10];	//initialization algorithm
char scrossover[10];		//crossover algorithm
char sselection[10];		//selection algorithm
char slsearch[10];			//local search algorithm
char srestart[10];			//restart algorithm
int inftype = 1;			//parameter for inf selection (1 or 2 or 3 or 4, 1 is default)
int nchilds;				//number of childs for the crossover (1 or 2, 2 is default)

//definitions of output variables declared in dep.h
int* gbest;					//global best so far
FitnessType fgbest;			//fitness of the global best so far
int nfesFoundAt;			//nfes when the global best has been found
unsigned long timeFoundAt;	//milliseconds when the global best has been found
int stageFoundAt;			//stage when the global best has been found
int nfes;					//number of evaluations performed
int ngen;					//number of generations performed
unsigned long execTime;		//for time save/resume scheme see save_resume.cpp
int nrestarts;				//number of restart performed
int nforcedrestarts;		//number of restart forced
int minStageLength;			//min length of an evolution stage (between two consecutive restarts)
int maxStageLength;			//max length of an evolution stage (between two consecutive restarts)
double avgStageLength;		//avg length of an evolution stage
int improvingStages;		//number of improving stages
int nfesls;					//evaluations spent in local searches
int nls;					//number of local searches performed
int nImprovingls;			//number of improving local searches
FitnessType totImprovingls;	//total delta of improving in local searches
bool gbestls;				//true if the gbest has been obtained in a local search
bool haveToSave;			//flag for saving execution state/memory
int improvingSteps;			//number of improvements so far
int lsImprovingSteps;		//number of improvements so far done by the local searches
double sfSuccAvg;			//f statistic
double sfSuccVar;			//f statistic
double sfSuccMin;			//f statistic
double sfSuccMax;			//f statistic
double crSuccAvg;			//cr statistic (needed for obxcr)
double crSuccVar;			//cr statistic (needed for obxcr)
double crSuccMin;			//cr statistic (needed for obxcr)
double crSuccMax;			//cr statistic (needed for obxcr)
int child1succ;				//two childs statistic
int child2succ;				//two childs statistic

//inner variables
static int** x;				//array of main individuals (Iliffe/display style)
static FitnessType* fx;		//array of main individuals' fitnesses (Iliffe/display style)
static int** ix;			//array of main individuals' inverses (Iliffe/display style)
static int** y1;			//array of childrens1 (Iliffe/display style)
static FitnessType* fy1;	//array of childrens1's fitnesses (Iliffe/display style)
static int** y2;			//array of childrens2 (Iliffe/display style)
static FitnessType* fy2;	//array of childrens1's fitnesses (Iliffe/display style)
static int* ps;				//population storage (Iliffe/display style) x-y1-y2-ix interleaved
static FitnessType* fs;		//fitness storage (Iliffe/display style) fx-fy1-fy2 interleaved
static double* sfx;			//main array of scale factors (Iliffe/display style)
static double* sfy;			//temp array of scale factors (Iliffe/display style)
static double* sfs;			//scale factor storage (Iliffe/display style) sfx-sfy interleaved
static double* crx;			//main array of crossover probabilities (Iliffe/display style) (needed by obxcr)
static double* cry;			//temp array of crossover probabilities (Iliffe/display style) (needed by obxcr)
static double* crs;			//crossover probability storage (Iliffe/display style) crx-cry interleaved
static int* tmpint;			//temporary memory
static bool sameFitness;	//flag indicating if fitnesses are all equal
static int lastRestart;		//generation when last restart happened (to count restart statistics)
static FitnessType fgbestAtStageStart;//fitness at the begin of the current stage (to count the number of improving stages)
static int permByteSize;	//permutation size in bytes (used for memcpy)
static int* lsHistory;		//local search history (only one perm)
static bool lsmode;			//flag indicating if local search is running or not
static unsigned long savedTime = 0;//for time save/resume scheme see save_resume.cpp
static int nfesWhenToForceRestart;//to manage forced restart
static int forcedRestartPeriod;//to manage forced restart
static int diameter;		//diameter of the generating set

//inner functions definitions
bool popEvolve();
inline void updateGbest(int* x, FitnessType fx);
inline bool termination();
void depSave();

//function pointers definitions
void (*popInit)(void);
void (*diffMutation)(int);
void (*crossover)(int);
void (*selection)(void);
void (*localSearch)(int*,FitnessType&);
bool (*popRestart)(void);
bool (*popForcedRestart)(void);

//import some code
#include "depCommon.cpp"
#include "depPopInit.cpp"
#include "depDiffMutation.cpp"
#include "depCrossover.cpp"
#include "depSelection.cpp"
#include "depLocalSearch.cpp"
#include "depPopRestart.cpp"

//import online printing functions
#ifdef ONLINEPRINT
#include "onlineprint.h"
#endif

//import save_resume code
#include "save_resume.cpp"



//set dep default parameters
void depDefaultParameters() {
//the manifest constants TFT/MAKESPAN/LOP/LOPCC are mutually exclusive
#ifdef TFT
	np = 100;
	finit = 0.5;
	crinit = 0.5; //anche se inutile dato che il crossover "standard" su pfsp non usa parametri
	alpha = 0.01;
	heu = 0; //sarebbe 1 ma metto 0 perche' devo passare il file euristica al main
	ls = B_LS;
	frfactor = 0.25;
	ffmin = 0.1; //as standard jde rule
	ffmax = 1.0; //as standard jde rule
	crmin = 0.0; //anche se non usato qui
	crmax = 1.0; //anche se non usato qui
	strcpy(sgenerators,"asw");
	strcpy(sinitialization,"randheu");
	strcpy(scrossover,"tpii");
	strcpy(sselection,"alpha");
	strcpy(slsearch,"vns4");
	strcpy(srestart,"randls");
	nchilds = 2;
#endif
#ifdef MAKESPAN
	np = 20;
	finit = 0.5;
	crinit = 0.5; //anche se inutile dato che il crossover "standard" su pfsp non usa parametri
	alpha = 0.01;
	heu = 0; //sarebbe 1 ma metto 0 perche' devo passare il file euristica al main
	ls = L_LS;
	frfactor = 0.25;
	ffmin = 0.1; //as standard jde rule
	ffmax = 1.0; //as standard jde rule
	crmin = 0.0; //anche se non usato qui
	crmax = 1.0; //anche se non usato qui
	strcpy(sgenerators,"asw");
	strcpy(sinitialization,"randheu");
	strcpy(scrossover,"tpii");
	strcpy(sselection,"alpha");
	strcpy(slsearch,"vns4");
	strcpy(srestart,"randls");
	nchilds = 2;
#endif
#ifdef LOP
	np = 50;
	finit = 0.5;
	crinit = 0.5;
	alpha = 0.0; //NOTA QUI!!!
	heu = 0; //non usata per lop
	ls = L_LS; //non c'e' differenza fra L_LS e B_ls con il restart shrandls default di lop
	frfactor = 0.25; //DA RIVEDERE!!!
	ffmin = 0.1; //as standard jde rule
	ffmax = 1.0; //as standard jde rule
	crmin = 0.0; //as standard jde rule
	crmax = 1.0; //as standard jde rule
	strcpy(sgenerators,"asw");
	strcpy(sinitialization,"randheu");
	strcpy(scrossover,"obxcr");
	strcpy(sselection,"crowding");
	strcpy(slsearch,"ins");
	strcpy(srestart,"shrandls");
	nchilds = 2;
#endif
#ifdef LOPCC
	np = 80;
	finit = 0.5;
	crinit = 0.5;
	alpha = 0.0; 					//non usato con crowding
	heu = 0;
	ls = L_LS;						//non c'e' differenza fra L_LS e B_ls con il restart shrandls default di lop e non usato con crowding
	frfactor = 1.00;				//COSI' NON FA MAI RESTART FORZATI!!!
	ffmin = 0.1;					//as standard jde rule
	ffmax = 1.2;
	crmin = 0.0;					//as standard jde rule
	crmax = 1.0;					//as standard jde rule
	strcpy(sgenerators,"exc");
	strcpy(sinitialization,"randheu");
	strcpy(scrossover,"obxcr");
	strcpy(sselection,"crowding");
	strcpy(slsearch,"ins");			//NON SUCCEDE MAI
	strcpy(srestart,"shrandls");	//NON SUCCEDE MAI
	nchilds = 1;
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



//allocate dep memory and decide function pointers
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
	fs = new FitnessType[3*np];
	fx = fs;
	fy1 = fx+np;
	fy2 = fy1+np;
	//main and temporary scale factor memory (if information based selection, we need more memory)
	sfs = new double[ strcmp(sselection,"inf")!=0 || strcmp(sselection,"fitsep")!=0 ? 2*np : 3*np ];
	sfx = sfs;
	sfy = sfx+np;
	//main and temporary crossover probability memory (if information based selection, we need more memory)
	crs = new double[ strcmp(sselection,"inf")!=0 || strcmp(sselection,"fitsep")!=0 ? 2*np : 3*np ];
	crx = crs;
	cry = crx+np;
	//set the permutation byte size
	permByteSize = sizeof(int)*n;
	//lsHistory
	lsHistory = new int[n];
	memset(lsHistory,0,permByteSize); //all zeros
	//init to false the save flag
	haveToSave = false;
	//initialize popInit function pointer
	if (strcmp(sinitialization,"randheu")==0)
		popInit = popInit_randheu;
	else {
		cerr << "ERROR: the initialization \"" << sinitialization << "\" is not one of: \"randheu\"" << endl;
		exit(EXIT_FAILURE);
	}
	//init diffMutation function pointer, diameter, temporary memory basing on generators
	if (strcmp(sgenerators,"asw")==0) {			//ADJACENT SWAPS - BUBBLE SORT
		diffMutation = diffMutationBS;
		diameter = n*(n-1)/2;
		tmpint = new int[n*(n-1)/2+2*n]; //one sorting sequence + two permutations (consider also the crossover)
	} else if (strcmp(sgenerators,"ins")==0) {	//INSERTIONS - INSERTION SORT
		diffMutation = diffMutationIS;
		diameter = n-1;
		tmpint = new int[6*n]; //one insertions sequence (2*n) + four permutations (consider also the crossover)
		if (ffmax>1.0)		//to address default configuration
			ffmax = 1.0;
		if (finit>1.0 || ffmin>1.0 || ffmax>1.0) { //check f params for insertions!
			cerr << "FINIT, MINF, MAXF HAS TO BE <= 1.0 FOR INSERTIONS!!!" << endl;
			exit(EXIT_FAILURE);
		}
	} else if (strcmp(sgenerators,"exc")==0) {	//EXCHANGES - SELECTION SORT
		diffMutation = diffMutationSS;
		diameter = n-1;
		tmpint = new int[n*n+3*n]; //1 matrix cycles + 1 exchanges seq. + 1 perms. (consider also the crossover)
	} else if (strcmp(sgenerators,"ins2")==0) {
		diffMutation = diffMutationIS2;		//INSERTIONS 2 - INSERTION SORT FASTER BUT WITH LESS ENTROPY
		diameter = n-1;
		tmpint = new int[5*n]; //one insertions sequence (2*n) + three permutations (consider also the crossover)
		if (ffmax>1.0)		//to address default configuration
			ffmax = 1.0;
		if (finit>1.0 || ffmin>1.0 || ffmax>1.0) { //check f params for insertions!
			cerr << "FINIT, MINF, MAXF HAS TO BE <= 1.0 FOR INSERTIONS!!!" << endl;
			exit(EXIT_FAILURE);
		}
	} else {									//error
		cerr << "ERROR: the generators \"" << sgenerators << "\" is not one of: \"asw,ins,exc,ins2\"" << endl;
		exit(EXIT_FAILURE);
	}
	//initialize crossover function pointer
	if (strcmp(scrossover,"tpii")==0) 			//TPII (DI AGA)
		crossover = crossover_tpii;
	else if (strcmp(scrossover,"obxcr")==0) 	//OBXCR (ORIGINALE DI LOP)
		crossover = crossover_obxcr;
	else if (strcmp(scrossover,"obx")==0) 		//OBX (RANDOM)
		crossover = crossover_obx;
	else if (strcmp(scrossover,"tpiicr")==0)	//TPIICR (TPII CON PARAMETETRO CR)
		crossover = crossover_tpiicr;
	else {
		cerr << "ERROR: the crossover \"" << scrossover << "\" is not one of \"tpii,tpiicr,obxcr,obx\"" << endl;
		exit(EXIT_FAILURE);
	}
	//initialize selection function pointer
	if (strcmp(sselection,"alpha")==0)
		selection = selection_alpha;
	else if (strcmp(sselection,"crowding")==0)
		selection = selection_crowding;
#ifdef LOP	//currently these two are only supported by LOP!!!
	else if (strcmp(sselection,"inf")==0) {
		selection = selection_informationBased;
		initLoglookup(n);	//init also the logarithms lookup table for this selection
	}  else if (strcmp(sselection,"fitsep")==0) 
		selection = selection_fitnessSeparation;
#endif
	else {
		cerr << "ERROR: the selection \"" << sselection << "\" is not one of: \"alpha,crowding,inf\"" << endl;
		exit(EXIT_FAILURE);
	}
	//initialize localSearch function pointer
	if (strcmp(slsearch,"vns4")==0)
		localSearch = localSearch_vns4;
	else if (strcmp(slsearch,"ins")==0)
		localSearch = localSearch_ins;
	else {
		cerr << "ERROR: the localSearch \"" << slsearch << "\" is not one of: \"vns4,ins\"" << endl;
		exit(EXIT_FAILURE);
	}
	//initialize restart function pointers
	if (strcmp(srestart,"randls")==0) {
		popRestart = popRestart_randls;
		popForcedRestart = popForcedRestart_randls;
	} else if (strcmp(srestart,"shrandls")==0) {
		popRestart = popRestart_shrandls;
		popForcedRestart = popForcedRestart_shrandls;
	} else {
		cerr << "ERROR: the restart \"" << srestart << "\" is not one of: \"randls,shrandls\"" << endl;
		exit(EXIT_FAILURE);
	}
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
	delete[] crs;
	delete[] tmpint;
	delete[] lsHistory;
	freeLoglookup();
}



//evolve population
bool popEvolve() {
	//differential mutation, crossover, evaluation, update best, check termination
	int i;
	for (i=0; i<np; i++) {
		//generate one mutant
		diffMutation(i);
		//generate two (or one) offspring
		crossover(i);
		//evaluate first offspring
		fy1[i] = eval(y1[i]);
		nfes++;
		updateGbest(y1[i],fy1[i]);
		//check termination
		if (termination())
			return false;
		//if there are 2 childs are allowed (and not one)
		if (nchilds==2) {
			//evaluate first offspring
			fy2[i] = eval(y2[i]);
			nfes++;
			updateGbest(y2[i],fy2[i]);
			//check termination
			if (termination())
				return false;
		}
		//population loop done
	}
	//selection and compute sameFitness
	selection();
	//increment generation counter
	ngen++;
	//termination has not been triggered so return true
	return true;
	//done
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
		stageFoundAt = nrestarts;
		gbestls = lsmode;
		if (nfesWhenToForceRestart<INT_MAX)
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
	unsigned int now = maxTime>0 || maxStagnTime>0 ? (unsigned int)getTimer() : 0;
#if defined(LOPCC)
	return nfes>=maxnfes || (maxTime>0 && now>=maxTime) || (maxStagnTime>0 && now-timeFoundAt>=maxStagnTime) || fgbest<targetFit+.0005;	//.0005 (and <) because the optima are with 3 decimal precision
#elif defined(MINIMIZATION)
	return nfes>=maxnfes || (maxTime>0 && now>=maxTime) || (maxStagnTime>0 && now-timeFoundAt>=maxStagnTime) || fgbest<=fitbound;
#else
	return nfes>=maxnfes || (maxTime>0 && now>=maxTime) || (maxStagnTime>0 && now-timeFoundAt>=maxStagnTime) || fgbest>=fitbound;
#endif
}

