#include "global.h"
#include "random.h"
#include "problem.h"
#include "dep.h"
#include "utils.h"
#include <cstdlib>   //exit, atoi, EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>  //cout, endl, cerr
#include <cstring>   //strcpy, strcmp, strlen, strcat, strcmp
#include <cstdio>    //sscanf, FILE, fopen, fclose, fprintf, sprintf, fscanf
#include <csignal>   //signal, SIGALRM
#include <unistd.h>  //alarm
#include <climits>   //INT_MAX
#include <cmath>	 //sqrt
using namespace std;


//these variables are declared extern in global.h (since have to be visible to dep.cpp)
unsigned int seed = 0;
char out[STRING_LENGTH];
char exe[STRING_LENGTH];
unsigned int saveSeconds = 0;
char saveFile[STRING_LENGTH];
char scriptFile[STRING_LENGTH];	//for grid version
unsigned int maxTime = 0;
unsigned int maxStagnTime = 0;
FitnessType targetFit;
int maxStagnNfes = 0;

//inner variables
static bool resume;


void usage();
void readArguments(int argc, char** argv);
void writeResults();
void defaultSaveFile();
void saveHandler(int signum);
void preResume(char* filename);
void defaultMaxnfes();


int main(int argc, char** argv) {
	//set default parameters for dep
	depDefaultParameters();
	//read arguments, set parameters and load problem instance
	readArguments(argc,argv);
	//check if maxTime was set, in the case set maxnfes to +inf
	if (maxTime>0)
		maxnfes = INT_MAX;
	//check if it's normal or resume execution
	if (!resume) { //normal execution
		//install the alarm for saving in the case save was set
		if (saveSeconds>0) {
			signal(SIGALRM,saveHandler);
			alarm(saveSeconds);
		}
		//random seed if zero
		if (!seed)
			seed = randSeed(); //from /dev/random
		//init random number generator
		initRand(seed);
		//set default save file if the case
		if (saveSeconds>0 && saveFile[0]=='\0')
			defaultSaveFile();
		//done
	} else { //resume execution
		//call preResume to read inputs and instance
		preResume(saveFile);
		//install the alarm for saving in the case save was set
		signal(SIGALRM,saveHandler);
		alarm(saveSeconds);
		//done
	}
	//allocate dep memory
	depAlloc();
	//some printing before dep execution
	cout << "STARTING DEP:" << endl;
	cout << "               exe = " << exe << endl;
	cout << " -        instance = " << instance << endl;
	cout << " -               n = " << n << endl;
#if defined(TFT) || defined(MAKESPAN)
	cout << " -               m = " << m << endl;
#endif
	cout << " -        fitbound = " << fitbound << endl;
	cout << " -    problem type = " <<
#ifdef MINIMIZATION
									   "minimization" << endl;
#else
									   "maximization" << endl;
#endif
	cout << " -             gen = " << sgenerators << endl;
	cout << " -            init = " << sinitialization << endl;
	cout << " -           cross = " << scrossover << endl;
	cout << " -             sel = " << sselection << endl;
	cout << " -         lsearch = " << slsearch << endl;
	cout << " -         restart = " << srestart << endl;
	cout << " -             out = " << out << endl;
	cout << " -         maxnfes = " << maxnfes << endl;
	cout << " -    maxStagnNfes = " << maxStagnNfes << endl;
	cout << " -         maxTime = " << maxTime << endl;
	cout << " -    maxStagnTime = " << maxStagnTime << endl;
	cout << " -       targetFit = " << targetFit << endl;
	cout << " -            seed = " << seed << endl;
	cout << " -              np = " << np << endl;
	cout << " -           finit = " << finit << endl;
	cout << " -            fmin = " << ffmin << endl;
	cout << " -            fmax = " << ffmax << endl;
	cout << " -          crinit = " << crinit << endl;
	cout << " -           crmin = " << crmin << endl;
	cout << " -           crmax = " << crmax << endl;
	cout << " -         nchilds = " << nchilds << endl;
	cout << " -           alpha = " << alpha << endl;
	cout << " -         inftype = " << inftype << endl;
	cout << " -             heu = " << heu << endl;
	cout << " -         memetic = " << memetic << endl;
	cout << " -              ls = " << ls << " ";
	switch (ls) {
		case N_LS:
			cout << "(N_LS)" << endl;
			break;
		case B_LS:
			cout << "(B_LS)" << endl;
			break;
		case L_LS:
			cout << "(L_LS)" << endl;
			break;
		default:
			cerr << "ERROR: ls has to be 0,1,2!!!" << endl;
			exit(EXIT_FAILURE);
	}
	cout << " -        frfactor = " << frfactor << endl;
	cout << " -     saveSeconds = " << saveSeconds;
	if (saveSeconds==0)
		cout << " (no save)" << endl;
	else {
		cout << endl << " -        saveFile = " << saveFile << endl;
		cout << " -      scriptFile = " << scriptFile << endl;
	}
	cout << "... DEP IS RUNNING ..." << endl;
	//execute dep basing on normal or resume execution mode
	if (!resume)
		dep(); //normal execution
	else
		depResume(saveFile); //resume execution
	//some printing after dep execution
	cout << "DEP DONE:" << endl;
	cout << " -          fgbest = " << fgbest << endl;
	cout << " -     nfesFoundAt = " << nfesFoundAt << endl;
	char stime[256];
	millis2str(timeFoundAt,stime);
	cout << " -     timeFoundAt = " << stime << endl;
	cout << " -    stageFoundAt = " << stageFoundAt << endl;
	cout << " -            nfes = " << nfes << endl;
	cout << " -            ngen = " << ngen << endl;
	cout << " -       nrestarts = " << nrestarts << endl;
	cout << " - nforcedrestarts = " << nforcedrestarts << endl;
	cout << " -  minStageLength = " << minStageLength << endl;
	cout << " -  maxStageLength = " << maxStageLength << endl;
	cout << " -  avgStageLength = " << avgStageLength << endl;
	cout << " - improvingStages = " << improvingStages << endl;
	cout << " -             nls = " << nls << endl;
	cout << " -          nfesls = " << nfesls << endl;
	cout << " -    nImprovingls = " << nImprovingls << endl;
	cout << " -  totImprovingls = " << totImprovingls << endl;
	cout << " -         gbestls = " << gbestls << endl;
	cout << " -  improvingSteps = " << improvingSteps << endl;
	cout << " -lsImprovingSteps = " << lsImprovingSteps << endl;
	cout << " -       sfSuccAvg = " << sfSuccAvg << endl;
	cout << " -       sfSuccStd = " << sqrt(sfSuccVar) << endl;
	cout << " -       sfSuccMin = " << sfSuccMin << endl;
	cout << " -       sfSuccMax = " << sfSuccMax << endl;
	cout << " -       crSuccAvg = " << crSuccAvg << endl;
	cout << " -       crSuccStd = " << sqrt(crSuccVar) << endl;
	cout << " -       crSuccMin = " << crSuccMin << endl;
	cout << " -       crSuccMax = " << crSuccMax << endl;
	cout << " -      child1succ = " << child1succ << endl;
	cout << " -      child2succ = " << child2succ << endl;
	millis2str(execTime,stime);
	cout << " -        execTime = " << stime << endl;
	//check gbest validity and fitness
	bool gbestValid = permValid(gbest,n);
	bool gbestFitness = eval(gbest)==fgbest;
	//write results in output file
	writeResults();
	//some final printing
	cout << "Full results are available in " << out << endl;
	//free dep memory
	depFree();
	//destroy instance memory
	destroyInstance();
	//check gbestValid/gbestFitness and print error
	if (!gbestValid)
		cerr << "ERROR: GBEST IS NOT A VALID PERMUTATION!!!" << endl;
	if (!gbestFitness) {
		cerr << "ERROR: GBEST FITNESS DOES NOT MATCH!!!" << endl;
		cerr << "PERM = "; printPerm(gbest,n,' ',true);
		cerr << "FIT = " << fgbest << endl;
	}
	//return success/failure basing on gbestValid/gbestFitness
	return gbestValid && gbestFitness ? EXIT_SUCCESS : EXIT_FAILURE;
	//done
}


void usage() {
	cout << "-------------" << endl;
	cout << "DEP USAGE:" << endl;
	cout << "(1) dep INSTANCE_FILE OUTPUT_FILE" << endl;
	cout << "    [--gen STR (asw,ins,exc,ins2) NOTE: ins2 are faster insertions with less entropy]" << endl;
	cout << "    [--init STR (randheu)]" << endl;
	cout << "    [--cross STR (tpii,tpiicr,obxcr,obx)]" << endl;
	cout << "    [--nchilds INT (1,2)]" << endl;
	cout << "    [--sel STR (alpha,crowding,inf,fitsep)]" << endl;
	cout << "    [--lsearch STR (vns4,ins)]" << endl;
	cout << "    [--restart STR (randls,shrandls)]" << endl;
	cout << "    [--seed UINT]" << endl;
	cout << "    [--maxnfes INT]" << endl;
	cout << "    [--maxstagnnfes INT]" << endl;
	cout << "    [--np INT]" << endl;
	cout << "    [--finit DOUBLE]" << endl;
	cout << "    [--fmin DOUBLE]" << endl; //minimum value for f
	cout << "    [--fmax DOUBLE]" << endl; //maximum value for f
	cout << "    [--crinit DOUBLE]" << endl;
	cout << "    [--crmin DOUBLE]" << endl;
	cout << "    [--crmax DOUBLE]" << endl;
	cout << "    [--f DOUBLE] (if you set this value, please don't set --fmin, --fmax, and --finit)" << endl;
	cout << "    [--cr DOUBLE] (if you set this value, please don't set --crmin, --crmax, and --crinit)" << endl;
	cout << "    [--alpha DOUBLE]" << endl;
	cout << "    [--inftype INT (1,2,3,4)]" << endl;
	cout << "    [--heu STRING(filename)]" << endl;
	cout << "    [--ls 0(N_LS)/1(B_LS)/2(L_LS)]" << endl;
	cout << "    [--frfactor DOUBLE]" << endl;
	cout << "    [--memetic]" << endl;
	cout << "    [--save UINT(seconds) and optionally STRING(filename) and optionally STRING(scriptfile)]" << endl;
	cout << "    [--maxtime UINT(milliseconds)]" << endl;
	cout << "    [--maxstagntime UINT(milliseconds)]" << endl;
	cout << "    [--targetfit" <<
#ifdef FIT_INT
									"INT]" << endl;
#else
									"DOUBLE]" << endl;
#endif
	cout << "(2) dep resume SAVING_FILE" << endl;
	cout << "-------------" << endl;
}


void readArguments(int argc, char** argv) {
	//check if enough arguments
	if (argc<3) {
		usage();
		exit(EXIT_FAILURE);
	}
	//init targetFit to +inf or -inf
#if defined(MINIMIZATION) && defined(FIT_REAL)
	targetFit = MINUS_INF;
#elif defined(MINIMIZATION) && defined(FIT_INT)
	targetFit = INT_MIN;
#elif defined(MAXIMIZATION) && defined(FIT_REAL)
	targetFit = PLUS_INF;
#elif defined(MAXIMIZATION) && defined(FIT_INT)
	targetFit = INT_MAX;
#endif
	//check if it's a normal or resume execution
	if (strcmp(argv[1],"resume")==0) {
		//if resume only set the flag to true and read saveFile
		resume = true;
		strcpy(saveFile,argv[2]);
	} else {
		//if normal, set flag to false, and do the rest
		resume = false;
		//read instance and save out filename and exe filename
		strcpy(exe,argv[0]);
		readInstance(argv[1]);
		strcpy(out,argv[2]);
		//init heu to 0 and memetic to false
		heu = 0;
		memetic = false;
		//init default maxnfes basing on instance size
		defaultMaxnfes();
		//read and set the optional arguments
		int i = 3;
		while (i<argc) {
			if (strcmp(argv[i],"--seed")==0) {
				sscanf(argv[i+1],"%u",&seed);
				i += 2;
			} else if (strcmp(argv[i],"--maxnfes")==0) {
				maxnfes = atoi(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--maxstagnnfes")==0) {
				maxStagnNfes = atoi(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--np")==0) {
				np = atoi(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--finit")==0) {
				finit = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--crinit")==0) {
				crinit = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--alpha")==0) {
				alpha = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--heu")==0) {
				readHeuristic(argv[i+1]);
				heu = 1; //otherwise it is 0
				i += 2;
			} else if (strcmp(argv[i],"--ls")==0) {
				ls = atoi(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--frfactor")==0) {
				frfactor = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--fmin")==0) {
				ffmin = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--fmax")==0) {
				ffmax = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--crmin")==0) {
				crmin = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--crmax")==0) {
				crmax = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--f")==0) {
				//subtle way to set a fixed F parameter
				finit = ffmin = ffmax = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--cr")==0) {
				//subtle way to set a fixed CR parameter
				crinit = crmin = crmax = atof(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--save")==0) {
				sscanf(argv[i+1],"%u",&saveSeconds);
				if (i+2<argc && argv[i+2][0]!='-') {
					strcpy(saveFile,argv[i+2]);
					if (i+3<argc && argv[i+3][0]!='-') {
						strcpy(scriptFile,argv[i+3]);
						i += 4;
					} else {
						strcpy(scriptFile,"none");
						i += 3;
					}
				} else {
					//default fileSave is set later if the case (have to known the seed)
					saveFile[0] = '\0';
					i += 2;
				}
			} else if (strcmp(argv[i],"--maxtime")==0) {
				sscanf(argv[i+1],"%u",&maxTime);
				i += 2;
			} else if (strcmp(argv[i],"--maxstagntime")==0) {
				sscanf(argv[i+1],"%u",&maxStagnTime);
				i += 2;
			} else if (strcmp(argv[i],"--gen")==0) {
				strcpy(sgenerators,argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--init")==0) {
				strcpy(sinitialization,argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--cross")==0) {
				strcpy(scrossover,argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--sel")==0) {
				strcpy(sselection,argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--lsearch")==0) {
				strcpy(slsearch,argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--restart")==0) {
				strcpy(srestart,argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--inftype")==0) {
				inftype = atoi(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--nchilds")==0) {
				nchilds = atoi(argv[i+1]);
				i += 2;
			} else if (strcmp(argv[i],"--targetfit")==0) {
#ifdef FIT_INT
				targetFit = atoi(argv[i+1]);
#else
				targetFit = atof(argv[i+1]);
#endif
				i += 2;
			} else if (strcmp(argv[i],"--memetic")==0) {
				memetic = true;
				i++;
			} else {
				cerr << "COMMAND LINE PARAMETERS WRONG!" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	//done
}


void writeResults() {
	//check if out file exists and if not write the csv header
	FILE* f = fopen(out,"r");
	if (!f) {
		f = fopen(out,"w");
#if defined(TFT) || defined(MAKESPAN)
		fprintf(f,"exe,instance,n,m,gen,init,cross,nchilds,sel,lsearch,restart,memetic,maxnfes,maxstagnnfes,maxtime,maxstagntime,targetfit,seed,np,finit,fmin,fmax,crinit,crmin,crmax,alpha,heu,ls,frfactor,inftype"); //input
#endif
#if defined(LOP) || defined(LOPCC)
		fprintf(f,"exe,instance,n,gen,init,cross,nchilds,sel,lsearch,restart,memetic,maxnfes,maxstagnnfes,maxtime,maxstagntime,targetfit,seed,np,finit,fmin,fmax,crinit,crmin,crmax,alpha,heu,ls,frfactor,inftype"); //input
#endif
		fprintf(f,",fgbest,nfesFoundAt,timeFoundAt,stageFoundAt,nfes,ngen,nrestarts,nforcedrestarts,improvingSteps,lsImprovingSteps,execTime,minStageLength,maxStageLength,avgStageLength,improvingStages,nls,nfesls,nImprovingls,totImprovingls,gbestls,sfSuccAvg,sfSuccStd,sfSuccMin,sfSuccMax,crSuccAvg,crSuccStd,crSuccMin,crSuccMax,child1succ,child2succ,gbest\n"); //output
	}
	fclose(f);
	//write this csv line in append
	char sgbest[PERMSTR_SIZE];
	char sfinit[32];
	char sfmin[32];
	char sfmax[32];
	char scrinit[32];
	char scrmin[32];
	char scrmax[32];
	char salpha[32];
	char sfrfactor[32];
	char savgStageLength[32];
	char ssfSuccAvg[32];
	char ssfSuccStd[32];
	char ssfSuccMax[32];
	char ssfSuccMin[32];
	char scrSuccAvg[32];
	char scrSuccStd[32];
	char scrSuccMax[32];
	char scrSuccMin[32];
	perm2str(gbest,n,sgbest);
	double2str(finit,sfinit);
	double2str(ffmin,sfmin);
	double2str(ffmax,sfmax);
	double2str(crinit,scrinit);
	double2str(crmin,scrmin);
	double2str(crmax,scrmax);
	double2str(alpha,salpha);
	double2str(frfactor,sfrfactor);
	double2str(avgStageLength,savgStageLength);
	double2str(sfSuccAvg,ssfSuccAvg);
	double sfSuccStd = sqrt(sfSuccVar);
	double2str(sfSuccStd,ssfSuccStd);
	double2str(sfSuccMax,ssfSuccMax);
	double2str(sfSuccMin,ssfSuccMin);
	double2str(crSuccAvg,scrSuccAvg);
	double crSuccStd = sqrt(crSuccVar);
	double2str(crSuccStd,scrSuccStd);
	double2str(crSuccMax,scrSuccMax);
	double2str(crSuccMin,scrSuccMin);
#ifdef FIT_REAL
	char sfgbest[64];
	char stotImprovingls[64];
	double2str(fgbest,sfgbest);
	double2str(totImprovingls,stotImprovingls);
#endif
	f = fopen(out,"a");
	//input print
#if defined(TFT) || defined(MAKESPAN)
	fprintf(f,"%s,%s,%d,%d,%s,%s,%s,%d,%s,%s,%s,%d,%d,%d,%d,%d,%d,%u,%d,%s,%s,%s,%s,%s,%s,%s,%d,%d,%s,%d",
				exe,instance,n,m,sgenerators,sinitialization,scrossover,nchilds,sselection,slsearch,srestart,memetic?1:0,maxnfes,maxStagnNfes,maxTime,maxStagnTime,targetFit,seed,np,sfinit,sfmin,sfmax,scrinit,scrmin,scrmax,salpha,heu,ls,sfrfactor,inftype);
#endif
#if defined(LOP)
	fprintf(f,"%s,%s,%d,%s,%s,%s,%d,%s,%s,%s,%d,%d,%d,%d,%d,%d,%u,%d,%s,%s,%s,%s,%s,%s,%s,%d,%d,%s,%d",
				exe,instance,n,sgenerators,sinitialization,scrossover,nchilds,sselection,slsearch,srestart,memetic?1:0,maxnfes,maxStagnNfes,maxTime,maxStagnTime,targetFit,seed,np,sfinit,sfmin,sfmax,scrinit,scrmin,scrmax,salpha,heu,ls,sfrfactor,inftype);
#endif
#if defined(LOPCC)
	fprintf(f,"%s,%s,%d,%s,%s,%s,%d,%s,%s,%s,%d,%d,%d,%d,%d,%.10lf,%u,%d,%s,%s,%s,%s,%s,%s,%s,%d,%d,%s,%d",
				exe,instance,n,sgenerators,sinitialization,scrossover,nchilds,sselection,slsearch,srestart,memetic?1:0,maxnfes,maxStagnNfes,maxTime,maxStagnTime,targetFit,seed,np,sfinit,sfmin,sfmax,scrinit,scrmin,scrmax,salpha,heu,ls,sfrfactor,inftype);
#endif
	//output print
#ifdef FIT_INT
	fprintf(f,",%d,%d,%lu,%d,%d,%d,%d,%d,%d,%d,%lu,%d,%d,%s,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%s\n",
			fgbest,nfesFoundAt,timeFoundAt,stageFoundAt,nfes,ngen,nrestarts,nforcedrestarts,
			improvingSteps,lsImprovingSteps,execTime,
			minStageLength,maxStageLength,savgStageLength,improvingStages,
			nls,nfesls,nImprovingls,totImprovingls,gbestls?1:0,ssfSuccAvg,ssfSuccStd,ssfSuccMin,ssfSuccMax,scrSuccAvg,scrSuccStd,scrSuccMin,scrSuccMax,child1succ,child2succ,sgbest);
#else
	fprintf(f,",%s,%d,%lu,%d,%d,%d,%d,%d,%d,%d,%lu,%d,%d,%s,%d,%d,%d,%d,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%s\n",
			sfgbest,nfesFoundAt,timeFoundAt,stageFoundAt,nfes,ngen,nrestarts,nforcedrestarts,
			improvingSteps,lsImprovingSteps,execTime,
			minStageLength,maxStageLength,savgStageLength,improvingStages,
			nls,nfesls,nImprovingls,stotImprovingls,gbestls?1:0,ssfSuccAvg,ssfSuccStd,ssfSuccMin,ssfSuccMax,scrSuccAvg,scrSuccStd,scrSuccMin,scrSuccMax,child1succ,child2succ,sgbest);
#endif
	fclose(f);
	//done
}


void defaultSaveFile() {
	//concatenate instance and seed in saveFile
	sprintf(saveFile,"%s_%u",instance,seed);
	//replace '/' and '.' with '_'
	for (unsigned int i=0; i<strlen(saveFile); i++)
		if (saveFile[i]=='/' || saveFile[i]=='.')
			saveFile[i] = '_';
	//cat .sav to saveFile
	strcat(saveFile,".sav");
	//done
}


void saveHandler(int signum) {
	//signum should be always SIGALRM, thus it is ignored
	//set the dep flag for saving (it is in dep.h/cpp)
	haveToSave = true;
	//set again the alarm (no need to install again the signal handler)
	alarm(saveSeconds);
	//done
}


void preResume(char* filename) {
	//pre-resume to load input variables and instance
	//some variables
	int nowarning; //nowarning is to avoid compiler complaints
	unsigned int j,b;
	unsigned char* bytes; //used to cast double and assume that sizeof(unsigned char)==1
	//open filename and check for errors
	FILE* fsav = fopen(filename,"r");
	if (!fsav) { //if error print and exit
		cerr << "WARNING: I was unable to open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	//read from the file
	//read saveSeconds and saveFile read first
	nowarning = fscanf(fsav,"%u",&saveSeconds);           //saveSeconds 1
	nowarning = fscanf(fsav,"%s",saveFile);               //saveFile 2
	nowarning = fscanf(fsav,"%s",scriptFile);             //scriptFile 2bis
	//read the rng
	readRandState(fsav);                                  //rng state 3
	//values of tft.h/cpp: read "instance", call readInstance, and read maxnfes
	nowarning = fscanf(fsav,"%s",instance);               //instance 4
	readInstance(instance);
	nowarning = fscanf(fsav,"%d",&maxnfes);               //maxnfes 5
	//other values of global.h (except saveSeconds and saveFile read at the begin
	nowarning = fscanf(fsav,"%u",&seed);                  //seed 6
	nowarning = fscanf(fsav,"%s",out);                    //out 7
	nowarning = fscanf(fsav,"%s",exe);                    //exe 8
	//values of dep.h (extern ones and only input) without call depAlloc()
	nowarning = fscanf(fsav,"%d",&np);                    //np 9
	bytes = (unsigned char*)&finit;                       //finit 10
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&ffmin;                       //fmin 10bis
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&ffmax;                       //fmax 10tris
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%s",sgenerators);			  //sgenerators 10quadris
	nowarning = fscanf(fsav,"%s",sinitialization);		  //sinitialization 10penta
	nowarning = fscanf(fsav,"%s",sselection);			  //sselection 10esa
	nowarning = fscanf(fsav,"%s",scrossover);			  //scrossover 10epta
	nowarning = fscanf(fsav,"%s",slsearch);				  //slsearch 10octa
	nowarning = fscanf(fsav,"%s",srestart);				  //srestart 10nine
	bytes = (unsigned char*)&crinit;                      //crinit 10ten
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crmin;                       //crmin 10eleven
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crmax;                       //crmax 10twelve
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&alpha;                       //alpha 11
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%d",&heu);                   //heu 12
	nowarning = fscanf(fsav,"%d",&ls);                    //ls 13
	bytes = (unsigned char*)&frfactor;                    //frfactor 14
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%d",&inftype);               //inftype 14bis
	nowarning = fscanf(fsav,"%d",&nchilds);				  //nchilds 14tris
	nowarning = fscanf(fsav,"%u",&maxStagnTime);		  //maxStagnTime 14quadris
#ifdef FIT_INT
	nowarning = fscanf(fsav,"%d",&targetFit);			  //targetFit 14penta
#else
	bytes = (unsigned char*)&targetFit;
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
#endif
	nowarning = fscanf(fsav,"%d",&maxStagnNfes);		  //maxStagnNfes 14esa
	nowarning = fscanf(fsav,"%d",&b);					  //memetic 14epta
	memetic = b==1;
	//close the file
	fclose(fsav);
	//statement to avoid compiler complaints
	nowarning++;
	//done
}


void defaultMaxnfes() {
	maxnfes = 100000; //default "stupid" value
#if defined(TFT) || defined(MAKESPAN)
	//for taillard instances
	if (n==20) {
		if (m==5)
			maxnfes = 182224100;
		else if (m==10)
			maxnfes = 224784800;
		else if (m==20)
			maxnfes = 256896400;
	} else if (n==50) {
		if (m==5)
			maxnfes = 220712150;
		else if (m==10)
			maxnfes = 256208100;
		else if (m==20)
			maxnfes = 275954150;
	} else if (n==100) {
		if (m==5)
			maxnfes = 235879800;
		else if (m==10)
			maxnfes = 266211000;
		else if (m==20)
			maxnfes = 283040000;
	} else if (n==200) {
		if (m==10)
			maxnfes = 272515500;
		else if (m==20)
			maxnfes = 287728850;
	} else if (n==500) {
		if (m==20)
			maxnfes = 260316750;
	}
#endif
#if defined(LOP)
	//Ceberio et al. use: n^2 * (1000, 5000, 10000)
	#define MAXNFES_MULTIPLIER 10000 //1000, 5000, 10000 (as in Ceberio et al)
	maxnfes = n*n*MAXNFES_MULTIPLIER;
#endif
#if defined(LOPCC)
	//do not use maxnfes but maxStagnNfes and maxtime
	maxnfes = INT_MAX;
	//maxStagnNfes = n*n*200;	//derived from preliminar experiments	//old!!!
	switch (n) {
		case 100:	maxStagnNfes = 10000000; maxTime = 3600000; break;	//10M nfes of stagnation + 1 hour total
		case 150:	maxStagnNfes = 10000000; maxTime = 3600000; break;	//10M nfes of stagnation + 1 hour total
		default:	break;												//no limit
	}
#endif
}

