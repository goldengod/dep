#ifndef DEP_H
#define DEP_H

#include "problem.h"			//to include the typedef for FitnessType

//constants
#define N_LS 0					//No local search
#define B_LS 1					//Baldwinian local search
#define L_LS 2					//Lamarckian local search

//input
extern int np;					//population size declaration
extern double finit;			//initial scale factor
extern double crinit;			//initial crossover probability (needed for obxcr)
extern double alpha;			//alpha value (needed for alpha selection)
extern int heu;					//0 or 1 for use or not heuristic in randheu initialization
extern int ls;					//0 (no ls), 1 (baldwinian ls), 2 (lamarckian ls) for using local search in some restarts
extern double frfactor;			//forced restart factor (0,1+) - restart forced if no gbest update for maxnfes*frfactor nfes
extern double ffmin;			//min value for scale factor f (used by jde rule) (fmin conflicts with cmath)
extern double ffmax;			//max value for scale factor f (used by jde rule) (fmax conflicts with cmath)
extern double crmin;			//min value for crossover probability cr (used by jde rule) (needed for obxcr)
extern double crmax;			//max value for crossover probability cr (used by jde rule) (needed for obxcr)
extern char sgenerators[10];	//generating set for differential mutation
extern char sinitialization[10];//initialization algorithm
extern char scrossover[10];		//crossover algorithm
extern char sselection[10];		//selection algorithm
extern char slsearch[10];		//local search algorithm
extern char srestart[10];		//restart algorithm
extern int inftype;				//parameter for inf selection (1 or 2 or 3 or 4)
extern int nchilds;				//number of childs for the crossover (1 or 2)
extern bool memetic;			//memetic flag

//output
extern int* gbest;				//global best
extern FitnessType fgbest;		//fitness of the global best
extern int nfesFoundAt;			//evaluatins when the global best has been found
extern unsigned long timeFoundAt;//milliseconds when the global best has been found
extern int stageFoundAt;		//stage when the global best has been found
extern int nfes;				//number of evaluations performed
extern int ngen;				//number of generations performed
extern unsigned long execTime;	//time in millis needed
extern int nrestarts;			//number of restarts
extern int nforcedrestarts;		//number of forced restarts
extern int minStageLength;		//restart statistic
extern int maxStageLength;		//restart statistic
extern double avgStageLength;	//restart statistic
extern int improvingStages;		//restart statistic
extern int nfesls;				//local search statistic
extern int nls;					//local search statistic
extern int nImprovingls;		//local search statistic
extern FitnessType totImprovingls;//local search statistic
extern bool gbestls;			//local search statistic
extern int improvingSteps;		//total number of improvements
extern int lsImprovingSteps;	//number of improvements due to local searches
extern double sfSuccAvg;		//f statistic
extern double sfSuccVar;		//f statistic
extern double sfSuccMin;		//f statistic
extern double sfSuccMax;		//f statistic
extern double crSuccAvg;		//cr statistic (needed for obxcr)
extern double crSuccVar;		//cr statistic (needed for obxcr)
extern double crSuccMin;		//cr statistic (needed for obxcr)
extern double crSuccMax;		//cr statistic (needed for obxcr)
extern int child1succ;			//two childs statistic
extern int child2succ;			//two childs statistic

//others
extern bool haveToSave;			//flag for saving execution state/memory declaration

//functions
void dep();						//main dep function
void depDefaultParameters();	//set default dep parameters
void depAlloc();				//allocate dep memory
void depFree();					//free dep memory
void depResume(char* filename);	//resume a saved execution of dep

#endif

