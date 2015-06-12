#ifndef DEP_H
#define DEP_H

#define N_LS 0 //No local search
#define B_LS 1 //Baldwinian local search
#define L_LS 2 //Lamarckian local search

//population size declaration
extern int np;

//initial scale factor declaration
extern double finit;

//alpha value declaration
extern double alpha;

//heu and ls value declaration
extern int heu; //0 or 1
extern int ls; //0 (no ls), 1 (baldwinian ls), 2 (lamarckian ls)

//forced restart factor (0,1+) - restart forced if no gbest update for maxnfes*frfactor nfes
extern double frfactor;

//global best declaration
extern int* gbest;
extern int fgbest;
extern int nfesFoundAt;
extern int stageFoundAt;

//number of evaluations and generation declaration
extern int nfes;
extern int ngen;

//execution time declaration
extern unsigned long execTime;

//number of restarts performed declaration
extern int nrestarts;
extern int nforcedrestarts;

//restarts statistics declaration
extern int minStageLength;
extern int maxStageLength;
extern double avgStageLength;
extern int improvingStages;

//local search statistics declaration
extern int nfesls;
extern int nls;
extern int nImprovingls;
extern int totImprovingls;
extern bool gbestls;

//flag for saving execution state/memory declaration
extern bool haveToSave;

//number of improving steps performed declaration
extern int improvingSteps;
extern int lsImprovingSteps;

//main dep function
void dep();

//set default dep parameters
void depDefaultParameters();

//allocate/free dep memory
void depAlloc();
void depFree();

//dep function to resume a saved execution
void depResume(char* filename);

#endif

