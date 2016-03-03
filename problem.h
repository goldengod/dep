#ifndef PROBLEM_H
#define PROBLEM_H


//BEGIN PERMUTATION FLOWSHOP PART
#if defined(TFT) || defined(MAKESPAN)
#define MINIMIZATION
#define FIT_INT
typedef int FitnessType;
extern int n;							//instance size
extern int m;							//instance size
extern int** ptime;						//ptime matrix (n x m)
//extern int best;						//best tft/makespan
#ifdef GFC
void computeGFC(int* x);				//compute GFC data structure
int evalGFC(int* x, int k);				//partial evaluation using GFC
int evalUpdateGFC(int* x, int k);		//partial evaluation using GFC + update of the GFC data structure
#endif
#endif
//ENDPERMUTATION FLOWSHOP PART

//BEGIN LOP PART
#ifdef LOP
#define MAXIMIZATION
#define FIT_INT
typedef int FitnessType;
extern int n;							//instance size
extern int** h;							//i/o matrix (n x n)
#endif
//END LOP PART

//BEGIN LOPCC PART
#ifdef LOPCC
#define MINIMIZATION
#define FIT_REAL
typedef double FitnessType;
extern int n;							//instance size
extern double** c;						//cost matrix
extern double* p;						//additional cost vector
#endif
//END LOPCC PART

//BEGIN COMMON PART
#define STRING_LENGTH 256
extern char instance[STRING_LENGTH];	//instance filename declaration
extern int maxnfes;						//maxnfes declaration
extern int* heup;						//heuristic permutation for the problem at hand (liu-reeves/neh/...)
extern FitnessType heufit;				//heuristic fitness for the problem at hand (liu-reeves/neh/...)
extern FitnessType fitbound;			//fitness bound (lower bound for minimization, upper bound for maximization) used for termination
void readInstance(char* filename);		//read the instance from file
void destroyInstance();					//destroy the instance memory
FitnessType eval(int* x);				//evaluate the permutation x (full evaluation)
void readHeuristic(char* filename);		//read heuristic file (to fill heufit and heup)
//END COMMON PART


#endif

