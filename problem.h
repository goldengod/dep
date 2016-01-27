#ifndef PROBLEM_H
#define PROBLEM_H


//BEGIN COMMON PART
#define STRING_LENGTH 256
extern char instance[STRING_LENGTH];	//instance filename declaration
extern int maxnfes;						//maxnfes declaration
extern int heufit;						//heuristic fitness for the problem at hand (liu-reeves/neh/...)
extern int* heup;						//heuristic permutation for the problem at hand (liu-reeves/neh/...)
extern int fitbound;					//fitness bound (lower bound for minimization, upper bound for maximization) used for termination
void readInstance(char* filename);		//read the instance from file
void destroyInstance();					//destroy the instance memory
int eval(int* x);						//evaluate the permutation x (full evaluation)
void readHeuristic(char* filename);		//read heuristic file (to fill heufit and heup)
//END COMMON PART


//BEGIN PERMUTATION FLOWSHOP PART
#if defined(TFT) || defined(MAKESPAN)
#define MINIMIZATION
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
extern int n;							//instance size
extern int** h;							//i/o matrix (n x n)
#endif
//END LOP PART


#endif

