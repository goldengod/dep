#ifndef PROBLEM_H
#define PROBLEM_H

#define STRING_LENGTH 256

//instance filename declaration
extern char instance[STRING_LENGTH];

//floshop n,m declaration
extern int n;
extern int m;

//ptime matrix (n x m) declaration
extern int** ptime;

//liu reeves fitness and permutation declaration
extern int heufit;
extern int* heup;

//maxnfes declaration
extern int maxnfes;

//best tft declaration
extern int best;

//read the instance from file
void readInstance(char* filename);

//destroy the instance memory
void destroyInstance();

//evaluate the permutation x (full evaluation)
int eval(int* x);

#ifdef GFC
//GFC functions
void computeGFC(int* x);
int evalGFC(int* x, int k);
int evalUpdateGFC(int* x, int k);
#endif

//read heuristic file (to fill heufit and heup)
void readHeuristic(char* filename);

#endif

