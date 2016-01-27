#ifndef TFT_H
#define TFT_H

#define STRING_LENGTH 256

//instance filename declaration
extern char instance[STRING_LENGTH];

//floshop n,m declaration
extern int n;
extern int m;

//ptime matrix (n x m) declaration
extern int** ptime;

//liu reeves fitness and permutation declaration
extern int lrtft;
extern int* lrp;

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

#endif

