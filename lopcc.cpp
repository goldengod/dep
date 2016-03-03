#include "problem.h"
#include "utils.h"
#include <cstdio>    //FILE, fopen, fclose, fscanf
#include <cstdlib>   //EXIT_FAILURE, exit
#include <iostream>  //cerr, endl
#include <cmath>     //sqrt
#include <cstring>   //strcpy
using namespace std;


#define STRING_LENGTH 256


char instance[STRING_LENGTH];			//instance filename definition (declared extern in problem.h)
int n;									//permutation length definition (declared extern in problem.h)
double** c;								//cost matrix interface definition (declared extern in problem.h)
double* p;								//additional cost vector definition (declared extern in problem.h)
double fitbound;						//upper bound definition (declared extern in problem.h)
double heufit = MINUS_INF;				//not used here for lopcc (but declared extern in problem.h)
int* heup = 0;							//not used here for lopcc (but declared extern in problem.h)
int maxnfes;							//maxnfes definition (declared extern in problem.h)

static double* cs;						//inner cost matrix + additional vector storage definition (Iliffe/display style)
static double* alphaCost;				//inner n-length array for LOPCC fitness computation


void readInstance(char* filename) {
	//setup instance name
	strcpy(instance,filename);
	//open the file, read how many lines there are, set n, close the file
	FILE* f = fopen(filename,"r");
	if (!f) {
		cerr << "ERROR: Unable to open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	n = nlines(f);
	fclose(f);
	//open again the file and check for errors
	f = fopen(filename,"r");
	if (!f) {
		cerr << "ERROR: Unable to open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	//allocate memory and setup the pointers (Iliffe/display style)
	int tot = n*n+n;	//cost matrix (n*n) + additional vector (n)
	cs = new double[tot];
	c = new double*[n];
	for (int i=0; i<n; i++)
		c[i] = &cs[i*n];
	p = cs+(n*n);
	alphaCost = new double[n];
	//read the matrix entries and check for errors
	for (int i=0; i<tot; i++) {
		if (fscanf(f,"%lf",&cs[i])!=1) {
			cerr << "ERROR: Unable to read " << i << "th (0-based) entry from " << filename << endl;
			exit(EXIT_FAILURE);
		}
	}
	//close the file
	fclose(f);
	//done
}


void destroyInstance() {
	//deallocate memory
	delete[] c;
	delete[] cs;
	delete[] alphaCost;
	//done
}


double eval(int* x) {
	//See "the linear ordering problem with cumulative costs" [Bertacco, Brunetta, Fischetti]
	//fit = SUM_i=0^n-1 alpha_{x_i}
	//alpha_{x_i} = p_{x_i} + SUM_j=i+1^n-1 ( c_{x_i,x_j} * alpha_{x_j} ) for i = n-1, n-2, ..., 0
	//some variables
	int i,j;
	double f;
	//compute last alphaCost and initialize fitness
	f = alphaCost[n-1] = p[x[n-1]];							//simplification alphaCost[n-1] instead of alphaCost[x[n-1]]
	//outer + inner summation
	for (i=n-2; i>=0; i--) {
		alphaCost[i] = p[x[i]];								//simplification alphaCost[i] instead of alphaCost[x[i]]
		for (j=i+1; j<n; j++)
			alphaCost[i] += c[x[i]][x[j]] * alphaCost[j];	//simplification alphaCost[i] ([j]) instead of alphaCost[x[i]] (x[j])
		f += alphaCost[i];									//simplification alphaCost[i] instead of alphaCost[x[i]]
	}
	//return fitness
	return f;
	//done
}


void readHeuristic(char* filename) {
	//NOT USED FOR LOPCC!!!
	cerr << "ERROR: WE DONT HAVE HEURISTIC IMPLEMENTED FOR LOPCC!!!" << endl;
	exit(EXIT_FAILURE);
}

