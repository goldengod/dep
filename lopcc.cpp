#include "problem.h"
#include <cstdio>    //FILE, fopen, fclose, fscanf
#include <cstdlib>   //EXIT_FAILURE, exit
#include <iostream>  //cerr, endl
#include <cmath>     //sqrt
#include <cstring>   //strcpy
#include <climits>   //INT_MAX
using namespace std;


#define STRING_LENGTH 256


//instance filename definition (declared extern in problem.h)
char instance[STRING_LENGTH];

//permutation length definition (declared extern in problem.h)
int n;

//i/o matrix interface definition (declared extern in problem.h)
int** h;

//upper bound definition (declared extern in problem.h)
int fitbound;

//not used here for lop (but declared extern in problem.h)
int heufit = INT_MAX;
int* heup = 0;

//maxnfes definition (declared extern in problem.h)
int maxnfes;

//inner i/o matrix storage definition (Iliffe/display style)
static int* hs;


void readInstance(char* filename) {
	//setup instance name
	strcpy(instance,filename);
	//open the file and check for errors
	FILE* f = fopen(filename,"r");
	if (!f) {
		cerr << "ERROR: Unable to open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	//read n and check for errors
	if (fscanf(f,"%d",&n)!=1) {
		cerr << "ERROR: Unable to read N from " << filename << endl;
		exit(EXIT_FAILURE);
	}
	//allocate memory and setup the pointers (Iliffe/display style)
	int n2 = n*n;
	hs = new int[n2];
	h = new int*[n];
	for (int i=0; i<n; i++)
		h[i] = &hs[i*n];
	//read the matrix entries and check for errors
	for (int i=0; i<n2; i++) {
		if (fscanf(f,"%d",&hs[i])!=1) {
			cerr << "ERROR: Unable to read " << i << "th (0-based) entry from " << filename << endl;
			exit(EXIT_FAILURE);
		}
	}
	//upper bound in the instance file are for LOP, so always INT_MAX here
	fitbound = INT_MAX;
	//close the file
	fclose(f);
	//done
}


void destroyInstance() {
	//deallocate memory
	delete[] h;
	delete[] hs;
	//done
}


int eval(int* x) {
	//TODO!!!
	return -1;



/*
	//sum up upper triangular part of h
	static int last = n-1;
	int fx = 0;
	int i,j,*r;
	for (i=0; i<last; i++) {
		r = h[x[i]];
		for (j=i+1; j<n; j++)
			fx += r[x[j]]; //the same of h[x[i]][x[j]] but faster
	}
	//return the sum
	return fx;
	//done
*/
}


void readHeuristic(char* filename) {
	//NOT USED FOR LOP!!!
	cerr << "ERROR: WE DONT HAVE HEURISTIC IMPLEMENTED FOR LOP!!!" << endl;
	exit(EXIT_FAILURE);
}


/*
double sparsity() {
	//count the 0 entries of h
	int s = 0;
	int n2 = n*n;
	for (int i=0; i<n2; i++)
		if (hs[i]==0)
			s++;
	//return s as percentual
	return 100*s/(double)n2;
	//done
}


double vc() {
	//compute avg of h entries
	int t = 0;
	int n2 = n*n;
	for (int i=0; i<n2; i++)
		t += hs[i];
	double avg = t/(double)n2;
	//compute std.dev. of h entries
	double std = 0.;
	for (int i=0; i<n2; i++)
		std += (hs[i]-avg)*(hs[i]-avg);
	std = sqrt(std/n2);
	//return the variation coefficient (std/avg)
	return std/avg;
	//done
}


double skewness() {
	//compute avg of h entries
	int t = 0;
	int n2 = n*n;
	for (int i=0; i<n2; i++)
		t += hs[i];
	double avg = t/(double)n2;
	//compute std.dev. of h entries
	double std = 0.;
	for (int i=0; i<n2; i++)
		std += (hs[i]-avg)*(hs[i]-avg);
	std = sqrt(std/n2);
	//compute Paerson gamma1 skewness (see http://en.wikipedia.org/wiki/Skewness)
	double gamma1 = 0.;
	for (int i=0; i<n2; i++)
		gamma1 += ((hs[i]-avg)/std) * ((hs[i]-avg)/std) * ((hs[i]-avg)/std);
	gamma1 /= n2;
	//return the skewness
	return gamma1;
	//done
}
*/

