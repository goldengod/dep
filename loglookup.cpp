#include "loglookup.h"
#include <cmath>

//lookup table
static double* loglookupTable = 0;

//factor to convert ln(x) to log_2(x)
static double ln2inv = 1./log(2);

//initialize
void initLoglookup(int max) {
	loglookupTable = new double[max]; //from 1 to max included
	loglookupTable[0] = 0.; //log(1) is 0
	for (int n=1; n<max; n++)
		loglookupTable[n] = log((double)(n+1));
}

//finalize
void freeLoglookup() {
	if (loglookupTable) {
		delete[] loglookupTable;
		loglookupTable = 0;
	}
}

//base 2
double log2lookup(int n) {
	return ln2inv * loglookupTable[n-1];
}

//base e
double loglookup(int n) {
	return loglookupTable[n-1];
}

