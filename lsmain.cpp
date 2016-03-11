#include "problem.h"
#include "utils.h"
#include "random.h"
#include <cstdlib>
#include <iostream>
#include <cstring>
using namespace std;


#define USAGE_STRING "Usage: ./ls (vns4|ins) INSTANCE P[0] P[1] ... P[n-1]"


#include "depLocalSearch.cpp"


int main(int argc, char** argv) {
	//check command line parameters
	if (argc<2) {
		cerr << USAGE_STRING << endl;
		return EXIT_FAILURE;
	}
	//read ls type
	char lstype[100];
	strcpy(lstype,argv[1]);
	//read instance
	readInstance(argv[2]);
	//read seed solution
	int x[n];
	for (int i=0; i<n; i++)
		x[i] = atoi(argv[i+3]);
	//set local search function pointer
	void (*localSearch)(int*,FitnessType&);
	if (strcmp(lstype,"vns4")==0)
		localSearch = localSearch_vns4;
	else if (strcmp(lstype,"ins")==0)
		localSearch = localSearch_ins;
	else {
		cerr << "LS TYPE IS UNKNOWN!!!" << endl;
		return EXIT_FAILURE;
	}
	//run the local search
//TODO!!!
	//return success
	return EXIT_SUCCESS;
	//done
}

