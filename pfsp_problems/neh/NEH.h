#ifndef NEH_H
#define NEH_H

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cfloat>
using namespace std;

typedef vector<int> vectint;

/* NEH Heuristics */

class NEH {
	vector< vector<int> > comp_time; 
	vector<int> sum_time;
	int nmach, njobs;
public:
	NEH(char* problem);
	int eval(vectint &p);
	int partial(vectint &,int kk);
	int sum_idle(vectint &);
};

#endif
