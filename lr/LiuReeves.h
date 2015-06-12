#ifndef LIUREEVES_H
#define LIUREEVES_H

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cfloat>
#include "tft.h"
using namespace std;

typedef vector<int> vectint;

/* Liu Reeves Heuristics */

class LiuReeves {
	vector< vector<double> > comp_time; 
	vector<double> comp_time_artif, artif_proc_time;
	int nmach, njobs;
	double weighted_idle_time(int k,int i);
	double artificial_tft(int k,int i,set<int> remaining);
	double sequence(int i,vectint &p);
	double index(int i,set<int> rem);
public:
	LiuReeves(char* problem);
	double eval(int x,vectint &p);
};

#endif
