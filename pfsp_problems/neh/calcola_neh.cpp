#include "makespan.h"
#include "NEH.h"
#include <iostream>

using namespace std;

int main(int argc,char **argv) {
	NEH h(argv[1]);
	vector<int> p(n);
	int v=h.eval(p);
	cout << "neh_fitness, neh_perm :" << endl;
	cout << "\t" << v << "\t";
	for(int i=0;i<n;i++) cout<<p[i]<<" ";
	cout << endl;
	return 0;
	/*
	cout << v << " ";
	for(int i=0;i<n;i++) cout << p[i] << " ";
	cout << endl;
	cout << h.sum_idle(p) << endl;
	*/
}


