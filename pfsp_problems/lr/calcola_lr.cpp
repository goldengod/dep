#include "tft.h"
#include "LiuReeves.h"
#include <iostream>

using namespace std;

int main(int argc,char **argv) {
	LiuReeves h(argv[1]);
	vector<int> p(n);
	int v=h.eval(n/m,p); //LR(n/m)
	cout << "lrnm fitness, lrnm perm :" << endl;
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


