#ifndef ONLINEPRINT_H
#define ONLINEPRINT_H

#include <iostream>  //cerr, endl
//#include <algorithm>//VALE!!!!!!!!!!!!!!!!!!!!!!!!
using namespace std;


inline void genPrint() {
	//print nothing
	
	/*
	static int* c = new int[n*n];
	memset(c,0,sizeof(int)*(n*n));
	for (int k=0; k<np; k++) {
		for (int i=0; i<n-1; i++)
			for (int j=i+1; j<n; j++)
				c[x[k][i]*n+x[k][j]]++;
	}
	*/
	
	/*
	double fitavg = 0.;
	for (int i=0; i<np; i++)
		fitavg += fx[i];
	fitavg /= np;
	
	int fit[np];
	memcpy(fit,fx,sizeof(int)*np);
	sort(fit,fit+np);
	
	cout << "fitness media = " << ((int)fitavg) << "\tfit 1st = " << fit[np-1] << "\t fit 2nd = " << fit[np-2] << "\t best sofar = " << fgbest << endl;
	*/
	
}


inline void bestPrint() {
	//print nfes,ngen,nrestarts,fgbest
	cerr << "nfes=" << nfes << "\tngen=" << ngen
		<< "\tnrestarts=" << nrestarts << "\tfgbest=" << fgbest
		<< "\tls=" << gbestls << endl;
}

#endif

