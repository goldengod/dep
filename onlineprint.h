#ifndef ONLINEPRINT_H
#define ONLINEPRINT_H

#include <iostream>  //cerr, endl
using namespace std;


inline void genPrint() {
   //print nothing
}


inline void bestPrint() {
   //print nfes,ngen,nrestarts,fgbest
   cerr << "nfes=" << nfes << "\tngen=" << ngen
        << "\tnrestarts=" << nrestarts << "\tfgbest=" << fgbest
        << "\tls=" << gbestls << endl;
}

#endif

