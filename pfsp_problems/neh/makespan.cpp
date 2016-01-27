#include "makespan.h"
#include "utils.h"
#include <cstdio>    //FILE, fopen, fclose, fscanf
#include <cstdlib>   //EXIT_FAILURE, exit
#include <iostream>  //cerr, endl
#include <cstring>   //strcpy
using namespace std;


//Branchless max between x and y, bits has to be sizeof(int)*8-1 (see http://aggregate.org/MAGIC)
#define FAST_MAX(x,y,bits) ((x)-((((x)-(y))>>(bits))&((x)-(y))))


//instance filename definition (declared extern in tft.h)
char instance[STRING_LENGTH];

//flowshop n,m definition (declared extern in tft.h)
int n;
int m;

//ptime matrix (n x m) interface definition (declared extern in tft.h)
int** ptime;

//liu reeves fitness and permutation definition (declared extern in tft.h)
int nehms;
int* nehp;

//maxnfes definition (declared extern in tft.h)
int maxnfes;

//best tft definition (declared extern in tft.h)
int best;

//inner ptimes matrix storage definition (Iliffe/display style)
static int* ptimes;

//inner aux memory to compute fitness
static int* aux;

//inner constants
const static int bits = sizeof(int)*8-1; //used in FAST_MAX


void readInstance(char* filename) {
   //allocate and setup instance name
   strcpy(instance,filename);
   //open the file and check for errors
   FILE* f = fopen(filename,"r");
   if (!f) {
      cerr << "ERROR: Unable to open " << filename << endl;
      exit(EXIT_FAILURE);
   }
   //discard a line (contains some text from ceberio) and check for eof error
   if (discardLine(f)==EOF)
      goto EOF_ERROR;
   //read n,m and discard the rest (it's useless) and check for eof error
   if (fscanf(f,"%d %d",&n,&m)!=2)
      goto EOF_ERROR;
   if (discardLine(f)==EOF)
      goto EOF_ERROR;
   //allocate memory and setup the pointers (Iliffe/display style)
   ptimes = new int[n*m];
   ptime = new int*[n];
   for (int i=0; i<n; i++)
      ptime[i] = ptimes+i*m;
   nehp = new int[n];
   aux = new int[m];
   //discard a line (contains some text from ceberio) and check for eof error
   if (discardLine(f)==EOF)
      goto EOF_ERROR;
   //read matrix entries and check for eof error
   //in the file it is m x n, but ptime is n x m to favor memory locality
   for (int i=0; i<m; i++)
      for (int j=0; j<n; j++)
         if (fscanf(f,"%d",&ptime[j][i])!=1)
            goto EOF_ERROR;
/*
   //finish to read last matrix line and check for eof error
   if (discardLine(f)==EOF)
      goto EOF_ERROR;
   //discard a line (it's text by me) and check for eof error
   if (discardLine(f)==EOF)
      goto EOF_ERROR;
   //read best and maxEvals, discard the rest of the line and check for eof error
   if (fscanf(f,"%d %d",&best,&maxnfes)!=2)
      goto EOF_ERROR;
   if (discardLine(f)==EOF)
      goto EOF_ERROR;
   //discard a line (it's text by me) and check for eof error
   if (discardLine(f)==EOF)
      goto EOF_ERROR;
   //read NEH makespan and permutation (by NEH heuristic) and check for eof error
   if (fscanf(f,"%d",&nehms)!=1)
      goto EOF_ERROR;
   for (int i=0; i<n; i++)
      if (fscanf(f,"%d",&nehp[i])!=1)
         goto EOF_ERROR;
*/
   //close the file
   fclose(f);
   //return
   return;
   //eof error handling code
   EOF_ERROR:
   cerr << "ERROR: " << filename << " has a bad format" << endl;
   exit(EXIT_FAILURE);
}


void destroyInstance() {
   //deallocate memory
   delete[] ptime;
   delete[] ptimes;
   delete[] nehp;
   delete[] aux;
   //done
}


/*
int eval(int* x) {
   //evaluate tft using "big" equation in ceberio's paper (using dynamic programming)
#ifdef MYDEBUG
   if (!permValid(x,n)) {
      cout<<"eval1"<<endl;
      exit(1);
   }
#endif
   static int lastm = m-1;
   int i,j;
   int job = x[0];
   aux[0] = ptime[job][0];          //part 1 of equation (2)
   for (j=1; j<m; j++)              //part 3 of equation (2)
      aux[j] = ptime[job][j] + aux[j-1];
   int ms = aux[lastm];
   for (i=1; i<n; i++) {
      job = x[i];
      aux[0] += ptime[job][0];      //part 2 of equation (2)
      for (j=1; j<m; j++)           //part 4 of equation (2)
         aux[j] = ptime[job][j] + FAST_MAX(aux[j],aux[j-1],bits);
      ms = FAST_MAX(ms,aux[lastm],bits);
   }
   return ms;
   //done
}
*/


int eval(int* x) {
   //evaluate tft using "big" equation in ceberio's paper (using dynamic programming)
#ifdef MYDEBUG
   if (!permValid(x,n)) {
      cout<<"eval1"<<endl;
      exit(1);
   }
#endif
   int i,j,ms,*pa,*pt; //next statement set also pa=aux and pt=ptime[x[0]]=ptimes+x[0]*m
   *(pa=aux) = *(pt=ptimes+(*x)*m); //part1 of eq. (2) -> aux[0]=ptime[x[0]][0]
   for (j=1; j<m; j++) {            //part3 of eq. (2) -> aux[j]=aux[j-1]+ptime[x[0]][j]
      ++pa;
      *pa = *(pa-1) + *(++pt);
   }
   ms = *pa;  //ms=aux[m-1]
   for (i=1; i<n; i++) { //next statement set also pa=aux and pt=ptime[x[i]]=ptimes+x[i]*m
      *(pa=aux) += *(pt=ptimes+(*(++x))*m);  //part2 of eq. (2) -> aux[0]+=ptime[x[i]][0]
      for (j=1; j<m; j++) {         //part4 of eq. (2) -> aux[j]=max(aux[j],aux[j-1])+ptime[x[i]][j]
         ++pa;
         *pa = FAST_MAX(*pa,*(pa-1),bits) + *(++pt);
      }
      ms = FAST_MAX(ms,*pa,bits); //ms=max(ms,aux[m-1])
   }
   return ms;
   //done
}

