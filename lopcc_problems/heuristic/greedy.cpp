#include <cstdlib>
#include <cstdio>
#include <iostream>
using namespace std;


//global instance variables
int n;		//instance size
double** c;	//cost matrix (n*n)
double* p;	//cost vector (n)


//useful functions
void readInstance(char* filename);
double peval(int* x, int k);


//main function
int main(int argc, char** argv) {
	//some local variables
	int i,j,k,nu,jstar,*x,*u;
	double f,fstar;
	//check for command line parameter
	if (argc<2) {
		cerr << "Please provide the instance filename as command line parameter!!!" << endl;
		return EXIT_FAILURE;
	}
	//read the instance
	readInstance(argv[1]);
	//allocate memory
	x = new int[n];
	u = new int[n];
	//initialize unused items (all)
	for (i=0; i<n; i++)
		u[i] = i;
	nu = n;
	//loop to fill all the positions (k) backward
	for (k=n-1; k>0; k--) {
		//scan all the unused item to find the best one (jstar)
		fstar = 1./0.;	//+inf since it is a minimization problem
		for (j=0; j<nu; j++) {
			x[k] = u[j];
			f = peval(x,k);
			if (f<fstar) {	//< since it is a minimization problem
				fstar = f;
				jstar = j;
			}
		}
		//add u[jstar] to x and remove from u
		x[k] = u[jstar];
		u[jstar] = u[--nu];
		//done	
	}
	//only one choice for the last (first) position
	x[0] = u[0];
	//compute the full fitness value
	f = peval(x,0);	//full fitness is partial evaluation till 0
	//print the objective value and the solution x
	cout << "fitness\tsolution" << endl;
	cout << f << "\t";
	for (i=0; i<n; i++)
		cout << x[i] << " ";
	cout << endl;
	//free memory
	delete[] x;
	delete[] u;
	delete[] p;
	for (i=0; i<n; i++)
		delete[] c[i];
	delete[] c;
	//return success
	return EXIT_SUCCESS;
	//done
}


//skip a line in a file
int discardLine(FILE* f) {
	int c;
	do {
		c = fgetc(f);
	} while (c!='\n' && c!=EOF);
	return c;
}


//count the number of lines in a file
int nlines(FILE* f) {
	int nl = 0;
	while (discardLine(f)!=EOF)
		nl++;
	return nl;
}


//read the instance global variables from file
void readInstance(char* filename) {
	//open the file, read how many lines there are, set n, close the file
	FILE* f = fopen(filename,"r");
	if (!f) {
		cerr << "ERROR: Unable to open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	n = nlines(f);
	fclose(f);
	//open again the file and check for errors
	f = fopen(filename,"r");
	if (!f) {
		cerr << "ERROR: Unable to open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	//allocate memory
	c = new double*[n];
	for (int i=0; i<n; i++)
		c[i] = new double[n];
	p = new double[n];
	//read the matrix entries and check for errors
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			if (fscanf(f,"%lf",&c[i][j])!=1) {
				cerr << "ERROR: Unable to read the matrix entry (" << i << "," << j << ")" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	//read the vector entries and check for errors
	for (int i=0; i<n; i++) {
		if (fscanf(f,"%lf",&p[i])!=1) {
			cerr << "ERROR: Unable to read the vector entry " << i << endl;
			exit(EXIT_FAILURE);
		}
	}
	//close the file
	fclose(f);
	//done
}


//partially evaluate a solution x in the chunk [k,n)
double peval(int* x, int k) {
	//fit = SUM_i=0^n-1 alpha_{x_i}
	//alpha_{x_i} = p_{x_i} + SUM_j=i+1^n-1 ( c_{x_i,x_j} * alpha_{x_j} ) for i = n-1, n-2, ..., 0
	//some variables
	int i,j;
	double f;
	//alpha cost array need to be allocated the first time (doesnt mind if I will not free it at the end)
	static double* a = 0;
	if (!a)
		a = new double[n];
	//compute last a and initialize fitness
	f = a[n-1] = p[x[n-1]];
	//outer + inner summation (the outer only till k because it is a partial evaluation)
	for (i=n-2; i>=k; i--) {
		a[i] = p[x[i]];
		for (j=i+1; j<n; j++)
			a[i] += c[x[i]][x[j]] * a[j];
		f += a[i];
	}
	//return the partial fitness
	return f;
	//done
}

