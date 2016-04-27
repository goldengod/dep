#include "problem.h"
#include "utils.h"
#include "random.h"
#include <cstdlib>   //EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>  //cerr, cout, endl
#include <iomanip>
using namespace std;


#define USAGE_STRING "Usage: ./test INSTANCE P[0] P[1] ... P[n-1]"


int main(int argc, char** argv) {
	
	if (argc<2) {
		cerr << USAGE_STRING << endl;
		return EXIT_FAILURE;
	}
	readInstance(argv[1]);
	int x[n];
	for (int i=0; i<n; i++)
		x[i] = atoi(argv[i+2]);
	cout << "fitness = " << fixed << setprecision(10) << eval(x) << endl;
	char s[PERMSTR_SIZE];
	perm2str(x,n,s);
	cout << "permutation = " << s << endl;

	return EXIT_SUCCESS;
}

