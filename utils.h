#ifndef UTILS_H
#define UTILS_H

#include <cstdio> //FILE

#define PERMSTR_SIZE 2048		//suggested size for a permutation string
#define PLUS_INF (1.0/0.0)		//plus infinity
#define MINUS_INF (-1.0/0.0)	//minus infinity

//convert the permutation p of length n to a string s
char* perm2str(int* p, int n, char* s=0);

//check if p of length n is a valid permutation
bool permValid(int* p, int n);

//convert the milliseconds millis to a time-string s (days, hours, minutes, seconds, milliseconds)
char* millis2str(unsigned long millis, char* s=0);

//convert the double n to a string s without trailing zeros
char* double2str(double n, char* s=0);

//discard a line from a file and return the last character ('\n' or EOF)
int discardLine(FILE* f);

//print a permutation to standard output
void printPerm(int* p, int n);

//print a permutation in python format
void printPermPython(int* p, int n);

//insertion sort algorithm that sort a basing on the values in v
void insertionSortWithValues(int* a, int l, int* v);

//check if an array is sorted basing on the values in v
bool isSortedWithValues(int* a, int l, int* v);

//check if an array is sorted
bool isSorted(int* a, int l);

//convert a string to a permutation
void str2perm(char* s, int* p, int n);

//read a permutation from the keyboard
void readPermFromKeyboard(int* x, int n);

//return the number of lines in the file f (it must be opened before and closed after)
int nlines(FILE* f);

#endif

