#ifndef UTILS_H
#define UTILS_H

#include <cstdio> //FILE

#define PERMSTR_SIZE 2048  //suggested size for a permutation string

//convert the permutation p of length n to a string s
char* perm2str(int* p, int n, char* s=0);

//check if p of length n is a valid permutation
bool permValid(int* p, int n);

//convert the milliseconds millis to a time-string s (days, hours, minutes, seconds, milliseconds)
char* millis2str(unsigned long millis, char* s=0);

//convert the double n to a string s without traling zeros
char* double2str(double n, char* s=0);

//discard a line from a file and return the last character ('\n' or EOF)
int discardLine(FILE* f);

//print a permutation to standard output
void printPerm(int* p, int n);

#endif

