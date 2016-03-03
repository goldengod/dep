#include "utils.h"
#include <cstdio>    //sprintf, FILE, fgetc, EOF, sscanf
#include <cstring>   //strlen, memset
#include <iostream>  //cout, endl, cin, cerr
#include <cstdlib>	 //exit, EXIT_FAILURE, atoi
using namespace std;



char* perm2str(int* p, int n, char* s) {
	if (!s)
		s = new char[PERMSTR_SIZE];
	sprintf(s,"%d",p[0]);
	char* ptr;
	for (int i=1; i<n; i++) {
		ptr = s+strlen(s);
		sprintf(ptr," %d",p[i]);
	}
	return s;
}


bool permValid(int* p, int n) {
	bool objects[n];
	memset(objects,0,sizeof(bool)*n);
	for (int i=0; i<n; i++) {
		if (p[i]<0 || p[i]>=n || objects[p[i]])
			return false;
		objects[p[i]] = true;
	}
	return true;
}


char* millis2str(unsigned long millis, char* s) {
	if (!s)
		s = new char[256]; //should be enough
	int ms = millis%1000;
	int sec = millis/1000;
	int min = sec/60;
	sec %= 60;
	int h = min/60;
	min %= 60;
	int d = h/24;
	h %= 24;
	if (d>0)
		sprintf(s,"%d days %d hours %d minutes %d.%d seconds",d,h,min,sec,ms);
	else if (h>0)
		sprintf(s,"%d hours %d minutes %d.%d seconds",h,min,sec,ms);
	else if (min>0)
		sprintf(s,"%d minutes %d.%d seconds",min,sec,ms);
	else
		sprintf(s,"%d.%d seconds",sec,ms);
	return s;
}


char* double2str(double n, char* s) {
	if (!s)
		s = new char[32]; //should be enough
	sprintf(s,"%.10lf",n); //10 decimal digits
	int i = strlen(s);
	do {
		--i;
	} while (s[i]=='0' && i>0);
	if (s[i]=='.')
	s[i] = '\0';
	else
	s[i+1] = '\0';
	return s;
}


int discardLine(FILE* f) {
	int c;
	do {
		c = fgetc(f);
	} while (c!='\n' && c!=EOF);
	return c;
}


void printPerm(int* p, int n) {
	for (int i=0; i<n; i++)
		cout << p[i] << (i<n-1?",":"");
	cout << endl;
}


void printPermPython(int* p, int n) {
	cout << "[ ";
	for (int i=0; i<n; i++)
		cout << p[i] << (i<n-1 ? "," : "") << " ";
	cout << "]" << endl;	
}


void insertionSortWithValues(int* a, int l, int* v) {
	//https://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Insertion_sort
	int i,j,t;
	for(i=1; i<l; i++) {
		t = a[i];
		for (j=i-1; j>=0 && v[a[j]]>v[t]; j--)
			a[j+1] = a[j];
		a[j+1] = t;
	}
}


bool isSortedWithValues(int* a, int l, int* v) {
	for (int i=1; i<l; i++)
		if (v[a[i-1]]>v[a[i]])
			return false;
	return true;
}


bool isSorted(int* a, int l) {
	for (int i=1; i<l; i++)
		if (a[i-1]>a[i])
			return false;
	return true;
}


void str2perm(char* s, int* p, int n) {
	for (int i=0; i<n; i++)
		if (sscanf(s,"%d",&p[i])!=1) {
			cerr << "ERROR IN STR2PERM!!!" << endl;
			exit(EXIT_FAILURE);
		}
}


void readPermFromKeyboard(int* x, int n) {
	char s[10];
	for (int i=0; i<n; i++) {
		cin >> s;
		x[i] = atoi(s);
	}
}


int nlines(FILE* f) {
	int nl = 0;
	while (discardLine(f)!=EOF)
		nl++;
	return nl;
}

