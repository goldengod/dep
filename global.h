#ifndef GLOBAL_H
#define GLOBAL_H

#include "problem.h"					//I need it for FitnessType

#define STRING_LENGTH 256

//these variables are defined in deptft.cpp and must be visible in dep.cpp
extern unsigned int seed;
extern char out[STRING_LENGTH];
extern char exe[STRING_LENGTH];
extern unsigned int saveSeconds;
extern char saveFile[STRING_LENGTH];
extern char scriptFile[STRING_LENGTH];	//for grid version
extern unsigned int maxTime;			//for stop after maxTime milliseconds
extern unsigned int maxStagnTime;		//for stop after maxStagnTime milliseconds of stagnation
extern FitnessType targetFit;			//for stop when target fitness is reached

#endif

