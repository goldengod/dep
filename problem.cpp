#include "problem.h"

//THE MANIFEST CONSTANTS "TFT" AND "MAKESPAN" ARE MUTUALLY EXCLUSIVE

#ifdef TFT
	#include "tft.cpp"
#endif

#ifdef MAKESPAN
	#include "makespan.cpp"
#endif

