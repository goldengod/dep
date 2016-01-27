#include "problem.h"

//THE MANIFEST CONSTANTS TFT/MAKESPAN/LOP ARE MUTUALLY EXCLUSIVE

#ifdef TFT
	#include "tft.cpp"
#endif

#ifdef MAKESPAN
	#include "makespan.cpp"
#endif

#ifdef LOP
	#include "lop.cpp"
#endif

