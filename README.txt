Please use the "Download ZIP" button of github page
---------------------------------------------------

to compile:
make

to run:
depms FSP_INSTANCE_FILE CSV_OUTPUT_FILE --heu NEH_FILE
deptft FSP_INSTANCE_FILE CSV_OUTPUT_FILE --heu LR_FILE
deplop LOP_INSTANCE_FILE CSV_OUTPUT_FILE
deplopcc LOPCC_INSTANCE_FILE CSV_OUTPUT_FILE

FSP_INSTANCE_FILE, NEH_FILE, LR_FILE are contained in the directory "pfsp_problems"
LOP_INSTANCE_FILE are contained in the directory "lop_problems"
LOPCC_INSTANCE_FILE are contained in the directory "lopcc_problems"

----------------------------------------------------

the directories "lr" and "neh" contain the software to generate LR_FILE and NEH_FILE given a FSP_INSTANCE_FILE

----------------------------------------------------

MANIFEST CONSTANTS FOR THE COMPILATION:

1) ONLINEPRINT: it enables printing-while-running code
2) MYDEBUG: to have some debug infos
3) GFC: enable/disable the use of the GFC data structure for fast fitness computation during the local search phase (ONLY VALID FOR PFSP PROBLEMS)

NOTE: You can enable the constants above by adding "-DONLINEPRINT" (et similia) at the row "CPPFLAGS" of the "Makefile"

----------------------------------------------------

VERSION HISTORY:

v1.0.0 - uploaded the 12 jun 2015 - first version
v1.0.1 - uploaded the 03 nov 2015 - differential mutation with insertions' generators added (dep.cpp.INSERTIONS)
v1.0.2 - uploaded the 16 nov 2015 - sorting sequence (instead of decomposition) on insertions (modified dep.cpp.INSERTIONS)
v1.0.3 - uploaded the 19 nov 2015 - differential mutation with exchanges' generators added (dep.cpp.EXCHANGES)
v1.0.4 - uploaded the 24 nov 2015 - limited/bounded version for randis,randss (dep.cpp.INSERTIONS,dep.cpp.EXCHANCES)
v1.0.5 - uploaded the 24 nov 2015 - bug fix = there was rounding instead of ceil on the diffMutations (dep.cpp.ORIGINAL,dep.cpp.INSERTIONS,dep.cpp.EXCHANGES,dep.cpp)
v1.0.6 - uploaded the 10 dec 2015 - more entropy on dep.cpp.INSERTIONS (previous version is now in dep.cpp.INSERTIONS_FASTER_LESS_ENTROPY)
v1.0.7 - uploaded the 15 dec 2015 - implemented the case F>1 for adjacent swaps (dep.cpp and dep.cpp.ORIGINAL), added the --minf,--maxf,--f parameter (currently maxf>1 does not work with other generators or other mutations), the output file also write minf,maxf, the save mechanism also saves minf,maxf
v1.0.8 - uploaded the 07 jan 2016 - bug fix for fixed f, insertions works only with F<=1, added python library for small experiments
v1.0.9 - uploaded the 13 jan 2016 - exchanges with F>1 implemented, reindentation everywhere, updated plib, minf/maxf renamed to fmin/fmax, bugfix on makefile
v1.1.0 - uploaded the 15 jan 2016 - different generators for diffMutation merged in a single file and using function pointers, write also maxtime on output, bugfixes on makefile (bad mutations files are no more valid for the moment)

v1.5.0 - uploaded the 27 jan 2016 - revisited everything, lop implemented, function pointers, more statistics, new crossovers, new selection, directory names changed, ...
v1.5.1 - uploaded the 01 mar 2016 - new selection methods introduced (some are not finalized or not working yet), FitnessType, lopcc sketched, new lopcc instances, ...
v1.5.2 - uploaded the 11 mar 2016 - bugfix on cr for obxcr, tpiicr crossover, 1 child possibility, maxStagnationTime, lopcc fully working, heuristics for lopcc, ...
v1.5.3 - uploaded the 14 mar 2016 - local searches implemented as standalone executable, only one getTimer in termination()

----------------------------------------------------







OLD CONTENT OF README NO MORE VALID FROM 15-01-2016!!!!!!!!!!!!!!!!!!!!!!!

OTHER MUTATIONS FOR TEST PURPOSES:

1) dep.cpp.ORIGINAL
it contains the code for the original differential mutation introduced in the paper (generators are adjacent swaps and randomized decomposition is performed by means of randomized bubble sort)
there is implementation for F>1 introduced on date 15 dec 2015

2) dep.cpp.PATHRELINKING 
it contains the code for the path relinking implementation

3) dep.cpp.RANDOMINDIVIDUAL 
it contains the code for the "random individual" mutation

4) dep.cpp.COMPLETELYRANDOM
it contains the code for the "completely random" mutation

5) dep.cpp.INSERTIONS
introduced on date 03-11-2015. Final version on date 10-12-2015. It is like (1) but it uses insertions as generators and randomized insertion sort
there is also the variant "dep.cpp.INSERTIONS_FASTER_LESS_ENTROPY" that it is a bit faster (though asymptotically the same") but it has "less entropy" on the produced decomposition
there is NO implementation for F>1 (lot of open problems for this case)

6) dep.cpp.EXCHANGES
introduced on date 19-11-2015. It is like (1) but it uses exchanges as generators and randomized (and generalized) selection sort
there is implementation for F>1 introduced on date 13 jan 2016

NOTE: In order to compile with the chosen mutation scheme, please copy the chosen file to dep.cpp and then compile
TODO: currently PATHRELINKING,RANDOMINDIVIDUAL,COMPLETELYRANDOM do not compile anymore for few errors due to some updates
