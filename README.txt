Please use the "Download ZIP" button of github page
---------------------------------------------------

to compile:
make

to run:
depms FSP_INSTANCE_FILE CSV_OUTPUT_FILE --heu NEH_FILE
deptft FSP_INSTANCE_FILE CSV_OUTPUT_FILE --heu LR_FILE

FSP_INSTANCE_FILE, NEH_FILE, LR_FILE are contained in the directory "test_problems"

----------------------------------------------------

the directories "lr" and "neh" contain the software to generate LR_FILE and NEH_FILE given a FSP_INSTANCE_FILE

----------------------------------------------------

MANIFEST CONSTANTS FOR THE COMPILATION:

1) ONLINEPRINT: it enables printing-while-running code
2) MYDEBUG: to have some debug infos
3) GFC: enable/disable the use of the GFC data structure for fast fitness computation during the local search phase

NOTE: You can enable the constants above by adding "-DONLINEPRINT" (et similia) at the row "CPPFLAGS" of the "Makefile"

----------------------------------------------------

OTHER MUTATIONS FOR TEST PURPOSES:

1) dep.cpp.ORIGINAL
it contains the code for the original differential mutation introduced in the paper (generators are adjacent swaps and randomized decomposition is performed by means of randomized bubble sort)

2) dep.cpp.PATHRELINKING 
it contains the code for the path relinking implementation

3) dep.cpp.RANDOMINDIVIDUAL 
it contains the code for the "random individual" mutation

4) dep.cpp.COMPLETELYRANDOM
it contains the code for the "completely random" mutation

5) dep.cpp.INSERTIONS
introduced on date 03-11-2015. Final version on date 10-12-2015. It is like (1) but it uses insertions as generators and randomized insertion sort
there is also the variant "dep.cpp.INSERTIONS_FASTER_LESS_ENTROPY" that it is a bit faster (though asymptotically the same") but it has "less entropy" on the produced decomposition

6) dep.cpp.EXCHANGES
introduced on date 19-11-2015. It is like (1) but it uses exchanges as generators and randomized (and generalized) selection sort

NOTE: In order to compile with the chosen mutation scheme, please copy the chosen file to dep.cpp and then compile

----------------------------------------------------

VERSION HISTORY:

v1.00 - uploaded the 12 jun 2015 - first version
v1.01 - uploaded the 03 nov 2015 - differential mutation with insertions' generators added (dep.cpp.INSERTIONS)
v1.02 - uploaded the 16 nov 2015 - sorting sequence (instead of decomposition) on insertions (modified dep.cpp.INSERTIONS)
v1.03 - uploaded the 19 nov 2015 - differential mutation with exchanges' generators added (dep.cpp.EXCHANGES)
v1.04 - uploaded the 24 nov 2015 - limited/bounded version for randis,randss (dep.cpp.INSERTIONS,dep.cpp.EXCHANCES)
v1.05 - uploaded the 24 nov 2015 - bug fix = there was rounding instead of ceil on the diffMutations (dep.cpp.ORIGINAL,dep.cpp.INSERTIONS,dep.cpp.EXCHANGES,dep.cpp)
v1.06 - uploaded the 10 dec 2015 - more entropy on dep.cpp.INSERTIONS (previous version is now in dep.cpp.INSERTIONS_FASTER_LESS_ENTROPY)


