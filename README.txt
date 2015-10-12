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

-------------------------------------------------------------------------------------------------------

OTHER MUTATIONS FOR TEST PURPOSES:

1) dep.cpp.ORIGINAL
it contains the code for the original differential mutation introduced in the paper

2) dep.cpp.PATHRELINKING 
it contains the code for the path relinking implementation

3) dep.cpp.RANDOMINDIVIDUAL 
it contains the code for the "random individual" mutation

4) dep.cpp.COMPLETELYRANDOM
it contains the code for the "completely random" mutation

NOTE: In order to compile with the chosen mutation scheme, please copy the chosen file to dep.cpp and then compile
