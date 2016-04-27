CPP = g++
CC = gcc
CPPFLAGS = -O2 -Wall -DNDEBUG -DGFC -DONLINEPRINT
CCFLAGS = -Wall
INCLUDES = -I"/usr/include"
LFLAGS = -L"/usr/lib"
LIBS = 

all: deptft depms deplop deplopcc testtft testms testlop testlopcc lstft lsms lslop lslopcc

clean:
	rm *.o deptft depms deplop deplopcc testtft testms testlop testlopcc lstft lsms lslop lslopcc *.pyc

cleantilde:
	rm *~

SFMT.o: SFMT.c SFMT.h SFMT-common.h SFMT-params.h SFMT-params19937.h
	$(CC) -O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing --param max-inline-insns-single=1800 -Wmissing-prototypes -Wall -std=c99 -DSFMT_MEXP=19937 -c SFMT.c

random.o: random.cpp random.h
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c random.cpp $(LFLAGS) $(LIBS)

timer.o: timer.cpp timer.h
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c timer.cpp $(LFLAGS) $(LIBS)

utils.o: utils.cpp utils.h
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c utils.cpp $(LFLAGS) $(LIBS)

loglookup.o: loglookup.cpp loglookup.h
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c loglookup.cpp $(LFLAGS) $(LIBS)

problem_tft.o: problem.cpp problem.h tft.cpp
	$(CPP) -DTFT $(CPPFLAGS) $(INCLUDES) -c problem.cpp $(LFLAGS) $(LIBS) -o problem_tft.o

problem_ms.o: problem.cpp problem.h makespan.cpp
	$(CPP) -DMAKESPAN $(CPPFLAGS) $(INCLUDES) -c problem.cpp $(LFLAGS) $(LIBS) -o problem_ms.o
	
problem_lop.o: problem.cpp problem.h lop.cpp
	$(CPP) -DLOP $(CPPFLAGS) $(INCLUDES) -c problem.cpp $(LFLAGS) $(LIBS) -o problem_lop.o

problem_lopcc.o: problem.cpp problem.h lopcc.cpp
	$(CPP) -DLOPCC $(CPPFLAGS) $(INCLUDES) -c problem.cpp $(LFLAGS) $(LIBS) -o problem_lopcc.o
	
DEPINCLUDES = depCommon.cpp depPopInit.cpp depDiffMutation.cpp depCrossover.cpp depSelection.cpp depLocalSearch.cpp depPopRestart.cpp

dep_tft.o: dep.cpp dep.h $(DEPINCLUDES) save_resume.cpp
	$(CPP) -DTFT $(CPPFLAGS) $(INCLUDES) -c dep.cpp $(LFLAGS) $(LIBS) -o dep_tft.o

dep_ms.o: dep.cpp dep.h $(DEPINCLUDES) save_resume.cpp
	$(CPP) -DMAKESPAN $(CPPFLAGS) $(INCLUDES) -c dep.cpp $(LFLAGS) $(LIBS) -o dep_ms.o
	
dep_lop.o: dep.cpp dep.h $(DEPINCLUDES) save_resume.cpp
	$(CPP) -DLOP $(CPPFLAGS) $(INCLUDES) -c dep.cpp $(LFLAGS) $(LIBS) -o dep_lop.o

dep_lopcc.o: dep.cpp dep.h $(DEPINCLUDES) save_resume.cpp
	$(CPP) -DLOPCC $(CPPFLAGS) $(INCLUDES) -c dep.cpp $(LFLAGS) $(LIBS) -o dep_lopcc.o

DEPTFT_OBJECTS = SFMT.o random.o timer.o utils.o loglookup.o problem_tft.o dep_tft.o

DEPMS_OBJECTS = SFMT.o random.o timer.o utils.o loglookup.o problem_ms.o dep_ms.o

DEPLOP_OBJECTS = SFMT.o random.o timer.o utils.o loglookup.o problem_lop.o dep_lop.o

DEPLOPCC_OBJECTS = SFMT.o random.o timer.o utils.o loglookup.o problem_lopcc.o dep_lopcc.o

deptft: depmain.cpp dep_tft.o $(DEPTFT_OBJECTS)
	$(CPP) -DTFT -static $(CPPFLAGS) $(INCLUDES) depmain.cpp $(DEPTFT_OBJECTS) $(LFLAGS) $(LIBS) -lm -o deptft

depms: depmain.cpp dep_ms.o $(DEPMS_OBJECTS)
	$(CPP) -DMAKESPAN -static $(CPPFLAGS) $(INCLUDES) depmain.cpp $(DEPMS_OBJECTS) $(LFLAGS) $(LIBS) -lm -o depms
	
deplop: depmain.cpp dep_lop.o $(DEPLOP_OBJECTS)
	$(CPP) -DLOP -static $(CPPFLAGS) $(INCLUDES) depmain.cpp $(DEPLOP_OBJECTS) $(LFLAGS) $(LIBS) -lm -o deplop
	
deplopcc: depmain.cpp dep_lopcc.o $(DEPLOPCC_OBJECTS)
	$(CPP) -DLOPCC -static $(CPPFLAGS) $(INCLUDES) depmain.cpp $(DEPLOPCC_OBJECTS) $(LFLAGS) $(LIBS) -lm -o deplopcc

testtft: test.cpp problem_tft.o utils.o random.o SFMT.o
	$(CPP) -DTFT -static $(CPPFLAGS) $(INCLUDES) test.cpp problem_tft.o utils.o random.o SFMT.o $(LFLAGS) $(LIBS) -o testtft

testms: test.cpp problem_ms.o utils.o random.o SFMT.o
	$(CPP) -DMAKESPAN -static $(CPPFLAGS) $(INCLUDES) test.cpp problem_ms.o utils.o random.o SFMT.o $(LFLAGS) $(LIBS) -o testms

testlop: test.cpp problem_lop.o utils.o random.o SFMT.o
	$(CPP) -DLOP -static $(CPPFLAGS) $(INCLUDES) test.cpp problem_lop.o utils.o random.o SFMT.o $(LFLAGS) $(LIBS) -o testlop

testlopcc: test.cpp problem_lopcc.o utils.o random.o SFMT.o
	$(CPP) -DLOPCC -static $(CPPFLAGS) $(INCLUDES) test.cpp problem_lopcc.o utils.o random.o SFMT.o $(LFLAGS) $(LIBS) -o testlopcc

lstft: lsmain.cpp depLocalSearch.cpp problem_tft.o utils.o random.o SFMT.o timer.o
	$(CPP) -DTFT -static $(CPPFLAGS) $(INCLUDES) problem_tft.o utils.o random.o SFMT.o timer.o lsmain.cpp $(LFLAGS) $(LIBS) -lm -o lstft

lsms: lsmain.cpp depLocalSearch.cpp problem_ms.o utils.o random.o SFMT.o timer.o
	$(CPP) -DMAKESPAN -static $(CPPFLAGS) $(INCLUDES) problem_ms.o utils.o random.o SFMT.o timer.o lsmain.cpp $(LFLAGS) $(LIBS) -lm -o lsms

lslop: lsmain.cpp depLocalSearch.cpp problem_lop.o utils.o random.o SFMT.o timer.o
	$(CPP) -DLOP -static $(CPPFLAGS) $(INCLUDES) problem_lop.o utils.o random.o SFMT.o timer.o lsmain.cpp $(LFLAGS) $(LIBS) -lm -o lslop

lslopcc: lsmain.cpp depLocalSearch.cpp problem_lopcc.o utils.o random.o SFMT.o timer.o
	$(CPP) -DLOPCC -static $(CPPFLAGS) $(INCLUDES) problem_lopcc.o utils.o random.o SFMT.o timer.o lsmain.cpp $(LFLAGS) $(LIBS) -lm -o lslopcc

