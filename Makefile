CPP = g++
CC = gcc
CPPFLAGS = -O2 -Wall -DNDEBUG -DGFC
CCFLAGS = -Wall
INCLUDES = -I"/usr/include"
LFLAGS = -L"/usr/lib"
LIBS = 

all: deptft depms testtft testms

clean:
	rm *.o deptft depms testtft testms *.pyc

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

problem_tft.o: problem.cpp problem.h tft.cpp
	$(CPP) -DTFT $(CPPFLAGS) $(INCLUDES) -c problem.cpp $(LFLAGS) $(LIBS) -o problem_tft.o

problem_ms.o: problem.cpp problem.h makespan.cpp
	$(CPP) -DMAKESPAN $(CPPFLAGS) $(INCLUDES) -c problem.cpp $(LFLAGS) $(LIBS) -o problem_ms.o

dep_tft.o: dep.cpp dep.h save_resume.cpp
	$(CPP) -DTFT $(CPPFLAGS) $(INCLUDES) -c dep.cpp $(LFLAGS) $(LIBS) -o dep_tft.o

dep_ms.o: dep.cpp dep.h save_resume.cpp
	$(CPP) -DMAKESPAN $(CPPFLAGS) $(INCLUDES) -c dep.cpp $(LFLAGS) $(LIBS) -o dep_ms.o

DEPTFT_OBJECTS = SFMT.o random.o timer.o utils.o problem_tft.o dep_tft.o

DEPMS_OBJECTS = SFMT.o random.o timer.o utils.o problem_ms.o dep_ms.o

deptft: depmain.cpp save_resume.cpp $(DEPTFT_OBJECTS)
	$(CPP) -static $(CPPFLAGS) $(INCLUDES) depmain.cpp $(DEPTFT_OBJECTS) $(LFLAGS) $(LIBS) -o deptft

depms: depmain.cpp save_resume.cpp $(DEPMS_OBJECTS)
	$(CPP) -static $(CPPFLAGS) $(INCLUDES) depmain.cpp $(DEPMS_OBJECTS) $(LFLAGS) $(LIBS) -o depms

testtft: test.cpp problem_tft.o utils.o random.o SFMT.o
	$(CPP) -static $(CPPFLAGS) $(INCLUDES) test.cpp problem_tft.o utils.o random.o SFMT.o $(LFLAGS) $(LIBS) -o testtft

testms: test.cpp problem_ms.o utils.o random.o SFMT.o
	$(CPP) -static $(CPPFLAGS) $(INCLUDES) test.cpp problem_ms.o utils.o random.o SFMT.o $(LFLAGS) $(LIBS) -o testms

