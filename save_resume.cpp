/*
EXECTIME SCHEME IN CASE OF SAVE/RESUME

at init (during variable definition in dep.cpp):
- savedTime=0

during save (depSave):
- tillNow=getTimer()
- save on file savedTime+tillNow

during load (depLoad):
- load savedTime from file

at end (dep and depResume):
- execTime=getTimer()+savedTime
*/


#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cstdio>


//function declarations
void depLoad(char* filename);

//save execution state/memory if haveToSave is  set and reset to false the flag
void depSave() {
	//some variables
	int i,k,*p;
	unsigned int j;
	unsigned char* bytes; //used to cast double and assume that sizeof(unsigned char)==1
	//time till now to manage execTime in case of save/resume
	unsigned long tillNow = getTimer();
	//reset the save flag to false
	haveToSave = false;
	//set again the alarm
	alarm(saveSeconds);
	//open saveFile and check for errors
	FILE* fsav = fopen(saveFile,"w");
	if (!fsav) { //if error print and return
		cerr << "WARNING: I was unable to open " << saveFile << endl;
		return;
	}
	//write on saveFile
	//save saveSeconds and saveFile before everything
	fprintf(fsav,"%u\n",saveSeconds);         //saveSeconds 1
	fprintf(fsav,"%s\n",saveFile);            //saveFile 2
	fprintf(fsav,"%s\n",scriptFile);          //scriptFile 2bis
	//save the rng
	writeRandState(fsav);                     //rng state 3
	//from tft.h/cpp save only "instance" and "maxnfes", other values are read using "readInstance"
	fprintf(fsav,"%s\n",instance);            //instance 4
	fprintf(fsav,"%d\n",maxnfes);             //maxnfes 5
	//other from global.h (except saveSeconds and saveFile written before)
	fprintf(fsav,"%u\n",seed);                //seed 6
	fprintf(fsav,"%s\n",out);                 //out 7
	fprintf(fsav,"%s\n",exe);                 //exe 8
	//values from dep.h (extern ones)
	fprintf(fsav,"%d\n",np);                  //np 9
	bytes = (unsigned char*)&finit;           //finit 10
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&ffmin;           //fmin 10bis
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&ffmax;           //fmax 10tris
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	fprintf(fsav,"%s\n",sgenerators);		  //sgenerators 10quadris
	fprintf(fsav,"%s\n",sinitialization);	  //sinitialization 10penta
	fprintf(fsav,"%s\n",sselection);          //sselection 10esa
	fprintf(fsav,"%s\n",scrossover);          //scrossover 10epta
	fprintf(fsav,"%s\n",slsearch);			  //slsearch 10octa
	fprintf(fsav,"%s\n",srestart);			  //srestart 10nine
	bytes = (unsigned char*)&crinit;          //crinit 10ten
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&crmin;           //crmin 10eleven
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&crmax;           //crmax 10twelve
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&alpha;           //alpha 11
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	fprintf(fsav,"%d\n",heu);                 //heu 12
	fprintf(fsav,"%d\n",ls);                  //ls 13
	bytes = (unsigned char*)&frfactor;        //frfactor 14
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	fprintf(fsav,"%d\n",inftype);			  //inftype 14bis
	fprintf(fsav,"%d\n",nchilds);			  //nchilds 14tris
	fprintf(fsav,"%u\n",maxStagnTime);		  //maxStagnTime 14quadris
#ifdef FIT_INT
	fprintf(fsav,"%d\n",targetFit);			  //targetFit 14penta
#else
	bytes = (unsigned char*)&targetFit;
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
#endif
	for (i=0; i<n; i++)                       //gbest 15
		fprintf(fsav,"%d ",gbest[i]);
	fprintf(fsav,"\n");
#ifdef FIT_INT
	fprintf(fsav,"%d\n",fgbest);              //fgbest 16
#else
	bytes = (unsigned char*)&fgbest;          //fgbest 16
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
#endif
	fprintf(fsav,"%d\n",nfesFoundAt);         //nfesFoundAt 17
	fprintf(fsav,"%lu\n",timeFoundAt);        //timeFoundAt 17b
	fprintf(fsav,"%d\n",stageFoundAt);        //stageFoundAt 18
	fprintf(fsav,"%d\n",nfes);                //nfes 19
	fprintf(fsav,"%d\n",ngen);                //ngen 20
	fprintf(fsav,"%lu\n",savedTime+tillNow);  //savedTime 21
	fprintf(fsav,"%d\n",nrestarts);           //nrestarts 22
	fprintf(fsav,"%d\n",nforcedrestarts);     //nforcedrestarts 23
	fprintf(fsav,"%d\n",minStageLength);      //minStageLength 24
	fprintf(fsav,"%d\n",maxStageLength);      //maxStageLength 25
	bytes = (unsigned char*)&avgStageLength;  //avgStageLength 26
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	fprintf(fsav,"%d\n",improvingStages);     //improvingStages 27
	fprintf(fsav,"%d\n",nfesls);              //nfesls 28
	fprintf(fsav,"%d\n",nls);                 //nls 29
	fprintf(fsav,"%d\n",nImprovingls);        //nImprovingls 30
#ifdef FIT_INT
	fprintf(fsav,"%d\n",totImprovingls);      //totImprovingls 31
#else
	bytes = (unsigned char*)&totImprovingls;  //totImprovingls 31
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
#endif
	fprintf(fsav,"%d\n",gbestls?1:0);         //gbestls 32
	//dep.cpp (static ones): ps and x-y1-y2-ix
	for (i=0; i<4*np*n; i++)                  //ps saved completely 33
		fprintf(fsav,"%d ",ps[i]);
	fprintf(fsav,"\n");
	for (i=0; i<np; i++) {                    //x-y1-y2-ix saved as a 4*np string of x|1|2|i 34
		for (k=0; k<4; k++) {
			p = ps + (i*4+k)*n;
			if (p==x[i])
				fprintf(fsav,"x");
			else if (p==y1[i])
				fprintf(fsav,"1");
			else if (p==y2[i])
				fprintf(fsav,"2");
			else //if (p==ix[i])
				fprintf(fsav,"i");
		}
	}
	fprintf(fsav,"\n");
	//dep.cpp (static ones): fs and no need to save fx-fy1-fy2
#ifdef FIT_INT
	for (i=0; i<3*np; i++)                    //fs saved completely 35
		fprintf(fsav,"%d ",fs[i]);
	fprintf(fsav,"\n");
#else
	for (i=0; i<3*np; i++) {                  //fs saved completely 35
		bytes = (unsigned char*)&(fs[i]);
		for (j=0; j<sizeof(double); j++)
			fprintf(fsav,"%u ",bytes[j]);
	}
	fprintf(fsav,"\n");
#endif
	//dep.cpp (static ones): sfs and no need to save sfx-sfy
	for (i=0; i<2*np; i++) {                  //sfs saved completely 36
		bytes = (unsigned char*)&(sfs[i]);
		for (j=0; j<sizeof(double); j++)
			fprintf(fsav,"%u ",bytes[j]);
	}
	fprintf(fsav,"\n");
	//dep.cpp (static ones): crs and no need to save crx-cry
	for (i=0; i<2*np; i++) {                  //crs saved completely 36bis
		bytes = (unsigned char*)&(crs[i]);
		for (j=0; j<sizeof(double); j++)
			fprintf(fsav,"%u ",bytes[j]);
	}
	fprintf(fsav,"\n");
	//dep.cpp other static values (no need to save tmpint)
	fprintf(fsav,"%d\n",sameFitness?1:0);     //sameFitness 37
	fprintf(fsav,"%d\n",lastRestart);         //lastRestart 38
#ifdef FIT_INT
	fprintf(fsav,"%d\n",fgbestAtStageStart);  //fgbestAtStageStart 39
#else
	bytes = (unsigned char*)&fgbestAtStageStart;//fgbestAtStageStart 39
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
#endif
	fprintf(fsav,"%d\n",permByteSize);        //permByteSize (maybe useless) 40
	for (i=0; i<n; i++)                       //lsHistory 41
		fprintf(fsav,"%d ",lsHistory[i]);
	fprintf(fsav,"\n");
	fprintf(fsav,"%d\n",lsmode?1:0);          //lsmode (useless since I cant save in lsmode) 42
	fprintf(fsav,"%d\n",nfesWhenToForceRestart);//nfesWhenToForceRestart 43
	fprintf(fsav,"%d\n",forcedRestartPeriod); //forcedRestartPeriod 44
	fprintf(fsav,"%d\n",diameter);            //diameter 44bis
	//dep.cpp/h improving steps
	fprintf(fsav,"%d\n",improvingSteps);      //improvingSteps 45
	fprintf(fsav,"%d\n",lsImprovingSteps);    //lsImprovingSteps 46
	//dep.cpp/h f,cr,child statistics
	bytes = (unsigned char*)&sfSuccAvg;       //sfSuccAvg 47
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&sfSuccVar;       //sfSuccVar 48
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&sfSuccMin;       //sfSuccMin 49
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&sfSuccMax;       //sfSuccMax 50
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&crSuccAvg;       //crSuccAvg 51
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&crSuccVar;       //crSuccVar 52
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&crSuccMin;       //crSuccMin 53
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	bytes = (unsigned char*)&crSuccMax;       //crSuccMax 54
	for (j=0; j<sizeof(double); j++)
		fprintf(fsav,"%u ",bytes[j]);
	fprintf(fsav,"\n");
	fprintf(fsav,"%d\n",child1succ);          //child1succ 55
	fprintf(fsav,"%d\n",child2succ);          //child2succ 56
	//close saveFile
	fclose(fsav);
	//call the script if it is not "none"
	if (strcmp(scriptFile,"none")!=0) {
		char cmd[STRING_LENGTH];
		sprintf(cmd,"./%s %s",scriptFile,saveFile);
		system(cmd); //this is a blocking call!!!
	}
	//write a state message on std error
	time_t now = time(0);
	char* snow = ctime(&now);
	snow[strlen(snow)-1] = '\0';
	cerr << "[" << snow << "] State saved on " << saveFile << endl;
	//done
}


//dep function to resume a saved execution
void depResume(char* filename) {
	//load memory/state from the file filename
	depLoad(filename);
	//evolve population as in "dep" but without "popInit"
	//init timer
	setTimer();
	//set a flag to distinguish termination inside popRestart (due to local search) or popEvolve
	bool lsTermination = false; //needed to correctly handle restart statistics
	//evolution loop
	do {
		//some printing
#ifdef ONLINEPRINT
		genPrint();
#endif
		//check and save dep execution state/memory
		if (haveToSave)
			depSave();
		//check and restart population
		if (popRestart() || popForcedRestart()) { //considering the short circuit
			//if true, maxnfes was exhausted during local search, so break evolution loop
			lsTermination = true;
			break;
		}
		//evolve population
	} while (popEvolve());
	//compute final execution time
	execTime = getTimer() + savedTime;
	//adjust ngen and restart statistics (but only if the termination doesnt happen in a local search)
	if (lsTermination) {
		ngen++;
		nrestarts++; //let restartStatistics() believes that there was a restart
		restartStatistics();
		nrestarts--; //undo the dummy restart
	}
	//some final printing
#ifdef ONLINEPRINT
	genPrint();
#endif
	//done
}


//load memory/state from the file filename
void depLoad(char* filename) {
	//some variables
	int i,k,*p,nowarning; //nowarning is to avoid compiler complaints
	unsigned int j,b;
	unsigned char* bytes; //used to cast double and assume that sizeof(unsigned char)==1
	char c;
	//open filename and check for errors
	FILE* fsav = fopen(filename,"r");
	if (!fsav) { //if error print and exit
		cerr << "WARNING: I was unable to open " << filename << endl;
		exit(EXIT_FAILURE);
	}
	//read from the file
	//read saveSeconds and saveFile read first - ALREADY READ BUT I HAVE TO SKIP
	nowarning = fscanf(fsav,"%u",&saveSeconds);           //saveSeconds 1
	nowarning = fscanf(fsav,"%s",saveFile);               //saveFile 2
	nowarning = fscanf(fsav,"%s",scriptFile);             //scriptFile 2bis
	//read the rng - ALREADY READ BUT I HAVE TO SKIP
	readRandState(fsav);                                  //rng state 3
	//values of tft.h/cpp: read "instance" and read maxnfes  - ALREADY READ BUT I HAVE TO SKIP
	nowarning = fscanf(fsav,"%s",instance);               //instance 4
	nowarning = fscanf(fsav,"%d",&maxnfes);               //maxnfes 5
	//other of global.h (except saveSeconds and saveFile read at the begin)  - ALREADY READ BUT I HAVE TO SKIP
	nowarning = fscanf(fsav,"%u",&seed);                  //seed 6
	nowarning = fscanf(fsav,"%s",out);                    //out 7
	nowarning = fscanf(fsav,"%s",exe);                    //exe 8
	//values of dep.h (extern ones but only input)  - ALREADY READ BUT I HAVE TO SKIP
	nowarning = fscanf(fsav,"%d",&np);                    //np 9
	bytes = (unsigned char*)&finit;                       //finit 10
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&ffmin;                       //fmin 10bis
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&ffmax;                       //fmax 10tris
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%s",sgenerators);			  //sgenerators 10quadris
	nowarning = fscanf(fsav,"%s",sinitialization);		  //sinitialization 10penta
	nowarning = fscanf(fsav,"%s",sselection);             //sselection 10esa
	nowarning = fscanf(fsav,"%s",scrossover);			  //scrossover 10epta
	nowarning = fscanf(fsav,"%s",slsearch);			  	  //slsearch 10octa
	nowarning = fscanf(fsav,"%s",srestart);               //srestart 10nine
	bytes = (unsigned char*)&crinit;                      //crinit 10ten
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crmin;                       //crmin 10eleven
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crmax;                       //crmax 10twelve
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&alpha;                       //alpha 11
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%d",&heu);                    //heu 12
	nowarning = fscanf(fsav,"%d",&ls);                    //ls 13
	bytes = (unsigned char*)&frfactor;                    //frfactor 14
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%d",&inftype);               //inftype 14bis
	nowarning = fscanf(fsav,"%d",&nchilds);				  //nchilds 14tris
	nowarning = fscanf(fsav,"%u",&maxStagnTime);		  //maxStagnTime 14quadris
#ifdef FIT_INT
	nowarning = fscanf(fsav,"%d",&targetFit);			  //targetFit 14penta
#else
	bytes = (unsigned char*)&targetFit;
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
#endif
	//values of dep.h (extern ones but only output)
	for (i=0; i<n; i++)                                   //gbest 15
		nowarning = fscanf(fsav,"%d ",&(gbest[i]));
#ifdef FIT_INT
	nowarning = fscanf(fsav,"%d",&fgbest);                //fgbest 16
#else
	bytes = (unsigned char*)&fgbest;                      //fgbest 16
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
#endif
	nowarning = fscanf(fsav,"%d",&nfesFoundAt);           //nfesFoundAt 17
	nowarning = fscanf(fsav,"%lu",&timeFoundAt);          //timeFoundAt 17b
	nowarning = fscanf(fsav,"%d",&stageFoundAt);          //stageFoundAt 18
	nowarning = fscanf(fsav,"%d",&nfes);                  //nfes 19
	nowarning = fscanf(fsav,"%d",&ngen);                  //ngen 20
	nowarning = fscanf(fsav,"%lu",&savedTime);            //savedTime 21
	nowarning = fscanf(fsav,"%d",&nrestarts);             //nrestarts 22
	nowarning = fscanf(fsav,"%d",&nforcedrestarts);       //nforcedrestarts 23
	nowarning = fscanf(fsav,"%d",&minStageLength);        //minStageLength 24
	nowarning = fscanf(fsav,"%d",&maxStageLength);        //maxStageLength 25
	bytes = (unsigned char*)&avgStageLength;              //avgStageLength 26
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%d",&improvingStages);       //improvingStages 27
	nowarning = fscanf(fsav,"%d",&nfesls);                //nfesls 28
	nowarning = fscanf(fsav,"%d",&nls);                   //nls 29
	nowarning = fscanf(fsav,"%d",&nImprovingls);          //nImprovingls 30
#ifdef FIT_INT
	nowarning = fscanf(fsav,"%d",&totImprovingls);        //totImprovingls 31
#else
	bytes = (unsigned char*)&totImprovingls;              //totImprovingls 31
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
#endif
	nowarning = fscanf(fsav,"%d",&i);                     //gbestls 32
	gbestls = i==0 ? false : true;
	//dep.cpp (static ones): ps and x-y1-y2-ix
	for (i=0; i<4*np*n; i++)                              //ps read completely 33
		nowarning = fscanf(fsav,"%d",ps+i);
	char x12i[4*np];                                      //x-y1-y2-ix read as a string 4*np of x|1|2|i 34
	nowarning = fscanf(fsav,"%s",x12i);
	for (i=0; i<np; i++) {
		for (k=0; k<4; k++) {
			c = x12i[i*4+k];
			p = ps+(i*4+k)*n;
			if (c=='x')
				x[i] = p;
			else if (c=='1')
				y1[i] = p;
			else if (c=='2')
				y2[i] = p;
			else //c=='i'
				ix[i] = p;
		}
	}
	//dep.cpp (static ones): fs and no need to touch fx-fy1-fy2
#ifdef FIT_INT
	for (i=0; i<3*np; i++)                                //fs read completely 35
		nowarning = fscanf(fsav,"%d",&(fs[i]));
#else
	for (i=0; i<3*np; i++) {                              //fs read completely 35
		bytes = (unsigned char*)&(fs[i]);
		for (j=0; j<sizeof(double); j++) {
			nowarning = fscanf(fsav,"%u",&b);
			bytes[j] = (unsigned char)b;
		}
	}
#endif
	//dep.cpp (static ones): sfs and no need to touch sfx-sfy
	for (i=0; i<2*np; i++) {                              //sfs read completely 36
		bytes = (unsigned char*)&(sfs[i]);
		for (j=0; j<sizeof(double); j++) {
			nowarning = fscanf(fsav,"%u",&b);
			bytes[j] = (unsigned char)b;
		}
	}
	//dep.cpp (static ones): crs and no need to touch crx-cry
	for (i=0; i<2*np; i++) {                              //crs read completely 36bis
		bytes = (unsigned char*)&(crs[i]);
		for (j=0; j<sizeof(double); j++) {
			nowarning = fscanf(fsav,"%u",&b);
			bytes[j] = (unsigned char)b;
		}
	}
	//dep.cpp other static values (no need to touch tmpint)
	nowarning = fscanf(fsav,"%d",&i);                     //sameFitness 37
	sameFitness = i==0 ? false : true;
	nowarning = fscanf(fsav,"%d",&lastRestart);           //lastRestart 38
#ifdef FIT_INT
	nowarning = fscanf(fsav,"%d",&fgbestAtStageStart);    //fgbestAtStageStart 39
#else
	bytes = (unsigned char*)&fgbestAtStageStart;          //fgbestAtStageStart 39
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
#endif
	nowarning = fscanf(fsav,"%d",&permByteSize);          //permByteSize (maybe useless) 40
	for (i=0; i<n; i++)                                   //lsHistory 41
		nowarning = fscanf(fsav,"%d",&(lsHistory[i]));
	nowarning = fscanf(fsav,"%d",&i);                     //lsmode (useless since I cant save in lsmode) 42
	lsmode = i==0 ? false : true;
	nowarning = fscanf(fsav,"%d",&nfesWhenToForceRestart);//nfesWhenToForceRestart 43
	nowarning = fscanf(fsav,"%d",&forcedRestartPeriod);   //forcedRestartPeriod 44
	nowarning = fscanf(fsav,"%d",&diameter);              //diameter 44bis
	//dep.cpp/h improving steps
	nowarning = fscanf(fsav,"%d",&improvingSteps);        //improvingSteps 45
	nowarning = fscanf(fsav,"%d",&lsImprovingSteps);      //lsImprovingSteps 46
	//dep.cpp/h f,cr,child statistics
	bytes = (unsigned char*)&sfSuccAvg;                   //sfSuccAvg 47
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&sfSuccVar;                   //sfSuccVar 48
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&sfSuccMin;                   //sfSuccMin 49
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&sfSuccMax;                   //sfSuccMax 50
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crSuccAvg;                   //crSuccAvg 51
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crSuccVar;                   //crSuccVar 52
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crSuccMin;                   //crSuccMin 53
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	bytes = (unsigned char*)&crSuccMax;                   //crSuccMax 54
	for (j=0; j<sizeof(double); j++) {
		nowarning = fscanf(fsav,"%u",&b);
		bytes[j] = (unsigned char)b;
	}
	nowarning = fscanf(fsav,"%d",&child1succ);            //child1succ 55
	nowarning = fscanf(fsav,"%d",&child2succ);            //child2succ 56
	//close the file
	fclose(fsav);
	//statement to avoid compiler complaints
	nowarning++;
	//done
}

