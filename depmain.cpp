#include "global.h"
#include "random.h"
#include "problem.h"
#include "dep.h"
#include "utils.h"
#include <cstdlib>   //exit, atoi, EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>  //cout, endl, cerr
#include <cstring>   //strcpy, strcmp, strlen, strcat, strcmp
#include <cstdio>    //sscanf, FILE, fopen, fclose, fprintf, sprintf, fscanf
#include <csignal>   //signal, SIGALRM
#include <unistd.h>  //alarm
#include <climits>   //INT_MAX
using namespace std;


//these variables are declared extern in global.h (since have to be visible to dep.cpp)
unsigned int seed = 0;
char out[STRING_LENGTH];
char exe[STRING_LENGTH];
unsigned int saveSeconds = 0;
char saveFile[STRING_LENGTH];
char scriptFile[STRING_LENGTH];	//for grid version
unsigned int maxTime = 0;

//inner variables
static bool resume;


void usage();
void readArguments(int argc, char** argv);
void writeResults();
void defaultSaveFile();
void saveHandler(int signum);
void preResume(char* filename);
void defaultMaxnfes();


int main(int argc, char** argv) {   
   //set default parameters for dep
   depDefaultParameters();
   //read arguments, set parameters and load pfsp-tft instance
   readArguments(argc,argv);
   //check if maxTime was set, in the case set maxnfes to +inf
   if (maxTime>0)
      maxnfes = INT_MAX;
   //check if it's normal or resume execution
   if (!resume) { //normal execution
      //install the alarm for saving in the case save was set
      if (saveSeconds>0) {
         signal(SIGALRM,saveHandler);
         alarm(saveSeconds);
      }
      //random seed if zero
      if (!seed)
         seed = randSeed(); //from /dev/random
      //init random number generator
      initRand(seed);
      //set default save file if the case
      if (saveSeconds>0 && saveFile[0]=='\0')
         defaultSaveFile();
      //done
   } else { //resume execution
      //call preResume to read inputs and instance
      preResume(saveFile);
      //install the alarm for saving in the case save was set
      signal(SIGALRM,saveHandler);
      alarm(saveSeconds);
      //done
   }
   //allocate dep memory
   depAlloc();
   //some printing before dep execution
   cout << "STARTING DEPTFT:" << endl;
   cout << "               exe = " << exe << endl;
   cout << " -        instance = " << instance << endl;
   cout << " -               n = " << n << endl;
   cout << " -               m = " << m << endl;
   cout << " -             out = " << out << endl;
   cout << " -         maxnfes = " << maxnfes << endl;
   cout << " -         maxTime = " << maxTime << endl;
   cout << " -            seed = " << seed << endl;
   cout << " -              np = " << np << endl;
   cout << " -           finit = " << finit << endl;
   cout << " -           alpha = " << alpha << endl;
   cout << " -             heu = " << heu << endl;
   cout << " -              ls = " << ls << " ";
   switch (ls) {
      case N_LS:
         cout << "(N_LS)" << endl;
         break;
      case B_LS:
         cout << "(B_LS)" << endl;
         break;
      case L_LS:
         cout << "(L_LS)" << endl;
         break;
      default:
         cerr << "ERROR: ls has to be 0,1,2!!!" << endl;
         exit(EXIT_FAILURE);
   }
   cout << " -        frfactor = " << frfactor << endl;
   cout << " -     saveSeconds = " << saveSeconds;
   if (saveSeconds==0)
      cout << " (no save)" << endl;
   else {
      cout << endl << " -        saveFile = " << saveFile << endl;
      cout << " -      scriptFile = " << scriptFile << endl;
   }
   cout << "... DEPTFT IS RUNNING ..." << endl;
   //execute dep basing on normal or resume execution mode
   if (!resume)
      dep(); //normal execution
   else
      depResume(saveFile); //resume execution
   //some printing after dep execution
   cout << "DEPTFT DONE:" << endl;
   cout << " -          fgbest = " << fgbest << endl;
   cout << " -     nfesFoundAt = " << nfesFoundAt << endl;
   cout << " -    stageFoundAt = " << stageFoundAt << endl;
   cout << " -            nfes = " << nfes << endl;
   cout << " -            ngen = " << ngen << endl;
   cout << " -       nrestarts = " << nrestarts << endl;
   cout << " - nforcedrestarts = " << nforcedrestarts << endl;
   cout << " -  minStageLength = " << minStageLength << endl;
   cout << " -  maxStageLength = " << maxStageLength << endl;
   cout << " -  avgStageLength = " << avgStageLength << endl;
   cout << " - improvingStages = " << improvingStages << endl;
   cout << " -             nls = " << nls << endl;
   cout << " -          nfesls = " << nfesls << endl;
   cout << " -    nImprovingls = " << nImprovingls << endl;
   cout << " -  totImprovingls = " << totImprovingls << endl;
   cout << " -         gbestls = " << gbestls << endl;
   cout << " -  improvingSteps = " << improvingSteps << endl;
   cout << " -lsImprovingSteps = " << lsImprovingSteps << endl;
   char stime[256];
   millis2str(execTime,stime);
   cout << " -        execTime = " << stime << endl;
   //check gbest validity and fitness
   bool gbestValid = permValid(gbest,n);
   bool gbestFitness = eval(gbest)==fgbest;
   //write results in output file
   writeResults();
   //some final printing
   cout << "Full results are available in " << out << endl;
   //free dep memory
   depFree();
   //destroy instance memory
   destroyInstance();
   //check gbestValid/gbestFitness and print error
   if (!gbestValid)
      cerr << "ERROR: GBEST IS NOT A VALID PERMUTATION!!!" << endl;
   if (!gbestFitness)
      cerr << "ERROR: GBEST FITNESS DOES NOT MATCH!!!" << endl;
   //return success/failure basing on gbestValid/gbestFitness
   return gbestValid && gbestFitness ? EXIT_SUCCESS : EXIT_FAILURE;
   //done
}


void usage() {
   cout << "-------------" << endl;
   cout << "DEP USAGE:" << endl;
   cout << "(1) dep INSTANCE_FILE OUTPUT_FILE" << endl;
   cout << "    [--seed UINT]" << endl;
   cout << "    [--maxnfes INT]" << endl;
   cout << "    [--np INT]" << endl;
   cout << "    [--finit DOUBLE]" << endl;
   cout << "    [--alpha DOUBLE]" << endl;
   cout << "    [--heu STRING(filename)]" << endl;
   cout << "    [--ls 0(N_LS)/1(B_LS)/2(L_LS)]" << endl;
   cout << "    [--frfactor DOUBLE]" << endl;
   cout << "    [--save UINT(seconds) and optionally STRING(filename) and optionally STRING(scriptfile)]" << endl;
   cout << "    [--maxtime UINT(milliseconds)]" << endl;
   cout << "(2) dep resume SAVING_FILE" << endl;
   cout << "-------------" << endl;
}


void readArguments(int argc, char** argv) {
   //check if enough arguments
   if (argc<3) {
      usage();
      exit(EXIT_FAILURE);
   }
   //check if it's a normal or resume execution
   if (strcmp(argv[1],"resume")==0) {
      //if resume only set the flag to true and read saveFile
      resume = true;
      strcpy(saveFile,argv[2]);
   } else {
      //if normal, set flag to false, and do the rest
      resume = false;
      //read instance and save out filename and exe filename
      strcpy(exe,argv[0]);
      readInstance(argv[1]);
      strcpy(out,argv[2]);
      //init heu to 0
      heu = 0;
      //init default maxnfes basing on nxm
      defaultMaxnfes();
      //read and set the optional arguments
      int i = 3;
      while (i<argc) {
         if (strcmp(argv[i],"--seed")==0) {
            sscanf(argv[i+1],"%u",&seed);
            i += 2;
         } else if (strcmp(argv[i],"--maxnfes")==0) {
            maxnfes = atoi(argv[i+1]);
            i += 2;
         } else if (strcmp(argv[i],"--np")==0) {
            np = atoi(argv[i+1]);
            i += 2;
         } else if (strcmp(argv[i],"--finit")==0) {
            finit = atof(argv[i+1]);
            i += 2;
         } else if (strcmp(argv[i],"--alpha")==0) {
            alpha = atof(argv[i+1]);
            i += 2;
         } else if (strcmp(argv[i],"--heu")==0) {
            readHeuristic(argv[i+1]);
            heu = 1; //otherwise it is 0
            i += 2;
         } else if (strcmp(argv[i],"--ls")==0) {
            ls = atoi(argv[i+1]);
            i += 2;
         } else if (strcmp(argv[i],"--frfactor")==0) {
            frfactor = atof(argv[i+1]);
            i += 2;
         } else if (strcmp(argv[i],"--save")==0) {
            sscanf(argv[i+1],"%u",&saveSeconds);
            if (i+2<argc && argv[i+2][0]!='-') {
               strcpy(saveFile,argv[i+2]);
               if (i+3<argc && argv[i+3][0]!='-') {
                  strcpy(scriptFile,argv[i+3]);
                  i += 4;
               } else {
                  strcpy(scriptFile,"none");
                  i += 3;
               }
            } else {
               //default fileSave is set later if the case (have to known the seed)
               saveFile[0] = '\0';
               i += 2;
            }
         } else if (strcmp(argv[i],"--maxtime")==0) {
            sscanf(argv[i+1],"%u",&maxTime);
            i += 2;
         }
      }
   }
   //done
}


void writeResults() {
   //check if out file exists and if not write the csv header
   FILE* f = fopen(out,"r");
   if (!f) {
      f = fopen(out,"w");
      fprintf(f,"exe,instance,n,m,maxnfes,seed,np,finit,alpha,heu,ls,frfactor"); //input
      fprintf(f,",fgbest,nfesFoundAt,stageFoundAt,nfes,ngen,nrestarts,nforcedrestarts,improvingSteps,lsImprovingSteps,execTime,minStageLength,maxStageLength,avgStageLength,improvingStages,nls,nfesls,nImprovingls,totImprovingls,gbestls,gbest\n"); //output
   }
   fclose(f);
   //write this csv line in append
   char sgbest[PERMSTR_SIZE];
   char sfinit[32];
   char salpha[32];
   char sfrfactor[32];
   char savgStageLength[32];
   perm2str(gbest,n,sgbest);
   double2str(finit,sfinit);
   double2str(alpha,salpha);
   double2str(frfactor,sfrfactor);
   double2str(avgStageLength,savgStageLength);
   f = fopen(out,"a");
   //input print
   fprintf(f,"%s,%s,%d,%d,%d,%u,%d,%s,%s,%d,%d,%s",
             exe,instance,n,m,maxnfes,seed,np,sfinit,salpha,heu,ls,sfrfactor);
   //output print
   fprintf(f,",%d,%d,%d,%d,%d,%d,%d,%d,%d,%lu,%d,%d,%s,%d,%d,%d,%d,%d,%d,%s\n",
           fgbest,nfesFoundAt,stageFoundAt,nfes,ngen,nrestarts,nforcedrestarts,
           improvingSteps,lsImprovingSteps,execTime,
           minStageLength,maxStageLength,savgStageLength,improvingStages,
           nls,nfesls,nImprovingls,totImprovingls,gbestls?1:0,sgbest);
   fclose(f);
   //done
}


void defaultSaveFile() {
   //concatenate instance and seed in saveFile
   sprintf(saveFile,"%s_%u",instance,seed);
   //replace '/' and '.' with '_'
   for (unsigned int i=0; i<strlen(saveFile); i++)
      if (saveFile[i]=='/' || saveFile[i]=='.')
         saveFile[i] = '_';
   //cat .sav to saveFile
   strcat(saveFile,".sav");
   //done
}


void saveHandler(int signum) {
   //signum should be always SIGALRM, thus it is ignored
   //set the dep flag for saving (it is in dep.h/cpp)
   haveToSave = true;
   //set again the alarm (no need to install again the signal handler)
   alarm(saveSeconds);
   //done
}


void preResume(char* filename) {
   //pre-resume to load input variables and instance
   //some variables
   int nowarning; //nowarning is to avoid compiler complaints
   unsigned int j,b;
   unsigned char* bytes; //used to cast double and assume that sizeof(unsigned char)==1
   //open filename and check for errors
   FILE* fsav = fopen(filename,"r");
   if (!fsav) { //if error print and exit
      cerr << "WARNING: I was unable to open " << filename << endl;
      exit(EXIT_FAILURE);
   }
   //read from the file
   //read saveSeconds and saveFile read first
   nowarning = fscanf(fsav,"%u",&saveSeconds);           //saveSeconds 1
   nowarning = fscanf(fsav,"%s",saveFile);               //saveFile 2
   nowarning = fscanf(fsav,"%s",scriptFile);             //scriptFile 2bis
   //read the rng
   readRandState(fsav);                                  //rng state 3
   //values of tft.h/cpp: read "instance", call readInstance, and read maxnfes
   nowarning = fscanf(fsav,"%s",instance);               //instance 4
   readInstance(instance);
   nowarning = fscanf(fsav,"%d",&maxnfes);               //maxnfes 5
   //other values of global.h (except saveSeconds and saveFile read at the begin
   nowarning = fscanf(fsav,"%u",&seed);                  //seed 6
   nowarning = fscanf(fsav,"%s",out);                    //out 7
   nowarning = fscanf(fsav,"%s",exe);                    //exe 8
   //values of dep.h (extern ones and only input) without call depAlloc()
   nowarning = fscanf(fsav,"%d",&np);                    //np 9
   bytes = (unsigned char*)&finit;                       //finit 10
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
   //close the file
   fclose(fsav);
   //statement to avoid compiler complaints
   nowarning++;
   //done
}


void defaultMaxnfes() {
   maxnfes = 1000; //default "stupid" value
   //for taillard instances
   if (n==20) {
      if (m==5)
         maxnfes = 182224100;
      else if (m==10)
         maxnfes = 224784800;
      else if (m==20)
         maxnfes = 256896400;
   } else if (n==50) {
      if (m==5)
         maxnfes = 220712150;
      else if (m==10)
         maxnfes = 256208100;
      else if (m==20)
         maxnfes = 275954150;
   } else if (n==100) {
      if (m==5)
         maxnfes = 235879800;
      else if (m==10)
         maxnfes = 266211000;
      else if (m==20)
         maxnfes = 283040000;
   } else if (n==200) {
      if (m==10)
         maxnfes = 272515500;
      else if (m==20)
         maxnfes = 287728850;
   } else if (n==500) {
      if (m==20)
         maxnfes = 260316750;
   }
}

