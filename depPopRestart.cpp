//RESTART PROCEDURES

//this line is needed to make the following function know that "localSearch" exists
extern void (*localSearch)(int*,FitnessType&);

//this is to manage the restart statistics
inline void restartStatistics() {
	//compute restart statistics (call after increasing nrestarts)
	int stageLength = ngen-lastRestart;
	if (stageLength<minStageLength)
		minStageLength = stageLength;
	if (stageLength>maxStageLength)
		maxStageLength = stageLength;
	avgStageLength += (stageLength-avgStageLength)/nrestarts; //moving average (see wikipedia)
	lastRestart = ngen;
#ifdef MINIMIZATION
	if (fgbest<fgbestAtStageStart)
#else
	if (fgbest>fgbestAtStageStart)
#endif
		improvingStages++;
	fgbestAtStageStart = fgbest;
	//done
}

//BEGIN RANDLS RESTART
//restart population if restart has been triggered
bool popRestart_randls() {
	//if all fitnesses are not the same do nothing
	if (!sameFitness)
		return false; //false since maxnfes was not exhausted
	//the first (it's a best) individual is keeped and the other are randomized
	for (int i=1; i<np; i++) {
		prand(n,x[i]);
#if defined(MINIMIZATION) && defined(FIT_INT)
		fx[i] = INT_MAX;
#elif defined(MINIMIZATION) && defined(FIT_REAL)
		fx[i] = PLUS_INF;
#elif defined(MAXIMIZATION) && defined(FIT_INT)
		fx[i] = INT_MIN;
#elif defined(MAXIMIZATION) && defined(FIT_REAL)
		fx[i] = MINUS_INF;
#endif
		sfx[i] = finit;
		crx[i] = crinit;
		//since there was replacement, update ix[i] by computing inverse
		int* inv = ix[i];
		int* t = x[i];
		for (int k=0; k<n; k++)
			inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
#ifdef MYDEBUG
		if (!permValid(x[i],n)) {
			cout<<"restart x["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(ix[i],n)) {
			cout<<"restart ix["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//anyway reset sfx[0] and crx[0]
	sfx[0] = finit;
	crx[0] = crinit;
	//local search on x[0]
	if (ls==B_LS) { //baldwin
		memcpy(tmpint,x[0],permByteSize);
		FitnessType ft = fx[0];
		localSearch(tmpint,ft);
	} else if (ls==L_LS) {//lamarck
		localSearch(x[0],fx[0]);
		for (int k=0; k<n; k++)
			ix[0][x[0][k]] = k;
	}
	//update nfesWhenToForceRestart, thus if normal restart don't do forcedrestart
	if (nfesWhenToForceRestart<INT_MAX)
		nfesWhenToForceRestart = nfes + forcedRestartPeriod;
	//increase restarts counter
	nrestarts++;
	//restart statistics
	restartStatistics();
	//return true/false if budget of nfes was exhausted in the local search
	return termination();
	//done
}

bool popForcedRestart_randls() {
	//check if restart has to be forced, return if not the case
	if (nfes<nfesWhenToForceRestart)
		return false; //false since maxnfes was not exhausted
	//find best individual
	int ibest = 0;
	int i;
	for (i=1; i<np; i++)
		if (fx[i]<fx[0])
			ibest = i;
	//restart population except ibest
	for (i=0; i<np; i++) {
		if (i!=ibest) {
			prand(n,x[i]);
#if defined(MINIMIZATION) && defined(FIT_INT)
			fx[i] = INT_MAX;
#elif defined(MINIMIZATION) && defined(FIT_REAL)
			fx[i] = PLUS_INF;
#elif defined(MAXIMIZATION) && defined(FIT_INT)
			fx[i] = INT_MIN;
#elif defined(MAXIMIZATION) && defined(FIT_REAL)
			fx[i] = MINUS_INF;
#endif
			//update ix
			int* inv = ix[i];
			int* t = x[i];
			for (int k=0; k<n; k++)
				inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
#ifdef MYDEBUG
			if (!permValid(x[i],n)) {
				cout<<"forcedRestart x["<<i<<"]"<<endl;
				exit(1);
			}
			if (!permValid(ix[i],n)) {
				cout<<"forcedRestart ix["<<i<<"]"<<endl;
				exit(1);
			}
#endif
		}
		sfx[i] = finit; //finit also for the best
		crx[i] = crinit; //crinit also for the best
	}
	//lamarckian local search on x[ibest]
	localSearch(x[ibest],fx[ibest]);
	for (int k=0; k<n; k++)
		ix[ibest][x[ibest][k]] = k;
	//update nfesWhenToForceRestart
	if (nfesWhenToForceRestart<INT_MAX)
		nfesWhenToForceRestart = nfes + forcedRestartPeriod;
	//update statistics
	nrestarts++;
	nforcedrestarts++;
	restartStatistics();
	//return true/false if budget of nfes was exhausted in the local search
	return termination();
	//done
}
//END RANDLS RESTART



//BEGIN SHAKE+HALFRAND+LS RESTART (LARMKIAN OR BALDWINIAN LOCAL SEARCH DO NOT APPLY HERE, ONLY N_LS AND ONE OF THE OTHERS)
inline void shake(int* x, int* y, int ns) {
   int i,j,t;
   if (x!=y)
      memcpy(x,y,sizeof(int)*n);
   for (i=0; i<ns; i++) {
      j = irand(n-1);
      t = x[j];
      x[j] = x[j+1];
      x[j+1] = t;
   }
}

inline void restart_shrandls(int lstype) {
	//local search on the global best so far and put it in x[0]
	int i;
	memcpy(x[0],gbest,sizeof(int)*n);
	if (lstype!=N_LS) //LAMARKIAN OR BALDWINIAN DOES NOT CHANGE NOTHING HERE!!!
		localSearch(x[0],fx[0]);
	//first half takes a shake of the improved global best
	for (i=0; i<np/2; i++) {
		shake(x[i],x[0],irand(n*(n-1)/2));
#if defined(MINIMIZATION) && defined(FIT_INT)
		fx[i] = INT_MAX;
#elif defined(MINIMIZATION) && defined(FIT_REAL)
		fx[i] = PLUS_INF;
#elif defined(MAXIMIZATION) && defined(FIT_INT)
		fx[i] = INT_MIN;
#elif defined(MAXIMIZATION) && defined(FIT_REAL)
		fx[i] = MINUS_INF;
#endif
		sfx[i] = finit;
		crx[i] = crinit;
	}
	//second half is completely randomized
	for (; i<np; i++) {
		prand(n,x[i]);
#if defined(MINIMIZATION) && defined(FIT_INT)
		fx[i] = INT_MAX;
#elif defined(MINIMIZATION) && defined(FIT_REAL)
		fx[i] = PLUS_INF;
#elif defined(MAXIMIZATION) && defined(FIT_INT)
		fx[i] = INT_MIN;
#elif defined(MAXIMIZATION) && defined(FIT_REAL)
		fx[i] = MINUS_INF;
#endif
		sfx[i] = finit;
		crx[i] = crinit;
	}
	//since there was replacement, update ix[i] by computing inverse
	for (i=0; i<np; i++) {
		int* inv = ix[i];
		int* t = x[i];
		for (int k=0; k<n; k++)
			inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
#ifdef MYDEBUG
		if (!permValid(x[i],n)) {
			cout<<"restart x["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(ix[i],n)) {
			cout<<"restart ix["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	  }
	  //done
}

bool popRestart_shrandls() {
	//if all fitnesses are not the same do nothing
	if (!sameFitness)
		return false; //false since maxnfes was not exhausted
	//do the restart
	restart_shrandls(ls); //using the parameter for the local search type
	//update nfesWhenToForceRestart, thus if normal restart don't do forcedrestart
	if (nfesWhenToForceRestart<INT_MAX)
		nfesWhenToForceRestart = nfes + forcedRestartPeriod;
	//increase restarts counter
	nrestarts++;
	//restart statistics
	restartStatistics();
	//return true/false if budget of nfes was exhausted in the local search
	return termination();
	//done
}

bool popForcedRestart_shrandls() {
	//check if restart has to be forced, return if not the case
	if (nfes<nfesWhenToForceRestart)
		return false; //false since maxnfes was not exhausted
	//do the restart
	restart_shrandls(L_LS); //always apply local search in case of forced restart (LARMKIAN OR BALDWINIAN DO NOT MIND)
	//update nfesWhenToForceRestart
	if (nfesWhenToForceRestart<INT_MAX)
		nfesWhenToForceRestart = nfes + forcedRestartPeriod;
	//update statistics
	nrestarts++;
	nforcedrestarts++;
	restartStatistics();
	//return true/false if budget of nfes was exhausted in the local search
	return termination();
	//done
}
//END SHAKE+HALFRAND+LS RESTART

