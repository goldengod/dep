//POPULATION INITIALIZATION PROCEDURES

//ALL RANDOM EXCEPT (OPTIONALLY) ONE BUILD USING AN HEURISTIC
void popInit_randheu() {
	//init nfes and fgbest
	nfes = 0;
#if defined(MINIMIZATION) && defined(FIT_INT)
	fgbest = INT_MAX;
#elif defined(MINIMIZATION) && defined(FIT_REAL)
	fgbest = PLUS_INF;
#elif defined(MAXIMIZATION) && defined(FIT_INT)
	fgbest = INT_MIN;
#elif defined(MAXIMIZATION) && defined(FIT_REAL)
	fgbest = MINUS_INF;
#endif
	//random population initialization apart one individual from heuristic
	for (int i=0; i<np; i++) {
		if (i==0 && heu==1)
			memcpy(x[0],heup,permByteSize); //heuristic permutation
		else
			prand(n,x[i]); //random
		fx[i] = eval(x[i]);
		nfes++;
		updateGbest(x[i],fx[i]);
		sfx[i] = finit;
		crx[i] = crinit;
		for (int k=0; k<n; k++) //ix[i] = inverse of x[i]
			ix[i][x[i][k]] = k;
#ifdef MYDEBUG
		if (!permValid(x[i],n)) {
			cout<<"popInit x["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(ix[i],n)) {
			cout<<"popInit ix["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//init ngen, nrestarts, samefitness
	ngen = 1; //initialization was the first generation
	nrestarts = 0;
	nforcedrestarts = 0;
	sameFitness = false;
	//init restart statistics variables
	lastRestart = 0;
	minStageLength = maxnfes; //it's impossible to be more
	maxStageLength = 0; //it's impossible to be less
	avgStageLength = 0.;
#if defined(MINIMIZATION) && defined(FIT_INT)
	fgbestAtStageStart = INT_MAX;	//for sure at stage end it will be better
#elif defined(MINIMIZATION) && defined(FIT_REAL)
	fgbestAtStageStart = PLUS_INF;	//for sure at stage end it will be better
#elif defined(MAXIMIZATION) && defined(FIT_INT)
	fgbestAtStageStart = INT_MIN;	//for sure at stage end it will be better
#elif defined(MAXIMIZATION) && defined(FIT_REAL)
	fgbestAtStageStart = MINUS_INF;	//for sure at stage end it will be better
#endif	
	improvingStages = 0;
	//init local search statistics variables
	nfesls = 0;
	nls = 0;
	nImprovingls = 0;
	totImprovingls = 0;	//ok both for FIT_INT and FIT_REAL
	//init lsmode to false
	lsmode = false;
	//init variables for forced restart
	nfesWhenToForceRestart = forcedRestartPeriod = maxnfes==INT_MAX ? INT_MAX : maxnfes*frfactor;
	//init improving steps
	improvingSteps = lsImprovingSteps = 0;
	//init f, cr, child statistics variables
	sfSuccAvg = 0.;
	sfSuccVar = 0.;
	sfSuccMin = 1.;
	sfSuccMax = 0.;
	crSuccAvg = 0.;
	crSuccVar = 0.;
	crSuccMin = 1.;
	crSuccMax = 0.;
	child1succ = 0;
	child2succ = 0;
	//memetic part
	if (memetic) {
		for (int i=0; i<np; i++)
			localSearch(x[i],fx[i]);
	}
	//done
}

