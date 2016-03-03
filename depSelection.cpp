//SELECTION PROCEDURES

//ALPHA SELECTION BEGIN
//some declarations
inline void selection_alpha_1o1(int i);

//main selection routine
void selection_alpha() {
	//compute sameFitness and perform np 1o1 selections
	sameFitness = true;
	for (int i=0; i<np; i++) {
		selection_alpha_1o1(i);
#ifdef FIT_INT
		sameFitness &= !(fx[i]-fx[0]); //it's the same of "sameFitness = sameFitness && fx[i]==fx[0]"
#else
		sameFitness = sameFitness && fx[i]==fx[0];
#endif
	}
	//done
}

//relative deviation of y wrt x (positive for minimization / negative for maximization)
#define REL_DEV(y,x) (((y)-(x))/(double)(x))

//selection of individual i
void selection_alpha_1o1(int i) {
	//variables
	int* t;
	//crisp pre-selection between y1[i]/y2[i] (exchange pointers without copy)
#ifdef MINIMIZATION
	if (fy2[i]<fy1[i]) { //tie favors y1[i]
#else
	if (fy2[i]>fy1[i]) { //tie favors y1[i]
#endif
		fy1[i] = fy2[i];
		t = y1[i];
		y1[i] = y2[i];
		y2[i] = t;
#ifdef MYDEBUG
		if (!permValid(y1[i],n)) {
			cout<<"selection1 y1["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(y2[i],n)) {
			cout<<"selection1 y2["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//alpha-selection between y1[i]/x[i] (exchange pointers without copy)
	//the manifest constants TFT/MAKESPAN/LOP/LOPCC are mutually exclusive
#if defined(TFT) || defined(LOPCC) //minimization
	if ( fy1[i]<fx[i] || urand()<alpha-REL_DEV(fy1[i],fx[i]) ) { //tie favors, in some way, x[i]
#endif
#ifdef MAKESPAN //minimization
	if ( fy1[i]<=fx[i] || urand()<alpha-REL_DEV(fy1[i],fx[i]) ) { //tie favors, in some way, y1[i]
#endif
#if defined(LOP) //maximization
	if ( fy1[i]>fx[i] || urand()<alpha+REL_DEV(fy1[i],fx[i]) ) { //tie favors, in some way, x[i]
#endif
		fx[i] = fy1[i];
		t = y1[i];
		y1[i] = x[i];
		x[i] = t; //t is now the new individual
		sfx[i] = sfy[i];  //jde rule
		crx[i] = cry[i];  //jde rule
		//since there was replacement, update ix[i] by computing the inverse
		int* inv = ix[i];
		for (int k=0; k<n; k++)
			inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
		//update cr,f,child statistics (for running var see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
		if (fy1[i]==fy2[i]) //NOT PRECISE BUT OK!!!
			child2succ++;
		else
			child1succ++;
		int oldavg = sfSuccAvg;
		sfSuccAvg = (sfSuccAvg*(child1succ+child2succ-1)+sfx[i])/(child1succ+child2succ);
		sfSuccVar = ((child1succ+child2succ-1)*sfSuccVar+(sfx[i]-oldavg)*(sfx[i]-sfSuccAvg))/(child1succ+child2succ);
		if (sfx[i]>sfSuccMax)
			sfSuccMax = sfx[i];
		if (sfx[i]<sfSuccMin)
			sfSuccMin = sfx[i];
		oldavg = crSuccAvg;
		crSuccAvg = (crSuccAvg*(child1succ+child2succ-1)+crx[i])/(child1succ+child2succ);
		crSuccVar = ((child1succ+child2succ-1)*crSuccVar+(crx[i]-oldavg)*(crx[i]-crSuccAvg))/(child1succ+child2succ);
		if (crx[i]>crSuccMax)
			crSuccMax = crx[i];
		if (crx[i]<crSuccMin)
			crSuccMin = crx[i];
		//some debug
#ifdef MYDEBUG
		if (!permValid(y1[i],n)) {
			cout<<"selection2 y1["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(x[i],n)) {
			cout<<"selection2 x["<<i<<"]"<<endl;
			exit(1);
		}
		if (!permValid(ix[i],n)) {
			cout<<"selection2 ix["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//done
}
//ALPHA SELECTION END



//BEGIN CROWDING SELECTION
inline int footrule_distance(int *p1, int *p2);

void selection_crowding() {
	//some variables
	static int* child_to_replace = new int[np];			//to contain the indexes of the child to put at position i in the population (or -1)	//NO FREE MA PAZIENZA!!!
	static FitnessType* new_fit = new FitnessType[np];	//to contain the fitness of the new population individual in the case of selection		//NO FREE MA PAZIENZA!!!
	int i,j,k,*t,i_closest,min_dist,dist,*inv;
	//initialize child_to_replace and new_fit arrays
	for (i=0; i<np; i++) {
		child_to_replace[i] = -1;		//-1 is an impossible value
		new_fit[i] = fx[i];				//take the old pop fitness
	}
	//scan all the childs
	for (i=0; i<np; i++) {
		//child1 becomes the best child
#ifdef MINIMIZATION
		if (fy2[i]<fy1[i]) { //tie favors y1[i]
#else
		if (fy2[i]>fy1[i]) { //tie favors y1[i]
#endif
			fy1[i] = fy2[i];
			t = y1[i];
			y1[i] = y2[i];
			y2[i] = t;
			
		}
		//finds the most similar population element to child1
		i_closest = 0;
		min_dist = footrule_distance(y1[i],x[0]); //or bubblesort_distance???
		for (j=1; j<np; j++) {
			dist = footrule_distance(y1[i],x[j]); //or bubblesort_distance???
			if (dist<min_dist) {
				i_closest = j;
				min_dist = dist;
			}
		}
		//update child_to_replace and new_fit
#if defined(TFT) || defined(LOPCC)	//minimization
		if (fy1[i]<new_fit[i_closest]) { //favors, in case of ties, the old population individual
#elif defined(MAKESPAN)				//minimization
		if (fy1[i]<=new_fit[i_closest]) { //favors, in case of ties, the selected child
#elif defined(LOP)					//maximization
		if (fy1[i]>new_fit[i_closest]) { //favors, in case of ties, the old population individual
#endif
			child_to_replace[i_closest] = i;
			new_fit[i] = fy1[i];
		}
		//goto next child
	}
	//copies the selected children on the population
	for (i=0; i<np; i++) {
		int j = child_to_replace[i];
		if (j!=-1) { //if -1 no replacement
			fx[i] = fy1[j];
			t = y1[j];
			y1[j] = x[i];
			x[i] = t; //t is now the new individual
			sfx[i] = sfy[j];  //jde rule
			crx[i] = cry[j];  //jde rule
			//since there was replacement, update ix[i] by computing the inverse
			inv = ix[i];
			for (k=0; k<n; k++)
				inv[*t++] = k; //remember t is x[i] ... (inv[*t++] is the same of inv[t[k]])
			//update cr,f,child statistics (for running var see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance)
			if (fy1[j]==fy2[j]) //NOT PRECISE BUT OK!!!
				child2succ++;
			else
				child1succ++;
			int oldavg = sfSuccAvg;
			sfSuccAvg = (sfSuccAvg*(child1succ+child2succ-1)+sfx[i])/(child1succ+child2succ);
			sfSuccVar = ((child1succ+child2succ-1)*sfSuccVar+(sfx[i]-oldavg)*(sfx[i]-sfSuccAvg))/(child1succ+child2succ);
			if (sfx[i]>sfSuccMax)
				sfSuccMax = sfx[i];
			if (sfx[i]<sfSuccMin)
				sfSuccMin = sfx[i];
			oldavg = crSuccAvg;
			crSuccAvg = (crSuccAvg*(child1succ+child2succ-1)+crx[i])/(child1succ+child2succ);
			crSuccVar = ((child1succ+child2succ-1)*crSuccVar+(crx[i]-oldavg)*(crx[i]-crSuccAvg))/(child1succ+child2succ);
			if (crx[i]>crSuccMax)
				crSuccMax = crx[i];
			if (crx[i]<crSuccMin)
				crSuccMin = crx[i];
		}
	}
	//compute sameFitness
	sameFitness = true;
	for (i=1; i<np; i++) {
		if (fx[i]!=fx[0]) {
			sameFitness = false;
			break;
		}
	}
	//done
}

int footrule_distance(int *p1, int *p2) {
	static int* p3 = tmpint; //size n, I can reuse tmpint memory
	int i,dist=0;
	for (i=0; i<n; i++)
		p3[p1[i]] = i;
	for (i=0; i<n; i++)
		dist += abs(i-p3[p2[i]]);
	return dist;
}
//END CROWDING SELECTION



//BEGIN ENTROPY SELECTION #1 (INFORMATION + FITNESS SEPARATELY)
inline void updatePrecCounters(int* x, int* c) {
	int i,j;
	for (i=0; i<n-1; i++)
		for (j=i+1; j<n; j++)
			c[x[i]*n+x[j]]++;				//increase c[ x[i], x[j] ]
}


inline double absInf(int* x, int* c) {
	//compute SUM_prec(x){log(c[prec(x)])}
	int i,j;
	double s = 0.;
	for (i=0; i<n-1; i++)
		for (j=i+1; j<n; j++)
			if (inftype==1 || inftype==3)
				s += loglookup(c[x[i]*n+x[j]]+1);	//accumulate log(c[x[i],x[j]] + 1)  the +1 is because I have to consider also him (and to avoid log0)
			else //inftype==2 || inftype==4
				s += h[x[i]][x[j]]*loglookup(c[x[i]*n+x[j]]+1);	//MAXIMIZATION!!!
	return s;
}


void selection_informationBased() {
	//print an error if minimization
#ifndef LOP
	cerr << "CURRENTLY THIS SELECTION SCHEME WORKS ONLY FOR LOP!!!"
	exit(EXIT_FAILURE);
#endif
	//some variables
	static int totPrec = n*(n-1)/2;			//total number of precedences in a permutation (n choose 2)
	static int np3 = 3*np;					//size of the candidate set
	static int*** pset = new int**[np3];	//USARE TMPINT!!!
	static int* fset = fs;					//fx-fy1-fy2 into fset (directly, since they are stored consecutively in fs)
	static int* c = new int[n*n];			//USARE TMPINT!!!
	static int* ch = new int[np3];			//USARE TMPINT!!!																		//PER STATISTICHE!!!
	int i,k,ifit,tf,*tp;
	double logFactor,fit,inf,fitMax,td;
	//(1) INITIALIZE VARIABLES
	//put x-y1-y2 population pointers into pset
	for (i=0; i<np; i++) {
		pset[i] = &x[i];
		pset[i+np] = &y1[i];
		pset[i+2*np] = &y2[i];
		ch[i] = 0;																													//PER STATISTICHE!!!
		ch[i+np] = 1;																												//PER STATISTICHE!!!
		ch[i+2*np] = 2;																												//PER STATISTICHE!!!
	}
	//initialize precedence counters to all zeros
	memset(c,0,sizeof(int)*(n*n));
	//copy sfy/cry for y2
	memcpy(sfs+2*np,sfs+np,sizeof(double)*np);
	memcpy(crs+2*np,crs+np,sizeof(double)*np);
	//(2) ITERATION 0
	//put fittest individual into the new population (elitism)
	ifit = 0;
	for (i=1; i<np3; i++)
		if (fset[i]>fset[ifit])				//MAXIMIZATION!!!
			ifit = i;
	if (ifit!=0) {							//update only if ifit is not already the current individual to add
		//genotype
		tp = *(pset[0]);
		*(pset[0]) = *(pset[ifit]);
		*(pset[ifit]) = tp;
		//true fitness
		tf = fset[0];
		fset[0] = fset[ifit];
		fset[ifit] = tf;
		//scale factor
		td = sfs[0];
		sfs[0] = sfs[ifit];
		sfs[ifit] = td;
		//crossover probability
		td = crs[0];
		crs[0] = crs[ifit];
		crs[ifit] = td;
		//child types																												//PER STATISTICHE!!!
		tf = ch[0];																													//PER STATISTICHE!!!
		ch[0] = ch[ifit];																											//PER STATISTICHE!!!
		ch[ifit] = tf;																												//PER STATISTICHE!!!
		//selected individual correctly inserted in the new population
	}
	//update precedence counters using the precedences of the inserted individual
	updatePrecCounters(*(pset[0]),c);
	//init sameFitness
	sameFitness = true;
	//(3) ITERATIONS FROM 1 TO NP-1
	//fill the new population till its size is np
	for (k=1; k<np; k++) {
		//useful later to compute information
		logFactor = totPrec*loglookup(k+1);	//the +1 is because I have to consider also him (and to avoid log0) (see above)
		//initialize maximum "modified fitness"
		fitMax = -1./0.; //-inf				//MAXIMIZATION!!!
		//scan the remaining individuals to compute their "modified fitness" using their information wrt to counters and find the fittest
		for (i=k; i<np3; i++) {
			//compute the information of the individual i	 =   (n choose 2)*log(k+1) - SUM_prec(x){log(c[prec(x)]+1)}		the +1 is because I have to consider also him (and to avoid log0) (see above)
			//inf = logFactor - absInf(pset[i],c);
			inf = (logFactor - absInf(*(pset[i]),c))/logFactor; //now it is in [0,1]
			//compute modified fitness of the individual i   =   fset[i]^ALPHA * inf^BETA	(ALPHA and BETA are temporarily set to 1)
			if (inftype==1)
				fit = fset[i] * inf;		//MAXIMIZATION!!!
			else if (inftype==2)
				fit = inf;
			else { //inftype is 3 or 4
				if (k<np/2) {
					fit = fset[i];
					for (int j=0; j<k; j++) {
						if (fit==fset[j]) {
							fit = -1;		//MAXIMIZATION
							break;
						}
					}
				} else
					fit = inftype==3 ? fset[i]*inf : inf;	//MAXIMIZATION
			}
			//update fitMax and ifit
			if (fit>fitMax) {				//MAXIMIZATION!!!
				fitMax = fit;
				ifit = i;
			}
			//goto next candidate
		}
		//put fittest individual (using modified fitness) into the new population
		if (ifit!=k) {						//update only if ifit is not already the current individual to add
			//genotype
			tp = *(pset[k]);
			*(pset[k]) = *(pset[ifit]);
			*(pset[ifit]) = tp;
			//true fitness
			tf = fset[k];
			fset[k] = fset[ifit];
			fset[ifit] = tf;
			//scale factor
			td = sfs[k];
			sfs[k] = sfs[ifit];
			sfs[ifit] = td;
			//crossover probability
			td = crs[k];
			crs[k] = crs[ifit];
			crs[ifit] = td;
			//child types																											//PER STATISTICHE!!!
			tf = ch[k];																												//PER STATISTICHE!!!
			ch[k] = ch[ifit];																										//PER STATISTICHE!!!
			ch[ifit] = tf;																											//PER STATISTICHE!!!
			//selected individual correctly inserted in the new population
		}
		//update sameFitness
		if (fx[k]!=fx[0])
			sameFitness = false;
		//update precedence counters using the precedences of the inserted individual
		updatePrecCounters(*(pset[k]),c);
		//goto next choice
	}
#ifdef MYDEBUG
	for (i=0; i<np; i++) {
		if (*(pset[i])!=x[i]) {
			cout << "PROBLEMA IN INF SELECTION 1!!!" << endl;
			exit(1);
		}
		if (fx[i]!=fset[i]) {
			cout << "PROBLEMA IN INF SELECTION 2!!!" << endl;
			exit(1);
		}
	}
#endif
	//(6) INVERSE COMPUTATION
	for (i=0; i<np; i++)
		for (k=0; k<n; k++)
			ix[i][x[i][k]] = k;
	//(7) UPDATE STATISTICS																											//PER STATISTICHE!!!
	for (i=0; i<np; i++) {																											//PER STATISTICHE!!!
		switch (ch[i]) {																											//PER STATISTICHE!!!
			case 1:																													//PER STATISTICHE!!!
				child1succ++;																										//PER STATISTICHE!!!
				break;																												//PER STATISTICHE!!!
			case 2:																													//PER STATISTICHE!!!
				child2succ++;																										//PER STATISTICHE!!!
				break;																												//PER STATISTICHE!!!
			default:																												//PER STATISTICHE!!!
				break;																												//PER STATISTICHE!!!
		}																															//PER STATISTICHE!!!
		if (ch[i]!=0) {																												//PER STATISTICHE!!!
			int oldavg = sfSuccAvg;																									//PER STATISTICHE!!!
			sfSuccAvg = (sfSuccAvg*(child1succ+child2succ-1)+sfx[i])/(child1succ+child2succ);										//PER STATISTICHE!!!
			sfSuccVar = ((child1succ+child2succ-1)*sfSuccVar+(sfx[i]-oldavg)*(sfx[i]-sfSuccAvg))/(child1succ+child2succ);			//PER STATISTICHE!!!
			if (sfx[i]>sfSuccMax)																									//PER STATISTICHE!!!
				sfSuccMax = sfx[i];																									//PER STATISTICHE!!!
			if (sfx[i]<sfSuccMin)																									//PER STATISTICHE!!!
				sfSuccMin = sfx[i];																									//PER STATISTICHE!!!
			oldavg = crSuccAvg;																										//PER STATISTICHE!!!
			crSuccAvg = (crSuccAvg*(child1succ+child2succ-1)+crx[i])/(child1succ+child2succ);										//PER STATISTICHE!!!
			crSuccVar = ((child1succ+child2succ-1)*crSuccVar+(crx[i]-oldavg)*(crx[i]-crSuccAvg))/(child1succ+child2succ);			//PER STATISTICHE!!!
			if (crx[i]>crSuccMax)																									//PER STATISTICHE!!!
				crSuccMax = crx[i];																									//PER STATISTICHE!!!
			if (crx[i]<crSuccMin)																									//PER STATISTICHE!!!
				crSuccMin = crx[i];																									//PER STATISTICHE!!!
		}																															//PER STATISTICHE!!!
	}																																//PER STATISTICHE!!!
	//done
}
//END ENTROPY SELECTION #1 (INFORMATION + FITNESS SEPARATELY)




//BEGIN FITNESS SEPARATION SELECTION
#define MAX_GAP 10000
#define MAX_INTERVAL 10000

void selection_fitnessSeparation() {
	//print an error if minimization
	cerr << "CURRENTLY THIS SELECTION SCHEME DOES NOT WORK!!!" << endl;		//TEMPORARY CODE!!!
	exit(EXIT_FAILURE);														//TEMPORARY CODE!!!
#ifdef MINIMIZATION
	cerr << "CURRENTLY THIS SELECTION SCHEME WORKS ONLY FOR MAXIMIZATION PROBLEMS!!!" << endl;
	exit(EXIT_FAILURE);
#endif
	//some variables
	static int np3 = 3*np;					//size of the candidate set
	static int*** pset = new int**[np3];	//USARE TMPINT!!!
	static int* gap = new int[np3];			//USARE TMPINT!!!
	static int* inter = new int[np3];		//USARE TMPINT!!!
	static int* ppos = new int[np3];		//USARE TMPINT!!!
	static int* when = new int[np3];		//USARE TMPINT!!!
	static int* fset = fs;					//fx-fy1-fy2 into fset (directly, since they are stored consecutively in fs)
	static int* ch = new int[np3];			//USARE TMPINT!!!																		//PER STATISTICHE!!!
	static bool first = true;
	int i,j,*tp,tf,f,g;
	double td;
	//copy sfy/cry for y2
	memcpy(sfs+2*np,sfs+np,sizeof(double)*np);
	memcpy(crs+2*np,crs+np,sizeof(double)*np);
	//check first execution and initialize
	if (first) {
		for (i=0; i<np; i++) {
			gap[i] = MAX_GAP;
			inter[i] = MAX_INTERVAL;
			when[i] = ngen + inter[i];
			ppos[i] = np;
		}
		first = false;
	}
	//set selection variables for the childs
	for (i=np; i<np3; i++) {
		gap[i] = MAX_GAP;
		inter[i] = MAX_INTERVAL;
		when[i] = ngen + inter[i];
		ppos[i] = np;
	}
	//put x-y1-y2 population pointers into pset
	for (i=0; i<np; i++) {
		pset[i] = &x[i];
		pset[i+np] = &y1[i];
		pset[i+2*np] = &y2[i];
		ch[i] = 0;																													//PER STATISTICHE!!!
		ch[i+np] = 1;																												//PER STATISTICHE!!!
		ch[i+2*np] = 2;	
	}
	//sort pset/fset (stable sort using insertion sort)
	for (i=1; i<np3; i++) {
		j = i;
		while (j>0 && fset[j-1]<fset[j]) {
			//begin swap
			//genotype
			tp = *(pset[j-1]);
			*(pset[j-1]) = *(pset[j]);
			*(pset[j]) = tp;
			//fitness
			tf = fset[j-1];
			fset[j-1] = fset[j];
			fset[j] = tf;
			//scale factor
			td = sfs[j-1];
			sfs[j-1] = sfs[j];
			sfs[j] = td;
			//crossover probability
			td = crs[j-1];
			crs[j-1] = crs[j];
			crs[j] = td;
			//gap
			tf = gap[j-1];
			gap[j-1] = gap[j];
			gap[j] = tf;
			//inter
			tf = inter[j-1];
			inter[j-1] = inter[j];
			inter[j] = tf;
			//when
			tf = when[j-1];
			when[j-1] = when[j];
			when[j] = tf;
			//ppos
			tf = ppos[j-1];
			ppos[j-1] = ppos[j];
			ppos[j] = tf;
			//child types																											//PER STATISTICHE!!!
			tf = ch[j-1];																											//PER STATISTICHE!!!
			ch[j-1] = ch[j];																										//PER STATISTICHE!!!
			ch[j] = tf;																												//PER STATISTICHE!!!
			//end swap
			j--;			
		}
	}
cout<<"fset = "; printPerm(fset,np3);//vale
cout<<"gap  = "; printPerm(gap,np3);//vale
cout<<"inter= "; printPerm(inter,np3);//vale
cout<<"when = "; printPerm(when,np3);//vale
cout<<"ppos = "; printPerm(ppos,np3);//vale
	//now select
	f = fset[0];	//first is taken
	g = gap[0];		//first is taken
	j = 1;
	for (i=1; i<np && j<np3; i++) {	//first is taken
		if (fset[j]<=f-g) {			//take the j-th individual
			if (i!=j) {
				//begin swap
				//genotype
				tp = *(pset[i]);
				*(pset[i]) = *(pset[j]);
				*(pset[j]) = tp;
				//fitness
				tf = fset[i];
				fset[i] = fset[j];
				fset[j] = tf;
				//scale factor
				td = sfs[i];
				sfs[i] = sfs[j];
				sfs[j] = td;
				//crossover probability
				td = crs[i];
				crs[i] = crs[j];
				crs[j] = td;
				//gap
				tf = gap[i];
				gap[i] = gap[j];
				gap[j] = tf;
				//inter
				tf = inter[i];
				inter[i] = inter[j];
				inter[j] = tf;
				//when
				tf = when[i];
				when[i] = when[j];
				when[j] = tf;
				//ppos
				tf = ppos[i];
				ppos[i] = ppos[j];
				ppos[j] = tf;
				//child types																											//PER STATISTICHE!!!
				tf = ch[i];																												//PER STATISTICHE!!!
				ch[i] = ch[j];																											//PER STATISTICHE!!!
				ch[j] = tf;																												//PER STATISTICHE!!!
				//end swap
			}
			f = fset[i];
			if (i!=ppos[i]) {
				gap[i] = MAX_GAP;
				inter[i] = MAX_INTERVAL;
				when[i] = ngen+inter[i];
				ppos[i] = i;
			} else if (when[i]==ngen) {
				gap[i] /= 2;
				inter[i] /= 2;
				when[i] = ngen+inter[i];
			}
			g = gap[i];
		}
		//move to next candidate
		j++;
	}
	//fill the remaining
	for (j=i; j<np; j++) {
		prand(n,x[j]);
		fset[j] = INT_MIN;
		sfs[j] = finit;
		crs[j] = crinit;
		gap[j] = MAX_GAP;
		inter[j] = MAX_INTERVAL;
		when[j] = ngen + inter[j];
		ppos[j] = np;
	}
	//inverse computation
	for (i=0; i<np; i++)
		for (j=0; j<n; j++)
			ix[i][x[i][j]] = j;
	//update statistics																												//PER STATISTICHE!!!
	for (i=0; i<np; i++) {																											//PER STATISTICHE!!!
		switch (ch[i]) {																											//PER STATISTICHE!!!
			case 1:																													//PER STATISTICHE!!!
				child1succ++;																										//PER STATISTICHE!!!
				break;																												//PER STATISTICHE!!!
			case 2:																													//PER STATISTICHE!!!
				child2succ++;																										//PER STATISTICHE!!!
				break;																												//PER STATISTICHE!!!
			default:																												//PER STATISTICHE!!!
				break;																												//PER STATISTICHE!!!
		}																															//PER STATISTICHE!!!
		if (ch[i]!=0) {																												//PER STATISTICHE!!!
			int oldavg = sfSuccAvg;																									//PER STATISTICHE!!!
			sfSuccAvg = (sfSuccAvg*(child1succ+child2succ-1)+sfx[i])/(child1succ+child2succ);										//PER STATISTICHE!!!
			sfSuccVar = ((child1succ+child2succ-1)*sfSuccVar+(sfx[i]-oldavg)*(sfx[i]-sfSuccAvg))/(child1succ+child2succ);			//PER STATISTICHE!!!
			if (sfx[i]>sfSuccMax)																									//PER STATISTICHE!!!
				sfSuccMax = sfx[i];																									//PER STATISTICHE!!!
			if (sfx[i]<sfSuccMin)																									//PER STATISTICHE!!!
				sfSuccMin = sfx[i];																									//PER STATISTICHE!!!
			oldavg = crSuccAvg;																										//PER STATISTICHE!!!
			crSuccAvg = (crSuccAvg*(child1succ+child2succ-1)+crx[i])/(child1succ+child2succ);										//PER STATISTICHE!!!
			crSuccVar = ((child1succ+child2succ-1)*crSuccVar+(crx[i]-oldavg)*(crx[i]-crSuccAvg))/(child1succ+child2succ);			//PER STATISTICHE!!!
			if (crx[i]>crSuccMax)																									//PER STATISTICHE!!!
				crSuccMax = crx[i];																									//PER STATISTICHE!!!
			if (crx[i]<crSuccMin)																									//PER STATISTICHE!!!
				crSuccMin = crx[i];																									//PER STATISTICHE!!!
		}																															//PER STATISTICHE!!!
	}																																//PER STATISTICHE!!!
	//same fitness
	sameFitness = true;
	for (i=1; i<np; i++)
		if (fx[0]!=fx[i]) {
			sameFitness = false;
			break;
		}
cout<<"fset = "; printPerm(fset,np3);//vale
cout<<"gap  = "; printPerm(gap,np3);//vale
cout<<"inter= "; printPerm(inter,np3);//vale
cout<<"when = "; printPerm(when,np3);//vale
cout<<"ppos = "; printPerm(ppos,np3);//vale
exit(1);//vale
	//done
}
//END FITNESS SEPARATION SELECTION






























