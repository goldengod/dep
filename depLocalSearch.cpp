//LOCAL SEARCH PROCEDURES

//BEGIN VNS4 CODE
#if (defined(TFT) || defined(MAKESPAN)) && defined(GFC)
//Branchless min between x and y, bits has to be sizeof(int)*8-1 (see http://aggregate.org/MAGIC)
#define FAST_MIN(x,y,bits) ((x)+((((y)-(x))>>(bits))&((y)-(x))))
#endif
#if (defined(TFT) || defined(MAKESPAN)) && defined(GFC)
bool intStep(int* x, int& fx, bool first);
#else
bool intStep(int* x, int& fx);
#endif
void insStep(int* x, int& fx);
#if (defined(TFT) || defined(MAKESPAN)) && defined(GFC)
inline void fwdins(int* x, int i, int j);
inline void bwdins(int* x, int i, int j);
inline void bwdins2(int* x, int i, int j);
#endif

//vns4 local search: full interchange ls (1st imp. style) + one step insert ls (best imp. style)
void localSearch_vns4(int* x, int& fx) {
	//if x is in the history, return
	int i;
	for (i=0; i<n; i++)
		if (lsHistory[i]!=x[i])
			break;
	if (i==n)
		return;
	//set lsmode
	lsmode = true;
	//save the fitness of the seed for statistics
	int fseed = fx;
	//perform full interchange + insertion step
	int fprev;
	do {
		fprev = fx;
#if (defined(TFT) || defined(MAKESPAN)) && defined(GFC)
		bool firstInt = true; //used the first time to call computeGFC
		while (intStep(x,fx,firstInt))   //interchange local search (using gfc)
			firstInt = false;
#else
		while (intStep(x,fx));           //interchange local search (without gfc)
#endif
		if (termination()) //break the do-while if termination
			break;
		insStep(x,fx);                   //insertion step
		if (termination()) //break the do-while if termination
			break;
#ifdef MINIMIZATION
	} while (fx<fprev);
#else
	} while (fx>fprev);
#endif
	//save x in the history
	memcpy(lsHistory,x,permByteSize);
	//update local search statistics
	nls++;
#ifdef MINIMIZATION
	if (fx<fseed) {
#else
	if (fx>fseed) {
#endif
		nImprovingls++;
#ifdef MINIMIZATION
		totImprovingls += fseed-fx;
#else
		totImprovingls += fx-fseed;
#endif
	}
	//unset lsmode
	lsmode = false;
	//done
}


#if (defined(TFT) || defined(MAKESPAN)) && defined(GFC)
//E' DI SICURO MINIMIZATION E NON MAXIMIZATION!!!
//VERSIONI DI INTSTEP E INSSTEP CHE USANO GFC!!!!!!!!!!!!!!
bool intStep(int* x, int& fx, bool first) {
	//GFC VERSION
	//one step of interchange local search
	//1st improvement style using a random permutation
	//variables
	const static int bits = sizeof(int)*8-1;
	static int* r = tmpint+n; //since tmpint is used to copy x in the case of B_LS
	static int nm1 = n-1;
	int i,j,ii,jj,t,fseed;
	//compute gfc data of x BUT ONLY IF FIRST application in vns4!!!
	if (first)
		computeGFC(x);
	//set fitness of the seed and get a random permutation
	fseed = fx;
	prand(n,r);
	//scan for interchange moves
	for (ii=0; ii<nm1; ii++) {
		for (jj=ii+1; jj<n; jj++) {
			//get i,j from the random permutation
			i = r[ii];
			j = r[jj];
			//interchange jobs at positions i and j in x
			t = x[i];
			x[i] = x[j];
			x[j] = t; //now t is equal to x[j]
			//evaluate x using gfc, update gbest and check termination
			fx = evalGFC(x,FAST_MIN(i,j,bits));
#ifdef MYDEBUG
			if (fx!=eval(x)) {
				cerr<<"intStep evaluation mismatch"<<endl;
				exit(1);
			}
#endif
			nfes++;
			nfesls++;
			updateGbest(x,fx);
			if (termination())
				return false;
			//if better UPDATE gfc data and return true, x and fx are modified
			if (fx<fseed) {
				evalUpdateGFC(x,FAST_MIN(i,j,bits));
				return true;
			}
			//if here, not better, so reset the interchange above by redoing it
			x[j] = x[i];
			x[i] = t; //remember that t was x[j]
			//done
		}
	}
	//if here, no return above, so fx has not been improved, thus reset fx and return false
	fx = fseed;
	return false;
	//done
}


void insStep(int* x, int& fx) {
	//GFC VERSION
	//one step of insertion local search
	//best improvement style, so no random permutation
	//(i,j) means "move job at pos i to pos j"
	//variables
	int i,j,fbest,bi,bj,ft;
	//no need for computing gfc data since intStep did it!!!
	//computeGFC(x)
	//init fbest to fx of the seed
	fbest = fx;
	//scan for insertion by doing rightmost insertions before
	for (i=n-1; i>=0; i--) {
		for (j=i+1; j<n; j++) { //do (i,j) and (j,i)
			//(i,j), since i<j for sure this is a forward insertion
			fwdins(x,i,j);
			//evaluate x using gfc, update gbest and check termination
			ft = evalGFC(x,i); //since i<j
#ifdef MYDEBUG
			if (ft!=eval(x)) {
				cerr<<"insStep evaluation mismatch"<<endl;
				exit(1);
			}
#endif
			nfes++;
			nfesls++;
			updateGbest(x,ft);
			if (termination())
				return;
			//update fbest,bi,bj if better ft is better
			if (ft<fbest) { //with <= I simulate < with the scan from 0 to n
				fbest = ft;
				bi = i;
				bj = j;
			}
			//before doing (j,i) check if it's case that (j,i)!=(i,j)
			if (i+1!=j) { //!= means <
				//perform two times (j,i): one to undo (i,j) and one to do (j,i)
				bwdins2(x,j,i);
				//evaluate x using gfc, update gbest and check termination
				ft = evalGFC(x,i); //since i<j
				nfes++;
				nfesls++;
				updateGbest(x,ft);
				if (termination())
					return;
				//update fbest,bi,bj if better ft is better
				if (ft<fbest) { //with <= I simulate < with the scan from 0 to n
					fbest = ft;
					bi = j; //inverted since the insertion is (j,i)
					bj = i; //inverted since the insertion is (j,i)
				}
				//undo (j,i) by doing (i,j)
				fwdins(x,i,j);
			} else //in the case that (j,i) was not performed, simply undo (i,j) by doing (j,i)
				bwdins(x,j,i);
			//done both (i,j) and (j,i)
		}
	}
	//if there was an improvement wrt the seed solution, update the solution by doing (bi,bj) insertion
	if (fbest<fx) {
		ins(x,bi,bj);
		fx = fbest;
	}
	//done
}

#else

//VERSIONI DI INTSTEP E INSSTEP CHE ***NON*** USANO GFC!!!!!!!!!!!!!!
bool intStep(int* x, int& fx) {
	//NO-GFC VERSION
	//one step of interchange local search
	//1st improvement style using a random permutation
	//variables
	static int* r = tmpint+n; //since tmpint is used to copy x in the case of B_LS
	static int nm1 = n-1;
	int i,j,ii,jj,t,fseed;
	//set fitness of the seed and get a random permutation
	fseed = fx;
	prand(n,r);
	//scan for interchange moves
	for (ii=0; ii<nm1; ii++) {
		for (jj=ii+1; jj<n; jj++) {
			//get i,j from the random permutation
			i = r[ii];
			j = r[jj];
			//interchange jobs at positions i and j in x
			t = x[i];
			x[i] = x[j];
			x[j] = t; //now t is equal to x[j]
			//evaluate x, update gbest and check termination
			fx = eval(x);
			nfes++;
			nfesls++;
			updateGbest(x,fx);
			if (termination())
				return false;
			//if better return true, x and fx are modified
#ifdef MINIMIZATION
			if (fx<fseed)
#else
			if (fx>fseed)
#endif
				return true;
			//if here, not better, so reset the interchange above by redoing it
			x[j] = x[i];
			x[i] = t; //remember that t was x[j]
			//done
		}
	}
	//if here, no return above, so fx has not been improved, thus reset fx and return false
	fx = fseed;
	return false;
	//done
}


void insStep(int* x, int& fx) {
	//NO-GFC VERSION
	//one step of insertion local search
	//best improvement style, so no random permutation
	//(i,j) means "move job at pos i to pos j"
	//variables
	int i,j,fbest,bi,bj,ft;
	//init fbest to fx of the seed
	fbest = fx;
	//scan for insertion moves in the same order of GFC version
	for (i=n-1; i>=0; i--) {
		for (j=i+1; j<n; j++) {
			ins(x,i,j);
			ft = eval(x);
			nfes++;
			nfesls++;
			updateGbest(x,ft);
			if (termination())
				return;
#ifdef MINIMIZATION
			if (ft<fbest) {
#else
			if (ft>fbest) {
#endif
				fbest = ft;
				bi = i;
				bj = j;
			}
			ins(x,j,i);
			if (i+1!=j) {
				ins(x,j,i);
				ft = eval(x);
				nfes++;
				nfesls++;
				updateGbest(x,ft);
				if (termination())
					return;
#ifdef MINIMIZATION
				if (ft<fbest) {
#else
				if (ft>fbest) {
#endif
					fbest = ft;
					bi = j; //note that j and i are reversed here!!!
					bj = i; //note that j and i are reversed here!!!
				}
				ins(x,i,j);
			}
		}
	}
	/*
	for (i=n-1; i>=0; i--) {
		for (j=n-1; j>=0; j--) {
			if (i!=j && i!=j+1) { //(i,i+1)==(i+1,i) so avoid to redo
				//perform the insertion (i,j) on x
				ins(x,i,j);
				//evaluate x, update gbest and check termination
				ft = eval(x);
				nfes++;
				nfesls++;
				updateGbest(x,ft);
				if (termination())
					return;
				//update fbest,bi,bj if better ft is better
#ifdef MINIMIZATION
				if (ft<fbest) {
#else
				if (ft>fbest) {
#endif
					fbest = ft;
					bi = i;
					bj = j;
				}
				//undo the insertion (i,j) on x by doing the insertion (j,i) on x
				ins(x,j,i);
			}
		}
	}
	*/
	//if there was an improvement wrt the seed solution, update the solution by doing (bi,bj) insertion
#ifdef MINIMIZATION
	if (fbest<fx) {
#else
	if (fbest>fx) {
#endif
		ins(x,bi,bj);
		fx = fbest;
	}
	//done
}
#endif

#if (defined(TFT) || defined(MAKESPAN)) && defined(GFC)
inline void fwdins(int* x, int i, int j) {
	//assume i<j and perform insertion (i,j) on x (this is a forward insertion)
	int t = x[i];
	memmove(x+i,x+i+1,sizeof(int)*(j-i));
	x[j] = t;
}


inline void bwdins(int* x, int i, int j) {
	//assume i>j and perform insertion (i,j) on x (this is a backward insertion)
	int t = x[i];
	memmove(x+j+1,x+j,sizeof(int)*(i-j));
	x[j] = t;
}


inline void bwdins2(int* x, int i, int j) {
	//assume i>j+1 and perform insertion (i,j) two times on x (these are two backward insertions)
	//note that, if i==j+1 then x remains unchanged but this case is not handled here
	int t1 = x[i-1];
	int t2 = x[i];
	memmove(x+j+2,x+j,sizeof(int)*(i-j-1));
	x[j] = t1;
	x[j+1] = t2;
}
#endif
//END VNS4 CODE



//BEGIN INS LOCAL SEARCH (used for lop)
//simple insert local search first improvement (no optimization)
void localSearch_ins(int* x, int& fx) {
	//if x is in the history, return
	int i;
	for (i=0; i<n; i++)
		if (lsHistory[i]!=x[i])
			break;
	if (i==n)
		return;
	//set lsmode
	lsmode = true;
	//save fseed
	int fseed = fx;
	//1st improvement insertion local search (using x to choose the indexes order)
	bool imp;
	int j,ii,jj,fy;
	do {
		imp = false;
		for (i=0; i<n && !imp; i++) { //!imp for the 1st improvement style
			for (j=0; j<n && !imp; j++) { //!imp for the 1st improvement style
				ii = x[i];
				jj = x[j];
				if (ii!=jj && jj!=ii-1) {
					ins(x,ii,jj);
					fy = eval(x);
#ifdef MYDEBUG
					if (!permValid(x,n)) {
						cout<<"localSearch_ins"<<endl;
						exit(1);
					}
#endif
					nfes++;
					nfesls++;
					updateGbest(x,fy);
					if (termination())
						return;
#ifdef MINIMIZATION
					if (fy<fx) {
#else
					if (fy>fx) {
#endif
						fx = fy;
						imp = true;
					} else
						ins(x,jj,ii);
				}
			}
		}
	} while (imp);
#ifdef MYDEBUG
	if (!permValid(x,n)) {
		cout<<"localSearch_ins"<<endl;
		exit(1);
	}
#endif
	//save x in the history
	memcpy(lsHistory,x,permByteSize);
	//update local search statistics
	nls++;
#ifdef MINIMIZATION
	if (fx<fseed) {
#else
	if (fx>fseed) {
#endif
		nImprovingls++;
#ifdef MINIMIZATION
		totImprovingls += fseed-fx;
#else
		totImprovingls += fx-fseed;
#endif
	}
	//unset lsmode
	lsmode = false;
	//done
}
//END INS LOCAL SEARCH

