//INCLUDE IT IN DEP.CPP

//MAIN FUNCTION DECLARATIONS
void diffMutationBS(int i);
void diffMutationIS(int i);
void diffMutationSS(int i);
void diffMutationIS2(int i);

//INLINE FUNCTION DECLARATIONS
inline void randbs(int* s, int& l, int* x, int* ainv, int nainv);
inline void randis(int* s, int& l, int* x, double f=1.0);
inline void randss(int* s, int& l, int* x, double f=1.0);
inline void randmergess(int* s, int& l, int* x, double f=1.0);
inline void randis2(int* s, int& l, int* x, double f=1.0);



//MAIN FUNCTION DEFINITION FOR ADJACENT SWAPS
void diffMutationBS(int i) {
	//DE/rand/1: y1[i] = x[r0] + F * (x[r1] - x[r2])
	//(0) initialize variables
	int r0,r1,r2,t,k,j,w,*p,*ixx,*pt,*pp,*p0;
	double tempDouble;                           //temp variable for ceil rounding
	static int* ss = tmpint;                     //use the temp memory (sorting sequence)
	static int* z = tmpint+(n*(n-1)/2);          //use the temp memory (just at right of ss)
	static int* ainv = z+n;                      //use the temp memory (just at right of z)
//int mydebug;//debug
	//(1) get 3 different indexes r0,r1,r2 different also from i
	threeRandIndicesDiffFrom(np,i,r0,r1,r2);
	//(2) compute scale factor using jde rule
	sfy[i] = urand()<.1 ? fmin+(fmax-fmin)*urandi() : sfx[i];
	//(3) check and handle the easy case F=1
	if (sfy[i]==1.) {
		//begin case F=1
		//compute y1[i] = x[r0] + (x[r1] - x[r2]) = x[r0] * x[r2]^-1 * x[r1]
		p = y1[i];
		p0 = x[r0];
		ixx = ix[r2];
		pp = x[r1];
		for (j=0; j<n; j++)
			p[j] = p0[ixx[pp[j]]];
#ifdef MYDEBUG
		int myp[n];
		for (j=0; j<n; j++)
			myp[j] = x[r0][ix[r2][x[r1][j]]];
		for (j=0; j<n; j++)
			if (y1[i][j]!=myp[j]) {
				cout << "Problem with the case F=1" << endl;
				exit(1);
			}
#endif
		//nothing else to do (this is an easy case) so return
		return;
		//end case F=1
	}
	//(4) distinguish the 2 cases F<1 and F>1
	if (sfy[i]<1.) {
		//begin case F<1: build z with its adjacent inversions, randbs on z, compute k, compute starting y
		//(5a) build z with its adjacent inversions ainv/k (k is #ainv)
		//z = x[r1]-x[r2] = x[r2]^-1*x[r1] = sort_seq(x[r1]^-1*x[r2]) (I already know x[r1]^-1=ix[r1])
		ixx = ix[r1];
		p = x[r2];
		*z = ixx[*p];           //first outside the for ... it's the same of "z[0] = ixx[p[0]]"
		k = 0;                  //init no. of adjacent inversions
		pt = ainv;              //due because ainv is static
		for (j=1; j<n; j++) {   //from the second to the last inside the for (so the if is ok)
			z[j] = ixx[*++p];    //"ixx[*++p]" is the same of "ixx[p[j]]"
			if (z[j-1]>z[j]) {
				*pt++ = j-1;      //considering also next line, it is the same of "ainv[k++]=j-1"
				k++;
			}
		}
		//(6a) ss,w = sort_sequence(z) using randbs sorting loop
		randbs(ss,w,z,ainv,k);
		//(7a) compute k, i.e. where to truncate ss, basing on w and sfy[i] and using ceil rounding
		k = tempDouble = sfy[i] * w;  //tempDouble needed to implement ceil rounding
		if (k<tempDouble)             //this is to implement ceil rounding
			k++;
		//(8a) compute starting y: y1[i] becomes a copy of x[r0]
		memcpy(y1[i],x[r0],permByteSize);
		//end case F<1
	} else {
		//begin case F>1: build z with its adjacent inversions, compute starting y, randbs on z, compute k
		//(5b) build z with its adjacent inversions ainv/k (k is #ainv)
		//DEC(R) = SORT(E->X) + SORT(X->R) = X + SORT(X->R) = X + SORT(R^-1*X) = X + SORT(R*X)
		//      thus the z to sort has to be R*(x[r1]-x[r2]) = R * x[r2]^-1 * x[r1]
		//      note: the multiplication at left by R can be simulated, i.e. (R*X)[i] = R[X[i]] = n-1-X[i]
		//(8b - anticipated) together compute also starting y, i.e. x[r0] * x[r2]^-1 * x[r1]
		ixx = ix[r2];
		p = x[r1];
		*z = n-1-ixx[*p];       //first outside the for ... it's the same of "z[0] = n-1-ixx[p[0]]"
		p0 = x[r0];             //for starting y (8b)
		pp = y1[i];             //for starting y (8b)
		*pp = p0[ixx[*p]];      //for starting y (8b), 1st outside for, equal to y1[i][0]=x[r0][ixx[x[r1][0]]]
		k = 0;                  //init no. of adjacent inversions
		pt = ainv;              //due because ainv is static
		for (j=1; j<n; j++) {   //from the second to the last inside the for (so the if is ok)
			z[j] = n-1-ixx[*++p];//"ixx[*++p]" is the same of "ixx[p[j]]"
			pp[j] = p0[ixx[*p]]; //for starting y (8b), equal to "y1[i][j] = x[r0][ixx[x[r1][j]]]"
			if (z[j-1]>z[j]) {
				*pt++ = j-1;      //considering also next line, it is the same of "ainv[k++]=j-1"
				k++;
			}
		}
		//(6b) ss,w = sort_sequence(z) using randbs sorting loop
		randbs(ss,w,z,ainv,k);
#ifdef MYDEBUG
		//I have to check that "x[r2]^-1 * x[r1] * ss = reverse_identity"
		int myp[n];
		for (j=0; j<n; j++)
			myp[j] = ix[r2][x[r1][j]];
		for (j=0; j<w; j++) {
			t = myp[ss[j]];
			myp[ss[j]] = myp[ss[j]+1];
			myp[ss[j]+1] = t;
		}
		for (j=0; j<n; j++)
			if (myp[j]!=n-1-j) {
				cout << "Problem with case F>1" << endl;
				exit(1);
			}
#endif
		//(8b) compute k, i.e. where to truncate ss, basing on w and sfy[i] and using ceil rounding
		//in this case (F>1), L=D-W, K = min{ceil(F*L),D}-L = min{ceil(F*(D-W)),D}-D+W
		k = tempDouble = sfy[i]*(diameter-w);  //tempDouble needed to implement ceil rounding
		if (k<tempDouble)                      //this is to implement ceil rounding
			k++;
		if (k>diameter)                        //this is the minimum part of the formula above
			k = diameter;
		k += w - diameter;                     //this is the last "-D+W" part
//mydebug=w;//debug
#ifdef MYDEBUG
		if (k<0 || k>w) {
			cout << "k is not ok in the case F>1" << endl;
			exit(1);
		}
#endif
		//end case F>1
	}
	//(9) apply the first k entries of ss to y1[i]
	p = y1[i];
	pt = ss;             //due beecause ss is a static variable
	for (j=0; j<k; j++) {
		w = *pt++;        //"*pt++" is the same of "ss[j]"
		t = p[w];         //this line and the next two implement the adjacent swap
		p[w] = p[w+1];
		p[w+1] = t;
	}
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"diffMutation y1["<<i<<"]"<<endl;
		exit(1);
	}
	//to check the following please uncomment the two lines with the variable "mydebug" above
	//check the correct number of generators
	//int myz[n];    //put x[r1]^-1 * x[r2] in myz
	//for (j=0; j<n; j++)
	//   myz[j] = ix[r1][x[r2][j]];
	//int myai[n];   //put the adj. inv. of myz in myai
	//int mynai = 0;
	//for (j=1; j<n; j++)
	//   if (myz[j-1]>myz[j])
	//      myai[mynai++] = j-1;
	//int myss[n*n]; //put the sort. seq. of myz in myss
	//int mylss;
	//randbs(myss,mylss,myz,myai,mynai);
	//int mygen;     //compute the no. of generators to check
	//double mytd;
	//mygen = mytd = sfy[i]*mylss;
	//if (mygen<mytd)
	//   mygen++;
	//if (sfy[i]<1.) { //check if everything ok - case f<1
	//   if (k!=mygen) {
	//      cout << "Number of generators mismatch when F<1" << endl;
	//      exit(1);
	//   }
	//} else {       //check if everything ok - case f>1
	//   //L=D-W, K = min{ceil(F*L),D}-L = min{ceil(F*(D-W)),D}-D+W
	//   if (mygen>diameter)
	//      mygen = diameter;
	//   int myl2 = mydebug; //length of 2nd part
	//   int myl1 = diameter-myl2; //length of 1st part
	//   if (myl1!=mylss) {
	//      cout << "Problem with lengths when F>1" << endl;
	//      exit(1);
	//   }
	//   if (mylss+k!=mygen) {
	//      cout << "Number of generators mismatch when F>1" << endl;
	//      exit(1);
	//   }
	//}
#endif
	//done
}



//INLINE FUNCTION DEFINITIONS FOR ADJACENT SWAPS
//s/l = sequence that sorts x, ainv/nainv = precomputed adjacent inversions of x
inline void randbs(int* s, int& l, int* x, int* ainv, int nainv) {
	int i,j,t;
	//initialize sequence length to 0
	l = 0;
	//sorting loop
	while (nainv>0) {
		//randomly select an adjacent inversion (i,j)=(i,i+1)
		t = irand(nainv);
		i = ainv[t];
		j = i+1;
		//remove (i,j), i.e. i, from the list of adjacent inversions, i.e. ainv
		ainv[t] = ainv[--nainv];
		//add previous inversion, i.e. (i-1,i)=(t,i), to ainv if the case
		t = i-1;
		if (i>0 && x[t]<x[i] && x[t]>x[j])
			ainv[nainv++] = t;
		//add next inversion, i.e. (i,i+1)=(j,t), to ainv if the case
		t = j+1;
		if (t<n && x[t]<x[i] && x[t]>x[j])
			ainv[nainv++] = j;
		//apply the swap (i,j) to x
		t = x[i];
		x[i] = x[j];
		x[j] = t;
		//save the swap in s (only i)
		s[l++] = i;
	}
	//done
}



//MAIN FUNCTION DEFINITION FOR INSERTIONS
void diffMutationIS(int i) {
	//DE/rand/1: y1[i] = x[r0] + F * (x[r1] - x[r2])
	//initialize variables
	int r0,r1,r2,*p,*pt,*ix1,lss,j,insi,insj;
	static int* ss = tmpint; //use the temp memory (sorting sequence - 2*n max length)
	static int* ix1dotx2 = tmpint+5*n; //use the last chunk of temp memory
	//get 3 different indexes r0,r1,r2 different also from i
	threeRandIndicesDiffFrom(np,i,r0,r1,r2);
	//compute scale factor and truncation bound using jde rule (need before randis for limit version)
	sfy[i] = urand()<.1 ? fmin+(fmax-fmin)*urandi() : sfx[i];
	//x[r1]-x[r2] = x[r2]^-1*x[r1] = sort_seq.(x[r1]^-1*x[r2]) in 2 steps (I already know x[r1]^-1=ix[r1]):
	//(1) ix1dotx2 = ix1 * x[r2]
	ix1 = ix[r1];
	p = x[r2];
	for (j=0; j<n; j++)
		ix1dotx2[j] = ix1[*p++]; //same as ix1[p[j]] = ix1[x[r2][j]]
#ifdef MYDEBUG
	if (!permValid(ix1dotx2,n)) {
		cout<<"diffMutation ix1dotx2"<<endl;
		exit(1);
	}
#endif
	//(2) ss,lss = randis(ix1dotx2) using limited version (sfy[i])
	randis(ss,lss,ix1dotx2,sfy[i]);
	//apply the inverse of the insertions in ss to x[r0] and put the result in y1[i]
	//all the insertions in ss, since we are using limited version of randis
	p = y1[i];
	memcpy(p,x[r0],permByteSize);
	pt = ss; //... since ss is a static variable
	for (j=0; j<lss; j++) { //ss is a sorting sequence, so "inversion happened automatically"
		insi = *pt++;
		insj = *pt++;
		ins(p,insi,insj);
	}
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"diffMutation y1["<<i<<"]"<<endl;
		exit(1);
	}
#endif
	//done
}


//INLINE FUNCTION DEFINITIONS FOR INSERTIONS
//s/l = sequence that sorts x using insertions (INS_k(i,j) is s[2*k]<->s[2*k+1])
//f is used to bound the insertions chain (it will be long l=ceil(f*original_l))
void randis(int* s, int& l, int* x, double f) {
	//variables
	int i,j,k,y,a,b,m,ql,ll,ul,ind,w;
	double tempDouble; //temp variable for ceil rounding
	//initialize working memory (note that s is tmpint and take 2*n space)
	static int* lis = tmpint+2*n; //just at the right of s   ("alive" together with s,u,su)
	static int* u = lis+n;        //just at the right of lis ("alive" together with s,lis,su)
	static int* su = u+n;         //just at the right of u   ("alive" together with s,u,lis)
	static int* t = tmpint;       //reuse part of s space    ("alive" together with lis,u,q but not s)
	static int* q = tmpint+n;     //reuse part of s space    ("alive" together with lis,u,t but not s)
	/*
	PHASE1) Compute the following data:
		- t = array s.t. t[x[i]] is the length of a lis ending with x[i]
		- ll = length of a lis
		- ind = index of a uniformly random element that ends a lis (with length k)
		COST: O(n*log(n))
	*/
	//initialize length of the priority queue and length of the lis to 0
	ql = ll = 0;
	//scan the array (left-to-right) in order to: set t[x[i]], update ll and ind
	for (i=0; i<n; i++) {
		//y is the current value
		y = x[i];
		//binary search to find the index a (==b) of the smallest value in q greater than y
		a = 0;
		b = ql;
		while (a<b) {
			m = (a+b)/2;
			if (y>q[m])
				a = m+1;
			else
				b = m;
		}
		//replace q[a] with y and increase ql if a==ql
		q[a] = y;
		if (a==ql)
			ql++;
		//set t[y]
		t[y] = a==0 ? 1 : 1+t[q[a-1]];
		//update ll (max length) and ind (random index with length ll using reservoir sampling)
		if (t[y]>ll) {
			ll = t[y];
			k = 1;            //init/reset reservoir sampling counter
			ind = i;
		} else if (t[y]==ll) {
			k++;              //update reservoir sampling counter
			if (irand(k)==0)  //same as rand01()<=1/k (reservoir sampling choice)
				ind = i;
		}
	} //end-for
	/*
	PHASE2) Compute the following data:
		- lis = array of length ll s.t. x[lis[i]] is the i-th _value_ of a RANDOM lis
		- u = array of length ul=N-ll s.t. x[u[i]] is _not_ a _value_ of the lis, u is sorted
		COST: O(nlogn) because we have to sort u at the end, otherwise O(n)
	*/
	//the indexes from ind+1 to N-1 go directly in u
	ul = n-ll;
	for (i=n-1; i>ind; i--) //are n-1-ind elements to insert
		u[--ul] = i;
	//set target length (m=ll-1) and last lis index (ind)
	m = ll - 1;
	lis[m] = ind;
	//set minimum length value for which an index entered in lis so far
	a = ll;
	//scan (right-to-left) from ind-1 to 0 in order to fill lis
	for (i=ind-1; i>=0; i--) {
		//y is the (current) length of a lis ending at x[i], thus i is a candidate for lis[y-1]
		y = t[x[i]];
		//current value x[i] is feasible iff:
		//(1) its length y (==t[x[i]]) is greater or equal to the target length, AND
		//(2) its length differs from lis length ll (other end-values excluded by res.sampl. of phase1), AND
		//(3) it is smaller than the lis value at position y (for sure already filled)
		if (y>=m && y<ll && x[i]<x[lis[y]] ) { //feasible
			//two cases to consider:
			//(1) normal: y==m (current length matches target length)
			//(2) backtracking: y>m (current length greater than target length)
			if (y==m) { //normal case
				//update min length value for which an index entered in lis so far
				//(this is the only place where new smaller length value can be discovered)
				//or, if not a new min, for sure there was a previous index in lis that has to move in u
				if (m<a)
					a = m;
				else
					u[--ul] = lis[m-1];
				//reset/set length observation counter to 1 (used for reservoir sampling) (no need to init)
				q[m] = 1;
				//set lis index and decrement target length (it is like reservoir sampling with probability 1)
				lis[--m] = i;
				//end of normal case
			} else { //backtracking case (y>m for sure, since y<m has been already discarded)
				//increment the counter of observed values with this length
				q[y]++;
				//using res.sampling: set lis/u index, update target length and u
				if (irand(q[y])==0) {   //same as rand01()<1/q[y]
					m = y-1;          //new target length
					u[--ul] = lis[m]; //old lis value goes in u
					lis[m] = i;       //replace lis[m] with i
				} else //if res.sampling say no, i goes directly in u
					u[--ul] = i;
				//end of backtracking case
			}
			//end of feasible branch
		} else { //unfeasible
			//i goes directly in u
			u[--ul] = i;
			//end of unfeasible branch
		}
		//proceed to next/left item
	}
	//sort u first part (till current ul) (the last part is already sorted for sure)
	//insertionSort(u,ul);
	//restore u length
	ul = n-ll;
	//sort u basing on the values in x
	insertionSortWithValues(u,ul,x);
#ifdef MYDEBUG
	if (!isSortedWithValues(u,ul,x)) {
		cout << "u is not sorted wrt x" << endl;
		exit(1);
	}
#endif
	/*
	PHASE3) The following part of code:
		- compute l (==N-ll), i.e., the length of a decomposition (sorting sequance) of x
		- compute a random sorting sequence of x in the output parameter s in a way that
			x^-1 = INS(s[0],s[1]) * ... * INS(s[2*(l-1)],s[2*(l-1)+1]) = PROD_{i=0}^{l-1} INS(s[2*i],s[2*i+1])
			x = INS(s[2*(l-1)+1],s[2*(l-1)] * ... * INS(s[1],s[0]) = PROD_{i=l-1}^{0} INS(s[2*i+1],s[2*i])
		- sort x
		COST: O(n^2)
	*/
	//compute target length for the sorting sequence (limit version)
//f=1.0;//debug
	l = tempDouble = f*ul;  //ceil rounding ok
	if (l<tempDouble)       //ceil rounding ok
		l++;                 //ceil rounding ok
	//compute successors of u in lis (su[i] such that lis[su[i]] is the first lis element > u[i])
	j = 0;
	for (i=0; i<ul; i++) {  //cost O(n) since ul+ll=n and both i and j advance of one step
		while (j<ll && x[lis[j]]<x[u[i]])
			j++;
		su[i] = x[lis[j]]>x[u[i]] ? j : ll; //because j cannot reach ll in the while above
	}
	//initialize s true length
	w = 0;
	//while target length not reached increase it and perform bookkeeping
	while (w<l) { //begin of sorting loop
		//select i from u using roulette wheel with "weighted reservoir sampling" technique
		m = 0;                                    //initialize sum of observed weights
		for (i=0; i<ul; i++) {
			b = su[i]==ll ? n : lis[su[i]];        //right bound of current item
			a = su[i]==0 ? 0 : lis[su[i]-1];       //left bound of current item
			y = a==b ? 1 : b-a;                    //weight for current item
			m += y;                                //update sum of observed weights
			if (irand(m)<y)                        //same as urand01()<weight/m
				k = i;                              //update reservoir (of size 1)
		}
		i = u[k]; //the chosen i with bounds a,b / Pos k in u / Succ in lis at pos su[k]
		//compute the left/right bounds for the allowed j
		b = su[k]==ll ? n : lis[su[k]];           //right bound of current i==u[k]
		a = su[k]==0 ? 0 : lis[su[k]-1];          //left bound of current i==u[k] 
		//choose a random index j from: (1) [a,a] if a==b, (2) [a,b) if i<a, (3) (a,b] if i>b
		j = a==b ? a : (a + (i<a?0:1) + irand(b-a));
		//append INS(i,j) in the sorting sequence s (and increase its length w)
		s[2*w] = i;
		s[2*w+1] = j;
		w++;
		//update x,lis,u,su by modifying their values after x=x*INS(i,j)
		//(LIS_UPDATE_PART1) insert index j at position su[k] in lis
		//    by also moving right lis elements from su[k] to ll and incrementing ll
		a = su[k];  //a is now the position in lis where index j has to be inserted
		if (a<ll)   //this if is needed only when i<j
			memmove(lis+a+1,lis+a,sizeof(int)*(ll-a));
		ll++;
		lis[a] = j;
		//(U/SU_UPDATE_PART1) remove item k from u/su by decreasing ul and move left the chunk [k+1,ul)
		ul--;
		if (k<ul) {
			memmove(u+k,u+k+1,sizeof(int)*(ul-k));
			memmove(su+k,su+k+1,sizeof(int)*(ul-k));
		}
		//(X_UPDATE_PART1) store in a temporary variable (y) x[i]
		y = x[i];
		//now 2 cases: i<j and i>j
		if (i<j) {
			//(LIS_UPDATE_PART2a) decrement lis[h] s.t. 0<=h<a and lis[h]>i (now lis is updated)
			for (; --a>=0 && lis[a]>i;)   //note that a still is su[k]
				lis[a]--;
			//(U_UPDATE_PART2a) decrement the u elements that are in (i,j] (now u is updated)
			for (a=0; a<ul; a++)
				if (u[a]>i && u[a]<=j)
					u[a]--;
			//(X_UPDATE_PART2a) move-left values x[i+1..j] because it is a forward insertion
			memmove(x+i,x+i+1,sizeof(int)*(j-i));
			//end case i<j
		} else { //i>j
			//(LIS_UPDATE_PART2b) increment lis[h] s.t. a<h<ll and lis[h]<i (now lis is updated)
			for (; ++a<ll && lis[a]<i;)   //note that a still is su[k]
				lis[a]++;
			//(U_UPDATE_PART2b) increment the u elements that are in [j,i) (now u is updated)
			for (a=0; a<ul; a++)
				if (u[a]>=j && u[a]<i)
					u[a]++;
			//(X_UPDATE_PART2b) move-right values x[j..i] because it is a backward insertion
			memmove(x+j+1,x+j,sizeof(int)*(i-j));
			//end case i>j
		}
		//(SU_UPDATE_PART2) increment su elements in the last chunk [k,ul) (now su is updated)
		for (a=k; a<ul; a++)
			su[a]++;
		//(X_UPDATE_PART3) finalize x*INS(i,j) by storing in x[j] the temporary value y (now x is updated)
		x[j] = y;
		//go to next sorting move
#ifdef MYDEBUG
		if (!isSortedWithValues(lis,ll,x)) {
			cout << "Problem with lis" << endl;
			exit(1);
		}
		if (!isSortedWithValues(u,ul,x)) {
			cout << "Problem with u" << endl;
			exit(1);
		}
		if (ll+ul!=n) {
			cout << "ll+ul not n" << endl;
			exit(1);
		}
		int myj = 0;
		int mysu[n];
		for (int myi=0; myi<ul; myi++) {
			while (myj<ll && x[lis[myj]]<x[u[myi]])
				myj++;
			mysu[myi] = x[lis[myj]]>x[u[myi]] ? myj : ll;
			int myb = mysu[myi]==ll ? n : lis[mysu[myi]];
			int mya = mysu[myi]==0 ? 0 : lis[mysu[myi]-1];
		}
		for (int myi=0; myi<ul; myi++)
			if (mysu[myi]!=su[myi]) {
				cout << "Problem with su" << endl;
				exit(1);
			}
#endif
	} //end of sorting loop
#ifdef MYDEBUG
	if (w!=l) {
		cout<<"w!=l"<<endl;
		exit(1);
	}
	//to check this, please uncomment the "f=1.0" above
	//if (!isSorted(x,n)) {
	//   cout << "x is not sorted" << endl;
	//   exit(1);
	//s}
#endif
	//done
}


//MAIN FUNCTION DEFINITION FOR EXCHANGES
void diffMutationSS(int i) { //EXCHANGES!!!
	//DE/rand/1: y1[i] = x[r0] + F * (x[r1] - x[r2])
	//initialize variables
	int r0,r1,r2,*p,*p0,*pp,*ixx,lss,j,exci,excj,t;
	static int* ss = tmpint;	//temp memory (max length 2*n)	//EXCHANGES
	static int* z = tmpint;		//tempmem (n) alive with ss so z DESTROYED ON RANDSS or RANDMERGESS //EXCHANGES!!!
	//(1) get 3 different indexes r0,r1,r2 different also from i
	threeRandIndicesDiffFrom(np,i,r0,r1,r2);
	//(2) compute scale factor and truncation bound using jde rule (need before randss for limit version)
	sfy[i] = urand()<.1 ? fmin+(fmax-fmin)*urandi() : sfx[i];
	//(3) check and handle the easy case F=1
	if (sfy[i]==1.) {
		//begin case F=1
		//compute y1[i] = x[r0] + (x[r1] - x[r2]) = x[r0] * x[r2]^-1 * x[r1]
		p = y1[i];
		p0 = x[r0];
		ixx = ix[r2];
		pp = x[r1];
		for (j=0; j<n; j++)
			p[j] = p0[ixx[pp[j]]];
#ifdef MYDEBUG
		int myp[n];
		for (j=0; j<n; j++)
			myp[j] = x[r0][ix[r2][x[r1][j]]];
		for (j=0; j<n; j++)
			if (y1[i][j]!=myp[j]) {
				cout << "Problem with the case F=1" << endl;
				exit(1);
			}
#endif
		//nothing else to do (this is an easy case) so return
		return;
		//end case F=1
	}
	//(4) distinguish the 2 cases F<1 and F>1
	if (sfy[i]<1.) {
		//begin case F<1: build z, randss on z, compute starting y
		//(5a) build z that will be input for randss (consider x[r1]-x[r2] = x[r2]^-1*x[r1] = sort_seq.(x[r1]^-1*x[r2])) thus z=x[r1]^-1*x[r2]
		ixx = ix[r1];
		p = x[r2];
		for (j=0; j<n; j++)
			z[j] = ixx[*p++]; //same as ixx[p[j]] == ixx[x[r2][j]]
		//(6a) ss,lss = randss(z) using limited version (sfy[i])
		randss(ss,lss,z,sfy[i]);
		//(7a) compute starting y, i.e., copy x[r0] on y1[i]
		memcpy(y1[i],x[r0],permByteSize);
		//end case F<1
	} else {
		//begin case F>1: build z, randss on z, compute starting y
		//(5b) build z that will be input for randmergess (consider I have to go from x[r1]-x[r2] = x[r2]^-1*x[r1] towards a 1cycle permutation) thus z=x[r2]^-1*x[r1]
		//(7b anticipated) compute starting y, i.e., assign to y1[i] the composition x[r0]*x[r2]^-1*x[r1] == x[r0]*z
		ixx = ix[r2];
		p = x[r1];
		p0 = x[r0];
		pp = y1[i];
		for (j=0; j<n; j++) {
			z[j] = ixx[*p++]; //same as ixx[p[j]] == ixx[x[r1][j]]
			pp[j] = p0[z[j]]; // same as y1[i][j] = x[r0][z[j]];
		}
		//(6b) ss,lss = randmergess(z) using limited version (sfy[i])
		randmergess(ss,lss,z,sfy[i]);
		//end case F>1
	}
	//(8) apply all the exchanges in ss to y1[i], since we are using limited version of randss
	p = y1[i];
	pp = ss; //... since ss is a static variable
	for (j=0; j<lss; j++) {
		exci = *pp++;
		excj = *pp++;
		t = p[exci];
		p[exci] = p[excj];
		p[excj] = t;
	}
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"diffMutation y1["<<i<<"]"<<endl;
		exit(1);
	}
#endif   
	//done	
}


//INLINE FUNCTION DEFINITIONS FOR EXCHANGES
//USED FOR F>1 (NO CHECK INSIDE)
//s/l = sequence that bring x towards a 1 cycle permutation (SWAP_k(i,j) is s[2*k]<->s[2*k+1])
//f is used to bound the insertions chain
void randmergess(int* s, int& l, int* x, double f) {
	//some variables
	int nc,nv,nexc,i,j,t,lim,r;
	double tempDouble;
	static bool* v = (bool*)(tmpint+n); //alive with x,c,lc (reuse part of s not used by x)
	static int* lc = tmpint+(2*n);      //alive with x,s,c,v (memory just after v)
	static int* c = tmpint+(3*n);       //alive with x,v,c,lc (last chunk of memory of size n^2)
	/* PHASE1 / Cycles decomposition of x (**ALL** cycles) / outputs will be:
		- nc is the TOTAL number of cycles of x
		- c is a (n*n)-matrix s.t. c[i] contains the items of i-th cycle of x
		- lc is a n-vector s.t. lc[i] is the length of c[i]
		- nexc is the (double of the) total number of exchanges for UNsorting x (it is SUM_{i=0}^{nc-1}(lc[i]*(n-lc[i])))
		- COST: O(n)
	*/
	//initialize v to all false and other counters to 0 (also the first cycle length)
	memset(v,0,sizeof(bool)*n); //all false
	nc = nexc = nv = i = lc[0] = 0;
	//scan all the items in the order given by the permutation (first open cycle is nc=0)
	for (;;) {
		//insert current index in the currently open cycle (and increase its length)
		c[nc*n+(lc[nc]++)] = i; //in matrix form = c[nc][lc[nc]++] = i;
		//set the visited flag of the i-th index and increase nv
		v[i] = true;
		nv++;
		//if all indexes have been visited, then close the cycle and exit the loop
		if (nv==n) {
			nexc += lc[nc]*(n-lc[nc]);
			nc++;
			break;
		}
		//try to move index by following the cycle
		i = x[i];
		//check if cycle has to be closed
		if (v[i]) { //same as if (i==c[nc*n]) { //matrix form i==c[nc][0]
			//close the cycle
			nexc += lc[nc]*(n-lc[nc]);
			nc++;
			//open next cycle by setting/resetting its length
			lc[nc] = 0;
			//move to next unvisited index (no need to move circularly, leftmost unv. index ensured)
			while (v[++i]);
		} //end-if
	} //end-infinite-for
#ifdef MYDEBUG
	int debugSum2 = 0;
	for (int debugi=0; debugi<nc; debugi++)
		debugSum2 += lc[debugi];
	if (debugSum2!=n) {
		cout << "sum of cycles lengths does not match n" << endl;
		exit(1);
	}
#endif
	/* PHASE2 / Sorting sequence of x towards a 1-cycle permutation / outputs will be:
		- l (==f*(n-nc)-(n-nc)), i.e. the length of a decomposition (sorting sequence) of x
		- s, i.e. a sequence of exchanges that moves x towards one of a 1 cycle permutation
		- does not touch x
		- COST: O(n^2)
	*/
	//compute lim: "n-nc" are the exchanges from E to X; thus I need "max{f*(n-nc),diameter}" exchanges; since "n-nc" are already in X the lim is "max{f*(n-nc),diameter} - (n-nc)"
	lim = tempDouble = f*(n-nc);	//main equation part
	if (lim<tempDouble)				//ceil rounding ok
		lim++;						//ceil rounding ok
	if (lim>diameter)				//max part
		lim = diameter;				//max part
	lim -= n-nc;					//discount part
#ifdef MYDEBUG //check if lim+(n-nc)<=diameter
if (lim+(n-nc)>diameter) {
	cout << "problem with lim in randmergess" << endl;
	exit(1);
}
#endif
	//initialize sequence length
	l = 0;
	//loop until lim exchanges have been done
	while (l<lim) {
#ifdef MYDEBUG //check if nexc is ok
int debugSum1 = 0;
for (int debugi=0; debugi<nc; debugi++)
	debugSum1 += lc[debugi]*(n-lc[debugi]);
if (nexc!=debugSum1) {
	cout << "problem with nexc in randmergess" << endl;
	exit(1);
}
#endif
		//1st roulette wheel to peak cycle i (each cycle "i" has probability "(lc[i]*(n-lc[i]))/nexc")
		r = irand(nexc);
		for (i=0; i<nc; i++) {
			t = lc[i]*(n-lc[i]);	//#exchanges for cycle i
			if (r<t)                //k is the chosen cycle
				break;
			r -= t;                 //update r to implement roulette wheel
		}
		//2nd roulette wheel to peak cycle j (each cycle "j" has weight "lc[j] if j!=i, 0 if j==i" and the sum of weights is "n-lc[i]")
		r = irand(n-lc[i]);
		for (j=0; j<nc; j++) {
			t = j!=i ? lc[j] : 0;	//#exchanges for cycle j
			if (r<t)				//h is the chosen cycle
				break;
			r -= t;					//update r to implement roulette wheel
		}
#ifdef MYDEBUG //cycle i and j must be different
if (i==j) {
	cout << "randmergess: cycles i and j are equal and this is not possible!!!" << endl;
	exit(1);
}
#endif
		//store the exchange indexes in s (by selecting uniformly random indexes from cycles i and j) and increase l (note that EXC(x,y)==EXC(y,x))
		s[2*l] = c[i*n+irand(lc[i])];
		s[2*l+1] = c[j*n+irand(lc[j])];
		l++;
		//update-part1 nexc for next iteration (subtract "lc[i]*(n-lc[i]) + lc[j]*(n-lc[j])")
		nexc -= lc[i]*(n-lc[i]) + lc[j]*(n-lc[j]);
		//make i the longest cycle between cycles i and j
		if (lc[i]<lc[j]) {
			t = i;
			i = j;
			j = t;
		}
		//merge cycles i and j by appending items in cycle j to cycle i, and increase the length of cycle i
		memcpy(c+(i*n+lc[i]),c+(j*n),sizeof(int)*lc[j]);
		lc[i] += lc[j];
		//copy last cycle in cycle j (so cycle j is removed) and decrease the number of cycles
		nc--;
		memcpy(c+(j*n),c+(nc*n),sizeof(int)*lc[nc]);
		lc[j] = lc[nc];
		//update-part2 nexc for next iteration (add the exchanges for the new cycle i, i.e., add "lc[i]*(n-lc[i])")
		nexc += lc[i]*(n-lc[i]);
		//goto next iteration
	} //end of exchanges loop
#ifdef MYDEBUG
	if (l!=lim) { //limited version
		cout << "in randmergess l does not match lim" << endl;
		exit(1);
	}
#endif
	//done
}


//USED FOR F<1 (NO CHECK INSIDE)
//s/l = sequence that sorts x using exchanges (SWAP_k(i,j) is s[2*k]<->s[2*k+1])
//f is used to bound the insertions chain (it will be long l=ceil(f*original_l))
void randss(int* s, int& l, int* x, double f) {
	//some variables
	int i,nc,vc,nv,nexc,r,e,p1,p2,t,lc1,lc2,lim; //lim used in limited version
	double tempDouble; //temp variable for ceil rounding
	static bool* v = (bool*)(tmpint+n); //alive with x,c,lc (reuse part of s not used by x)
	static int* lc = tmpint+(2*n);      //alive with x,s,c,v (memory just after v)
	static int* c = tmpint+(3*n);       //alive with x,v,c,lc (last chunk of memory of size n^2)
	/* PHASE1 / Cycles decomposition of x (store only cycles **GREATER THAN 1**) / outputs will be:
		- nc is the TOTAL number of cycles of x
		- vc is the number of cycles of x GREATER THAN 1 (so the length of c)
		- c is a (n*n)-matrix s.t. c[i] contains the items of i-th cycle (g.t. 1) of x
		- lc is a n-vector s.t. lc[i] is the length of c[i]
		- nexc is the total number of exchanges for sorting x (it is SUM_{i=0}^{vc-1}(lc[i] choose 2))
		- COST: O(n)
	*/
	//initialize v to all false and other counters to 0 (also the first cycle length)
	memset(v,0,sizeof(bool)*n); //all false
	nc = vc = nexc = nv = i = lc[0] = 0;
	//scan all the items in the order given by the permutation (first open cycle is nc=0)
	for (;;) {
		//insert current index in the currently open cycle (and increase its length)
		c[vc*n+(lc[vc]++)] = i; //in matrix form = c[vc][lc[vc]++] = i;
		//set the visited flag of the i-th index and increase nv
		v[i] = true;
		nv++;
		//if all indexes have been visited, then close the cycle and exit the loop
		if (nv==n) {
			if (lc[vc]>1) { //update nexc and increase vc iff its length is greater than 1
				nexc += lc[vc]*(lc[vc]-1)/2; //(lc[vc] choose 2)
				vc++;
			}
			nc++;
			break;
		}
		//try to move index by following the cycle
		i = x[i];
		//check if cycle has to be closed
		if (v[i]) { //same as if (i==c[vc*n]) { //matrix form i==c[vc][0]
			//close the cycle
			if (lc[vc]>1) {
				nexc += lc[vc]*(lc[vc]-1)/2; //(lc[vc] choose 2)
				vc++;
			}
			nc++;
			//open next cycle by setting/resetting its length
			lc[vc] = 0;
			//move to next unvisited index (no need to move circularly, leftmost unv. index ensured)
			while (v[++i]);
		} //end-if
	} //end-infinite-for
#ifdef MYDEBUG
	int debugSum2 = 0;
	for (int debugi=0; debugi<vc; debugi++)
		debugSum2 += lc[debugi];
	debugSum2 += (nc-vc);
	if (debugSum2!=n) {
		cout << "sum of cycles lengths does not match n" << endl;
		exit(1);
	}
#endif
	/* PHASE2 / Sorting sequence of x / outputs will be:
		- l (==n-nc), i.e. the length of a decomposition (sorting sequence) of x (NO MORE BECAUSE OF LIMITED VERSION)
		- s, i.e. a sorting sequence of x such that:
			x^-1 = EXC(s[0],s[1]) * ... * EXC(s[2*(l-1)],s[2*(l-1)+1]) = PROD_{i=0}^{l-1} EXC(s[2*i],s[2*i+1])
			x = EXC(s[2*(l-1)],s[2*(l-1)+1]) * ... * EXC(s[0],s[1]) = PROD_{i=l-1}^{0} EXC(s[2*i],s[2*i+1])
		- does not sort x
		- COST: O(n^2)
	*/
	//initialize decomposition length to 0 (it will be n-nc at the end, but in lim. version f*(n-nc))
	l = 0;
	//compute limit for limited version
	//lim = f*(n-nc) + .5; //ceil rounding //(THERE WAS A BUG, it is rounding and not ceil rounding)
	lim = tempDouble = f*(n-nc);  //ceil rounding ok
	if (lim<tempDouble)           //ceil rounding ok
		lim++;                    //ceil rounding ok
	//loop until there are n cycles of length 1
	//while (vc>0) { //same as "while (l<n-nc)"  //version without limit
	while (l<lim) {                              //version with limit
#ifdef MYDEBUG
		int debugSum = 0;
		for (int debugi=0; debugi<vc; debugi++)
			debugSum += lc[debugi]*(lc[debugi]-1)/2;
		if (debugSum!=nexc) {
			cout << "sum of (lc[i] choose 2) does not match nexc" << endl;
			exit(1);
		}
#endif
		//select a random cycle i using roulette wheel where lc[i]s are the slice sizes
		r = irand(nexc);
		for (i=0; i<vc; i++) {
			e = lc[i]*(lc[i]-1)/2;  //#exchanges in cycle i
			if (r<e)                //i is the chosen cycle
				break;
			r -= e;                 //update r to implement roulette wheel
		}
		//decrease nexc by the exchanges of cycle i
		nexc -= lc[i]*(lc[i]-1)/2; //(lc[i] choose 2)
		//select two random indexes in cycle i
		twoRandIndices(lc[i],p1,p2);
		//make p1<p2 (note that EXC(x,y)==EXC(y,x))
		if (p1>p2) {
			t = p1;
			p1 = p2;
			p2 = t;
		}
		//store the exchange indexes in s and increase l (note that EXC(x,y)==EXC(y,x))
		s[2*l] = c[n*i+p1];
		s[2*l+1] = c[n*i+p2];
		l++;
		//two new cycles: C'=[0,p1]+(p2,lc[i]) and C''=(p1,p2]
		//compute length of C' (lc1) and C'' (lc2)
		lc1 = p1-p2+lc[i];   //len[0,p1] + len(p2,lc[i]) = p1+1 + lc[i]-p2-1 = p1-p2+lc[i]
		lc2 = p2-p1;         //len(p1,p2]
		//4 cases:
		//(1) lc1>1 and lc2>1: update cycle i with C' and add cycle C'' at the end
		//(2) lc1>1 and lc2=1: update cycle i with C'
		//(3) lc1=1 and lc2>1: update cycle i with C''
		//(4) lc1=1 and lc2=1: remove cycle i by moving last cycle in position i
		if (lc1>1 && lc2>1) { //case 1
			//add cycle C'' at the last position and increment vc
			memcpy(c+n*vc,c+n*i+p1+1,sizeof(int)*lc2);
			lc[vc] = lc2;
			vc++;
			//update cycle i with C' by moving the slice (p2,lc[i]) (if not empty) after p1
			t = lc[i]-p2-1;   //len(p2,lc[i])
			if (t>0)
				memmove(c+n*i+p1+1,c+n*i+p2+1,sizeof(int)*t);
			lc[i] = lc1;
			//increase nexc by the exchanges of C' and C''
			nexc += (lc1*(lc1-1)+lc2*(lc2-1))/2; //(lc1 choose 2)+(lc2 choose 2)
			//case 1 done
		} else if (lc1>1 && lc2==1) { //case 2
			//update cycle i with C' by moving the slice (p2,lc[i]) (if not empty) after p1
			t = lc[i]-p2-1;   //len(p2,lc[i])
			if (t>0)
				memmove(c+n*i+p1+1,c+n*i+p2+1,sizeof(int)*t);
			lc[i] = lc1;
			//increase nexc by the exchanges in C'
			nexc += lc1*(lc1-1)/2; //(lc1 choose 2)
			//case 2 done
		} else if (lc2>1) { //case 3
			//update cycle i with C'' by moving the slice (p1,p2] at position 0
			memmove(c+n*i,c+n*i+p1+1,sizeof(int)*lc2);
			lc[i] = lc2;
			//increase nexc by the exchanges in C''
			nexc += lc2*(lc2-1)/2; //(lc2 choose 2)
			//case 3 done
		} else { //case 4
			//move last cycle (vc-1)th in position i and decrease vc (if cycle i is not the last)
			vc--;
			if (i!=vc) {
				memcpy(c+n*i,c+n*vc,sizeof(int)*lc[vc]);
				lc[i] = lc[vc];
			}
			//case 4 done
		} //end-cases
	} //end-while
#ifdef MYDEBUG
	if (l!=lim) { //limited version (prima era l!=n-nc)
		cout << "in randss l does not match lim" << endl;
		exit(1);
	}
#endif
	//done
}


//MAIN FUNCTION DEFINITION FOR INSERTIONS 2 (FASTER BUT WITH LESS ENTROPY)
void diffMutationIS2(int i) { //INSERTIONS!!!
	//DE/rand/1: y1[i] = x[r0] + F * (x[r1] - x[r2])
	//initialize variables
	int r0,r1,r2,*p,*pt,*ix1,lss,j,insi,insj;
	static int* ss = tmpint; //use the temp memory (sorting sequence - 2*n max length)
	static int* ix1dotx2 = tmpint+4*n; //use the last chunk of temp memory
	//get 3 different indexes r0,r1,r2 different also from i
	threeRandIndicesDiffFrom(np,i,r0,r1,r2);
	//compute scale factor and truncation bound using jde rule (need before randis2 for limit version)
	sfy[i] = urand()<.1 ? fmin+(fmax-fmin)*urandi() : sfx[i];
	//x[r1]-x[r2] = x[r2]^-1*x[r1] = sort_seq.(x[r1]^-1*x[r2]) in 2 steps (I already know x[r1]^-1=ix[r1]):
	//(1) ix1dotx2 = ix1 * x[r2]
	ix1 = ix[r1];
	p = x[r2];
	for (j=0; j<n; j++)
		ix1dotx2[j] = ix1[*p++]; //same as ix1[p[j]] = ix1[x[r2][j]]
#ifdef MYDEBUG
	if (!permValid(ix1dotx2,n)) {
		cout<<"diffMutation ix1dotx2"<<endl;
		exit(1);
	}
#endif
	//(2) ss,lss = randis2(ix1dotx2) using limited version (sfy[i])
	randis2(ss,lss,ix1dotx2,sfy[i]);
	//apply the inverse of the insertions in ss to x[r0] and put the result in y1[i]
	//all the insertions in ss, since we are using limited version of randis2
	p = y1[i];
	memcpy(p,x[r0],permByteSize);
	pt = ss; //... since ss is a static variable
	for (j=0; j<lss; j++) { //ss is a sorting sequence, so "inversion happened automatically"
		insi = *pt++;
		insj = *pt++;
		ins(p,insi,insj);
	}
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"diffMutation y1["<<i<<"]"<<endl;
		exit(1);
	}
#endif
	//done
}


//INLINE FUNCTION DEFINITIONS FOR INSERTIONS 2 (FASTER BUT WITH LESS ENTROPY)
//s/l = sequence that sorts x using insertions (INS_k(i,j) is s[2*k]<->s[2*k+1])
//f is used to bound the insertions chain (it will be long l=ceil(f*original_l))
void randis2(int* s, int& l, int* x, double f) {
	//variables
	int i,j,k,y,a,b,m,ql,ll,ul,jmin,jmax,ind,w;
	double tempDouble; //temp variable for ceil rounding
	//initialize working memory (note that s is tmpint and take 2*n space)
	static int* lis = tmpint+2*n; //just at the right of s   ("alive" together with s and u)
	static int* u = lis+n;        //just at the right of lis ("alive" together with s and lis)
	static int* t = tmpint;       //reuse part of s space    ("alive" together with lis,u,q but not s)
	static int* q = tmpint+n;     //reuse part of s space    ("alive" together with lis,u,t but not s)
	/*
	PHASE1) Compute the following data:
		- t = array s.t. t[x[i]] is the length of a lis ending with x[i]
		- ll = length of a lis
		- ind = index of a uniformly random element that ends a lis (with length k)
		COST: O(n*log(n))
	*/
	//initialize length of the priority queue and length of the lis to 0
	ql = ll = 0;
	//scan the array (left-to-right) in order to: set t[x[i]], update ll and ind
	for (i=0; i<n; i++) {
		//y is the current value
		y = x[i];
		//binary search to find the index a (==b) of the smallest value in q greater than y
		a = 0;
		b = ql;
		while (a<b) {
			m = (a+b)/2;
			if (y>q[m])
				a = m+1;
			else
				b = m;
		}
		//replace q[a] with y and increase ql if a==ql
		q[a] = y;
		if (a==ql)
			ql++;
		//set t[y]
		t[y] = a==0 ? 1 : 1+t[q[a-1]];
		//update ll (max length) and ind (random index with length ll using reservoir sampling)
		if (t[y]>ll) {
			ll = t[y];
			k = 1;            //init/reset reservoir sampling counter
			ind = i;
		} else if (t[y]==ll) {
			k++;              //update reservoir sampling counter
			if (irand(k)==0)  //same as rand01()<=1/k (reservoir sampling choice)
				ind = i;
		}
	} //end-for
	/*
	PHASE2) Compute the following data:
		- lis = array of length ll s.t. x[lis[i]] is the i-th _value_ of a RANDOM lis
		- u = array of length ul=N-ll s.t. x[u[i]] is _not_ a _value_ of the lis
		COST: O(n)
	*/
	//the indexes from ind+1 to N-1 go directly in u
	ul = 0;
	for (i=ind+1; i<n; i++)
		u[ul++] = i;
	//set target length (m=ll-1) and last lis index (ind)
	m = ll - 1;
	lis[m] = ind;
	//set minimum length value for which an index entered in lis so far
	a = ll;
	//scan (right-to-left) from ind-1 to 0 in order to fill lis
	for (i=ind-1; i>=0; i--) {
		//y is the (current) length of a lis ending at x[i], thus i is a candidate for lis[y-1]
		y = t[x[i]];
		//current value x[i] is feasible iff:
		//(1) its length y (==t[x[i]]) is greater or equal to the target length, AND
		//(2) its length differs from lis length ll (other end-values excluded by res.sampl. of phase1), AND
		//(3) it is smaller than the lis value at position y (for sure already filled)
		if (y>=m && y<ll && x[i]<x[lis[y]] ) { //feasible
			//two cases to consider:
			//(1) normal: y==m (current length matches target length)
			//(2) backtracking: y>m (current length greater than target length)
			if (y==m) { //normal case
				//update min length value for which an index entered in lis so far
				//(this is the only place where new smaller length value can be discovered)
				//or, if not a new min, for sure there was a previous index in lis that has to move in u
				if (m<a)
					a = m;
				else
					u[ul++] = lis[m-1];
				//reset/set length observation counter to 1 (used for reservoir sampling) (no need to init)
				q[m] = 1;
				//set lis index and decrement target length (it is like reservoir sampling with probability 1)
				lis[--m] = i;
				//end of normal case
			} else { //backtracking case (y>m for sure, since y<m has been already discarded)
				//increment the counter of observed values with this length
				q[y]++;
				//using res.sampling: set lis/u index, update target length and u
				if (irand(q[y])==0) {   //same as rand01()<1/q[y]
					m = y-1;          //new target length
					u[ul++] = lis[m]; //old lis value goes in u
					lis[m] = i;       //replace lis[m] with i
				} else //if res.sampling say no, i goes directly in u
					u[ul++] = i;
				//end of backtracking case
			}
			//end of feasible branch
		} else { //unfeasible
			//i goes directly in u
			u[ul++] = i;
			//end of unfeasible branch
		}
		//proceed to next/left item
	}
	/*
	PHASE3) The following part of code:
		- compute l (==N-ll), i.e., the length of a decomposition (sorting sequance) of x
		- compute a random sorting sequence of x in the output parameter s in a way that
			x^-1 = INS(s[0],s[1]) * ... * INS(s[2*(l-1)],s[2*(l-1)+1]) = PROD_{i=0}^{l-1} INS(s[2*i],s[2*i+1])
			x = INS(s[2*(l-1)+1],s[2*(l-1)] * ... * INS(s[1],s[0]) = PROD_{i=l-1}^{0} INS(s[2*i+1],s[2*i])
		- sort x
		COST: O(n^2)
	*/
	//decomposition d has length equal to initial u length
	//l = ul;      //version without limit
	//l = f*ul + .5; //version with limit (THERE WAS A BUG, it is rounding and not ceil rounding)
	l = tempDouble = f*ul;  //ceil rounding ok
	if (l<tempDouble)       //ceil rounding ok
		l++;                 //ceil rounding ok
	//initialize s true length
	w = 0;
	//while the lis length is less than N increase it and perform bookkeeping
	//while (ll<n) { //same as "while (w<l)"  //version without limit
	while (w<l) {                             //version with limit
		//select and remove a random index i from u (note that ul is decremented)
		k = irand(ul);
		i = u[k];
		ul--;
		if (k!=ul)   //if it is not the last, then swap it with the last
			u[k] = u[ul];
		//binary search to find the position a (==b) of the successor value of x[i] in the lis
		a = 0;
		b = ll;
		while (a<b) {
			m = (a+b)/2;
			if (x[i]>x[lis[m]])
				a = m+1;
			else
				b = m;
		}
		//succ (pred) value of x[i] is x[lis[a]] (x[lis[a-1]]), thus we take their indexes
		jmin = a==0 ? 0 : lis[a-1];
		jmax = a<ll ? lis[a] : n-1;
		//choose a random index j from [jmin,jmax) for i<jmin or from (jmin,jmax] for i>jmax
		if (jmin==jmax)
			j = jmin;
		else
			j = jmin + (i<jmin?0:1) + irand(jmax-jmin);
		//append INS(i,j) in the sorting sequence s (and increase its length)
		s[2*w] = i;
		s[2*w+1] = j;
		w++;
		//update x,lis,u by reflecting their values after x=x*INS(i,j):
		//insert index j at position a in lis (thus, move right lis elements from a to ll and increment ll)
		if (a<ll) //this if is needed only when i<j
			memmove(lis+a+1,lis+a,sizeof(int)*(ll-a));
		ll++;
		lis[a] = j;
		//store in a temporary variable (y) x[i]
		y = x[i];
		//2 cases: i<j and i>j
		if (i<j) {
			//decrement lis[k] s.t. 0<=k<a and lis[k]>i (now lis is updated)
			for (; --a>=0 && lis[a]>i;)
				lis[a]--;
			//move-left values x[i+1..j] because it is a forward insertion
			memmove(x+i,x+i+1,sizeof(int)*(j-i));
			//decrement u[k] s.t. i<u[k]<=j (now u is updated)
			for (k=0; k<ul; k++)
				if (u[k]>i && u[k]<=j)
					u[k]--;
			//end of i<j branch
		} else { //i>j
			//increment lis[k] s.t. a<k<ll and lis[k]<i (now lis is updated)
			for (; ++a<ll && lis[a]<i;)
				lis[a]++;
			//move-right values x[j..i] because it is a backward insertion
			memmove(x+j+1,x+j,sizeof(int)*(i-j));
			//increment u[k] s.t. j<=u[k]<i (now u is updated)
			for (k=0; k<ul; k++)
				if (u[k]<i && u[k]>=j)
					u[k]++;
			//end of i>j branch
		}
		//finalize x*INS(i,j) by storing in x[j] the temporary value y (now x is updated)
		x[j] = y;
		//proceed to next insertion
	}
#ifdef MYDEBUG
	if (w!=l) {
		cout<<"w!=l"<<endl;
		exit(1);
	}
#endif
	//done
}

