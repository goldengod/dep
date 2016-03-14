//CROSSOVER PROCEDURES

//BEGIN CROSSOVER TPII
//crossover TPII
inline void tpii(int* ch, int* fa, int* mo, int* ifa, int c1, int c2) {
	//child takes [c1,c2] from father, the other from mother using the order in mother
	//ifa is the inverse of fa
	//variables
	int k,t;
	//father part (between c1 and c2, both included)
	memcpy(ch+c1,fa+c1,sizeof(int)*(c2-c1+1));
	//left part (before c1) from mother (with j initialization to zero)
	//j = 0; //... needed by classical impl commented
	for (k=0; k<c1; k++) {
		while (ifa[(t=*mo++)]>=c1 && ifa[t]<=c2); //... same of what is commented below
		//do {
		//   t = *mo++; //... same as t = mo[j++];
		//} while (ifa[t]>=c1 && ifa[t]<=c2);
		*ch++ = t; //... same as ch[k] = t;
	}
	//right part (after c2) from mother (j was already set to the right position)
	ch += c2-c1; //... now ch is &ch_original[c2]
	for (k=c2+1; k<n; k++) {
		while (ifa[(t=*mo++)]>=c1 && ifa[t]<=c2); //... same of what is commented below
		//do {
		//   t = *mo++; //... same as t = mo[j++];
		//} while (ifa[t]>=c1 && ifa[t]<=c2);
		*++ch = t; //... same as ch[k] = t; //note that ch before the for was set to &ch_original[c2]
	}
	//done
}


//crossover of individual
void crossover_tpii(int i) {
	//TWO POINT CROSSOVER OF AGA: it produces two sons y1[i], y2[i] from the parents x[i], y1[i]
	//a son takes [c1,c2] from father, the other from mother using the order in mother
	//variables
	static int* tch = tmpint;     //temporary child for offspring 2
	static int* imut = tmpint+n;  //mutant inverse for offspring 2
	int c1,c2,t,*p;
	//two cut points c1<=c2 (also equal cutpoints)
	c1 = irand(n);
	c2 = irand(n);
	if (c1>c2) {
		t = c1;
		c1 = c2;
		c2 = t;
	}
	//(1) x[i] is father, y1[i] is mother, y2[i] becomes son
	if (nchilds==2) {	//build y2 only if 2 childs are allowed
		tpii(y2[i],x[i],y1[i],ix[i],c1,c2);
#ifdef MYDEBUG
		if (!permValid(y2[i],n)) {
			cout<<"crossover y2["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//(2) y1[i] is father, x[i] is mother, y1[i] becomes son (computing father inverse and using a temp array)
	p = y1[i];
	for (t=0; t<n; t++)
		imut[*p++] = t; //... same as imut[p[t]] = t;
#ifdef MYDEBUG
	if (!permValid(imut,n)) {
		cout<<"crossover imut"<<endl;
		exit(1);
	}
#endif
	tpii(tch,y1[i],x[i],imut,c1,c2);
	memcpy(y1[i],tch,permByteSize);
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"crossover y1["<<i<<"]"<<endl;
		exit(1);
	}
#endif
	//done
}
//END CROSSOVER TPII


//BEGIN CROSSOVER TPII-CR
//crossover of individual
void crossover_tpiicr(int i) {
	//TWO POINT CROSSOVER OF AGA: it produces two sons y1[i], y2[i] from the parents x[i], y1[i]
	//a son takes [c1,c2] from father, the other from mother using the order in mother
	//variables
	static int* tch = tmpint;     //temporary child for offspring 2
	static int* imut = tmpint+n;  //mutant inverse for offspring 2
	int c1,c2,t,*p,length;
	double cr;
	//compute crossover probability using jde rule (THERE WAS A BUG: NO NEED TO HALVE CR IN ANY CASE!!!)
	cr = cry[i] = urand()<.1 ? crmin+urandi()*(crmax-crmin) : crx[i];
	//two cut points c1<=c2 (also equal cutpoints) basing on crossover probability
	/* COSI' ERA "ESPONENZIALE"!!!
	length = 0;				//(1) init the length of the common part [c1,c1+length]
	for (j=1; j<n; j++) {	//(2) no more than n-1
		if (urand()<cr)		//(3a) increase length if the guess is less than cr
			length++;
		else
			break;			//(3b) at the first guess greater or equal to cr exit from the for
	}
	c1 = irand(n-length);	//(4) c1 is in [0,n-1-length]
	c2 = c1+length;			//(5) c2 is c1 + length
	*/
	length = n * cr;
	c1 = irand(n-length);
	c2 = c1 + length;
	//(1) x[i] is father, y1[i] is mother, y2[i] becomes son
	if (nchilds==2) {	//build y2 only if 2 childs are allowed
		tpii(y2[i],x[i],y1[i],ix[i],c1,c2);
#ifdef MYDEBUG
		if (!permValid(y2[i],n)) {
			cout<<"crossover y2["<<i<<"]"<<endl;
			exit(1);
		}
#endif
	}
	//(2) y1[i] is father, x[i] is mother, y1[i] becomes son (computing father inverse and using a temp array)
	p = y1[i];
	for (t=0; t<n; t++)
		imut[*p++] = t; //... same as imut[p[t]] = t;
#ifdef MYDEBUG
	if (!permValid(imut,n)) {
		cout<<"crossover imut"<<endl;
		exit(1);
	}
#endif
	tpii(tch,y1[i],x[i],imut,c1,c2);
	memcpy(y1[i],tch,permByteSize);
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout<<"crossover y1["<<i<<"]"<<endl;
		exit(1);
	}
#endif
	//done
}
//END CROSSOVER TPII-CR


//BEGIN CROSSOVER OBX-CR
inline void obx(int* ch, int* fa, int* mo, int* ifa, bool* mask) {
	//obx crossover for lop using the preset mask
	int j,k;
	k = 0;
	for (j=0; j<n; j++) {
		if (mask[j])
			ch[j] = fa[j];
		else {
			while (mask[ifa[mo[k]]]) //mo[k] = k-th object of mo, ifa[mo[k]] = index in fa of the k-th object of mo, mask[ifa[mo[k]]] = boolean mask for the index
				k++;
			ch[j] = mo[k++];
		}
	}
	/*for (j=0; j<n; j++) {
		if (mask[j])
			ch[j] = fa[j];
	}
	k = 0;
	for (j=0; j<n; j++) {
		o = mo[j];
		if (!mask[ifa[o]]) {
			while (mask[k])
				k++;
			ch[k++] = o;
		}
	}*/
	//done
}
 
void crossover_obxcr(int i) {
	//obx crossover with double child and in average
	//variables
	int j,*p;
	double cr;
	static int* tch = tmpint;					//temporary child for offspring 2
	static int* imut = tmpint+n;				//mutant inverse for offspring 2
	static bool* mask = (bool*)(tmpint+(2*n));	//mask memory
	//compute crossover probability using jde rule (THERE WAS A BUG: NO NEED TO HALVE CR IN ANY CASE!!!)
	cr = cry[i] = urand()<.1 ? crmin+urandi()*(crmax-crmin) : crx[i];
	//assign bool mask basing on cr and compute inverse of mutant
	p = y1[i];
	for (j=0; j<n; j++) {
		mask[j] = urand()<cr;
		imut[*p++] = j;
	}
	//father x[i], mother y1[i], child y2[i], inverse of father ix[i]
	if (nchilds==2)		//build y2 only if 2 childs are allowed
		obx(y2[i],x[i],y1[i],ix[i],mask);
	//father y1[i], mother x[i], child y1[i] but use tch, inverse of father imut / then copy tch in y1[i]
	obx(tch,y1[i],x[i],imut,mask);
	memcpy(y1[i],tch,permByteSize);
	//some debug
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout << "PROBLEMA IN OBXCR CROSSOVER!!!" << endl;
		exit(1);
	}
	if (nchilds==2 && !permValid(y2[i],n)) {
		cout << "PROBLEMA IN OBXCR CROSSOVER!!!" << endl;
		exit(1);
	}
#endif
	//done
}
//END CROSSOVER OBX-CR


//BEGIN CROSSOVER OBX
void crossover_obx(int i) {
	//obx crossover with double child and random (COME SOPRA CAMBIA SOLO IL CALCOLO DEL VETTORE MASK)
	//variables
	int j,*p,r,tmp,r2;
	double cr;
	static int* tch = tmpint;					//temporary child for offspring 2
	static int* imut = tmpint+n;				//mutant inverse for offspring 2
	static bool* mask = (bool*)(tmpint+(2*n));	//mask memory
	static int* rp = tch;						//for the random permutation
	//assign bool mask (first items from a random perm basing on a random number in [1,n-1]) (simultaneously compute inverse of mutant) (1 n-1 to avoid 0 and n because it is symmetric)
	r = irand(n-1)+1; //[0,n-1) + 1 = [1,n) = [1,n-1]
	p = y1[i];
	for (j=0; j<n; j++) {
		rp[j] = j; //identity
		mask[j] = false; //all false
		imut[*p++] = j; //inverse of mutant
	}
	for (j=0; j<r; j++) { //for r times, choose a random index on the right part of rp and keep the used indexes in the left part of rp, always assign the current rp index (rp[j])
		r2 = irand(n-j)+j;
		tmp = rp[r2];
		rp[r2] = rp[j];
		rp[j] = tmp;
		mask[tmp] = true; //tmp is rp[r2]
	}
	//father x[i], mother y1[i], child y2[i], inverse of father ix[i]
	if (nchilds==2)		//build y2 only if 2 childs are allowed
		obx(y2[i],x[i],y1[i],ix[i],mask);
	//father y1[i], mother x[i], child y1[i] but use tch, inverse of father imut / then copy tch in y1[i]
	obx(tch,y1[i],x[i],imut,mask);
	memcpy(y1[i],tch,permByteSize);
	//some debug
#ifdef MYDEBUG
	if (!permValid(y1[i],n)) {
		cout << "PROBLEMA IN OBXCR CROSSOVER!!!" << endl;
		exit(1);
	}
	if (nchilds==2 && !permValid(y2[i],n)) {
		cout << "PROBLEMA IN OBXCR CROSSOVER!!!" << endl;
		exit(1);
	}
#endif
	//done
}
//END CROSSOVER OBX

