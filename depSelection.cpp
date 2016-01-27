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
		sameFitness &= !(fx[i]-fx[0]); //it's the same of "sameFitness = sameFitness && fx[i]==fx[0]"
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
	//the manifest constants TFT/MAKESPAN/LOP are mutually exclusive
#ifdef TFT //minimization
	if ( fy1[i]<fx[i] || urand()<alpha-REL_DEV(fy1[i],fx[i]) ) { //tie favors, in some way, x[i]
#endif
#ifdef MAKESPAN //minimization
	if ( fy1[i]<=fx[i] || urand()<alpha-REL_DEV(fy1[i],fx[i]) ) { //tie favors, in some way, y1[i]
#endif
#ifdef LOP //maximization
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
inline int footrule_distance(int *p1,int *p2);

void selection_crowding() {
	//some variables
	static int* child_to_replace = new int[np];	//to contain the indexes of the child to put at position i in the population (or -1)	//NO FREE MA PAZIENZA!!!
	static int* new_fit = new int[np];			//to contain the fitness of the new population individual in the case of selection		//NO FREE MA PAZIENZA!!!
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
		min_dist = footrule_distance(y1[i],x[0]); //or bubblesort_distance
		for (j=1; j<np; j++) {
			dist = footrule_distance(y1[i],x[j]); //or bubblesort_distance
			if (dist<min_dist) {
				i_closest = j;
				min_dist = dist;
			}
		}
		//update child_to_replace and new_fit
#ifdef TFT //minimization
		if (fy1[i]<new_fit[i_closest]) { //favors, in case of ties, the old population individual
#endif
#ifdef MAKESPAN //minimization
		if (fy1[i]<=new_fit[i_closest]) { //favors, in case of ties, the selected child
#endif
#ifdef LOP //maximization
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

int footrule_distance(int *p1,int *p2) {
	static int* p3 = tmpint; //size n, I can reuse tmpint memory
	int i,dist=0;
	for (i=0; i<n; i++)
		p3[p1[i]] = i;
	for (i=0; i<n; i++)
		dist += abs(i-p3[p2[i]]);
	return dist;
}
//END CROWDING SELECTION

