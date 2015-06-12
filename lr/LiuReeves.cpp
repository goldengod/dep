#include "LiuReeves.h"

LiuReeves::LiuReeves(char* problem) {
   readInstance(problem);
	nmach=m;
	njobs=n;
	comp_time.resize(nmach);
	for(int i=0;i<nmach;i++)
		comp_time[i].resize(njobs);
	comp_time_artif.resize(nmach);
	artif_proc_time.resize(nmach);
}



double LiuReeves::sequence(int first,vectint &p) {
	int j, last=first;
	set<int> rem;
	double csum;
	p[0]=first;
	for(j=0;j<njobs;j++) if(j!=first)
		rem.insert(j);
	comp_time[0][0]=ptime[first][0];
	for(j=1;j<nmach;j++)
		comp_time[j][0]=comp_time[j-1][0]+ptime[first][j];
	csum=comp_time[nmach-1][0];
	for(int k=1;k<njobs;k++) {
		double ix_min=DBL_MAX; int imin;
		if(k<njobs-1) {
			for(set<int>::iterator rr=rem.begin();rr!=rem.end();++rr) {
				int r=*rr;
				set<int> rem2=rem;
				rem2.erase(r);
				double ix=index(r,rem2);
				if(ix<ix_min) {
					ix_min=ix; imin=r;
				}
			}
			last=imin; }
		else {
			last=*(rem.begin());
		}	
		p[k]=last;
		rem.erase(last);
		comp_time[0][k]=comp_time[0][k-1]+ptime[last][0];
		for(j=1;j<nmach;j++)
			comp_time[j][k]=max(comp_time[j][k-1],comp_time[j-1][k])+ptime[last][j];
		csum+=comp_time[nmach-1][k];
	}
	return csum;
}
	
	

struct my_pair { int j; double f; };

bool orderby(const my_pair &p1,const my_pair &p2) { return p1.f<p2.f; }

double LiuReeves::eval(int x,vectint &pmin) {
	int i,j;
	set<int> remaining;
	for(i=0;i<njobs;i++) remaining.insert(i);
	vector<my_pair> jf(njobs);
	for(i=0;i<njobs;i++) {
		remaining.erase(i);
		jf[i].j=i;
		jf[i].f=index(i,remaining);
		//cerr << "jobs " << i << " f=" << jf[i].f << endl;
		remaining.insert(i);
	}
	vectint p(njobs);
	sort(jf.begin(),jf.end(),orderby);
	double cmin=DBL_MAX;
	for(j=0;j<x;j++) {
		int i=jf[j].j;
		//cerr << j+1 << "-th job " << i << " f=" << jf[j].f << endl;
		double c=sequence(i,p);
		//cerr << " csum=" << c << endl;
		if(c<cmin) {
			cmin=c;
			pmin=p;
		}
	}
	return cmin;
}


double LiuReeves::index(int i,set<int> remaining) {
	int k=njobs-remaining.size()-1;
	double it=weighted_idle_time(k,i), at=artificial_tft(k,i,remaining);
	//if(k==0) cout << "job " << i << " " << it << " " << at << endl;
	double ix=(njobs-k-2)*it+at;
	//cout << "job " << i << " at " << k << " ix=" << ix << endl;
	return ix;
}

double LiuReeves::weighted_idle_time(int k,int i) {
	double it=0;
	comp_time[0][k]=(k>0 ? comp_time[0][k-1] : 0)+ptime[i][0];
	for(int j=1;j<nmach;j++) {
		double prev=k>0 ? comp_time[j][k-1] : 0;
		comp_time[j][k]=max(comp_time[j-1][k],prev)+ptime[i][j];
		double w=(double)nmach/((j+1)+k*(nmach-(j+1))/(njobs-2));
		if(comp_time[j-1][k]>prev)
			it+=w*(comp_time[j-1][k]-prev);
	}
	return it;
}

double LiuReeves::artificial_tft(int k,int i,set<int> rem) {
	for(int j=0;j<nmach;j++) {
		double pt=0;
		for(set<int>::iterator iq=rem.begin();iq!=rem.end();++iq) {
			int q=*iq;
			pt+=ptime[q][j];
		}
		artif_proc_time[j]=pt/(njobs-k-1);
	}
// maybe it's useless, because it was already computed in method weighted_idle_time
	comp_time[0][k]=(k>0 ? comp_time[0][k-1] : 0)+ptime[i][0];
	for(int j=1;j<nmach;j++) {
		double prev=k>0 ? comp_time[j][k-1] : 0;
		comp_time[j][k]=max(comp_time[j-1][k],prev)+ptime[i][j];
	}
	comp_time_artif[0]=comp_time[0][k]+artif_proc_time[0];
	for(int j=1;j<nmach;j++)
		comp_time_artif[j]=max(comp_time_artif[j-1],comp_time[j][k])+artif_proc_time[j];
	return comp_time[nmach-1][k]+comp_time_artif[nmach-1];
}


