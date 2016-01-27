#include "NEH.h"
#include "makespan.h"

NEH::NEH(char* problem) {
	readInstance(problem);
	nmach=m;
	njobs=n;
	comp_time.resize(nmach);
	for(int i=0;i<nmach;i++)
		comp_time[i].resize(njobs);
	sum_time.resize(njobs);
	for(int i=0;i<njobs;i++) {
		int sum=0;
		for(int j=0;j<nmach;j++) sum+=ptime[i][j];
		sum_time[i]=sum;
	}
}

int NEH::partial(vectint &p,int kk) {
	int f=p[0], i, j, ms;
	comp_time[0][f]=ptime[f][0];
	/*cout << " evaluating [";
	for(i=0;i<kk;i++) cout << p[i] << " ";*/
	for(i=1;i<nmach;i++) 
		comp_time[i][f]=comp_time[i-1][f]+ptime[f][i];
	ms=comp_time[nmach-1][f];
	for(j=1;j<kk;j++) {
		int k=p[j], h=p[j-1];
		comp_time[0][k]=comp_time[0][h]+ptime[k][0];
		for(i=1;i<nmach;i++) {
			comp_time[i][k]=max(comp_time[i-1][k],comp_time[i][h])+ptime[k][i];
		}
		ms=max(ms,comp_time[nmach-1][k]);
	}
	//cout << "] " << ms << endl;
	return ms;
}	


struct my_pair { int j; int f; };

bool orderby(const my_pair &p1,const my_pair &p2) { return p1.f>p2.f; }

int NEH::eval(vectint &pmin) {
	int i,j;
	vector<my_pair> arr(njobs);
	for(i=0;i<njobs;i++) {
		arr[i].j=i;
		arr[i].f=sum_time[i]; 
	}
	sort(arr.begin(),arr.end(),orderby);
	//for(i=0;i<njobs;i++) cout << arr[i].j << "/" << arr[i].f << " "; cout << endl;
	pmin[0]=arr[0].j; pmin[1]=arr[1].j;
	int v1=partial(pmin,2);
	pmin[0]=arr[1].j; pmin[1]=arr[0].j;
	int v2=partial(pmin,2);
	if(v1<v2) {
		pmin[0]=arr[0].j; pmin[1]=arr[1].j; }
	else {
		pmin[0]=arr[1].j; pmin[1]=arr[0].j;
	}
	int v_best, i_best;
	for(i=2;i<njobs;i++) {
		//for(j=0;j<i;j++) cout << pmin[j] << " ";
		int x=arr[i].j;
		//cout << " -- " << x << endl;
		pmin[i]=x;
		for(j=i;j>=0;j--) {
			int v=partial(pmin,i+1);
			if(j==i || v<v_best) {
				v_best=v;
				i_best=j; }
			if(j>0) { int temp=pmin[j-1]; pmin[j-1]=pmin[j]; pmin[j]=temp; }
		}
		//cout << v_best << " " << i_best << endl;
		for(j=0;j<i_best;j++) pmin[j]=pmin[j+1];
		pmin[i_best]=x;
	}
	return v_best;
}

int NEH::sum_idle(vectint& p) {
	partial(p,njobs);
	int sum=0;
	for(int i=1;i<njobs;i++) 
		for(int j=1;j<nmach;j++) 
			if(comp_time[j][i-1]>comp_time[j-1][i])
				sum+=(comp_time[j][i-1]-comp_time[j-1][i]);
	return sum;
}

