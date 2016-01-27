//COMMON FUNCTIONS FOR DEP.CPP

inline void ins(int* x, int i, int j) {
	//perform insertion (i,j) on x assuming nothing on i and j
	//(i,j) means "move job at pos i to pos j"
	int t = x[i];
	if (i<j) //forward
		memmove(x+i,x+i+1,sizeof(int)*(j-i));
	else //backward
		memmove(x+j+1,x+j,sizeof(int)*(i-j));
	x[j] = t;
}

