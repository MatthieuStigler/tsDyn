/*Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2006-01-23 20:07:32 +0100 (lun, 23 gen 2006) $*/
#include "tseriesChaos.h"

#define BIN(x, id, eps) ( (int) ( (x)[(id)]/(eps)) % BOX )

void fill_boxes(double *, int, int, int, double, int **, int *);

/*
Init of boxSearch struct.
x: time series (scaled between 0 and 1)
m, d: embedding dimension and time delay
blength: total number of points in the embedding space
eps: neighborhood size
*/
boxSearch init_boxSearch(double *x, int m, int d, int blength, 
	double eps){
	boxSearch res;
	int i;

	res.series = x;
	res.m = m;
	res.d = d;
	res.blength = blength;
	res.eps = eps;

	res.jpntr = (int*) R_alloc(blength, sizeof(int));
	res.jh = (int**) R_alloc(BOX, sizeof(int*));
	for(i=0; i<BOX; i++)
		res.jh[i] = (int*) R_alloc(BOX, sizeof(int));
	fill_boxes(x, m, d, blength, eps, res.jh, res.jpntr);

	return res;
}

/*
Find all neighbours of a specific point.
bs: Properly initialized boxSearch struct
t: theiler window
n: considered point
founds: indexes of found points
dists: euclidean distances of found points
nfounds: total number of neighbours found
*/
void find_nearests(boxSearch bs, int t, int n, int *founds, double *dists, int *nfounds) {
	int i, j, id, currx, start, end;
	int xbox, ybox;
	double tmpd;
	int md;
/**/
	double *x, eps;
	int m, d, blength;
	int **jh, *jpntr;
/**/

	x = bs.series;
	m = bs.m;
	d = bs.d;
	eps = bs.eps;
	blength = bs.blength;
	jh = bs.jh;
	jpntr = bs.jpntr;

	xbox = BIN(x,n,eps);		/*Grid position     */
	ybox = BIN(x, n+(m-1)*d, eps);  /*  of current point*/
	md = m*d;
	(*nfounds)=0;
	eps = sqr(eps);

	for(i=xbox-1; i<(xbox+2); i++) { /*Scan surrounding boxes*/
	if ((i<0)||(i>=BOX)) continue;   /*     ...              */
	for(j=ybox-1; j<(ybox+2); j++) { /*     ...             */
		if((j<0)||(j>=BOX)) continue; /*...              */
		start= jh[i][j];         /*Starting candidate index*/
		end  = ((j+1)==BOX) ? (blength-1): jh[i][j+1];/*Ending candidate index*/
		for(id=start; id<end; id++) { /*For all candidates..*/
			currx = jpntr[id];    /*current candidate*/
			if (abs(currx-n)<=t) continue; /*test for theiler window*/
			DIST2EPS(x, n, currx, md, d, eps, tmpd); /*compute real distance*/
			if (tmpd<eps) { /*If a real neighbour...*/
				dists[(*nfounds)++] = sqrt(tmpd); /*store distance*/
				founds[(*nfounds)-1]= currx; /*store point*/
			} /*end if a real neighbour*/
		} /*End for all candidates*/
	} /*End scan ...           */
	} /* ... surrounding boxes */
}

/*
Fill boxes for fast neighborhood search.
x: time series (scaled between 0 and 1)
m, d: embedding dimension and time delay
blength: total number of points in the embedding space
eps: neighborhood size
jh: a BOX*BOX array, i.e., the boxes histogram
jpntr: the boxes pointers
*/
void fill_boxes(double *x, int m, int d, int blength, double eps, int **jh, int *jpntr) {
	int i, j, binx, biny, offset;
	int next;

	offset = (m-1)*d;
	for (i=0; i<BOX; i++)
		for (j=0; j<BOX; j++)
			jh[i][j]=0;
	for (i=0; i<blength; i++) {
		binx = BIN(x, i, eps);
		biny = BIN(x, i+offset, eps);
		jh[binx][biny] ++;
	}
	for (i=0; i<(BOX*BOX-1); i++) {
		next = i+1;
		jh[(next - (next % BOX)) / BOX ][next % BOX] += 
			jh[(i - (i % BOX) ) / BOX ][i % BOX];
	}

	for (i=0; i<blength; i++) {
		binx = BIN(x, i, eps);
		biny = BIN(x, i+offset, eps);
		jpntr[--jh[binx][biny]] = i;
	}
}
