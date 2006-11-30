/*Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2006-01-23 20:07:32 +0100 (lun, 23 gen 2006) $*/
#ifndef tseriesChaos_h
#define tseriesChaos_h

#include <math.h>
#include <R.h>

#define BOX 100
#define sqr(a) (a)*(a)
#define MIN(a,b) (a)<(b) ? (a) : (b)
#define MAX(a,b) (a)>(b) ? (a) : (b)

/*Converts [a][b] matrix notation to vector-type single index notation,
as used in R.
a: row number; b: col. number; nr: total number of rows*/
#define INDEX(a, b, nr) ( (b)*(nr) + (a))

/*Copy matrix 'mat' contents in vector 'vec'*/
extern int MVCONV_i, MVCONV_j;
#define MATRIX2VEC(mat, vec, nr, nc) \
	for(MVCONV_i=0; MVCONV_i<(nr); MVCONV_i++) \
		for(MVCONV_j=0; MVCONV_j<(nc); MVCONV_j++) \
			(vec)[INDEX(MVCONV_i, MVCONV_j, (nr))] = (mat)[MVCONV_i][MVCONV_j];

/*Copy vector 'vec' contents in matrix 'mat'*/
#define VEC2MATRIX(vec, mat, nr, nc) \
	for(MVCONV_i=0; MVCONV_i<(nr); MVCONV_i++) \
		for(MVCONV_j=0; MVCONV_j<(nc); MVCONV_j++) \
			(mat)[MVCONV_i][MVCONV_j] = (vec)[INDEX(MVCONV_i, MVCONV_j, (nr))];

extern int DIST2_i;
/*Squared euclidean distance between points 'a' and 'b',
in the time delay embedding space of time series 'x', with time delay
'd', embedding dimension*time-delay 'md'.
Result in 'out'.
*/
#define DIST2(x, a, b, md, d, out) \
	(out)=0.0; \
	for(DIST2_i=0; (DIST2_i<(md)); DIST2_i+=(d)) \
		(out)+=sqr(x[(a)+DIST2_i] - x[(b)+DIST2_i]);

/*Same as DIST2, but stops computation if result exceeds 'eps'
Result in 'out'.
*/
#define DIST2EPS(x, a, b, md, d, eps, out) \
	(out)=0.0; \
	for(DIST2_i=0; (DIST2_i<(md)) && ((out)<(eps)); DIST2_i+=(d)) \
		(out)+=sqr(x[(a)+DIST2_i] - x[(b)+DIST2_i]);

typedef struct {
	double *series;
	int m, d;
	int blength;
	double eps;
	int **jh;
	int *jpntr;
} boxSearch;

extern boxSearch init_boxSearch(double *, int, int, int, double);
extern void find_nearests(boxSearch, int, int, int *, double *, int *);

#endif
