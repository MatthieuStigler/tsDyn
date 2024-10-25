/*
## Copyright (C) 2006  Antonio, Fabio Di Narzo. Last Modified: $Date: 2006-03-12 17:50:40 +0100 (dom, 12 mar 2006) $
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.
*/
/*Code ispired by ll-ar routine in the TISEAN package, by Rainer Hegger*/
#include <math.h>
#include <R.h>
#include "tseriesChaos.h"

void F77_CALL(dqrls)(double*, int*, int*, double*, int*, double*, double*, double*, double*, int*, int*, double*, double*);

/* Local linear autoregressive fit.
in_series: time series 0-1 scaled
in_length: time series length
in_m: embedding dimension
in_d: time delay
in_steps: steps ahead to forecast
in_tol: tolerance for OLS fitting
epsSeq: sequence of eps values
NEPS: length of epsSeq
trace: should debugging infos be printed? (0, 1, or >1)
out_res: local linear fit error as function of neighbourhood size
out_nok: number of 'good' points as function of neighborhood size
out_avfound: average number of neighbours for each 'good' point
*/
void llarinc(double *in_series, int *in_length, int *in_m, int *in_d, int *in_steps, double *in_tol, 
	  double *epsSeq, int *NEPS, int *trace, double *out_res, int *out_nok, double *out_avfound) {
	int i,j;
	double eps, *dists;
	double *series;
	int length, blength, m, d, steps;
	int *founds, nfounds;
	boxSearch bs;

	int pfound;
	double hav, avfound, hrms, err;
	double *X, *Y, *B, *QRAUX, *WORK, frc;
	int *PIVOT, curp, ti, tj, tp, ny, rank;

	/*INIT variables*/
	series = in_series;
	length = *in_length;
	m = *in_m;
	d = *in_d;
	steps = *in_steps;
	blength = length - (m-1)*d - steps;
	dists = (double*) R_alloc(blength, sizeof(double));
	founds = (int*) R_alloc(blength, sizeof(int));
	X = (double *) R_alloc(blength*(m+1), sizeof(double));
	Y = (double *) R_alloc(blength, sizeof(double));
	B = (double *) R_alloc(m+1, sizeof(double));
	PIVOT = (int *) R_alloc(m+1, sizeof(int));
	QRAUX = (double *) R_alloc(m+1, sizeof(double));
	WORK = (double *) R_alloc(2*(m+1), sizeof(double));
	ny = 1;
	tp = m+1;
	/*END INIT variables*/
	for(i = 0; i<(*NEPS); i++) {	/*for each neighborhood size*/
		eps = epsSeq[i];
		err=avfound=hrms=hav=0.0;
		pfound=0;
		bs = init_boxSearch(series, m, d, blength, eps);
		for(j=0; j<blength; j++) { /*for each point*/ 
			find_nearests(bs, steps, j, founds, dists, &nfounds);
			if((*trace)>1) Rprintf("j=%d n=%d\n", j, nfounds);
			if(nfounds > (2*(m+1))) {	/*if enough neighbours*/
				for(ti=0; ti<nfounds; ti++) {
					curp = founds[ti];
					Y[ti] = series[curp + (m-1)*d + steps];
					X[INDEX(ti, 0, nfounds)] = 1;	/*intercept*/
					for(tj=1; tj<tp; tj++)
						X[INDEX(ti, tj, nfounds)] = series[curp + (m-tj)*d];
				}
				for(tj=0; tj< tp; tj++) 
					PIVOT[tj] = tj;
				F77_CALL(dqrls)(X, &nfounds, &tp, Y, &ny, in_tol, B, Y, Y, &rank, PIVOT, QRAUX, WORK);
				/*see lm.fit R code in stats package as a clean example of using dqrls*/
				frc = B[PIVOT[0]];
				for(tj = 1; tj<rank; tj++)
					frc += series[j + ((m-PIVOT[tj])*d)] * B[PIVOT[tj]];
				err += sqr(frc - series[j + (m-1)*d + steps]);
				pfound++;
				avfound = avfound + nfounds;
				hrms += sqr(series[j+(m-1)*d+steps]);
				hav += series[j+(m-1)*d+steps];
			} /*end if enough neighbours*/
		} /*end for each point*/ 
		if(*trace) Rprintf("eps = %f\t n = %d\n", eps, pfound);
		out_nok[i] = pfound;
		if(pfound > 1) {
			avfound /= (double) pfound;
			hav /= (double) pfound;
			hrms = sqrt(fabs(hrms/(pfound-1.0) - sqr(hav)*((double)pfound/(pfound-1.0))));	/*standard deviation*/
			err = sqrt(err/(double)pfound)/hrms;	/*normalized forecasting error*/
			out_res[i] = err;
			out_avfound[i] = avfound;
		}
	} /*end for each neighborhood size*/
}

/* Local linear autoregressive fit (fitted values).
in_series: time series 0-1 scaled
in_length: time series length
in_m: embedding dimension
in_d: time delay
in_steps: steps ahead to forecast
in_tol: tolerance for OLS fitting
eps: neighborhood window
trace: should debugging infos be printed? (0, 1, or >1)
out_fitted: local linear fits
out_nOK: number of neighbours for each point
*/
void fittedllar(double *in_series, int *in_length, int *in_m, int *in_d, int *in_steps, double *in_tol, 
	  double *in_eps, int *in_trace, double *out_fitted, int *out_nOK) {
	int j;
	double eps, *dists;
	double *series;
	int length, blength, m, d, steps;
	int *founds, nfounds;
	boxSearch bs;

	double *X, *Y, *B, *QRAUX, *WORK, frc;
	int *PIVOT, curp, ti, tj, tp, ny, rank;

	/*INIT variables*/
	series = in_series;
	length = *in_length;
	m = *in_m;
	d = *in_d;
	steps = *in_steps;
	blength = length - (m-1)*d - steps;
	dists = (double*) R_alloc(blength, sizeof(double));
	founds = (int*) R_alloc(blength, sizeof(int));
	X = (double *) R_alloc(blength*(m+1), sizeof(double));
	Y = (double *) R_alloc(blength, sizeof(double));
	B = (double *) R_alloc(m+1, sizeof(double));
	PIVOT = (int *) R_alloc(m+1, sizeof(int));
	QRAUX = (double *) R_alloc(m+1, sizeof(double));
	WORK = (double *) R_alloc(2*(m+1), sizeof(double));
	ny = 1;
	tp = m+1;
	/*END INIT variables*/
	eps = *in_eps;
	bs = init_boxSearch(series, m, d, blength, eps);
	for(j=0; j<blength; j++) { /*for each point*/ 
	  R_CheckUserInterrupt();
	  find_nearests(bs, steps, j, founds, dists, &nfounds);
	  out_nOK[j] = nfounds;
	  if((*in_trace)>1) 
	    Rprintf("j=%d n=%d\n", j, nfounds);
	  if(nfounds > (2*(m+1))) {	/*if enough neighbours*/
	    for(ti=0; ti<nfounds; ti++) {
	      curp = founds[ti];
	      Y[ti] = series[curp + (m-1)*d + steps];
	      X[INDEX(ti, 0, nfounds)] = 1;	/*intercept*/
	      for(tj=1; tj<tp; tj++)
		X[INDEX(ti, tj, nfounds)] = series[curp + (m-tj)*d];
	    }
	    for(tj=0; tj< tp; tj++) 
	      PIVOT[tj] = tj;
	    F77_CALL(dqrls)(X, &nfounds, &tp, Y, &ny, in_tol, B, Y, Y, &rank, PIVOT, QRAUX, WORK);
	    /*see lm.fit R code in stats package as a clean example of using dqrls*/
	    frc = B[PIVOT[0]];
	    for(tj = 1; tj<rank; tj++)
	      frc += series[j + ((m-PIVOT[tj])*d)] * B[PIVOT[tj]];
	    out_fitted[j] = frc;
	  } /*end if enough neighbours*/
	} /*end for each point*/ 
}
