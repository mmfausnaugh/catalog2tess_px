/* mrq_fpg1.c
 * Alan M. Levine
 * June 6, 2017
 * Heritage: mrqmin.c and mrqcof.c - both from nrc (Numerical Recipes in C)
 */

/* Levenberg-Marquardt method for focal plane geometry, etc.
 *
 * Break functions into logically related pieces.
 * Indices start at 0 rather than 1.
 * Common variables are external rather than internal.
 * Use doubles rather than floats.
 * etc.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat_ra3.h"
#include "fpg1.h"
#include "mrq_fpg1.h"

#define NHIST 5

/* a[][] and lista[][] have 'ma' x 1 components of which 'mfit' are allowed to vary.
 * 
 * Instead of varying the order of the parameters, lista[i] == 1 for tthe
 * parameters to be determined and lista[i] == 0 otherwise.
 */

double **x, **y, **xsig, **ysig;
double **a, **dxda, **dyda, **covar, **alpha, **beta, **talpha, **tbeta;
double **oneda, **atry;
double chisq, ochisq, alamda;
int **lista, **indexa;
static int ndata1, ma1, mfit1;
double chisq_hist[NHIST];

/******************************************************************************/
// Call this immediately after setting lista[], which should be done after 
// calling mrq_define().

void set_indexa()
{
  int j, n;

  fprintf(stderr,"set_indexa(): beginning...\n");

  for(j=0;j<ma1;++j) {
    indexa[j][0] = 0;
  }

  n = 0;
  for(j=0;j<ma1;++j) {
    if(lista[j][0] == 1) {
      indexa[n][0] = j;
      ++n;
    }
  }
  if(n != mfit1) {
    // Error
  }
}

/******************************************************************************/
// Set values in lista[] outside of this file, after this function is called.

void mrq_define(int ndata, int ma, int mfit)
{
  fprintf(stderr,"mrq_define(): beginning...\n");

  x = matrix(ndata,1);
  y = matrix(ndata,1);
  xsig = matrix(ndata,1);
  ysig = matrix(ndata,1);

  a = matrix(ma,1);
  atry = matrix(ma,1);
  dxda = matrix(ma,1);
  dyda = matrix(ma,1);
  lista = imatrix(ma,1);
  indexa = imatrix(ma,1);

  oneda = matrix(mfit,1);
  covar = matrix(mfit,mfit);
  alpha = matrix(mfit,mfit);
  talpha = matrix(mfit,mfit);
  beta = matrix(mfit,1);
  tbeta = matrix(mfit,1);

  ndata1 = ndata;
  ma1 = ma;
  mfit1 = mfit;
}

/******************************************************************************/

void mrq_undefine()
{
  free_matrix(x,ndata1);
  free_matrix(y,ndata1);
  free_matrix(xsig,ndata1);
  free_matrix(ysig,ndata1);

  free_matrix(a,ma1);
  free_matrix(atry,ma1);
  free_matrix(dxda,ma1);
  free_matrix(dyda,ma1);
  free_imatrix(lista,ma1);
  free_imatrix(indexa,ma1);

  free_matrix(oneda,mfit1);
  free_matrix(covar,mfit1);
  free_matrix(alpha,mfit1);
  free_matrix(talpha,mfit1);
  free_matrix(beta,mfit1);
  free_matrix(tbeta,mfit1);
}

/******************************************************************************/

void print_index_lists()
{
  int i;

  fprintf(stderr,"print_index_lists(): beginning...\n");

  fprintf(stderr,"ndata1,ma1,mfit1= %d %d %d\n",ndata1,ma1,mfit1);
  for(i=0;i<ma1;++i) {
    fprintf(stderr,"i,lista,indexa = %d %d %d\n",i,lista[i][0],indexa[i][0]);
  }
}

/******************************************************************************/
// FIX mrqcof refs

void mrqcof_fpg(int iarr)
{
  int i, j, k, iccd, ng, ino;
  double xypix[2], wtx, wty, xsig2i, ysig2i, dx, dy, dchi;

  fprintf(stderr,"mrqcof_fpg(): beginning...\n");

  for(j=0;j<mfit1;j++) {
    for(k=0;k<mfit1;k++) {
      alpha[j][k] = 0.0;
    }
    beta[j][0] = 0.0;
  }
  chisq=0.0;
  fill_fpg_param_arrays(iarr);
  ng = 0;
  for(i=0;i<ndata1;i++) {
    iccd = one_star_fit_inputs(i,xypix);
    ino = get_ccdno(i);
    // fprintf(stderr,"mrqcof_fpg(): i,iccd,ino = %d %d %d\n",i,iccd,ino);
    if( (iccd >= 0) && (iccd == ino) ) {
      xsig2i = 1.0/(xsig[i][0]*xsig[i][0]);
      ysig2i = 1.0/(ysig[i][0]*ysig[i][0]);
      dx = x[i][0] - xypix[0];
      dy = y[i][0] - xypix[1];
      for(j=0;j<mfit1;j++) {
	wtx = dxda[indexa[j][0]][0]*xsig2i;
	wty = dyda[indexa[j][0]][0]*ysig2i;
	if(i==0) {
	  fprintf(stderr,"j,indexa,wtx,wty,xsig2i,ysig2i,xsig,ysig= %d %d %f %f    %f %f    %f %f\n",
		  j,indexa[j][0],wtx,wty,xsig2i,ysig2i,xsig[i][0],ysig[i][0]);
	}
	for (k=0;k<mfit1;k++) {
	  alpha[j][k] += (wtx*dxda[indexa[k][0]][0]) + (wty*dyda[indexa[k][0]][0]);
	}
	beta[j][0] += (dx*wtx) + (dy*wty);
      }
      ++ng;
      dchi = (dx*dx*xsig2i) + (dy*dy*ysig2i);
      chisq += dchi;
      fprintf(stderr,"i,dchi= %d %f\n",i,dchi);
    }
    else {
      fprintf(stderr,"i=%d; ccd no. mismatch: iccd, ccdno = %d %d\n",i,iccd,get_ccdno(i));
    }
  }
  fprintf(stderr,"ng = %d; chi-sq = %f\n",ng,chisq);
  chisq_hist_put(chisq);
}

/******************************************************************************/
// Call this function before the others below it (mrqmin_core() and mrq_finish())

void mrqmin_init()
{
  fprintf(stderr,"mrqmin_init(): beginning...\n");

  if (alamda < 0.0) {
    alamda = 0.001;
    mrqcof_fpg(0);
    ochisq = chisq;
  }
  else {
    // Error: print notice and exit
    fprintf(stderr,"mrqmin_init(): Error - alamda = %f is not < 0.0; exiting ...\n",
	    alamda);
    exit(-10);
  }
	
}

/******************************************************************************/

void chisq_hist_init()
{
  int i;

  fprintf(stderr,"chisq_hist_init(): beginning...\n");

  for(i=0;i<NHIST;++i) {
    chisq_hist[i] = -999.0;
  }
}

/******************************************************************************/

void chisq_hist_put(double csq)
{
  int i;

  fprintf(stderr,"chisq_hist_put(): beginning...\n");

  for(i=NHIST-2;i>=0;--i) {
    chisq_hist[i+1] = chisq_hist[i];
  }
  chisq_hist[0] = csq;
}

/******************************************************************************/
/* Check for convergence
 *
 * Return value = 0 => not converged
 * Return value = 1 =>  converged
 */

int mrq_test_converge(double delchi, int ncon)
{
  int icon;

  fprintf(stderr,"mrq_test_converge(): beginning...\n");

  icon = 0;
  if(fabs(chisq_hist[0] - chisq_hist[ncon]) <= delchi) {
    icon = 1;
    alamda = 0.0;
  }
  fprintf(stderr,"mrq_test_converge(): delchi,ncon,chisq_hist[0],chisq_hist[ncon],icon =\n");
  fprintf(stderr,"      %f %d %f %f %d\n",delchi,ncon,chisq_hist[0],chisq_hist[ncon],icon);

  return(icon);
}

/******************************************************************************/
// Call this function upon determination that the fitting should be done.

void mrq_finish()
{
  int j, k, ig;

  fprintf(stderr,"mrq_finish(): beginning...\n");

  if(alamda == 0.0) {
    for(j=0;j<mfit1;j++) {
      for(k=0;k<mfit1;k++) {
	covar[j][k]=alpha[j][k];
      }
      covar[j][j] = alpha[j][j]*(1.0+alamda);
      oneda[j][0] = beta[j][0];
    }
    ig = gaussj(covar,mfit1,oneda,1);
    if(ig == 1) {
      fprintf(stderr,"mrq_finish(): gaussj() received a singular matrix.\n");
      fprintf(stderr,"mrq_finish(): exiting... \n");
      exit(-13);
    }
  }
  else {
    fprintf(stderr,"mrq_finish(): Error - alamda = %f is not == 0.0; exiting ...\n",
	    alamda);
    exit(-11);
  }

}

/******************************************************************************/

void mrqmin_core()
{
  int j, k, ig;

  fprintf(stderr,"mrqmin_core(): beginning...\n");

  if (alamda > 0.0) {
    for(j=0;j<mfit1;j++) {
      for(k=0;k<mfit1;k++) {
	covar[j][k]=alpha[j][k];
	talpha[j][k]=alpha[j][k];
      }
      covar[j][j]=alpha[j][j]*(1.0 + alamda);
      oneda[j][0]=beta[j][0];
      tbeta[j][0]=beta[j][0];
    }
    ig = gaussj(covar,mfit1,oneda,1);
    if(ig == 1) {
      fprintf(stderr,"mrqmin_core(): gaussj() received a singular matrix.\n");
      // PRINT DIAGNOSTICS IF CODE COMES HERE IN TESTS
      fprintf(stderr,"mrqmin_core(): exiting... \n");
      exit(-14);
    }
    for(j=0;j<ma1;j++) {
      atry[j][0] = a[j][0];
    }
    for(j=0;j<mfit1;j++) {
      atry[indexa[j][0]][0] = a[indexa[j][0]][0] + oneda[j][0];
    }
    mrqcof_fpg(1);
    fprintf(stderr,"mrqmin_core(): chisq,ochisq,alamda= %f %f %f\n",chisq,ochisq,alamda);
    if (chisq < ochisq) {
      alamda *= 0.3;
      ochisq = chisq;
      for(j=0;j<mfit1;j++) {
	a[indexa[j][0]][0] = atry[indexa[j][0]][0];
      }
    } 
    else {
      alamda *= 3.3;
      chisq = ochisq;
      for(j=0;j<mfit1;j++) {
	for(k=0;k<mfit1;k++) {
	  alpha[j][k] = talpha[j][k];
	}
	beta[j][0] = tbeta[j][0];
      }
    }
    fprintf(stderr,"mrqmin_core(): chisq,ochisq,alamda= %f %f %f\n",chisq,ochisq,alamda);
    return;
  }
  else {
    // Error
    fprintf(stderr,"mrqmin_core(): Error - alamda = %f is not > 0.0; exiting ...\n",
	    alamda);
    exit(-12);
  }
}


/******************************************************************************/
