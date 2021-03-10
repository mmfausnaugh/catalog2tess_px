/* mrq_fpg4.h
 * Alan M. Levine
 * June 7, 2017
 *
 *
 * Header file for mrq_fpg4.c
 */

extern void set_indexa();
extern void mrq_define(int ndata, int ma, int mfit);
extern void mrq_undefine();
extern void print_index_lists();
extern void mrqcof_fpg(int iarr);
extern void mrqmin_init();
extern void chisq_hist_init();
extern void chisq_hist_put(double csq);
extern int mrq_test_converge(double delchi, int ncon);
extern void mrq_finish();
extern void mrqmin_core();

extern double **x, **y, **xsig, **ysig, **xmod, **ymod;
extern double **a, **dxda, **dyda, **covar, **alpha, **beta;
extern double **oneda, **atry, **damax;
extern double chisq, ochisq;
extern double alamda, alam0, alamdec, alaminc;
extern int **lista, **indexa;
