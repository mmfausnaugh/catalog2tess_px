/* astron_rb.h
 * Alan M. Levine
 * June 10, 2004 - from astron_a.h
 *
 */

extern void sunpos(double mjd, double possun[2], double rsun[3],
		   double vearth[3]);
extern void precess_radec(double mjd_from, double radec_f[2], double mjd_to,
			double radec_t[2]);
extern void precess_mat(double mjd, double pmat[3][3]);
extern void celgal(double ra, double dec, double *glong, double *glat);
extern void celgalmat(double rmat[3][3]);
