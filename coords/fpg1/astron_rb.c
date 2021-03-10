/* astron_rb.c
 * Alan M. Levine
 * June 10, 2004 - from astron_a.c
 *
 * Various astronomical functions
 *
 * June 10, 2004 - redo precession routines
 */ 

#include <stdio.h>
#include <math.h>
#include "vec.h"
#include "mat_ra.h"
#include "astron_rb.h"

#define pi 3.1415926536
#define r2d 57.29578
#define d2r 1.7453292e-2
#define MjdOff 2400000.5
#define AU_DAY 1731.5      /* 1 AU/day = 1731.5 km/s */

/******************************************************************************/
/* .c adaptation of SunPos.C (Wei Cui)
 * Description : it finds the position of the sun (RA and DEC) at a given time
 *               (in modified Julian date).
 *
 * Modified - to make ansi-style - Alan M. Levine, March 26, 1998
 *     - to add unit vector to argument list - AML, May 1, 1998 (sunpos_b.c)
 *     - to compute Earth's orbital velocity (rather approximate)
 *       and to change to J2000 coordinates  - AML, August 26, 1998 (sunpos_c.c)
 *     - to put rsun[3] in AU - AML, August 14, 2000
 */

void sunpos(double mjd, double possun[2], double rsun[3], double vearth[3])
{
  double n, l, g, eps, lambda, ra, dec, r;
  int ilamb;

  n = mjd + MjdOff - 2451545.0;
  l = 280.466 + 0.9856474 * n;
  g = 357.528 + 0.9856003 * n;
  eps = 23.44 - 0.0000004 * n;
  lambda = l + 1.915 * sin(g*d2r) + 0.02 * sin(2.0*g*d2r);

/* A tiny correction term for conversion to epoch J2000 */

  lambda -= 1.3969713 * n/36525.0;

  ilamb = (int) (lambda/360.0);
  lambda -= 360.0*ilamb;
  if (lambda < 0.0) lambda += 360.0;

  r = 1.00014 - 0.01671*cos(g*d2r) - 0.00014*cos(2.0*g*d2r); /* AU */

  ra = cos(eps*d2r) * tan(lambda*d2r);
  ra = atan(ra);
  if ((lambda > 90) && (lambda < 270))
    ra += pi;
  if (lambda >= 270)
    ra += 2.0 * pi;
  dec = sin(eps*d2r) * sin(lambda*d2r);
  dec = asin(dec);

  possun[0] = ra * r2d;
  possun[1] = dec * r2d;
 
  /* xyz vector (AU) available if desired */
  rsun[0] = r*cos(dec)*cos(ra);
  rsun[1] = r*cos(dec)*sin(ra);
  rsun[2] = r*sin(dec);

  /* Velocity vector for circular orbit approximation - good to a few percent */
  vearth[0] = AU_DAY*0.0172*sin(d2r*lambda);
  vearth[1] = AU_DAY*(-0.0158)*cos(d2r*lambda);
  vearth[2] = AU_DAY*(-0.0068)*cos(d2r*lambda);
}

/******************************************************************************/
/* Precession routine
 *          - precess equatorial coordinates given for a coordinate system for
 *             date mjd_from (MJD) to a coordinate system for date mjd_to.
 *          - To get J2000 coordinates, mjd_to = 51544.5
 *          - radec[2] is ra, dec (degrees)
 *
 */

void precess_radec(double mjd_from, double radec_f[2], double mjd_to, double radec_t[2])
{
  double p[3][3], pinv[3][3], pfull[3][3];
  double rf[3], rt[3];

  dsphcr(d2r*radec_f[0],d2r*radec_f[1],rf);
  precess_mat(mjd_from,p);
  trans(p,pinv);
  precess_mat(mjd_to,p);
  matmat(p,pinv,pfull);
  matvec(pfull,rf,rt);
  dnorm(rt);
  dcrsph(rt,&radec_t[0],&radec_t[1]);
  radec_t[0] /= d2r;
  radec_t[1] /= d2r; 
}

/******************************************************************************/
/* Precess matrix calculation
 *          - from Astronomical Almanac, Year 2000, p. B18
 *            Formulas have been truncated to give ~1 arcsec accuracy for
 *            < 100 yrs from J2000 epoch.
 */

#define MJD_T0 51544.5       /* 2000 Jan. 1.5 */

void precess_mat(double mjd, double pmat[3][3])
{
  double tcap;
  double zeta, z, theta;
  double coszeta, sinzeta, cosz, sinz, costh, sinth;

  tcap = (mjd - MJD_T0)/36525.0;

  zeta = d2r*(0.6406161*tcap + 0.0000839*tcap*tcap + 0.000005*tcap*tcap*tcap);
  z = d2r*(0.6406161*tcap + 0.0003041*tcap*tcap + 0.0000051*tcap*tcap*tcap);
  theta = d2r*(0.5567530*tcap - 0.0001185*tcap*tcap - 0.0000116*tcap*tcap*tcap);

  coszeta = cos(zeta);
  sinzeta = sin(zeta);
  cosz = cos(z);
  sinz = sin(z);
  costh = cos(theta);
  sinth = sin(theta);

  pmat[0][0] = cosz*costh*coszeta - sinz*sinzeta;
  pmat[0][1] = -cosz*costh*sinzeta - sinz*coszeta;
  pmat[0][2] = -cosz*sinth;
  pmat[1][0] = sinz*costh*coszeta + cosz*sinzeta;
  pmat[1][1] = -sinz*costh*sinzeta + cosz*coszeta;
  pmat[1][2] = -sinz*sinth;
  pmat[2][0] = sinth*coszeta;
  pmat[2][1] = -sinth*sinzeta;
  pmat[2][2] = costh;
}

/******************************************************************************/
/* Convert equatorial celestial coordinates (J2000) to Galactic coordinates
 *  (l_II, b_II)
 *
 * Angles should be in degrees
 */

#define MJD_B1950 33281.923   /* B1950 = JD 2433282.423 */

void celgal(double ra, double dec, double *glong, double *glat)
{
  double rcel[3], rgal[3], rmat[3][3];
  double radec2000[2], radec1950[2];
  double glo, gla;

  radec2000[0] = ra;
  radec2000[1] = dec;
  precess_radec(MJD_T0,radec2000,MJD_B1950,radec1950);
  dsphcr(d2r*radec1950[0],d2r*radec1950[1],rcel);
  celgalmat(rmat);
  matvec(rmat,rcel,rgal);
  dcrsph(rgal,&glo,&gla);
  *glong = glo/d2r;
  *glat = gla/d2r;

}

/******************************************************************************/
/* This function computes the direction cosine matrix for an SSC field of
 * view specified by the RA and Dec of the center of the FOV and the position
 * angle of the long direction (all angles in degrees).
 *
 * The returned matrix multiplies celestial-coordinate vectors to give
 * vectors w/ SSC coordinates.
 *
 */

void celgalmat(double rmat[3][3])
{
  double euler[3];

  euler[0] = 282.25*d2r;
  euler[1] = 62.6*d2r;
  euler[2] = -33.0*d2r;
  eulerm313(euler,rmat);
  return;
}
