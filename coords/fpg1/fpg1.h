/* fpg1.h
 * Alan M. Levine
 * June 7, 2017
 *
 * Header file for fpg1.c
 */

extern void sky_to_sc_mat();
extern void sc_to_cam_mat(double euler[3]);
extern void prmat(double rm[3][3], FILE *fp);
extern void optics_fp(double lng_deg, double lat_deg, double xyfp[2]);
extern void get_ra_dec_roll(double rm[3][3], double angs[3]);
extern void xyrotate(double angle_deg, double xyin[2], double xyout[2]);
extern int mm_to_pix(double xy[2], double ccdpx[2], double fitpx[2]);
extern void set_fits_pix_lims();
extern int star_in_border(int iccd, double row, double col, double pxbor);
extern int get_ccdno(int istar);
extern void from_fits_pix_t(double xy[2], double fitp[2]);
extern void setup_mats(int ipr);
extern int star_pixels(int i, double ccdpx[2], double fitpx[2], 
			int icode, int ideriv);
extern void all_star_pixels();
extern int one_star_fit_inputs(int istar, double fitpx[2]);
extern int code_param(int itype, int iicam, int iccd, int inum);
extern void decode_param(int icode, int idec[4]);
extern void fill_fpg_param_arrays(int iarr);
extern void fill_mrq_param_arrays();
extern void read_fpg_pars(FILE *fpin, FILE *fpo, int iicam);
extern void print_fpg_pars_deriv(FILE *fpo);
extern void delta_for_deriv(double dxy_mm, int npixccd);
extern void get_perp(double s[3], double v[3], double vp[3]);
extern void aber_star_pos(double v[3], double ra, double dec, double *raap,
			  double *decap);
extern void read_fpg_stars(FILE *fpin, FILE *fpd, double pxbor);
extern void print_fpg_star(int i, FILE *fpo, FILE *fpd, double sx, double sy);
extern void fill_mrq_stars();

extern double ra_sc, dec_sc, roll_sc, radecroll[3], sc_vel[3];
extern double dtor;
