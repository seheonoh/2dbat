#ifdef __INTEL_COMPILER             // if the MultiNest library was compiled with ifort
       #define NESTRUN nested_mp_nestrun_
#elif defined __GNUC__              // if the MultiNest library was compiled with gfortran
       #define NESTRUN __nested_MOD_nestrun
#else
       #error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libnest3.a | grep nestrun) & edit example_eggbox_C++/eggbox.cc
#endif

#ifndef MULTINEST_H
#define MULTINEST_H

#include <string.h>

/***************************************** C Interface to MultiNest **************************************************/
/* Uniform priors of ring parameters */
typedef struct {
    // Process 
    char wdir[1000]; // workding directory
    // einasto output file
    char txtfile_multinest_output[1000];
    int loop_check;
    int n_freeParams; // number of free parametres for either TR or ISO fits
    double e_sigma; // sigma of error distribution in km/s
    double e_sigma_fitted;
    double e_sigma_fitted_t;
    // Parameters for multinest
    double maxLogLike[999];
    double logZ[999];
    double logZerr[999];
    double maxLogLikeF;
    double logZF;
    double logZerrF;

    // VF fields
    int nax1;
    int nax2;
    int decim_x0;
    int decim_y0;
    double decimX;
    double decimY;
    double decimX_einasto_halofit;
    double decimY_einasto_halofit;
    double decimX_trfit;
    double decimY_trfit;
    float pixelScale; // arcsec/pixel

    // simple ellipse fit parametres
    double N_all_pixels_in_a_ring;
    double ellipse_xpos_boxfiltered;
    double ellipse_ypos_boxfiltered;
    double ellipse_pa_boxfiltered;
    double pa_EllipseFit_e; // GIPSY ellipse fit error
    double ellipse_incl_boxfiltered;
    double incl_EllipseFit_e; // GIPSY ellipse fit error
    double ellipse_semi_mx_boxfiltered;
    double xpos1_from_ellipsefit;
    double xpos2_from_ellipsefit;
    double ypos1_from_ellipsefit;
    double ypos2_from_ellipsefit;
    double vsys1_from_vlosfit;
    double vsys2_from_vlosfit;
    double perimeter;

// 39
    // set rings
    char final_fit; // trfit
    char fullFit; // einasto fit
    int Nrings;
    int Nrings_intp;
    int Nrings_to_semi_mx;
    int N_reliable_rings;
    int tilted_ring[4024*4024][2]; // (i, j) in a given ring
    double Npoints_in_tilted_ring;
    double total_Npoints_allRings;
    double ring_radius[999];
    double ring_intp[999];
    double ring_s;
    double ring_e;
    double ring_w;
    double ring_s_for_einasto1Dfit;
    double ring_e_for_einasto1Dfit;
    double ring_w_for_einasto1Dfit;
    double Npoints_in_tilted_ring_decim0;
    double npoints_inaring[999999];
    double npoints_inaring_decim0[999999];
    // radius
    double rGalaxyPlane_pixel;
    double rGalaxyPlane_pixel_max; // for R normalisation
    // weight
    double free_angle; // free angle with which velocities are discarded
    double wpow; // weighting order for cosine
    double rwpow; // radial weighting order
// 64

    // XPOS
    char xpos_fix; // "Y" (fix) or "N" (free)
    double xposF; // finally derived XPOS
    double xposF_e; // finally derived XPOS error
    double xpos0[999]; // currently derived XPOS
    double xpos[999];
    double xpos_e[999];
    double xpos_temp[999];
    double xpos_mode;
    double xpos_std;
    double xpos1; // x1
    double xpos2; // x2
    double xposF_EinastoFit; // finally derived xpos using ISO profile
    double xposF_EinastoFit_e; // finally derivex xpos_e using ISO profile
    double xposF_EinastoFit_t; // currently derived xpos using ISO profile
    double xposF_EinastoFit_e_t; // currently derivex xpos_e using ISO profile
// 79
    // YPOS
    char ypos_fix; // "Y" (fix) or "N" (free)
    double yposF; // finally derived XPOS
    double yposF_e; // finally derived XPOS error
    double ypos0[999]; // currently derived XPOS
    double ypos[999];
    double ypos_e[999];
    double ypos_temp[999];
    double ypos_mode;
    double ypos_std;
    double ypos1; // y1
    double ypos2; // y2
    double yposF_EinastoFit; // finally derived ypos using ISO profile
    double yposF_EinastoFit_e; // finally derivex ypos_e using ISO profile
    double yposF_EinastoFit_t; // currently derived ypos using ISO profile
    double yposF_EinastoFit_e_t; // currently derivex ypos_e using ISO profile
//94
    // VSYS
    char vsys_fix; // "Y" (fix) or "N" (free)
    double vsysF; // finally derived XPOS
    double vsysF_e; // finally derived XPOS error
    double vsys0[999]; // currently derived XPOS
    double vsys[999];
    double vsys_e[999];
    double vsys_temp[999];
    double vsys_mode;
    double vsys_std;
    double vsys1; // y1
    double vsys2; // y2
    double vsysF_EinastoFit; // finally derived vsys using ISO profile
    double vsysF_EinastoFit_e; // finally derivex vsys_e using ISO profile
    double vsysF_EinastoFit_t; // currently derived vsys using ISO profile
    double vsysF_EinastoFit_e_t; // currently derivex vsys_e using ISO profile
//109
    // PA
    char pa_fix; // "Y" (fix) or "N" (free)
    double paF; // finally derived PA
    double paF_e; // finally derived PA error
    double pa0[999]; // currently derived PA
    double pa[999];
    double pa_e[999];
    double pa_temp[999];
    double pa_temp_e[999];
    double pa1; // x1
    double pa2; // x2
    double pa1_for_TRfit; // x1
    double pa2_for_TRfit; // x2
    double _p1_tr[999];
    double _p2_tr[999];
    double PA_MAX_in_degree; // maximum PA value used for normalisation : 360
    double paF_EinastoFit; // finally derivex xpos using Burkert profile
    double paF_EinastoFit_e; // finally derivex xpos using Burkert profile
    // bspline
    char pa_function[13]; // PA function: bspline
    int n_coeffs_bspline_pa;
    // PA bspline nbreak : 0: constant
    int pa_nbreak_bspline;
    // PA bspline order : linear=1, quadrature=2, cubic=3
    int pa_order_bspline;
    double _p_bs[99];  // spline coefficients of PA
    double _p_bs_e[99];  // spline coefficients of PA
    double _p_bs_t[99]; // 
    double _p_bs_e_t[99]; // 
    double bspline1pa[99];  // s1: uniform prior for spline coefficients of PA
    double bspline2pa[99];  // s1: uniform prior for spline coefficients of PA
    double _p_bs_tr[99];
    double _bspline_pa_hist_sigma;
// 138
    // INCL
    char incl_fix; // "Y" (fix) or "N" (free)
    double inclF; // finally derived INCL
    double inclF_e; // finally derived INCL error
    double incl0[999]; // currently derived INCL
    double incl[999];
    double incl_e[999];
    double incl_temp[999];
    double incl_temp_e[999];
    double incl1; // x1
    double incl2; // x2
    double incl1_for_TRfit; // x1
    double incl2_for_TRfit; // x2
    double _i1_tr[999];
    double _i2_tr[999];
    double INCL_MAX_in_degree; // maximum INCL value used for normalisation : 90
    double inclF_EinastoFit; // finally derivex xpos using Burkert profile
    double inclF_EinastoFit_e; // finally derivex xpos using Burkert profile
// 155
    // bspline
    char incl_function[13]; // INCL function: bspline
    int n_coeffs_bspline_incl;
    // INCL bspline nbreak : 0: constant
    int incl_nbreak_bspline;
    // INCL bspline order : linear=1, quadrature=2, cubic=3
    int incl_order_bspline;
    double _i_bs[99];  // spline coefficients of INCL
    double _i_bs_e[99];  // spline coefficients of INCL
    double _i_bs_t[99]; // 
    double _i_bs_e_t[99]; // 
    double bspline1incl[99];  // s1: uniform prior for spline coefficients of INCL
    double bspline2incl[99];  // s1: uniform prior for spline coefficients of INCL
    double _i_bs_tr[99];
    double _bspline_incl_hist_sigma;
// 167
    // VROT
    char vrot_fix; // "Y" (fix) or "N" (free)
    double vrotF; // finally derived VROT
    double vrotF_e; // finally derived VROT
    double vrot1;
    double vrot2;
    double vrot0[999]; 
    double vrot[999];
    double vrot_e[999];
    double vrot_rec[999];
    double vrot_e_rec[999];
    double vrot_app[999];
    double vrot_e_app[999];
    double vrot_temp[999];
    double vrot_temp_e[999];
    double vrot_intp[999];
    double vrot_e_intp[999];
// 183k
    // VRAD
    char vrad_fix;
    double vrad1;
    double vrad2;
    double vradF;
    double vradF_e;
    double vrad_max;
    double vrad[999];
    double vrad_e[999];
    double vrad_rec[999];
    double vrad_rec_e[999];
    double vrad_app[999];
    double vrad_app_e[999];
    char vrad_function[13];
    int n_coeffs_bspline_vrad;
    int vrad_nbreak_bspline;
    int vrad_order_bspline;
    double bspline1vrad[99];
    double bspline2vrad[99];
    double _vr_bs[99];
    double _vr_bs_e[99];
    double _vr_bs_t[99];
    double _vr_bs_e_t[99];
    double _vr_bs_tr[99];
    double vrad_temp[999];
    double vrad_temp_e[999];
    double vrad_einastofit_bs[999];
    double vrad_einastofit_bs_e[999];
    double _bspline_vrad_hist_sigma;
// 211
    // parameters for Gaussian profile fit
    int n_gauss; // Number of gaussian
    double g_param[99]; // best fit parameters for gaussian function
    double g01;
    double g02;
    double gA1;
    double gA2;
    double gS1;
    double gS2;
    double gX1;
    double gX2;
// 221
    // line-of-sight velocities
    double LOS_hist_Gfit_V0;
    double LOS_hist_Gfit_sigma;
    double LOS_vel_hist_rbm;
    double LOS_vel_hist_std;
    double vlos_lower;
    double vlos_upper;
    double vlos_lower_limit;
    double vlos_upper_limit;
    double vel_geo_app_side;
    double vel_geo_rec_side;
    double dispersion_VLOS;
// 232
    // sigma_factor
    char sigma_factor_fix;
    double sigma_factor;
    double sigma_factor_e;
    double sigma_factor1;
    double sigma_factor2;
    double mean_mom4;
    double std_mom4;
    double sigma_factor_mode;
    double sigma_factor_std;
    double e_sigma_student_TR[999];
    double scale_factor_const_vlose_w;
    double scale_factor_var_vlose_w;
    double scale_factor_mom0weightedvf_e_to_mom2;
    double hist_mean_vf_e_mom0_weighted_scaled;
    double hist_std_vf_e_mom0_weighted_scaled;
    double einastofit_BIC;
// 248
    // median box filter
    int box_x;
    int box_y;

    // Einasto profile
    // n
    char _n_fix;
    double _n;
    double _ne;
    double _n1;
    double _n2;
    // r_2
    char r_2_fix;
    double r_2;
    double r_2e;
    double r_21;
    double r_22;
    // rho_2
    char rho_2_fix;
    double rho_2;
    double rho_2e;
    double rho_21;
    double rho_22;
// 265
    double _n_t;
    double _r_2_t;
    double _rho_2_t;

    double _ne_t;
    double _r_2e_t;
    double _rho_2e_t;

    double einasto1D_ne;
    double einasto1D_r2e;
    double einasto1D_rho2e;

    // Einasto profile : temporary  
    double _n1_t;
    double _n2_t;
    double r_21_t;
    double r_22_t;
    double rho_21_t;
    double rho_22_t;
// 280
    // Einasto halo params limits
    double Ein_n_min;
    double Ein_n_max;
    double Ein_r_2_min;
    double Ein_r_2_max;
    double Ein_rho_2_min;
    double Ein_rho_2_max;

    int n_hist_post;
    double hist_x[9999]; // hist_x of posterior
    double hist_y[9999]; // hist_y of posterior
    double hist_ye[9999]; // hist_ye of posterior
// 290

    // being added...
    double decimX_einasto_halofit_d;
    double decimY_einasto_halofit_d;
    int use_allPixels;  
// 293

    double vrot_einasto_error[999];
    double vrot_asymmetry_error[999];
    double vrot_dispersion_error[999];
// 296

    double _nu_studenT;
// 297

    double e_sigma_tr;
    double e_sigma_e_tr;
// 299

    double vf_e_user;
    double vdisp_user;
// 301
    int _nfilter;
// 302

} TR_ringParameters;

extern void NESTRUN(int *, int *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *, 
char *, int *, int *, int *, int *, int *, int *, double *, int *, void (*Loglike)(double *, int *, int *, 
double *, TR_ringParameters *), void (*dumper)(int *, int *, int *, double **, double **, double **, double *, 
double *, double *, double *, TR_ringParameters *), TR_ringParameters *context);

void run(int IS, int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, 
int maxModes, int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, 
int initMPI, double logZero, int maxiter, void (*LogLike)(double *, int *, int *, double *, TR_ringParameters *), 
void (*dumper)(int *, int *, int *, double **, double **, double **, double *, double *, double *, double *, TR_ringParameters *), 
TR_ringParameters *context);


/***********************************************************************************************************************/

#endif // MULTINEST_H

