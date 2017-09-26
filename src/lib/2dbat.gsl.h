#ifndef __2DBAT_GSL_H__
#define __2DBAT_GSL_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "2dbat.cfitsio.h"
#include "2dbat.multinest.h"
#include "2dbat.trfit.h"
#include "2dbat.sort.h"
#include "2dbat.einastofit.h"
#include "2dbat.global_params.h"
#include "2dbat.2dmaps.h"
#include "2dbat.ellipsefit.h"
#include "2dbat.etc.h"
#include "2dbat.gfit.h"
#include "2dbat.memory.h"
#include "2dbat.mpi.h"
#include "2dbat.priors.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>

// 2DBAT user defined functions
// GSL related

// GSL histogram
void robust_mean_std(double *input, int n, double *robust_mean, double *robust_std);
void robust_mean_std_e(double *input, int n, double *robust_mean, double *robust_std);
void robust_mean_std_histogram(double *input, int n, double *robust_mean, double *robust_std);
void robust_mean_std_histogram_ac(double *input, double *input_err, int n, double *robust_mean_ac, double *robust_std_ac);
int sshist(double *x, int n_data);

// GSL nonlinear solver: approximation
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver(int i, int j, TR_ringParameters *TRparam, double x_lo, double x_hi, int max_iter, char *solver);
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver_given_paincl(int i, int j, double _xpos, double _ypos, double _pa, double _incl, TR_ringParameters *TRparam, double x_lo, double x_hi, int max_iter, char *solver);
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation(double x, void *params);
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_given_paincl(double x, void *params);

// GSL bspline
int gsl_nonlinearfit_Bspline_filtering(char *param, TR_ringParameters *TRparam);
int bsplinefit_set_unipriors(char *param, TR_ringParameters *TRparam);

// derive r_ij used for computing weight
double r_ij_pa_incl_bspline_W(int i, int j, TR_ringParameters *TRparam);
double r_ij_pa_incl_given_W(int i, int j, double _xpos, double _ypos, double _pa, double _incl, TR_ringParameters *TRparam);

// calculations
double getmedian_g(float *calcarray, int len);
void sortra_c(float *x, int *n);
void get_rbm(float *array, int arraySize, float filter_lower, float filter_upper, float *rbm_out, float *std_out);
void get_rbm_double(double *array, int arraySize, double filter_lower, double filter_upper, double *rbm_out, double *std_out);
double getMedian(double* array, size_t arraySize);
void get_wMean_STD(double* array, double* array_w, int arraySize, double *mean_array, double *std_array);

// Nonlinear-equation solver: Brent method
double zero_rc_one(double a, double b, double machep, double t, double f(double XPOS, double YPOS, double i, double j, double *_p, double *_i, double k, double alpha, double n, double rGalaxyPlane_pixel, double rGalaxyPlane_pixel_max, TR_ringParameters *TRparam), double XPOS, double YPOS, double i, double j, double *_p, double *_i, double k, double alpha, double n, double rGalaxyPlane_pixel_max, TR_ringParameters *TRparam);
void zero_rc(double a, double b, double t, double *arg, int *status, double value);
double r8_abs(double x);
double r8_epsilon(void);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// --- End of line

#endif


