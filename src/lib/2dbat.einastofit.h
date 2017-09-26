#ifndef __2DBAT_EINASTOFIT_H__
#define __2DBAT_EINASTOFIT_H__

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
#include "2dbat.gsl.h"
#include "2dbat.sort.h"
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
// Einasto halo fit related

// vEinasto multinest 1D fit
void v_einasto_1d_multinest4p_remove_outliers(TR_ringParameters *TRparam);
//
void v_einasto_1d_multinest4p(multinest_paramters *multinest_param, TR_ringParameters *TRparam);
//
void loglikelihood_vEinasto4p_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
//
void dumper_vEinasto4p_student(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam);
//
double vEinasto3p(double *Cube, double r);
//
void spline_intp(TR_ringParameters *TRparam);
//
double spline_intp_ringparam(TR_ringParameters *TRparam, char *param, double ring);

// vEinasto multinest 2p fit
void vEinasto_multinest2p(multinest_paramters *multinest_param, TR_ringParameters *TRparam);
//
void loglikelihood_vEinasto2p(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
//
void dumper_vEinasto2p(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam);
//
double vEinasto2p(double *Cube, double r, double r_2);
//
void estimate_vrot_error(multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, char *finalfit, char final_fit);
//
// Einasto multinest 2D fit
void Einasto_haloFits_multinest(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *sersic_function, char sersicpart, char *_n_Einasto, char _n_fix, char *r_2_Einasto, char r_2_fix, char *rho_2_Einasto, char rho_2_fix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int rank, char *mt_outputfile);
//
void loglikelihood_einasto_halofit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
//
int einasto_halomodel(double *Cube, int i, int j, TR_ringParameters *TRparam, double *Vmodel_Einasto);
//
void einasto_halofit_multinest_student(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *_n_Einasto, char _n_fix, char *r_2_Einasto, char r_2_fix, char *rho_2_Einasto, char rho_2_fix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int rank, char *mt_outputfile);
//
void loglikelihood_einasto_halofit_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
//
void dumper_einasto_halofit(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam);
void set_nfree_params_einasto_halofit_multinest_student(char xposfix, char yposfix, char vsysfix, char pafix, char inclfix, char _n_fix, char r_2_fix, char rho_2_fix, char vradfix, char sigmafactorfix, TR_ringParameters *TRparam);
//
double radius_galaxy_plane(double radius, double x, double y, double pa, double incl);
//
void init_einastohalo_vrot_intp(TR_ringParameters *TRparam);
//
double pdiGamma_pdn_gsl_numerical_integral(double a, double b); // from a to b
//
double read_einasto_posteriors_and_calculate_cov(TR_ringParameters *TRparam, char *param1, char *param2);
//
double lower_inc_gamma_variant(double x, void * params);

// --- End of line

#endif
