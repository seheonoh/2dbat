#ifndef __2DBAT_GFIT_H__
#define __2DBAT_GFIT_H__

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
#include "2dbat.einastofit.h"
#include "2dbat.global_params.h"
#include "2dbat.2dmaps.h"
#include "2dbat.ellipsefit.h"
#include "2dbat.etc.h"
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
// Gaussian fit related

// Gfit multinest 1D fit
void Gfit_multinest(multinest_paramters *multinest_param, TR_ringParameters *TRparam);
void loglikelihood_gfit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
void dumper_Gfits(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam);
double gauss_function(double *Cube, int n_gauss, double x);

// gsl Gaussian fit
int gsl_nonlinearfit_gauss(TR_ringParameters *TRparam, int nn, double g0_init, double gA_init, double gX_init, double gS_init, double *g0, double *gA, double *gX, double *gS, int *result);
int gauss_f(const gsl_vector * x, void *fourier_inputdata, gsl_vector * f);
int gauss_df(const gsl_vector * x, void *fourier_inputdata, gsl_matrix * J);
int gauss_fdf(const gsl_vector * x, void *fourier_inputdata, gsl_vector * f, gsl_matrix * J);

// --- End of line

#endif

