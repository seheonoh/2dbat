#ifndef __2DBAT_PRIORS_H__
#define __2DBAT_PRIORS_H__

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
#include "2dbat.gfit.h"
#include "2dbat.memory.h"
#include "2dbat.mpi.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>

//++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++
// USER DEFINED FUNCTIONS
//++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++

// multinest Gaussian + uniform priors
double dierfc(double y);
double JeffreysPrior(double r, double x1, double x2);
double gaussian_prior(double r, double mu, double sigma);
double gaussian_prior_paincl_skew(double r, double mu, double sigma_N, double sigma_W);
double gaussian_prior_nrrho_skew(double r, double mu, double sigma_N, double sigma_W);
double gaussian_prior_esigma_skew(double r, double mu, double mu_max, double sigma_N, double sigma_W);
double uniform_priors(double r, double x1, double x2);
double DeltaFunctionPrior(double r, double x1);
double loguniform_priors(double r, double x1, double x2);
double logNormalPrior(double r, double a, double sigma);



// --- End of line

#endif


