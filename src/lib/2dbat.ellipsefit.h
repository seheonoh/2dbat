#ifndef __2DBAT_ELLIPSEFIT_H__
#define __2DBAT_ELLIPSEFIT_H__

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

// ellipse fit
double ellipse_function_rect(double a, double b, double x);
double ellipse_function_polar(double a, double e, double theta);

// Ellipse fitting and set initial uniform priors
void do_ellipseFit_update_uniformPriors_ringParam(TR_ringParameters *TRparam);

// gipsy ellipse fit1
int ellipse1_c(double *n ,      /* number of points */
                    float *x ,      /* X coordinates */
                    float *y ,      /* Y coordinates */
                    float *p );     /* ellipse parameters */

// gipsy ellipse fit2
int ellipse2_c(double *n ,      /* number of coordinates */
                    float *x ,      /* X coordinates */
                    float *y ,      /* Y coordinates */
                    float *p ,      /* ellipse parameters */
                    float *e ,      /* errors in ellipse parms. */
                    int *m);        /* fixed/free mask */

int invmat( double matrix[PARS][PARS], int nfree );

double inimat( double   s[PARS][PARS] ,
                        double  rl[PARS] ,
                        float   *x ,
                        float   *y ,
                        int n ,
                        float   *p ,
                        float   *e ,
                        int ip[PARS] ,
                        int nfree );

int inivec( double  s[PARS][PARS] ,
                        double  s1[PARS][PARS] ,
                        double  rl[PARS] ,
                        double  labda ,
                        double  *q ,
                        float   *p ,
                        float   *e ,
                        float   *x ,
                        float   *y ,
                        int n ,
                        int ip[PARS] ,
                        int nfree );

// --- End of line

#endif

