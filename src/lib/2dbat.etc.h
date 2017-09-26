#ifndef __2DBAT_ETC_h__
#define __2DBAT_ETC_h__

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
// ETC

// count the number of lines in a file
int countlines(FILE *ifp1);

// read multinest output file
void add_MAP_line_multinest_output(char *multinest_output, TR_ringParameters *TRparam);


// update params
void update_Einasto_params_priors(TR_ringParameters *TRparam, float Nsigma);
void update_multinest_params_dirtyfit(multinest_paramters *multinest_param);
void update_multinest_params_fullfit(multinest_paramters *multinest_param);
void update_vrot_prior(TR_ringParameters *TRparam);

int write_trfit_results(TR_ringParameters *TRparam, char *bayesFit_outputfile);
// read user input 
void read_user_input_init_params(TR_ringParameters *TRparam, multinest_paramters *multinest_param, filename_2dbat *fname, int argc, char *argv[]);
void usage_2dbat();
//
void print_priors_info(char *mode, TR_ringParameters *TRparam, int n_node, int iter_find_optimal_priors);

// --- End of line

#endif
