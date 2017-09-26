#ifndef __2DBAT_TRFIT_H__
#define __2DBAT_TRFIT_H__

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
#include "2dbat.priors.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>

// 2DBAT user defined functions
// Tilted-ring fit related

// 2DBAT global variables : defined in 2dbat.global_params.h
TR_ringParameters *TRparam;
multinest_paramters *multinest_param;
filename_2dbat *fname;


// set rings
void set_Nrings(TR_ringParameters *TRparam, double *Nall_ring, double *Navail_ring);

// print TR fit results
void print_trfit_results(TR_ringParameters *TRparam);

// check ring outlier
int ringOutliered(TR_ringParameters *TRparam, int ringN, double xposL, double xposU, double xpos_eL, double xpos_eU,
                           double yposL, double yposU, double ypos_eL, double ypos_eU,
                           double vsysL, double vsysU, double vsys_eL, double vsys_eU,
                           double paL, double paU, double pa_eL, double pa_eU,
                           double inclL, double inclU, double incl_eL, double incl_eU,
                           double vrotL, double vrotU, double vrot_eL, double vrot_eU,
                           double vradL, double vradU, double vrad_eL, double vrad_eU);

// TR fitting multinest
void loglikelihood_trfit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
//
void loglikelihood_trfit_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
//
int TiltedRingModel(double *Cube, int i, int j, TR_ringParameters *TRparam, double *Vmodel_TR); // Update TRmodels
//
void dumper_TRfits(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam);
//
void trfit_multinest_ellipsefit_rings_student(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *vrot, char vrotfix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, int find_symmetric_VF, char *finalfit, char final_fit);
//
void set_nfree_params_trfit_multinest_ellipsefit_rings_student(char xposfix, char yposfix, char vsysfix, char pafix, char inclfix, char vrotfix, char vradfix, char sigmafactorfix, TR_ringParameters *TRparam);
//
void define_tiltedRing(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side);
//
void define_tiltedRing_ellipse_rings(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side);
//
void define_tiltedRing_ellipse_rings_extend_rings(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side);
//
void find_Navail_Nall_pixels(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side, double *Nall_ring, double *Navail_ring);
//
void trfit_multinest_trfit_rings_normal(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *vrot, char vrotfix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, char *finalfit, char final_fit);
//
void trfit_multinest_trfit_rings_student(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *vrot, char vrotfix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, char *finalfit, char final_fit);
//
void set_nfree_params_trfit_multinest_trfit_rings_student(TR_ringParameters *TRparam);
void set_nfree_params_trfit_multinest_trfit_rings_student_test(TR_ringParameters *TRparam);
//
//void define_area_tofit(double ring0, double ring1, int decimX, int decimY, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side);

void define_area_tofit(double ring0, double ring1, int decimX, int decimY, multinest_paramters *multinest_param, TR_ringParameters *TRparam,  int side);
//
void define_area_tofit_set_geo_weight(double ring0, double ring1, int decimX, int decimY, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, int histogram);
//
void set_navail_nall_pixels(TR_ringParameters *TRparam, \
                            double *Nall_ring, double *Navail_ring, \
                            char *xpos, char xposfix,   \
                            char *ypos, char yposfix,   \
                            char *vsys, char vsysfix,   \
                            char *pa, char pafix,   \
                            char *incl, char inclfix,   \
                            char *_n, char _nfix,   \
                            char *_r_2, char _r_2fix,   \
                            char *_rho_2, char _rho_2fix,   \
                            char *fullfit, char fullfit_fix);



// --- End of line

#endif


