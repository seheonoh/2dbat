#ifndef __2DBAT_MPI_H__
#define __2DBAT_MPI_H__

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
#include "2dbat.priors.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>

// 2DBAT user defined functions
// MPI related

// sending or receiving MPI data 
void send_mpi_data(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi,  MPI_Status *mpistatus);
//
void receive_mpiData(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi, velocity_field *HI_VF_boxFiltered, velocity_field *HI_VF_boxFiltered_sigma, velocity_field *HI_VF_boxFiltered_sigma_e_norm, velocity_field *HI_VF_mom2, velocity_field *HI_VF_boxFiltered_decim0, velocity_field *HI_VF, MPI_Status *mpistatus);
//
void receive_mpiData_TRparam(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi, velocity_field *HI_VF_boxFiltered, velocity_field *HI_VF_boxFiltered_sigma, velocity_field *HI_VF_boxFiltered_sigma_e_norm, velocity_field *HI_VF_mom2, velocity_field *HI_VF_boxFiltered_decim0, velocity_field *HI_VF, MPI_Status *mpistatus);
//
MPI_Datatype allocate_mpi_dataset(TR_ringParameters *TRparam, MPI_Datatype *type, int *blocklen);


// --- End of line

#endif
