#ifndef __2DBAT_GLOBAL_PARAMS_H__
#define __2DBAT_GLOBAL_PARAMS_H__

// GIPSY ellipse fit parametres
#define cosd( x )   cos( x * 0.017453293 )
#define sind( x )   sin( x * 0.017453293 )
#define fmake( f, c )   {f.a=c;f.l=sizeof(c);}
#define SWAP(a,b) {float swap=(a);(a)=(b);(b)=swap;}
#define T   1000                    /* maximum iteration */
#define F   0.0174532925            /* conversion factor */
#define FAC 10.0                /* mixing * factor */
#define LAB 0.01                /* start mixing parameter */
#define LABMAX  1.0E+10             /* max. mixing parameter */
#define LABMIN  1.0E-10             /* min. mixing parameter */
#define TOL 0.0001              /* tolerance */
#define PARS 5

struct fourier_inputdata {
    int poly_order;
    size_t n;
    double *ring_radius;
    double *y;
    double *sigma;
    double ring_s;
    double ring_w;
    double rGalaxyPlane_pixel_max;
};

struct rGalaxyPlane_pixel_TR_nonlinearEquation_params {
    char pa_function[13];
    char incl_function[13];
    double PA_MAX_in_degree;
    double INCL_MAX_in_degree;
    double ring_radius[9999];
    int N_reliable_rings;
    double pa_temp[9999];
    double pa_temp_e[9999];
    double incl_temp[9999];
    double incl_temp_e[9999];
    double rGalaxyPlane_pixel_max;

    int PA_fourier_order;
    int INCL_fourier_order;
    int pa_order;
    int incl_order;
    int SersicPoly_order;
    int i;
    int j;
    double xpos;
    double ypos;
    double pa;
    double incl;
    //fourier
    double _p0;
    double _p_a[99];
    double _p_b[99];
    double _p_w;
    double _i0;
    double _i_a[99];
    double _i_b[99];
    double _i_w;
    // polynomial
    double _p[99];
    double _i[99];
    // bspline
    int n_coeffs_bspline_pa;
    int pa_nbreak_bspline;
    int n_coeffs_bspline_incl;
    int incl_nbreak_bspline;
    double _p_bs[99];
    double _i_bs[99];
    // sersic
    double kapa;
    double alpha;
    double sersic_n;
    double rmax;
};

struct rGalaxyPlane_pixel_TR_nonlinearEquation_Poly_params {
    int pa_order;
    int SersicPoly_order;
    int i;
    int j;
    double xpos;
    double ypos;
    double _p[99];
    double _i[99];
    double kapa;
    double alpha;
    double sersic_n;
    double rmax;
};


typedef struct {
    float a;
    float b;
    float e;
    float xpos;
    float ypos;
} Ellipse_Parameter;

typedef struct {
    float XPOS;
    float YPOS;
    float VSYS;
    float PA;
    float INCL;
    float radius;
    float width;
} ring_parameter;

typedef struct {
    char txtfile_multinest_output[1000];
    char fitsfile_2Dinput_VF[1000];
    char fitsfile_2Dinput_VF_e[1000];

    char fitsfile_mom0[1000];
    char fitsfile_mom2[1000];
    char pa_function_user[13];
    char incl_function_user[13];
    char vrad_function_user[13];

    char bayesFit_outputfile[1000];
    char bayesFit_outputfile_allfree[1000];
    char fitsfile_2Dinput_VF_boxfiltered[1000];
    char fitsfile_2Dinput_VF_boxfiltered_sigma[1000];
    char fitsfile_2Dinput_VF_boxfiltered_sigma_ew[1000];
    char fitsfile_HI_VF_sigma_geo_radial_angle_w[1000];
    char fitsfile_2Dinput_VF_boxfiltered_decim0[1000];
    char fitsfile_2Dinput_VF_boxfiltered_decim_user[1000];
    char fitsfile_trfit_model[1000];
    char fitsfile_einasto_modelvf[1000];
    char fitsfile_einasto_halofit_geo_radial_w[1000];
    char fitsfile_trfit_w[1000];
    char fitsfile_fract_Navail_Nall[1000];
    char fitsfile_res_input_minus_trfit[1000];
    char fitsfile_res_input_minus_einastofit[1000];
    char fitsfile_res_trfit_minus_einastofit[1000];
} filename_2dbat;


typedef struct {
    /* set the MultiNest sampling parameters */
    int is;                 // IS mode?
    int ceff;                   // run in constant efficiency mode?
    int mmodal;                 // do mode separation?
    int nlive;
    int nlive_einasto_halofit;              // number of live points

    double efr;             // set the required efficiency
    double tol;             // tol, defines the stopping criteria
    int maxiter;                // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 

    int ndims;                  // dimensionality (no. of free parameters)
    int nPar;                   // total no. of parameters including free & derived parameters
    int nClsPar;    
    int updInt;             // after how many iterations feedback is required & the output files should be updated
                            // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
    double Ztol;                // all the modes with logZ < Ztol are ignored
    int maxModes;               // expected max no. of modes (used only for memory allocation)
    int pWrap[9999];                // which parameters to have periodic boundary conditions?

    char root[500];     // root for output files
    int seed;                   // random no. generator seed, if < 0 then take the seed from system clock
    int fb;                 // need feedback on standard output?
    int resume;                 // resume from a previous job?
    int outfile;                // write output files?
    int initMPI;                // initialize MPI routines?, relevant only if compiling with MPI
                            // set it to F if you want your main program to handle MPI initialization
    double logZero;         // points with loglike < logZero will be ignored by MultiNest
                            // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

    int _is;
    int _ceff;
    int _nlive;
    int _nlive_einasto_halofit;
    int _fb;
    int _outfile;
    int _maxiter;
    double _efr;
    double _tol;

    int _is_d;
    int _ceff_d;
    int _nlive_d;
    int _nlive_einasto_halofit_d;
    int _fb_d;
    int _outfile_d;
    int _maxiter_d;
    double _efr_d;
    double _tol_d;
    
} multinest_paramters;


typedef struct {
    float **data;
} velocity_field;



/* global parameters: 2D velocity field array */
/*
typedef struct {
    velocity_field HI_VF[1]; // array for 2D vf : global
    velocity_field HI_VF_mom0[1]; // array for 2D mom0 : global
    velocity_field HI_VF_mom2[1]; // array for 2D mom2 : global
    velocity_field HI_VF_median[1]; // array for 2D vf : global
    velocity_field HI_VF_res[1]; // array for 2D vf : global
    velocity_field HI_VF_sn[1]; // array for 2D vf : global
    velocity_field HI_VF_sn_median[1]; // array for 2D vf : global
    velocity_field HI_VF_sn_res[1]; // array for 2D vf : global
    velocity_field HI_VF_sigma[1]; // array for 2D vf : global
    velocity_field HI_VF_boxFiltered_vlos_ew[1]; // array for 2D vf : global
    velocity_field HI_VF_sigma_res[1]; // array for 2D vf : global
    velocity_field HI_VF_boxFiltered[1]; // array for 2D vf : global
    velocity_field HI_VF_boxFiltered_sigma[1]; // array for 2D vf : global
    velocity_field HI_VF_boxFiltered_sigma_e_norm[1]; // array for 2D vf : global
    velocity_field HI_VF_boxFiltered_SN[1]; // array for 2D vf : global
    velocity_field HI_VF_boxFiltered_decim0[1]; // array for 2D vf : global
    velocity_field HI_VF_boxFiltered_decim_user[1]; // array for 2D vf : global
    velocity_field HI_VF_geo_radial_angle_w[1]; // array for 2D vf : global
    velocity_field HI_VF_weight_TRfit[1]; // array for 2D vf : global
    velocity_field HI_VF_fract_navail_nall[1]; // array for 2D vf : global
    velocity_field HI_VF_tr_model[1]; // array for 2D vf : global
    velocity_field HI_VF_einasto_halomodel[1]; // array for 2D vf : global
    velocity_field HI_VF_sigma_geo_flux_weighted[1]; // array for 2D vf : global
    velocity_field HI_VF_res_input_minus_trfit[1]; // array for 2D vf : global
    velocity_field HI_VF_res_input_minus_einastofit[1]; // array for 2D vf : global
    velocity_field HI_VF_res_trfit_minus_einastofit[1]; // array for 2D vf : global
} twod_maps_struct;
*/

// global variables : use extern to avoid multiple definition error
extern TR_ringParameters *TRparam; // connection line to multinest
extern multinest_paramters *multinest_param;
extern filename_2dbat *fname;
//extern twod_maps_struct *twod_maps;

// global 2d arrays : accessible to log_likelihood function in multinest
extern velocity_field HI_VF[1]; // array for 2D vf : global
extern velocity_field HI_VF_mom0[1]; // array for 2D mom0 : global
extern velocity_field HI_VF_mom2[1]; // array for 2D mom2 : global
extern velocity_field HI_VF_median[1]; // array for 2D vf : global
extern velocity_field HI_VF_res[1]; // array for 2D vf : global
extern velocity_field HI_VF_sn[1]; // array for 2D vf : global
extern velocity_field HI_VF_sn_median[1]; // array for 2D vf : global
extern velocity_field HI_VF_sn_res[1]; // array for 2D vf : global
extern velocity_field HI_VF_sigma[1]; // array for 2D vf : global
extern velocity_field HI_VF_boxFiltered_vlos_ew[1]; // array for 2D vf : global
extern velocity_field HI_VF_sigma_res[1]; // array for 2D vf : global
extern velocity_field HI_VF_boxFiltered[1]; // array for 2D vf : global
extern velocity_field HI_VF_boxFiltered_sigma[1]; // array for 2D vf : global
extern velocity_field HI_VF_boxFiltered_sigma_e_norm[1]; // array for 2D vf : global
extern velocity_field HI_VF_boxFiltered_SN[1]; // array for 2D vf : global
extern velocity_field HI_VF_boxFiltered_decim0[1]; // array for 2D vf : global
extern velocity_field HI_VF_boxFiltered_decim_user[1]; // array for 2D vf : global
extern velocity_field HI_VF_geo_radial_angle_w[1]; // array for 2D vf : global
extern velocity_field HI_VF_weight_TRfit[1]; // array for 2D vf : global
extern velocity_field HI_VF_fract_navail_nall[1]; // array for 2D vf : global
extern velocity_field HI_VF_tr_model[1]; // array for 2D vf : global
extern velocity_field HI_VF_einasto_halomodel[1]; // array for 2D vf : global
extern velocity_field HI_VF_sigma_geo_flux_weighted[1]; // array for 2D vf : global
extern velocity_field HI_VF_res_input_minus_trfit[1]; // array for 2D vf : global
extern velocity_field HI_VF_res_input_minus_einastofit[1]; // array for 2D vf : global
extern velocity_field HI_VF_res_trfit_minus_einastofit[1]; // array for 2D vf : global
extern velocity_field HI_VF_temp[1]; // array for 2D vf : global

// global parameter for _n of einasto halo which is used for deriving integral of a variant incomplete gamma function in lower_inc_gamma_variant()
extern double _n_einasto_global;
extern struct rGalaxyPlane_pixel_TR_nonlinearEquation_params *p;

// MPI datatype: see 2dbat.mpi_datatype.h
extern MPI_Datatype type[302];
extern int blocklen[302]; 

// --- End of line

#endif


