
#include "2dbat.einastofit.h"

// 2DBAT user defined functions
// Einasto halo fit related

double _n_einasto_global;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void v_einasto_1d_multinest4p_remove_outliers(TR_ringParameters *TRparam)
{
    int i=0, j=0;

    double chi2 = 0.0;
    double logsum_errors = 0.;
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double npoints = 0.0;
    double errorbar_norm_ratio;

    double vEinasto_model;
    double r_pc, D=1; // assuming the galaxy locates at 1 Mpc

    double *res = malloc(sizeof(double) * TRparam[0].Nrings_intp);
    double mean_res, std_res;
    double Cube[4];

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. priors for coefficients: x1 ~ x2

    // parameter 1: alpha
    Cube[0] = TRparam[0]._n;

    // parameter 2: r_2
    Cube[1] = TRparam[0].r_2;

    // parameter 3: rho_2
    Cube[2] = TRparam[0].rho_2;

    for(i=0; i<TRparam[0].Nrings_intp; i++)
    {
        if(TRparam[0].ring_intp[i] <= 0)
            TRparam[0].ring_intp[i] = 0.1;
        if(TRparam[0].vrot_e_intp[i] <= 0)
            TRparam[0].vrot_e_intp[i] = 1;

        r_pc = TRparam[0].ring_intp[i] * TRparam[0].pixelScale * D * pow(10, 3) / 206.265; // pix to pc
        vEinasto_model = vEinasto3p(Cube, r_pc);

        res[i] = TRparam[0].vrot_intp[i] - vEinasto_model;
	//printf("%f %f %f\n", TRparam[0].ring_intp[i], TRparam[0].vrot_intp[i], vEinasto_model);
    }
    robust_mean_std(res, TRparam[0].Nrings_intp, &mean_res, &std_res);
    robust_mean_std_histogram_ac(res, TRparam[0].vrot_e_intp, TRparam[0].Nrings_intp, &mean_res, &std_res);
    //printf("%f %f %e %f\n", TRparam[0]._n, TRparam[0].r_2, TRparam[0].rho_2, TRparam[0].e_sigma_fitted);

    for(i=0; i<TRparam[0].Nrings_intp; i++)
    {
        if(TRparam[0].ring_intp[i] <= 0)
            TRparam[0].ring_intp[i] = 0.1;
        if(TRparam[0].vrot_e_intp[i] <= 0)
            TRparam[0].vrot_e_intp[i] = 1;

        r_pc = TRparam[0].ring_intp[i] * TRparam[0].pixelScale * D * pow(10, 3) / 206.265; // pix to pc
        vEinasto_model = vEinasto3p(Cube, r_pc);

        if(fabs(TRparam[0].vrot_intp[i] - vEinasto_model) > mean_res + 5*std_res)
        {
            //TRparam[0].vrot_intp[i] = vEinasto_model; // replace outliers with models
        }
    }
    free(res);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void v_einasto_1d_multinest4p(multinest_paramters *multinest_param, TR_ringParameters *TRparam)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0;

    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    double efr, tol, Ztol, logZero;

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive;
    //if(nlive > 50) nlive = 50; // for a quick 1D einasto fit
    nlive = 100;
    efr = 0.1;
    tol = 0.1;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;

    strcpy(root, "v_einasto_1d_3p.");
    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    ndims = 4;
    nPar = 4;
    nClsPar = 4;

    int pWrap[ndims];
    for(i=0; i<ndims; i++)
    {
        pWrap[i] = multinest_param[0].pWrap[i];
    }

    /* Calling multinest */
    MPI_Barrier(MPI_COMM_WORLD);
    run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_vEinasto4p_student, dumper_vEinasto4p_student, TRparam);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_vEinasto4p_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0;

    double chi2 = 0.0;
    double logsum_errors = 0.;
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double npoints = 0.0;
    double e_sigma;
    double _r, _r1, _r2, _mu, _sigma;

    double vEinasto_model;
    double r_pc, D=1; // assuming the galaxy locates at 1 Mpc

    // for student-t distribution parameters
    double logsum_sigma2=0, logsum_student_chi=0, _nu, log_Likelihood_studenT;

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. priors for coefficients: x1 ~ x2

    // parameter 1: n
    _r = Cube[0];
    _r1 = TRparam[0]._n1_t;
    _r2 = TRparam[0]._n2_t;
    _mu = TRparam[0]._n_t;
    _sigma = _mu / 3.;
    if(TRparam[0].n_hist_post = 999)
        Cube[0] = fabs(gaussian_prior(_r, _mu, _sigma));
    else
        Cube[0] = uniform_priors(_r, _r1, _r2);

    //Cube[0] = loguniform_priors(_r, _r1, _r2);
    if(Cube[0] == 0) Cube[0] = _mu/20.0;

    // parameter 2: r_2
    _r = Cube[1];
    _r1 = TRparam[0].r_21_t;
    _r2 = TRparam[0].r_22_t;
    _mu = TRparam[0]._r_2_t;
    _sigma = _mu / 3.;

    if(TRparam[0].n_hist_post = 999)
        Cube[1] = fabs(gaussian_prior(_r, _mu, _sigma));
    else
        Cube[1] = loguniform_priors(_r, _r1, _r2);
    if(Cube[1] == 0) Cube[1] = _mu/20.0;

    // parameter 3: rho_2
    _r = Cube[2];
    _r1 = TRparam[0].rho_21_t;
    _r2 = TRparam[0].rho_22_t;
    _mu = TRparam[0]._rho_2_t;
    _sigma = _mu / 3.;
    if(TRparam[0].n_hist_post = 999)
        Cube[2] = fabs(gaussian_prior(_r, _mu, _sigma));
    else
        Cube[2] = loguniform_priors(_r, _r1, _r2);
    if(Cube[2] == 0) Cube[2] = _mu/20.0;

    // parameter 4: e_sigma
    _r = Cube[3];
    _mu = TRparam[0].e_sigma_fitted_t;
    _sigma = _mu / 3.;
    Cube[3] = fabs(gaussian_prior(_r, _mu, _sigma));
    if(Cube[3] == 0) Cube[3] = _mu/20.0;
    e_sigma = Cube[3];

    // 8. calculate log likelihood value
    // _nu of student-T distribution : _nu = 30 for normal : _nu = 1 recommended for best removing the outliers
    _nu = 3; // heavy tail
    logsum_student_chi = 0;
    logsum_errors = 0.;
    npoints = TRparam[0].Nrings_intp;
    //npoints = TRparam[0].Nrings;

    for(i=0; i<TRparam[0].Nrings_intp; i++)
    {
        if(TRparam[0].ring_intp[i] <= 0)
            TRparam[0].ring_intp[i] = 0.1;
        if(TRparam[0].vrot_e_intp[i] <= 0)
            TRparam[0].vrot_e_intp[i] = 1;

        r_pc = TRparam[0].ring_intp[i] * TRparam[0].pixelScale * D * pow(10, 3) / 206.265; // pix to pc
        vEinasto_model = vEinasto3p(Cube, r_pc);

        logsum_errors += 0.5*log(e_sigma*e_sigma);
        logsum_student_chi += ((1+_nu)/2.0)*log(1.0+pow((TRparam[0].vrot_intp[i] - vEinasto_model)/e_sigma, 2)/(_nu-2));
        //logsum_student_chi += ((1+_nu)/2.0)*log(1.0+pow((TRparam[0].vrot[i] - vEinasto_model)/e_sigma, 2)/(_nu-2));
    }
    log_Likelihood_studenT = npoints*(log(gsl_sf_gamma((_nu+1)/2.0)) - log(sqrt(M_PI*(_nu-2))*gsl_sf_gamma(_nu/2.0))) - logsum_errors - logsum_student_chi;
    *lnew = log_Likelihood_studenT;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dumper_vEinasto4p_student(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam)
{
    // convert the 2D Fortran arrays to C arrays
    int i, j;

    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    double postdist[*nSamples][*nPar + 2];
    for( i = 0; i < *nPar + 2; i++ )
    {
        for( j = 0; j < *nSamples; j++ )
        {
            postdist[j][i] = posterior[0][i * (*nSamples) + j];
        }

    }

    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
    double pLivePts[*nlive][*nPar + 1];
    for( i = 0; i < *nPar + 1; i++ )
    {
        for( j = 0; j < *nlive; j++ )
        {
            pLivePts[j][i] = physLive[0][i * (*nlive) + j];
        }
    }

    // Best fits
    TRparam[0]._n = paramConstr[0][*nPar*2+0];
    TRparam[0].r_2 = paramConstr[0][*nPar*2+1];
    TRparam[0].rho_2 = paramConstr[0][*nPar*2+2];
    TRparam[0].e_sigma_fitted = paramConstr[0][*nPar*2+3];

    // standard deviation of the parameters
    TRparam[0]._ne = paramConstr[0][*nPar*1+0];
    if(TRparam[0]._ne_t < TRparam[0]._ne)
    {
        TRparam[0]._ne_t = TRparam[0]._ne; // for setting up the efficient gap! to keep the error is large
        TRparam[0].einasto1D_ne = TRparam[0]._ne; // for setting up the efficient gap! to keep the error is large
    }

    TRparam[0].r_2e = paramConstr[0][*nPar*1+1];
    if(TRparam[0]._r_2e_t < TRparam[0].r_2e)
    {
        TRparam[0]._r_2e_t = TRparam[0].r_2e; 
        TRparam[0].einasto1D_r2e = TRparam[0].r_2e; 
    }

    TRparam[0].rho_2e = paramConstr[0][*nPar*1+2];
    if(TRparam[0]._rho_2e_t < TRparam[0].rho_2e)
    {
        TRparam[0]._rho_2e_t = TRparam[0].rho_2e; 
        TRparam[0].einasto1D_rho2e = TRparam[0].rho_2e; 
    }

    TRparam[0].maxLogLikeF = *maxLogLike;
    //TRparam[0].logZF = *logZ;
    //TRparam[0].logZerrF = *logZerr;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double vEinasto3p(double *Cube, double r)
{
    int i=0;
    double G = 4.302 * pow(10, -3); // pc * Mo^-1 * (km/s)^2
    double e = 2.71828;
    double _n, r_2, rho_2;
    double iGamma, _x;
    double vEinasto_model;

    // Cube[0] : alpha;
    // Cube[1] : r_2;
    // Cube[2] : rho_2;
    _n = Cube[0];
    r_2 = Cube[1]; // in pc
    rho_2 = Cube[2]; // in Mo/pc^3

    // incomplete Gamma function
    _x = r/r_2;

    iGamma = gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);
    vEinasto_model = sqrt((4.0*M_PI*G*_n*rho_2*pow(r_2, 3)/r)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) * iGamma));

    return vEinasto_model;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void spline_intp(TR_ringParameters *TRparam)
{
    int i=0;
    double xi, yi;

    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, TRparam[0].Nrings);
        gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].vrot_temp, TRparam[0].Nrings);

        for (xi=TRparam[0].ring_radius[0]; xi<=TRparam[0].ring_radius[TRparam[0].Nrings-1]; xi+=0.3*TRparam[0].ring_radius[0])
        {
            TRparam[0].ring_intp[i] = xi;
            yi = gsl_spline_eval(spline, xi, acc);
            //if(i > 0 && yi < 0.5*TRparam[0].vrot_intp[i-1]) // to avoid any weired VROT value
            //{
            //    TRparam[0].vrot_intp[i] = TRparam[0].vrot_intp[i-1];
            //}
            //else
            //{
            //    TRparam[0].vrot_intp[i] = yi;
            //}

            TRparam[0].vrot_intp[i] = yi;
            TRparam[0].vrot_e_intp[i] = 0.05*TRparam[0].vrot[(int)TRparam[0].Nrings/2];
            i++;
        }
        TRparam[0].Nrings_intp = i;
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);
    }
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double spline_intp_ringparam(TR_ringParameters *TRparam, char *param, double ring)
{
    int i=0;
    double intp_ringparam;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, TRparam[0].Nrings+1);
    double radius[9999];
    double xpos[9999];
    double ypos[9999];
    double vsys[9999];
    double pa[9999];
    double incl[9999];
    double vrot[9999];
    double vrad[9999];

    for(i=0; i<TRparam[0].Nrings+1; i++)
    {
        if(i==0)
        {
            radius[i] = 0;
            xpos[i] = TRparam[0].xpos0[0];
            ypos[i] = TRparam[0].ypos0[0];
            vsys[i] = TRparam[0].vsys0[0];
            pa[i] = TRparam[0].pa0[0];
            incl[i] = TRparam[0].incl0[0];
            vrot[i] = 0;
            vrad[i] = TRparam[0].vrad[0];
        }
        else
        {
            radius[i] = TRparam[0].ring_radius[i-1];
            xpos[i] = TRparam[0].xpos0[i-1];
            ypos[i] = TRparam[0].ypos0[i-1];
            vsys[i] = TRparam[0].vsys0[i-1];
            pa[i] = TRparam[0].pa0[i-1];
            incl[i] = TRparam[0].incl0[i-1];
            vrot[i] = TRparam[0].vrot0[i-1];
            vrad[i] = TRparam[0].vrad[i-1];
        }
    }

    if(strcmp(param, "XPOS") == 0)
    {
        //gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].xpos0, TRparam[0].Nrings);
        gsl_spline_init(spline, radius, xpos, TRparam[0].Nrings+1);
        intp_ringparam = gsl_spline_eval(spline, ring, acc);
    }
    else if(strcmp(param, "YPOS") == 0)
    {
        //gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].ypos0, TRparam[0].Nrings);
        gsl_spline_init(spline, radius, ypos, TRparam[0].Nrings+1);
        intp_ringparam = gsl_spline_eval(spline, ring, acc);
    }
    else if(strcmp(param, "VSYS") == 0)
    {
        //gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].vsys0, TRparam[0].Nrings);
        gsl_spline_init(spline, radius, vsys, TRparam[0].Nrings+1);
        intp_ringparam = gsl_spline_eval(spline, ring, acc);
    }
    else if(strcmp(param, "PA") == 0)
    {
        //gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].pa0, TRparam[0].Nrings);
        gsl_spline_init(spline, radius, pa, TRparam[0].Nrings+1);
        intp_ringparam = gsl_spline_eval(spline, ring, acc);
    }
    else if(strcmp(param, "INCL") == 0)
    {
        //gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl0, TRparam[0].Nrings);
        gsl_spline_init(spline, radius, incl, TRparam[0].Nrings+1);
        intp_ringparam = gsl_spline_eval(spline, ring, acc);
    }
    else if(strcmp(param, "VROT") == 0)
    {
        //gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].vrot0, TRparam[0].Nrings);
        gsl_spline_init(spline, radius, vrot, TRparam[0].Nrings+1);
        intp_ringparam = gsl_spline_eval(spline, ring, acc);
    }
    else if(strcmp(param, "VRAD") == 0)
    {
        //gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].vrad, TRparam[0].Nrings); // vrad[]
        gsl_spline_init(spline, radius, vrad, TRparam[0].Nrings+1); // vrad[]
        intp_ringparam = gsl_spline_eval(spline, ring, acc);
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return intp_ringparam;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void vEinasto_multinest2p(multinest_paramters *multinest_param, TR_ringParameters *TRparam)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0;

    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    double efr, tol, Ztol, logZero;

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive;
    if(nlive > 50) nlive = 500; // for a quick 1D einasto fit
    efr = 0.1;
    tol = 0.1;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;

    strcpy(root, "v_einasto_2d_2p.");
    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    ndims = 2;
    nPar = 2;
    nClsPar = 2;

    int pWrap[ndims];
    for(i=0; i<ndims; i++)
    {
        pWrap[i] = multinest_param[0].pWrap[i];
    }

    /* Calling multinest */
    run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_vEinasto2p, dumper_vEinasto2p, TRparam);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_vEinasto2p(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0;

    double chi2 = 0.0;
    double logsum_errors = 0.;
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double npoints = 0.0;
    double errorbar_norm_ratio;

    double vEinasto_model;
    double _n1, _n2, r_21, r_22, rho_21, rho_22;

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. priors for coefficients: x1 ~ x2

    // parameter 1: alpha
    _n1 = TRparam[0]._n1;
    _n2 = TRparam[0]._n2;
    Cube[0] = _n1 + Cube[0]*(_n2-_n1);

    // parameter 2: r_2
    //r_21 = TRparam[0].r_21;
    //r_22 = TRparam[0].r_22;
    //Cube[1] = r_21 + Cube[1]*(r_22-r_21);

    // parameter 3: rho_2
    rho_21 = TRparam[0].rho_21;
    rho_22 = TRparam[0].rho_22;
    Cube[1] = rho_21 + Cube[1]*(rho_22-rho_21);


    npoints = TRparam[0].Nrings_intp;
    GLLhood0 = -(npoints/2.0)*log(2.0*M_PI);
    chi2 = 0.;
    logsum_errors = 0.;

    for(i=0; i<TRparam[0].Nrings_intp; i++)
    {
        if(TRparam[0].ring_intp[i] <= 0)
            TRparam[0].ring_intp[i] = 0.1;
        if(TRparam[0].vrot_e_intp[i] <= 0)
            TRparam[0].vrot_e_intp[i] = 1;

        //TRparam[0].vrot_e_intp[i] = 0.1;

        vEinasto_model = vEinasto2p(Cube, TRparam[0].ring_intp[i], TRparam[0].r_2);
        logsum_errors += log(TRparam[0].vrot_e_intp[i]);
        chi2 += pow(((TRparam[0].vrot_intp[i] - vEinasto_model)/TRparam[0].vrot_e_intp[i]), 2);
    }
    slhood = GLLhood0 - logsum_errors - chi2/2.0;
    *lnew = slhood;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dumper_vEinasto2p(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam)
{
    // convert the 2D Fortran arrays to C arrays
    int i, j;

    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    double postdist[*nSamples][*nPar + 2];
    for( i = 0; i < *nPar + 2; i++ )
    {
        for( j = 0; j < *nSamples; j++ )
        {
            postdist[j][i] = posterior[0][i * (*nSamples) + j];
        }

    }

    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
    double pLivePts[*nlive][*nPar + 1];
    for( i = 0; i < *nPar + 1; i++ )
    {
        for( j = 0; j < *nlive; j++ )
        {
            pLivePts[j][i] = physLive[0][i * (*nlive) + j];
        }
    }

    // Best fits
    TRparam[0]._n = paramConstr[0][*nPar*2+0];
    //TRparam[0].r_2 = paramConstr[0][*nPar*2+0];
    TRparam[0].rho_2 = paramConstr[0][*nPar*2+1];

    // standard deviation of the parameters
    TRparam[0]._ne = paramConstr[0][*nPar*1+0];
    //TRparam[0].r_2e = paramConstr[0][*nPar*1+0];
    TRparam[0].rho_2e = paramConstr[0][*nPar*1+1];

    TRparam[0].maxLogLikeF = *maxLogLike;
    //TRparam[0].logZF = *logZ;
    //TRparam[0].logZerrF = *logZerr;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double vEinasto2p(double *Cube, double r, double r_2)
{
    int i=0;
    double G = 4.302;
    double e = 2.71828;
    double _n, rho_2;
    double iGamma, _x;
    double vEinasto_model;

    // Cube[0] : alpha;
    // Cube[1] : r_2;
    // Cube[2] : rho_2;
    //_n = Cube[0];
    _n = Cube[0];
    rho_2 = Cube[1];

    // incomplete Gamma function
    _x = r/r_2;

    iGamma = gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);
    vEinasto_model = sqrt((4.0*M_PI*G*_n*rho_2*pow(r_2, 3)/r)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) * iGamma));

    return vEinasto_model;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Einasto_haloFits_multinest(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *sersic_function, char sersicpart, char *_n_Einasto, char _n_fix, char *r_2_Einasto, char r_2_fix, char *rho_2_Einasto, char rho_2_fix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int rank, char *mt_outputfile)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0, j=0, i0=0, j0=0, string_length;
    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    double efr, tol, Ztol, logZero;
    double total_Npoints_allrings=0.;

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive_einasto_halofit;
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;
    nClsPar = multinest_param[0].nClsPar;
    ndims = multinest_param[0].ndims;

    //strncpy(root, multinest_param[0].root, strlen(multinest_param[0].root));
    strcpy(root, mt_outputfile);

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    set_nfree_params_einasto_halofit_multinest_student(xposfix, yposfix, vsysfix, pafix, inclfix, _n_fix, r_2_fix, rho_2_fix, vradfix, sigmafactorfix, TRparam);
    ndims = TRparam[0].n_freeParams;
    nPar = TRparam[0].n_freeParams;
    nClsPar = TRparam[0].n_freeParams;
    int pWrap[ndims];
    for(i = 0; i < ndims; i++)
    {
        pWrap[i] = multinest_param[0].pWrap[i];
    }

    /* Calling multinest */ 
    MPI_Barrier(MPI_COMM_WORLD);
    run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_einasto_halofit, dumper_einasto_halofit, TRparam);    
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_einasto_halofit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0, ii=0, i0=0, j0=0;
    int pa_n=0, pa_nF=0, pa_nS=0;
    int incl_n=0, incl_nF=0, incl_nS=0;
    int vrad_n=0, vrad_nF=0, vrad_nS=0;
    int n_freeParams = 0;
    int n_freeParams_uptoPA=0;
    int n_freeParams_uptoINCL=0;
    int w_n;

    double sigma_factor, v_sigma_weighted;
    double chi2 = 0.0;
    double ri, ro;
    double logsum_errors = 0., mean_errors=0, sum_errors=0;  
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double NobsPoints = 0.0, NobsPoints_available=0;
    double Vmodel_Einasto = 0.;
    double residual_VF_outliers_filter=0;
    double first_ring_in_pixels;
    float rbm_SN, rbm_sigma, std_SN, std_sigma;
    double error;
    double mom4;
    // log-normal distribution
    double logNormal_median, logNormal_sigma;
    double mu_log, sigma_log;

    // gaussian prior
    double _r, _mu, _sigma;
    double g_wing;

    // uniform prior
    double _r1, _r2; 


    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. XPOS prior : xpo1 ~ xpo2
    if(TRparam[0].xpos_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0].xpos_mode;
        _sigma = 1*TRparam[0].xpos_std;

        _r1 = TRparam[0].xpos1;
        _r2 = TRparam[0].xpos2;

        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        if(Cube[n_freeParams] < 0) Cube[n_freeParams] = 1E-5;
        n_freeParams++;
    }
    // 2. YPOS prior : ypo1 ~ ypo2
    if(TRparam[0].ypos_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0].ypos_mode;
        _sigma = 1*TRparam[0].ypos_std;

        _r1 = TRparam[0].ypos1;
        _r2 = TRparam[0].ypos2;

        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        if(Cube[n_freeParams] < 0) Cube[n_freeParams] = 1E-5;
        n_freeParams++;
    }
    // 3. VSYS prior : vsys1 ~ vsys2
    if(TRparam[0].vsys_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0].vsys_mode;
        _sigma = 1*TRparam[0].vsys_std;

        _r1 = TRparam[0].vsys1;
        _r2 = TRparam[0].vsys2;

        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        n_freeParams++;
    }
    // 4. PA prior :
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _r = Cube[n_freeParams];
                _mu = TRparam[0]._p_bs_t[pa_nS]; // mean: based on the first bspline fit to the first TRfit results

                if(_mu > 1.0-_mu) // skewed leftward 
                {
                    g_wing = _mu; // take the wider wing
                }
                else // skewed leftward
                {
                    g_wing = 1.0 - _mu;
                }
                _sigma = g_wing / 5.0;

                //_sigma = 1.0*TRparam[0]._p_bs_e_t[pa_nS]; // bspline_err : based on the first bspline fit to the first TRfit results
                TRparam[0]._p_bs_e_t[pa_nS] = _sigma;

                //_r1 = TRparam[0].bspline1pa[pa_nS];
                //_r2 = TRparam[0].bspline2pa[pa_nS];
                _r1 = _mu - 5*_sigma;
                _r2 = _mu + 5*_sigma;

                if(TRparam[0].fullFit == 'Y')
                {
                    //Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
                }
                else if(TRparam[0].fullFit == 'N')
                {
                    //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
                }

                //Cube[n_freeParams] = uniform_priors(_r, 1E-1, 0.9);
                Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
                //Cube[n_freeParams] = uniform_priors(_r, 1E-2, 1-1E-2);
                //Cube[n_freeParams] = JeffreysPrior(_r, 1E-1, 0.9);

                if(Cube[n_freeParams] <= 0) Cube[n_freeParams] = 1E-2;
                if(Cube[n_freeParams] >= 1) Cube[n_freeParams] = 1.0-1E-2;

                n_freeParams++;
            }
        }
    }
    n_freeParams_uptoINCL = n_freeParams;   

    // 5. INCL prior :
    if(TRparam[0].incl_fix == 'T')
    {
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _r = Cube[n_freeParams];
                _mu = TRparam[0]._i_bs_t[incl_nS]; // mean: based on the first bspline fit to the first TRfit results

                if(_mu > 1.0-_mu) // skewed leftward 
                {
                    g_wing = _mu; // take the wider wing
                }
                else // skewed rightward
                {
                    g_wing = 1.0 - _mu;
                }

                // inclination weight
                //_sigma = sin(_mu*M_PI/180.)*g_wing / 5.0;
                _sigma = g_wing / 5.0;


                TRparam[0]._i_bs_e_t[incl_nS] = _sigma;

                //_sigma = 1.0*TRparam[0]._i_bs_e_t[incl_nS]; // bspline_err : based on the first bspline fit to the first TRfit results

                _r1 = _mu - 5*_sigma;
                _r2 = _mu + 5*_sigma;
                if(_r1 <= 0) _r1 = 1E-3;
                if(_r2 >= 1) _r2 = 0.999;

                if(TRparam[0].fullFit == 'Y')
                {
                //    Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
                }
                else if(TRparam[0].fullFit == 'N')
                {
                 //   Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
                }

                //Cube[n_freeParams] = uniform_priors(_r, 1E-2, 1-1E-2);
                //Cube[n_freeParams] = JeffreysPrior(_r, 1E-1, 0.9);
                Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
                //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);

                if(Cube[n_freeParams] <= 0) Cube[n_freeParams] = 1E-5;
                if(Cube[n_freeParams] >= 1) Cube[n_freeParams] = 1.0-1E-5;

                n_freeParams++;
            }
        }
    }

    // 6. Einasto rotation curve : _n : _n1 ~ _n2
    if(TRparam[0]._n_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0]._n_t;
        _sigma = _mu / 3.0;

        _r1 = TRparam[0]._n1;
        _r2 = TRparam[0]._n2;

        _r1 = 0.1*_mu;
        _r2 = 3*_mu;

        //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);
        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
        //Cube[n_freeParams] = logNormalPrior(_r, mu_log, sigma_log);
        //Cube[n_freeParams] = DeltaFunctionPrior(_r, _mu);

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        if(Cube[n_freeParams] <= 0) Cube[n_freeParams] = _mu/3.0;
        //if(Cube[n_freeParams] >= _r2) Cube[n_freeParams] = _r2;

        n_freeParams++;
    }

    // 7. Einasto rotation curve : r_2 : r_21 ~ r_22
    if(TRparam[0].r_2_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0]._r_2_t;
        _sigma = _mu/3.;
        //_sigma = 0.1*TRparam[0]._r_2_t;

        _r1 = TRparam[0].r_21;
        _r2 = TRparam[0].r_22;

        //Cube[n_freeParams] = DeltaFunctionPrior(_r, _mu);
        //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);
        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        //mu_log = _mu;
        //sigma_log = 0.5;
        //Cube[n_freeParams] = logNormalPrior(_r, mu_log, sigma_log);

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        //if(Cube[n_freeParams] < _r1) Cube[n_freeParams] = _r1;
        //if(Cube[n_freeParams] > _r2) Cube[n_freeParams] = _r2;

        if(Cube[n_freeParams] <= 0) Cube[n_freeParams] = _mu/3.0;
        //if(Cube[n_freeParams] >= _r2) Cube[n_freeParams] = _r2;

        n_freeParams++;
    }

    // 8. Einasto rotation curve : rho_2 : rho_21 ~ rho_22
    if(TRparam[0].rho_2_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = 1*TRparam[0]._rho_2_t;
        _sigma = _mu / 3.;

        _r1 = TRparam[0].rho_21;
        _r2 = TRparam[0].rho_22;


        //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        //Cube[n_freeParams] = DeltaFunctionPrior(_r, _mu);
        //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);
        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        //mu_log = _mu;
        //sigma_log = 0.5;
        //Cube[n_freeParams] = logNormalPrior(_r, mu_log, sigma_log);

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        //if(Cube[n_freeParams] < _r1) Cube[n_freeParams] = _r1;
        //if(Cube[n_freeParams] > _r2) Cube[n_freeParams] = _r2;

        if(Cube[n_freeParams] <= 0) Cube[n_freeParams] = _mu/3.0;
        //if(Cube[n_freeParams] >= _r2) Cube[n_freeParams] = _r2;

        n_freeParams++;
    }


    // 9. VRAD prior :
    if(TRparam[0].vrad_fix == 'T')
    {
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
        {
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _r = Cube[n_freeParams];
                _mu = TRparam[0]._vr_bs_t[vrad_nS]; // mean: based on the first bspline fit to the first TRfit results
                _sigma = fabs(_mu);
                //_sigma = 1.0*TRparam[0]._vr_bs_e_t[vrad_nS]; // bspline_err : based on the first bspline fit to the first TRfit results

                //_r1 = TRparam[0].bspline1vrad[vrad_nS];
                //_r2 = TRparam[0].bspline2vrad[vrad_nS];
                //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
                Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

                //if(Cube[n_freeParams] < 0) Cube[n_freeParams] = 1E-5;
                //if(Cube[n_freeParams] > 1) Cube[n_freeParams] = 1.0-1E-5;
                n_freeParams++;
            }
        }
    }

    // 10. sigma_factor prior : sigma_factor1 ~ sigma_factor2
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0].sigma_factor_mode;
        _sigma = 1*TRparam[0].sigma_factor_std;

        _r1 = TRparam[0].sigma_factor1;
        _r2 = TRparam[0].sigma_factor2;

        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        sigma_factor = Cube[n_freeParams];
        n_freeParams++;
    }
    else
        sigma_factor = TRparam[0].sigma_factor; // use the best fit value from Vmodel_Einasto fitting

    // 11. e_sigma_fitted : 
    if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0)
    {
        _r = Cube[n_freeParams];
        _r1 = 1E-3;
        _r2 = 1E1;

        Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);
        v_sigma_weighted = Cube[n_freeParams];
        n_freeParams++;
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma != 0)
        v_sigma_weighted = TRparam[0].e_sigma;


    // 8. calculate log likelihood value
    chi2 = 0.;
    logsum_errors = 0.;
    sum_errors = 0.;
    mean_errors = 0.;
    NobsPoints_available = 0.;

    //v_sigma_weighted = 0.5;
  
    for(i=0; i<TRparam[0].Npoints_in_tilted_ring; i++)
    {
        i0 = TRparam[0].tilted_ring[i][0];
        j0 = TRparam[0].tilted_ring[i][1];

        einasto_halomodel(Cube, i0, j0, TRparam, &Vmodel_Einasto); // Update Vmodel_Einasto
        if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma > 0 && HI_VF_geo_radial_angle_w[0].data[j0][i0] > 0) // constant vlos_e mode
        {   
            // geometry (perimeter + cos(theta)) weighted constant vlos_e
            // HI_VF_geo_radial_angle_w[0].data[j0][i0] is weighted by perimeter + cos(theta) and the normalised (0 ~ 1) 
            // geometrical weight
            v_sigma_weighted = TRparam[0].e_sigma * TRparam[0].scale_factor_const_vlose_w / HI_VF_geo_radial_angle_w[0].data[j0][i0];
            // no geometrical weight
            //v_sigma_weighted = TRparam[0].e_sigma;
            logsum_errors += log(v_sigma_weighted);
                
            chi2 += pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/v_sigma_weighted, 2);
            NobsPoints_available += 1;
        }
        else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && HI_VF_geo_radial_angle_w[0].data[j0][i0] > 0) // constant vlos_e mode : partial fit
        {   
            // geometry (perimeter + cos(theta)) weighted constant vlos_e
            // HI_VF_geo_radial_angle_w[0].data[j0][i0] is weighted by perimeter + cos(theta) and the normalised (0 ~ 1) 
            if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0 \
            && !isinf(HI_VF_mom0[0].data[j0][i0]) \
            && !isnan(HI_VF_mom0[0].data[j0][i0]) \
            && HI_VF_mom0[0].data[j0][i0] > 0.0 && HI_VF_fract_navail_nall[0].data[j0][i0] > 0.1) // no blank
            {
                if(HI_VF_mom0[0].data[j0][i0] > 1E5*TRparam[0].mean_mom4)
                {
                    mom4 = (TRparam[0].mean_mom4);
                }
                else
                {
                    mom4 = (HI_VF_mom0[0].data[j0][i0]);
                }

                //chi2 += pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/v_sigma_weighted, 2)*HI_VF_geo_radial_angle_w[0].data[j0][i0];

                //v_sigma_weighted = (TRparam[0].scale_factor_var_vlose_w * HI_VF_boxFiltered_sigma[0].data[j0][i0]);
                //printf("%f\n", v_sigma_weighted);
                //logsum_errors += log(v_sigma_weighted)*HI_VF_geo_radial_angle_w[0].data[j0][i0];
                //chi2 += pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/v_sigma_weighted, 2)*HI_VF_geo_radial_angle_w[0].data[j0][i0];
                //NobsPoints_available += HI_VF_geo_radial_angle_w[0].data[j0][i0];

                //logsum_errors += log(v_sigma_weighted/HI_VF_geo_radial_angle_w[0].data[j0][i0]);
                //NobsPoints_available += HI_VF_geo_radial_angle_w[0].data[j0][i0];
                //chi2 += pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/(v_sigma_weighted/HI_VF_geo_radial_angle_w[0].data[j0][i0]), 2);

                logsum_errors += log(v_sigma_weighted);
                NobsPoints_available += 1;
                chi2 += pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/v_sigma_weighted, 2);
            }
        }
        else if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0 && HI_VF_geo_radial_angle_w[0].data[j0][i0] > 0) // fitted vlos_e mode : full fit
        {
            if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0 \
            && !isinf(HI_VF_mom0[0].data[j0][i0]) \
            && !isnan(HI_VF_mom0[0].data[j0][i0]) \
            && HI_VF_mom0[0].data[j0][i0] > 0.0) // no blank
            {
                //v_sigma_weighted = HI_VF_boxFiltered_sigma[0].data[j0][i0] * TRparam[0].scale_factor_var_vlose_w / HI_VF_geo_radial_angle_w[0].data[j0][i0] / HI_VF_mom0[0].data[j0][i0];

                v_sigma_weighted = HI_VF_boxFiltered_sigma[0].data[j0][i0] * TRparam[0].scale_factor_var_vlose_w / HI_VF_geo_radial_angle_w[0].data[j0][i0] / HI_VF_mom0[0].data[j0][i0];

                // put sigma_factor which is derived from the Bayesian fit
                v_sigma_weighted = sigma_factor*v_sigma_weighted;
                logsum_errors += log(v_sigma_weighted);
                chi2 += pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/v_sigma_weighted, 2);
                NobsPoints_available += 1;
            }
        }
    }

    GLLhood0 = -(NobsPoints_available/2.0)*log(2.0*M_PI);
    //logsum_errors  NobsPoints_available*log(v_sigma_weighted);
    slhood = GLLhood0 - logsum_errors - chi2/2.0;
    *lnew = slhood;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void einasto_halofit_multinest_student(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *_n_Einasto, char _n_fix, char *r_2_Einasto, char r_2_fix, char *rho_2_Einasto, char rho_2_fix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int rank, char *mt_outputfile)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0, j=0, i0=0, j0=0, string_length;
    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    double efr, tol, Ztol, logZero;
    double total_Npoints_allrings=0.;

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive_einasto_halofit;
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;
    nClsPar = multinest_param[0].nClsPar;
    ndims = multinest_param[0].ndims;

    //strncpy(root, multinest_param[0].root, strlen(multinest_param[0].root));
    strcpy(root, mt_outputfile);

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    set_nfree_params_einasto_halofit_multinest_student(xposfix, yposfix, vsysfix, pafix, inclfix, _n_fix, r_2_fix, rho_2_fix, vradfix, sigmafactorfix, TRparam);
    ndims = TRparam[0].n_freeParams;
    nPar = TRparam[0].n_freeParams;
    nClsPar = TRparam[0].n_freeParams;
    int pWrap[ndims];
    for(i = 0; i < ndims; i++)
    {
        pWrap[i] = multinest_param[0].pWrap[i];
    }

    /* Calling multinest */ 
    MPI_Barrier(MPI_COMM_WORLD);
    run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_einasto_halofit_student, dumper_einasto_halofit, TRparam);    
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_einasto_halofit_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0, i0=0, j0=0;
    int pa_nS=0, incl_nS=0, vrad_nS=0;
    int n_freeParams = 0;
    int n_freeParams_uptoPA=0;
    int n_freeParams_uptoINCL=0;

    double sigma_factor, e_sigma;
    double chi2 = 0.0;
    double logsum_errors = 0., sum_errors=0;  
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double NobsPoints_available=0;
    double Vmodel_Einasto = 0.;
    double error;
    double mom4;
    // log-normal distribution
    //double logNormal_median, logNormal_sigma;
    //double mu_log, sigma_log;

    // gaussian prior
    double _r, _mu, _mu_max, _sigma, _sigma_N, _sigma_W;
    double g_wing;

    // gaussian_skew prior
    double _n_max, _r_max, _rho_max;

    // uniform prior
    double _r1, _r2; 

    // for student-t distribution parameters
    double logsum_sigma2=0, logsum_student_chi=0, _nu, log_Likelihood_studenT; 

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. XPOS prior : xpo1 ~ xpo2
    if(TRparam[0].xpos_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0].xposF_EinastoFit_t;
        _sigma = TRparam[0].ellipse_semi_mx_boxfiltered / 10.;
        Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        if(Cube[n_freeParams] >= TRparam[0].nax1) Cube[n_freeParams] = TRparam[0].nax1/2.0;
        n_freeParams++;

        /*
        _r1 = TRparam[0].xpos1;
        _r2 = TRparam[0].xpos2;
        if(TRparam[0].opt_priors_xpos != 1) // if an efficient prior is not found try uniform prior
        {
            Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        }
        else if(TRparam[0].opt_priors_xpos == 1) // if an efficient prior is found try gaussian prior
        {
            //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
            Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
        }
        if(Cube[n_freeParams] < 0) Cube[n_freeParams] = 1E-5;
        n_freeParams++;
        */
    }
    // 2. YPOS prior : ypo1 ~ ypo2
    if(TRparam[0].ypos_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0].yposF_EinastoFit_t;
        _sigma = TRparam[0].ellipse_semi_mx_boxfiltered / 10.;
        Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        if(Cube[n_freeParams] >= TRparam[0].nax2) Cube[n_freeParams] = TRparam[0].nax2/2.0;
        n_freeParams++;

        /*
        _r = Cube[n_freeParams];
        _mu = TRparam[0].ypos_mode;
        _sigma = 1*TRparam[0].ypos_std;
        _r1 = TRparam[0].ypos1;
        _r2 = TRparam[0].ypos2;
        if(TRparam[0].opt_priors_ypos != 1) // if an efficient prior is not found try uniform prior
        {
            Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        }
        else if(TRparam[0].opt_priors_ypos == 1) // if an efficient prior is found try gaussian prior
        {
            //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
            Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
        }
        if(Cube[n_freeParams] < 0) Cube[n_freeParams] = 1E-5;
        n_freeParams++;
        */
    }
    // 3. VSYS prior : vsys1 ~ vsys2
    if(TRparam[0].vsys_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0].vsysF_EinastoFit_t;
        _sigma = TRparam[0].LOS_hist_Gfit_sigma/1.0;
        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
        n_freeParams++;

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        //if(Cube[n_freeParams] <= TRparam[0].vlos_lower_limit) Cube[n_freeParams] = TRparam[0].vlos_lower_limit;
        //if(Cube[n_freeParams] >= TRparam[0].vlos_upper_limit) Cube[n_freeParams] = TRparam[0].vlos_upper_limit;

        /*
        _r = Cube[n_freeParams];
        _mu = TRparam[0].vsys_mode;
        _sigma = 1*TRparam[0].vsys_std;

        _r1 = TRparam[0].vsys1;
        _r2 = TRparam[0].vsys2;

        if(TRparam[0].opt_priors_vsys != 1) // if an efficient prior is not found try uniform prior
        {
            Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        }
        else if(TRparam[0].opt_priors_vsys == 1) // if an efficient prior is found try gaussian prior
        {
            //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
            Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
        }

        if(Cube[n_freeParams] < 0) Cube[n_freeParams] = 1E-5;
        n_freeParams++;
        */
    }
    // 4. PA prior :
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
/*
                _r = Cube[n_freeParams];
                // update priors based on either the tr or the dirty fits 
                _mu = TRparam[0]._p_bs_t[pa_nS]; // mean: based on the first bspline fit to the first TRfit results
                if(_mu < 1.0-_mu) // skewed leftward 
                {
                    g_wing = _mu; // take the narrower wing
                    _sigma = g_wing / 3.0;
                    _sigma_N = g_wing / 5.0; // narrower sigma
                    _sigma_W = (1.0-g_wing) / 5.0; // wider sigma
                }
                else
                {
                    g_wing = 1.0 - _mu; // take the narrower wing
                    _sigma = g_wing / 3.0;
                    _sigma_N = g_wing / 5.0; // narrower sigma
                    _sigma_W = (1.0-g_wing) / 5.0; // wider sigma
                }
                TRparam[0]._p_bs_e_t[pa_nS] = _sigma;
                Cube[n_freeParams] = gaussian_prior_paincl_skew(_r, _mu, _sigma_N, _sigma_W);
*/

                // update priors based on either the tr or the dirty fits 
                _r = Cube[n_freeParams];
                _mu = TRparam[0]._p_bs_t[pa_nS]; // mean: based on the first bspline fit to the first TRfit results
                if(_mu < 0.5) // take the wider wing
                {
                    _sigma = (1-_mu) / 5.0;
                }
                else
                {
                    _sigma = (_mu) / 5.0;
                }
                
                TRparam[0]._p_bs_e_t[pa_nS] = _sigma;
                Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

                //Cube[n_freeParams] = uniform_priors(_r, 1E-1, 0.9);
                //Cube[n_freeParams] = uniform_priors(_r, 1E-2, 1-1E-2);
                //Cube[n_freeParams] = JeffreysPrior(_r, 1E-2, 0.99);
                if(Cube[n_freeParams] > 1) Cube[n_freeParams] = Cube[n_freeParams] - (int)(Cube[n_freeParams]/1)*1;
                if(Cube[n_freeParams] < 0) Cube[n_freeParams] = Cube[n_freeParams] + (1-(int)(Cube[n_freeParams]/1))*1;

                n_freeParams++;
            }
        }
    }
    n_freeParams_uptoINCL = n_freeParams;   

    // 5. INCL prior :
    if(TRparam[0].incl_fix == 'T')
    {
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
/*
                _r = Cube[n_freeParams];

                // update priors based on either the tr or the dirty fits 
                _mu = TRparam[0]._i_bs_t[incl_nS]; // mean: based on the first bspline fit to the first TRfit results

                if(_mu < 1.0-_mu) // skewed leftward 
                {
                    g_wing = _mu; // take the narrower wing
                    _sigma = g_wing / 3.0;
                    _sigma_N = g_wing / 5.0; // narrower sigma
                    _sigma_W = (1.0-g_wing) / 5.0; // wider sigma
                }
                else
                {
                    g_wing = 1.0 - _mu; // take the narrower wing
                    _sigma = g_wing / 3.0;
                    _sigma_N = g_wing / 5.0; // narrower sigma
                    _sigma_W = (1.0-g_wing) / 5.0; // wider sigma
                }
                TRparam[0]._i_bs_e_t[incl_nS] = _sigma;
                Cube[n_freeParams] = gaussian_prior_paincl_skew(_r, _mu, _sigma_N, _sigma_W);
*/

                //_sigma = 1.0*TRparam[0]._i_bs_e_t[incl_nS]; // bspline_err : based on the first bspline fit to the first TRfit results

                //_r1 = _mu - 5*_sigma;
                //_r2 = _mu + 5*_sigma;
                //if(_r1 <= 0) _r1 = 1E-3;
                //if(_r2 >= 1) _r2 = 0.999;

                //if(TRparam[0].fullFit == 'Y')
                //{
                //    Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
                //}
                //else if(TRparam[0].fullFit == 'N')
                //{
                 //   Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
                //}

                //Cube[n_freeParams] = uniform_priors(_r, 1E-2, 1-1E-2);
                //Cube[n_freeParams] = JeffreysPrior(_r, 1E-2, 0.99);
                //Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
                //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);

                // update priors based on either the tr or the dirty fits 
                _r = Cube[n_freeParams];
                _mu = TRparam[0]._i_bs_t[incl_nS]; // mean: based on the first bspline fit to the first TRfit results
                if(_mu < 0.5) // take the wider wing
                {
                    _sigma = (1-_mu) / 5.0;
                }
                else
                {
                    _sigma = (_mu) / 5.0;
                }
                TRparam[0]._i_bs_e_t[incl_nS] = _sigma;
                Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

                if(Cube[n_freeParams] > 1 && Cube[n_freeParams] <= 2) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/2)+1)*2 - Cube[n_freeParams];
                if(Cube[n_freeParams] > 2 && Cube[n_freeParams] <= 3) Cube[n_freeParams] = 1 - (((int)(Cube[n_freeParams]/3)+1)*3 - Cube[n_freeParams]);
                if(Cube[n_freeParams] > 3 && Cube[n_freeParams] <= 4) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/4)+1)*4 - Cube[n_freeParams];

                if(Cube[n_freeParams] < 0 && Cube[n_freeParams] >= -1) Cube[n_freeParams] = fabs(Cube[n_freeParams]);
                if(Cube[n_freeParams] < -1 && Cube[n_freeParams] >= -2) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/2)+1)*2 + Cube[n_freeParams];

                if(Cube[n_freeParams] < -2 && Cube[n_freeParams] >= -3) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/3)+1)*3 + Cube[n_freeParams];
                n_freeParams++;
            }
        }
    }

    // 6. Einasto rotation curve : _n : _n1 ~ _n2
    if(TRparam[0]._n_fix == 'T')
    {
        _r = Cube[n_freeParams];

        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0]._n_t; // mean: based on the first bspline fit to the first TRfit results
        _sigma = _mu / 3.0;
        //g_wing = _mu; 

        //_n_max = _mu*4;
        //_sigma_N = g_wing / 3.0; // narrower sigma
        //_sigma_W = (_n_max-g_wing) / 3.0; // wider sigma

        //_r1 = 0.1*_mu;
        //_r2 = 2*_mu;

        //Cube[n_freeParams] = gaussian_prior_nrrho_skew(_r, _mu, _sigma_N, _sigma_W);
        //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);
        Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));
        //Cube[n_freeParams] = logNormalPrior(_r, mu_log, sigma_log);
        //Cube[n_freeParams] = DeltaFunctionPrior(_r, _mu);

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        if(Cube[n_freeParams] == 0) Cube[n_freeParams] = _mu/20.0;
        //if(Cube[n_freeParams] >= _r2) Cube[n_freeParams] = _r2;

        n_freeParams++;
    }

    // 7. Einasto rotation curve : r_2 : r_21 ~ r_22
    if(TRparam[0].r_2_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0]._r_2_t; // mean: based on the first bspline fit to the first TRfit results
        _sigma = _mu / 3.0;
        //_sigma = 10000;
        //g_wing = _mu; 

        //_r_max = _mu * 4.0;
        //_sigma_N = g_wing / 3.0; // narrower sigma
        //_sigma_N = 4009.8; // mean_THINGS (17 galaxies in Chemin et al. 2012) / 3.0 in pc
        //_sigma_W = (_r_max-g_wing) / 3.0; // wider sigma
        //_sigma_W = _sigma_N; // use the same sigma

        //_r1 = TRparam[0].r_21;
        //_r2 = TRparam[0].r_22;

        //Cube[n_freeParams] = gaussian_prior_nrrho_skew(_r, _mu, _sigma_N, _sigma_W);
        //Cube[n_freeParams] = DeltaFunctionPrior(_r, _mu);
        //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);
        Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));

        //mu_log = _mu;
        //sigma_log = 0.5;
        //Cube[n_freeParams] = logNormalPrior(_r, mu_log, sigma_log);

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        //if(Cube[n_freeParams] < _r1) Cube[n_freeParams] = _r1;
        //if(Cube[n_freeParams] > _r2) Cube[n_freeParams] = _r2;

        if(Cube[n_freeParams] == 0) Cube[n_freeParams] = _mu/20.0;
        //if(Cube[n_freeParams] >= _r2) Cube[n_freeParams] = _r2;

        n_freeParams++;
    }

    // 8. Einasto rotation curve : rho_2 : rho_21 ~ rho_22
    if(TRparam[0].rho_2_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0]._rho_2_t; // mean: based on the first bspline fit to the first TRfit results
        _sigma = _mu / 3.0;
        //g_wing = _mu; 

        //_r1 = TRparam[0].rho_21;
        //_r2 = TRparam[0].rho_22;

        //Cube[n_freeParams] = gaussian_prior_nrrho_skew(_r, _mu, _sigma_N, _sigma_W);
        //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        //Cube[n_freeParams] = DeltaFunctionPrior(_r, _mu);
        //Cube[n_freeParams] = loguniform_priors(_r, _r1, _r2);
        Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));

        //mu_log = _mu;
        //sigma_log = 0.5;
        //Cube[n_freeParams] = logNormalPrior(_r, mu_log, sigma_log);

        //Cube[n_freeParams] = JeffreysPrior(_r, _r1, _r2);
        //if(Cube[n_freeParams] < _r1) Cube[n_freeParams] = _r1;
        //if(Cube[n_freeParams] > _r2) Cube[n_freeParams] = _r2;

        if(Cube[n_freeParams] == 0) Cube[n_freeParams] = _mu/20.0;
        //if(Cube[n_freeParams] >= _r2) Cube[n_freeParams] = _r2;

        n_freeParams++;
    }

    // 9. VRAD prior :
    if(TRparam[0].vrad_fix == 'T')
    {
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
        {
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _r = Cube[n_freeParams];

                // update priors based on either the tr or the dirty fits 
                _mu = TRparam[0]._vr_bs_t[vrad_nS]; // mean: based on the first bspline fit to the first TRfit results
                _sigma = fabs(_mu)/3.0;
                Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

                n_freeParams++;
            }
        }
        else if(strcmp(TRparam[0].vrad_function, "const") == 0)
        {
            Cube[n_freeParams] = TRparam[0].vrad1 + Cube[n_freeParams]*(TRparam[0].vrad2-TRparam[0].vrad1);
            n_freeParams++;
        }
    }

    // 10. sigma_factor prior : sigma_factor1 ~ sigma_factor2
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0].sigma_factor_mode;
        _sigma = 1*TRparam[0].sigma_factor_std;

        _r1 = TRparam[0].sigma_factor1;
        _r2 = TRparam[0].sigma_factor2;

        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        sigma_factor = Cube[n_freeParams];
        n_freeParams++;
    }
    else
        sigma_factor = TRparam[0].sigma_factor; // use the best fit value from Vmodel_Einasto fitting

    // 11. e_sigma_fitted : 
    //if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0)
    if(TRparam[0].e_sigma >= 0)
    {
        _r = Cube[n_freeParams];

        if(TRparam[0].fullFit == 'F') // full einasto fit : update priors based on the dirty fit
        {
        // update priors based on either the dirty fits 
            _mu = TRparam[0].e_sigma_fitted_t;
            _sigma = _mu/3.;
            Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));
            if(Cube[n_freeParams] == 0) Cube[n_freeParams] = _mu/20.0;
        }
        else // dirty einasto fit
        {
        // update priors based on either the tr fits 
	/*
            _mu = TRparam[0].e_sigma_student_TR[(int)(TRparam[0].Nrings/2)]; // take the one at an intermediate ring
            _mu_max = 10*_mu;
            _sigma_N = _mu / 3.;
            _sigma_W = (_mu_max-_mu) / 5.; 
            Cube[n_freeParams] = gaussian_prior_esigma_skew(_r, _mu, _mu_max, _sigma_N, _sigma_W);
	*/
            _mu = TRparam[0].e_sigma_student_TR[(int)(TRparam[0].Nrings/2)]; // take the one at an intermediate ring
            _sigma = _mu/3.;
            Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));
            if(Cube[n_freeParams] == 0) Cube[n_freeParams] = _mu/20.0;
        }

        e_sigma = Cube[n_freeParams];
        n_freeParams++;
    }

    // 8. calculate log likelihood value
    chi2 = 0.;
    logsum_errors = 0.;
    sum_errors = 0.;
    NobsPoints_available = 0.;

    // _nu of student-T distribution : _nu = 30 for normal : _nu = 1 recommended for best removing the outliers
    _nu = TRparam[0]._nu_studenT;
    logsum_student_chi = 0;
    logsum_sigma2 = 0;
  
    for(i=0; i<TRparam[0].Npoints_in_tilted_ring; i++)
    {
        i0 = TRparam[0].tilted_ring[i][0];
        j0 = TRparam[0].tilted_ring[i][1];

        einasto_halomodel(Cube, i0, j0, TRparam, &Vmodel_Einasto); // Update Vmodel_Einasto

        //if(TRparam[0].sigma_factor_fix == 'F' && HI_VF_geo_radial_angle_w[0].data[j0][i0] > 0) // constant vlos_e mode : partial fit
        if(TRparam[0].e_sigma >= 0 && HI_VF_geo_radial_angle_w[0].data[j0][i0] > 0) // constant vlos_e mode : partial fit
        {   
            if(HI_VF_fract_navail_nall[0].data[j0][i0] > 0.1) // no blank
            {
                // geometrical weight
                NobsPoints_available += 1*HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] * HI_VF_geo_radial_angle_w[0].data[j0][i0];
                logsum_errors += 0.5*log(e_sigma*e_sigma) * HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] * HI_VF_geo_radial_angle_w[0].data[j0][i0];
                logsum_student_chi += ((1+_nu)/2.0)*log(1.0+pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/e_sigma, 2)/(_nu-2)) * HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] * HI_VF_geo_radial_angle_w[0].data[j0][i0];

                // geometrical weight
                //logsum_errors += 0.5*log(M_PI*_nu*e_sigma*e_sigma) * HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] * HI_VF_geo_radial_angle_w[0].data[j0][i0];
                //NobsPoints_available += 1*HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] * HI_VF_geo_radial_angle_w[0].data[j0][i0];
                //logsum_student_chi += ((1+_nu)/2.0)*log(1.0+pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/e_sigma, 2)/_nu) * HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] * HI_VF_geo_radial_angle_w[0].data[j0][i0];

                // no geometrical weight
                //logsum_errors += 0.5*log(M_PI*_nu*e_sigma*e_sigma) * HI_VF_boxFiltered_vlos_ew[0].data[j0][i0];
                //NobsPoints_available += 1*HI_VF_boxFiltered_vlos_ew[0].data[j0][i0];
                //logsum_student_chi += ((1+_nu)/2.0)*log(1.0+pow((HI_VF_boxFiltered[0].data[j0][i0] - Vmodel_Einasto)/e_sigma, 2)/_nu) * HI_VF_boxFiltered_vlos_ew[0].data[j0][i0];
            }
        }
    }

//log_Likelihood_studenT = NobsPoints_available*(log(gsl_sf_gamma((_nu+1)/2.0)) - log(gsl_sf_gamma(_nu/2.0))) - logsum_errors - logsum_student_chi;
    log_Likelihood_studenT = NobsPoints_available*(log(gsl_sf_gamma((_nu+1)/2.0)) - log(sqrt(M_PI*(_nu-2))*gsl_sf_gamma(_nu/2.0))) - logsum_errors - logsum_student_chi;
    *lnew = log_Likelihood_studenT;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dumper_einasto_halofit(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam)
{
    // convert the 2D Fortran arrays to C arrays
    int i=0, j=0;
    int pa_n=0, pa_nF=0, pa_nS=0;
    int incl_n=0, incl_nF=0, incl_nS=0;
    int vrad_n=0, vrad_nF=0, vrad_nS=0;

    int n_freeParams=0;

    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    double postdist[*nSamples][*nPar + 2];

    for( i = 0; i < *nPar + 2; i++ )
    {
        for( j = 0; j < *nSamples; j++ )
        {
            postdist[j][i] = posterior[0][i * (*nSamples) + j];
        }

    }
    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
    double pLivePts[*nlive][*nPar + 1];
    for( i = 0; i < *nPar + 1; i++ )
    {
        for( j = 0; j < *nlive; j++ )
        {
            pLivePts[j][i] = physLive[0][i * (*nlive) + j];
        }
    }

    // 1. XPOS
    if(TRparam[0].xpos_fix == 'T')
    {
        /* save the best fit for xpos */
        TRparam[0].xposF_EinastoFit = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].xposF_EinastoFit_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    // 2. YPOS
    if(TRparam[0].ypos_fix == 'T')
    {

        /* save the best fit for ypos */
        TRparam[0].yposF_EinastoFit = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].yposF_EinastoFit_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    // 3. VSYS
    if(TRparam[0].vsys_fix == 'T')
    {
        /* save the best fit for vsys */
        TRparam[0].vsysF_EinastoFit = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].vsysF_EinastoFit_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    // 4. PA : bspline
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++) // MAP
            {
                // best fit
                //TRparam[0]._p_bs[pa_nS] = paramConstr[0][*nPar*3+n_freeParams];
                //TRparam[0]._p_bs_e[pa_nS] = paramConstr[0][*nPar+n_freeParams];

                // MAP
                TRparam[0]._p_bs[pa_nS] = paramConstr[0][*nPar*0+n_freeParams];
                TRparam[0]._p_bs_e[pa_nS] = paramConstr[0][*nPar+n_freeParams];
                n_freeParams++;
            }
        }
    }

    // 5. INCL : 
    if(TRparam[0].incl_fix == 'T')
    {
        if(strcmp(TRparam[0].incl_function, "bspline") == 0) // MAP
        {
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                // best fit
                //TRparam[0]._i_bs[incl_nS] = paramConstr[0][*nPar*3+n_freeParams];
                //TRparam[0]._i_bs_e[incl_nS] = paramConstr[0][*nPar+n_freeParams];

                // MAP
                TRparam[0]._i_bs[incl_nS] = paramConstr[0][*nPar*0+n_freeParams];
                TRparam[0]._i_bs_e[incl_nS] = paramConstr[0][*nPar+n_freeParams];
                n_freeParams++;
            }
        }
    }

    // 6. Einasto _n
    if(TRparam[0]._n_fix == 'T') // MAP
    {
        /* save the best fit for Einasto profile _n */
        // best fit
        //TRparam[0]._n = paramConstr[0][*nPar*3+n_freeParams];
        //TRparam[0]._ne = paramConstr[0][*nPar+n_freeParams];

        // MAP
        TRparam[0]._n = paramConstr[0][*nPar*0+n_freeParams];
        TRparam[0]._ne = paramConstr[0][*nPar+n_freeParams];

        n_freeParams++;
    }

    // 7. Einasto r_2
    if(TRparam[0].r_2_fix == 'T') // MAP
    {
        /* save the best fit for Einasto profile r_2 */
        // best fit
        //TRparam[0].r_2 = paramConstr[0][*nPar*0+n_freeParams];
        //TRparam[0].r_2e = paramConstr[0][*nPar+n_freeParams];

        // MAP
        TRparam[0].r_2 = paramConstr[0][*nPar*0+n_freeParams];
        TRparam[0].r_2e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    // 8. Einasto rho_2
    if(TRparam[0].rho_2_fix == 'T') // MAP
    {
        /* save the best fit for Einasto profile rho_2 */
        // best fit
        //TRparam[0].rho_2 = paramConstr[0][*nPar*0+n_freeParams];
        //TRparam[0].rho_2e = paramConstr[0][*nPar+n_freeParams];

        // MAP
        TRparam[0].rho_2 = paramConstr[0][*nPar*0+n_freeParams];
        TRparam[0].rho_2e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    // 7. Einasto r_2 : from the rho_2 - r_2 relation
    //if(TRparam[0].r_2_fix == 'F')
    //{
    //  /* save the best fit for Einasto profile r_2 */
    //  TRparam[0].r_2 = pow(10, (-log10(TRparam[0].rho_2) - 0.81)/1.61) * pow(10, 3); // in pc
    //  TRparam[0].r_2e = pow(10, (-log10(TRparam[0].rho_2e) - 0.81)/1.61) * pow(10, 3); // in pc
    //}

    // 9. VRAD : bspline
    if(TRparam[0].vrad_fix == 'T')
    {
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
        {
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                TRparam[0]._vr_bs[vrad_nS] = paramConstr[0][*nPar*0+n_freeParams];
                TRparam[0]._vr_bs_e[vrad_nS] = paramConstr[0][*nPar+n_freeParams];
                n_freeParams++;
            }
        }
    }

    // 8. sigma_factor
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)
    {
        /* save the best fit for sigma_factor */
        TRparam[0].sigma_factor = paramConstr[0][*nPar*0+n_freeParams];
        TRparam[0].sigma_factor_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }
    else
        TRparam[0].sigma_factor = 1.0;

    // 9. sigma_factor
    if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma >= 0)
    {
        /* save the best fit */
        TRparam[0].e_sigma_fitted = paramConstr[0][*nPar*0+n_freeParams];
        TRparam[0].sigma_factor_e = paramConstr[0][*nPar+n_freeParams];
        // save MAP temoporarily: this is used for updating the priors of full einasto fit later
        TRparam[0].e_sigma_fitted_t = paramConstr[0][*nPar*0+n_freeParams];
        n_freeParams++;
    }

    /* Save current maximum loglikelihood value, evidence and error values */
    // Due to unknown reasons, it's possible to save these three values at the same time. only one variable is allowed...
    TRparam[0].maxLogLikeF = *maxLogLike;
    //TRparam[0].logZF = *logZ;
    //TRparam[0].logZerrF = *logZerr;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void set_nfree_params_einasto_halofit_multinest_student(char xposfix, char yposfix, char vsysfix, char pafix, char inclfix, char _n_fix, char r_2_fix, char rho_2_fix, char vradfix, char sigmafactorfix, TR_ringParameters *TRparam)
{
    int pa_n=0, pa_nF=0, pa_nS=0;
    int incl_n=0, incl_nF=0, incl_nS=0;
    int vrad_n=0, vrad_nF=0, vrad_nS=0;
    int n_freeParams = 0;

    // 1. XPOS
    TRparam[0].xpos_fix = xposfix; 
    if(TRparam[0].xpos_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }

    // 2. YPOS
    TRparam[0].ypos_fix = yposfix;
    if(TRparam[0].ypos_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }

    // 3. VSYS
    TRparam[0].vsys_fix = vsysfix;
    if(TRparam[0].vsys_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }

    // 4. PA
    TRparam[0].pa_fix = pafix;
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
                n_freeParams++;
        }
        else if(strcmp(TRparam[0].vrad_function, "const") == 0)
        {
            n_freeParams++;
        }
    }   
    
    // 5. INCL
    TRparam[0].incl_fix = inclfix;
    //TRparam[0].sersic_part = sersicpart;

    if(TRparam[0].incl_fix == 'T')
    {
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
                n_freeParams++;
        }
        else if(strcmp(TRparam[0].vrad_function, "const") == 0)
        {
            n_freeParams++;
        }
    }   

    // 6. _n
    TRparam[0]._n_fix = _n_fix;
    if(TRparam[0]._n_fix == 'T')
    {
        n_freeParams++;
    }

    // 7. r_2
    TRparam[0].r_2_fix = r_2_fix;
    if(TRparam[0].r_2_fix == 'T')
    {
        n_freeParams++;
    }

    // 8. rho_2
    TRparam[0].rho_2_fix = rho_2_fix;
    if(TRparam[0].rho_2_fix == 'T')
    {
        n_freeParams++;
    }

    // 9. VRAD
    TRparam[0].vrad_fix = vradfix;
    if(TRparam[0].vrad_fix == 'T')
    {
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
        {
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++) // bspline
                n_freeParams++;
        }
        else if(strcmp(TRparam[0].vrad_function, "const") == 0)
        {
            n_freeParams++;
        }
    }

    // 9. sigma_factor
    //if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)
    if(TRparam[0].e_sigma >= 0)
    {
        n_freeParams++;
    }

    // save the total number of n_freeParams for multinest
    TRparam[0].n_freeParams = n_freeParams;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double radius_galaxy_plane(double radius, double x, double y, double pa, double incl)
{
    double xr, yr;
    xr = ( -x * sind( pa ) + y * cosd( pa ) );
    yr = ( -x * cosd( pa ) - y * sind( pa ) ) / cosd( incl );
    return( radius - sqrt( xr * xr + yr * yr ) );
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void init_einastohalo_vrot_intp(TR_ringParameters *TRparam)
{
    int i0=0;

    double rbm_vrot_half1, std_vrot_half1, rbm_vrot_half2, std_vrot_half2, rbm_res, std_res;
    double _a, del_x, del_y;

    if(TRparam[0].Nrings < 8)
        spline_intp(TRparam);
    else
    {
        TRparam[0].Nrings_intp = TRparam[0].Nrings;
        for(i0=0; i0<TRparam[0].Nrings; i0++)
        {
            TRparam[0].ring_intp[i0] = TRparam[0].ring_radius[i0];
            TRparam[0].vrot_intp[i0] = TRparam[0].vrot_temp[i0];

            //if(i0 > 0 && TRparam[0].vrot_intp[i0] < 0.5*TRparam[0].vrot_intp[i0-1]) // to avoid any weired VROT value
            //{
            //    TRparam[0].vrot_intp[i0] = TRparam[0].vrot_intp[i0-1];
            //}
            TRparam[0].vrot_e_intp[i0] = TRparam[0].vrot_temp_e[i0];
        }
    }

    double *vrot_half1 = malloc(sizeof(double) * (int)TRparam[0].Nrings_intp/2);
    double *vrot_e_half1 = malloc(sizeof(double) * (int)TRparam[0].Nrings_intp/2);
    double *vrot_half2 = malloc(sizeof(double) * (TRparam[0].Nrings_intp - (int)TRparam[0].Nrings_intp/2));
    double *vrot_e_half2 = malloc(sizeof(double) * (TRparam[0].Nrings_intp - (int)TRparam[0].Nrings_intp/2));
    double *res = malloc(sizeof(double) * (int)TRparam[0].Nrings_intp);
    double *res_e = malloc(sizeof(double) * (int)TRparam[0].Nrings_intp);

    // first half of vrot_intp
    for(i0=0; i0<(int)TRparam[0].Nrings_intp/2; i0++)
    {
        vrot_half1[i0] = TRparam[0].vrot_intp[i0];
        vrot_e_half1[i0] = TRparam[0].vrot_e_intp[i0];
    }

    // second half of vrot_intp
    for(i0=0; i0<(TRparam[0].Nrings_intp-(int)TRparam[0].Nrings_intp/2); i0++)
    {
        vrot_half2[i0] = TRparam[0].vrot_intp[i0+(int)TRparam[0].Nrings_intp/2];
        vrot_e_half2[i0] = TRparam[0].vrot_e_intp[i0+(int)TRparam[0].Nrings_intp/2];
    }

    //robust_mean_std_histogram_ac(vrot_half1, vrot_e_half1, (int)TRparam[0].Nrings_intp/2, &rbm_vrot_half1, &std_vrot_half1);
    //robust_mean_std_histogram_ac(vrot_half2, vrot_e_half2, (TRparam[0].Nrings_intp-(int)TRparam[0].Nrings_intp/2), &rbm_vrot_half2, &std_vrot_half2);

    robust_mean_std(vrot_half1, (int)TRparam[0].Nrings_intp/2, &rbm_vrot_half1, &std_vrot_half1);
    robust_mean_std(vrot_half2, (TRparam[0].Nrings_intp-(int)TRparam[0].Nrings_intp/2), &rbm_vrot_half2, &std_vrot_half2);

    if(rbm_vrot_half1 <= 0) rbm_vrot_half1 = rbm_vrot_half2/2.0;
    if(rbm_vrot_half2 <= 0) rbm_vrot_half2 = rbm_vrot_half1/2.0;

    del_x = -(TRparam[0].Nrings_intp/2.);
    del_y = rbm_vrot_half1 - rbm_vrot_half2;
    _a = del_y/del_x;

    for(i0=0; i0<(int)TRparam[0].Nrings_intp; i0++)
    {
        res[i0] = TRparam[0].vrot_intp[i0] - (_a*(i0 - TRparam[0].Nrings_intp/4.) + rbm_vrot_half1);
    }

    robust_mean_std_histogram_ac(res, TRparam[0].vrot_e_intp, TRparam[0].Nrings_intp, &rbm_res, &std_res);

    for(i0=0; i0<(int)TRparam[0].Nrings_intp; i0++)
    {
        if(fabs(res[i0]) > 1.5*std_res)
        {
            TRparam[0].vrot_intp[i0] = _a*(i0 - TRparam[0].Nrings_intp/4.) + rbm_vrot_half1;
        }

        if(TRparam[0].vrot_intp[i0] <= 0)
        {
            if(rbm_vrot_half1 < rbm_vrot_half2)
                TRparam[0].vrot_intp[i0] = rbm_vrot_half1;
            else
                TRparam[0].vrot_intp[i0] = rbm_vrot_half2;
        }
    }

    // initialise Einasto priors for the first run
    // n
    TRparam[0]._n_t = 2;
    TRparam[0]._n1_t = 1E-3;
    TRparam[0]._n2_t = 1E1;
    TRparam[0].Ein_n_min = TRparam[0]._n1_t;
    TRparam[0].Ein_n_max = TRparam[0]._n2_t;

    // r_2
    TRparam[0]._r_2_t = 1*TRparam[0].rGalaxyPlane_pixel_max * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // pix to pc assuming 1 Mpc;
    TRparam[0].r_21_t = 1E-6*TRparam[0].rGalaxyPlane_pixel_max * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // pix to pc assuming 1 Mpc;; // in pc
    TRparam[0].r_22_t = 1E6*TRparam[0].rGalaxyPlane_pixel_max * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // pix to pc assuming 1 Mpc;
    TRparam[0].r_21_t = 1E-10; // in pc
    TRparam[0].r_22_t = 1E10; // in pc
    TRparam[0].Ein_r_2_min = TRparam[0].r_21_t;
    TRparam[0].Ein_r_2_max = TRparam[0].r_22_t;

    // rho_2
    TRparam[0]._rho_2_t = 1E1;
    //TRparam[0].rho_21_t = pow(2.718281828459045, -1.61*log(TRparam[0].r_22_t));
    //TRparam[0].rho_22_t = pow(2.718281828459045, -1.61*log(TRparam[0].r_21_t));
    TRparam[0].rho_21_t = 1E-10;
    TRparam[0].rho_22_t = 1E10;
    TRparam[0].Ein_rho_2_min = TRparam[0].rho_21_t;
    TRparam[0].Ein_rho_2_max = TRparam[0].rho_22_t;

    // e_sigma
    TRparam[0].e_sigma_fitted_t = 5; // in km/s

    free(vrot_half1);
    free(vrot_e_half1);
    free(vrot_half2);
    free(vrot_e_half2);
    free(res);
    free(res_e);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int einasto_halomodel(double *Cube, int i, int j, TR_ringParameters *TRparam, double *Vmodel_Einasto) // Update TRmodels
{
    int pa_n=0, pa_nF=0, pa_nS=0;
    int bspline_order;
    int incl_n=0, incl_nF=0, incl_nS=0;
    int vrad_n=0, vrad_nF=0, vrad_nS=0;
    int  n=0, n_freeParams=0;
    int pa_order_, incl_order_, root_n;
    double cos_theta, sin_theta;
    double XPOS, YPOS, VSYS, PA, INCL, VROT, VRAD, GG;
    double _n, r_2, rho_2, _x, e, iGamma;
    double _p[999], _i[999], kapa, alpha, sersic_n;
    double _p_bs[99];
    double _i_bs[99];
    double _vr_bs[99];
    double pixelScale=0.0;
    double x, rGalaxyPlane_pixel, rGalaxyPlane_arcsec, rGalaxyPlane_pc, perimeter, perimeter_max, a, b;
    double costheta, sintheta; 
    GG = 4.302 * pow(10, -3); // pc Msol^-1 (km/s)^2

    // parameters for nonlinear equation solver
    double rGalaxyPlane_pixel_max;
    double machep;
    double machine_epsilon;
    double theta, costh, r, ri, ro, ring0, ring1, Npoints_in_tiltedRing;
    double fx_at_starting_point=0;
    double fx_at_end_point=0;
    double _r, r_pix;
    double x_hi_init;
    double incl_intp;
    double _ring_temp, _value_temp, _e_temp, r_pc;

    int ii, jj, i0, j0;


    // B-spline variables
    int ncoeffs;
    int nbreak;
    double xi, yi, yerr, Bj;
    double chisq, Rsq, dof, tss;
    gsl_bspline_workspace *_bw_pa, *_bw_incl, *_bw_vrad;
    gsl_vector *_B_pa, *_B_incl, *_B_vrad;
    gsl_vector *_c_pa, *_c_incl, *_c_vrad, *_w_pa, *_w_incl, *_w_vrad;
    gsl_vector *_x_pa, *_x_incl, *_x_vrad, *_y_pa, *_y_incl, *_y_vrad;
    gsl_matrix *_X_pa, *_X_incl, *_X_vrad, *_cov_pa, *_cov_incl, *_cov_vrad;
    gsl_multifit_linear_workspace *_mw_pa, *_mw_incl, *_mw_vrad;

    //n = TRparam[0].N_reliable_rings;
    n = TRparam[0].Nrings;
    // dynamic 1D array
    double *x_dat = malloc(sizeof(double) * n);
    double *y_dat = malloc(sizeof(double) * n);
    double *e_dat = malloc(sizeof(double) * n);

    if(strcmp(TRparam[0].pa_function, "bspline") == 0)
    {
        bspline_order = TRparam[0].pa_order_bspline;
        nbreak = TRparam[0].pa_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_pa = ncoeffs;

        /* allocate a cubic bspline workspace (k = 4) */
        _bw_pa = gsl_bspline_alloc(bspline_order+1, nbreak);
        _B_pa = gsl_vector_alloc(ncoeffs);
        _x_pa = gsl_vector_alloc(n);
        _y_pa = gsl_vector_alloc(n);
        _X_pa = gsl_matrix_alloc(n, ncoeffs);
        _c_pa = gsl_vector_alloc(ncoeffs);
        _w_pa = gsl_vector_alloc(n);
        _cov_pa = gsl_matrix_alloc(ncoeffs, ncoeffs);
        _mw_pa = gsl_multifit_linear_alloc(n, ncoeffs);
    }

    if(strcmp(TRparam[0].incl_function, "bspline") == 0)
    {
        bspline_order = TRparam[0].incl_order_bspline;
        nbreak = TRparam[0].incl_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_incl = ncoeffs;

        /* allocate a cubic bspline workspace (k = 4) */
        _bw_incl = gsl_bspline_alloc(bspline_order+1, nbreak);
        _B_incl = gsl_vector_alloc(ncoeffs);
        _x_incl = gsl_vector_alloc(n);
        _y_incl = gsl_vector_alloc(n);
        _X_incl = gsl_matrix_alloc(n, ncoeffs);
        _c_incl = gsl_vector_alloc(ncoeffs);
        _w_incl = gsl_vector_alloc(n);
        _cov_incl = gsl_matrix_alloc(ncoeffs, ncoeffs);
        _mw_incl = gsl_multifit_linear_alloc(n, ncoeffs);
    }

    if(strcmp(TRparam[0].vrad_function, "bspline") == 0 && TRparam[0].vrad_nbreak_bspline > 1 && TRparam[0].vrad_order_bspline >= 0)
    {
        bspline_order = TRparam[0].vrad_order_bspline;
        nbreak = TRparam[0].vrad_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_vrad = ncoeffs;

        /* allocate a cubic bspline workspace (k = 4) */
        _bw_vrad = gsl_bspline_alloc(bspline_order+1, nbreak);
        _B_vrad = gsl_vector_alloc(ncoeffs);
        _x_vrad = gsl_vector_alloc(n);
        _y_vrad = gsl_vector_alloc(n);
        _X_vrad = gsl_matrix_alloc(n, ncoeffs);
        _c_vrad = gsl_vector_alloc(ncoeffs);
        _w_vrad = gsl_vector_alloc(n);
        _cov_vrad = gsl_matrix_alloc(ncoeffs, ncoeffs);
        _mw_vrad = gsl_multifit_linear_alloc(n, ncoeffs);
    }

    pixelScale = TRparam[0].pixelScale; // arcsec/pixel
    rGalaxyPlane_pixel_max = TRparam[0].rGalaxyPlane_pixel_max;

    // XPOS
    if(TRparam[0].xpos_fix == 'T')
    {
        XPOS = Cube[n_freeParams]; // xpos fit
        TRparam[0].xposF = XPOS;
        n_freeParams++;
    }
    else
    {
        XPOS = TRparam[0].xposF;
    }

    // YPOS
    if(TRparam[0].ypos_fix == 'T')
    {
        YPOS = Cube[n_freeParams]; // ypos fit
        TRparam[0].yposF = YPOS;
        n_freeParams++;
    }
    else
    {   
        YPOS = TRparam[0].yposF;
    }

    // VSYS
    if(TRparam[0].vsys_fix == 'T')
    {
        VSYS = Cube[n_freeParams]; // vsys fit
        n_freeParams++;
    }
    else
    {
        VSYS = TRparam[0].vsysF;
    }

    // PA
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _p_bs[pa_nS] = Cube[pa_nS+n_freeParams];
                TRparam[0]._p_bs[pa_nS] = _p_bs[pa_nS];
            }
            n_freeParams += TRparam[0].n_coeffs_bspline_pa;
        }
    }

    if(TRparam[0].incl_fix == 'T')
    {   
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _i_bs[incl_nS] = Cube[incl_nS+n_freeParams];
                TRparam[0]._i_bs[incl_nS] = _i_bs[incl_nS];
            }
            n_freeParams += TRparam[0].n_coeffs_bspline_incl;
        }
    }

    // Einasto _n
    if(TRparam[0]._n_fix == 'T')
    {
        _n = Cube[n_freeParams];
        n_freeParams++;
    }
    else
        _n = TRparam[0]._n;

    // Einasto r_2 
    if(TRparam[0].r_2_fix == 'T')
    {
        r_2 = Cube[n_freeParams];
        n_freeParams++;
    }
    else
        r_2 = TRparam[0].r_2;

    // Einasto rho_2 
    if(TRparam[0].rho_2_fix == 'T')
    {
        rho_2 = Cube[n_freeParams];
        n_freeParams++;
    }
    else
        rho_2 = TRparam[0].rho_2;


    // update r_2 with the value from the rho_2 - r_2 relation derived using THINGS (Chemin et al. 2012)
    // log(rho_2) = -1.61 log(r_2) - 0.81 (using a diet Salpeter IMF)
    // Einasto r_2 
    //if(TRparam[0].r_2_fix == 'F')
    //{
    //    r_2 = pow(10, (-log10(rho_2) - 0.81)/1.61) * pow(10, 3); // in pc
    //}

    // VRAD
    if(TRparam[0].vrad_fix == 'T')
    {   
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0 && TRparam[0].vrad_nbreak_bspline > 1 && TRparam[0].vrad_order_bspline >= 0)
        {
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _vr_bs[vrad_nS] = Cube[vrad_nS+n_freeParams];
                TRparam[0]._vr_bs[vrad_nS] = _vr_bs[vrad_nS];
            }
            n_freeParams += TRparam[0].n_coeffs_bspline_vrad;
        }
    }

    // sigma_factor
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)
    {
        n_freeParams++;
    }

    // I. at starting point PA & INCL:
    if(strcmp(TRparam[0].pa_function, "bspline") == 0) // B-spline
    {
        // at starting point
        r_pix = 0;
        for (i0=0; i0<n; i0++)
        {
            _ring_temp = TRparam[0].ring_radius[i0]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].pa_temp[i0]/TRparam[0].PA_MAX_in_degree; // extrapolate inward
            _e_temp = TRparam[0].pa_temp_e[i0]/TRparam[0].PA_MAX_in_degree; // extrapolate inward

            x_dat[i0] = _ring_temp;
            y_dat[i0] = _value_temp;
            e_dat[i0] = _e_temp;

            gsl_vector_set (_x_pa, i0, x_dat[i0]);
            gsl_vector_set (_y_pa, i0, y_dat[i0]);
            gsl_vector_set (_w_pa, i0, 1/e_dat[i0]);
        }

        /* use uniform breakpoints on [0, 1] */
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_pa);

        /* construct the fit matrix _X */
        for (i0=0; i0<n; i0++)
        {
            xi = gsl_vector_get(_x_pa, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_pa, _bw_pa);

            /* fill in row i of _X */
            for (j0=0; j0<TRparam[0].n_coeffs_bspline_pa; j0++)
            {
                //Bj = gsl_vector_get(_B, j);
                Bj = 0.01;
                gsl_matrix_set(_X_pa, i0, j0, Bj);
            }
        }

        /* construct the fit matrix _cov */
        for (i0=0; i0<TRparam[0].n_coeffs_bspline_pa; i0++)
        {
            xi = gsl_vector_get(_x_pa, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_pa, _bw_pa);

            /* fill in row i of X */
            for (j0=0; j0<TRparam[0].n_coeffs_bspline_pa; j0++)
            {
                // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                //Bj = gsl_vector_get(_B, _j);
                Bj = 0.01;
                gsl_matrix_set(_cov_pa, i0, j0, Bj);
            }
        }


        // load the spline efficients generated by MCMC
        for(i0=0; i0<TRparam[0].n_coeffs_bspline_pa; i0++)
        {
            _c_pa->data[i0] = TRparam[0]._p_bs[i0];
        }
            
        xi = r_pix;
        if(xi < x_dat[0]) xi = x_dat[0];
        if(xi > x_dat[n-1]) xi = x_dat[n-1];

        gsl_bspline_eval(xi, _B_pa, _bw_pa);
        gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
        PA = yi;
        PA = (PA*TRparam[0].PA_MAX_in_degree)*M_PI/180.; // in radian
    }

    if(strcmp(TRparam[0].incl_function, "bspline") == 0) // B-spline
    {
        // at starting point
        r_pix = 0;
        for (i0=0; i0<n; i0++)
        {
            _ring_temp = TRparam[0].ring_radius[i0]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].incl_temp[i0]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward
            _e_temp = TRparam[0].incl_temp_e[i0]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward

            x_dat[i0] = _ring_temp;
            y_dat[i0] = _value_temp;
            e_dat[i0] = _e_temp;

            gsl_vector_set (_x_incl, i0, x_dat[i0]);
            gsl_vector_set (_y_incl, i0, y_dat[i0]);
            gsl_vector_set (_w_incl, i0, 1/e_dat[i0]);
        }

        /* use uniform breakpoints on [0, 1] */
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_incl);

        /* construct the fit matrix _X */
        for (i0=0; i0<n; i0++)
        {
            xi = gsl_vector_get(_x_incl, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_incl, _bw_incl);

            /* fill in row i of _X */
            for (j0=0; j0<TRparam[0].n_coeffs_bspline_incl; j0++)
            {
                //Bj = gsl_vector_get(_B, j);
                Bj = 0.01;
                gsl_matrix_set(_X_incl, i0, j0, Bj);
            }
        }


        /* construct the fit matrix _cov */
        for (i0=0; i0<TRparam[0].n_coeffs_bspline_incl; i0++)
        {
            xi = gsl_vector_get(_x_incl, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_incl, _bw_incl);

            /* fill in row i of X */
            for (j0=0; j0<TRparam[0].n_coeffs_bspline_incl; j0++)
            {
                // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                //Bj = gsl_vector_get(_B, _j);
                Bj = 0.01;
                gsl_matrix_set(_cov_incl, i0, j0, Bj);
            }
        }


        // load the spline efficients generated by MCMC
        for(i0=0; i0<TRparam[0].n_coeffs_bspline_incl; i0++)
        {
            _c_incl->data[i0] = TRparam[0]._i_bs[i0];
        }
            
        xi = r_pix;
        if(xi < x_dat[0]) xi = x_dat[0];
        if(xi > x_dat[n-1]) xi = x_dat[n-1];

        gsl_bspline_eval(xi, _B_incl, _bw_incl);
        gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
        INCL = yi;
        INCL = (INCL*TRparam[0].INCL_MAX_in_degree)*M_PI/180.; // in radian
    }
    // calculate fx at _r=0 with PA + INCL given MCMC priors 
    fx_at_starting_point = sqrt(pow((-((double)i-XPOS)*sin(PA) + ((double)j-YPOS)*cos(PA)), 2) + pow(((((double)i-XPOS)*cos(PA) + ((double)j-YPOS)*sin(PA))/cos(INCL)), 2)) - r_pix;


    // II. at end point:
    root_n=0;
    while(1)
    {
        root_n++;
        if(strcmp(TRparam[0].pa_function, "bspline") == 0) // b-spline
        {
            for (i0=0; i0<n; i0++)
            {
                _ring_temp = TRparam[0].ring_radius[i0]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
                _value_temp = TRparam[0].pa_temp[i0]/TRparam[0].PA_MAX_in_degree; // extrapolate inward
                _e_temp = TRparam[0].pa_temp_e[i0]/TRparam[0].PA_MAX_in_degree; // extrapolate inward

                x_dat[i0] = _ring_temp;
                y_dat[i0] = _value_temp;
                e_dat[i0] = _e_temp;

                gsl_vector_set (_x_pa, i0, x_dat[i0]);
                gsl_vector_set (_y_pa, i0, y_dat[i0]);
                gsl_vector_set (_w_pa, i0, 1/e_dat[i0]);
            }

            /* use uniform breakpoints on [0, 1] */
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_pa);
            /* construct the fit matrix _X */
            for (i0=0; i0<n; i0++)
            {
                xi = gsl_vector_get(_x_pa, i0);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_pa, _bw_pa);

                /* fill in row i of _X */
                for (j0=0; j0<TRparam[0].n_coeffs_bspline_pa; j0++)
                {
                    Bj = gsl_vector_get(_B_pa, j0);
                    gsl_matrix_set(_X_pa, i0, j0, Bj);
                }
            }


            /* construct the fit matrix _cov */
            for (i0=0; i0<TRparam[0].n_coeffs_bspline_pa; i0++)
            {
                xi = gsl_vector_get(_x_pa, i0);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_pa, _bw_pa);

                /* fill in row i of X */
                for (j0=0; j0<TRparam[0].n_coeffs_bspline_pa; j0++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B_pa, _j_pa);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_pa, i0, j0, Bj);
                }
            }

            // load the spline efficients generated by MCMC
            for(i0=0; i0<TRparam[0].n_coeffs_bspline_pa; i0++)
            {
                _c_pa->data[i0] = TRparam[0]._p_bs[i0];
            }
                
            // increase r_pix in steps of 0.1*r_max
            r_pix = 0.05*rGalaxyPlane_pixel_max*root_n;
            xi = r_pix / rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];

            gsl_bspline_eval(xi, _B_pa, _bw_pa);
            gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            PA = yi;
            PA = (PA*TRparam[0].PA_MAX_in_degree)*M_PI/180.; // in radian
        }


        if(strcmp(TRparam[0].incl_function, "bspline") == 0) // b-spline
        {
            for (i0=0; i0<n; i0++)
            {
                _ring_temp = TRparam[0].ring_radius[i0]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
                _value_temp = TRparam[0].incl_temp[i0]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward
                _e_temp = TRparam[0].incl_temp_e[i0]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward

                x_dat[i0] = _ring_temp;
                y_dat[i0] = _value_temp;
                e_dat[i0] = _e_temp;

                gsl_vector_set (_x_incl, i0, x_dat[i0]);
                gsl_vector_set (_y_incl, i0, y_dat[i0]);
                gsl_vector_set (_w_incl, i0, 1/e_dat[i0]);
            }

            /* use uniform breakpoints on [0, 1] */
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_incl);

            /* construct the fit matrix _X */
            for (i0=0; i0<n; i0++)
            {
                xi = gsl_vector_get(_x_incl, i0);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_incl, _bw_incl);

                /* fill in row i of _X */
                for (j0=0; j0<TRparam[0].n_coeffs_bspline_incl; j0++)
                {
                    Bj = gsl_vector_get(_B_incl, j0);
                    gsl_matrix_set(_X_incl, i0, j0, Bj);
                }
            }

            /* construct the fit matrix _cov */
            for (i0=0; i0<TRparam[0].n_coeffs_bspline_incl; i0++)
            {
                xi = gsl_vector_get(_x_incl, i0);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_incl, _bw_incl);

                /* fill in row i of X */
                for (j0=0; j0<TRparam[0].n_coeffs_bspline_incl; j0++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B, _j);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_incl, i0, j0, Bj);
                }
            }

            // load the spline efficients generated by MCMC
            for(i0=0; i0<TRparam[0].n_coeffs_bspline_incl; i0++)
            {
                _c_incl->data[i0] = TRparam[0]._i_bs[i0];
            }
                
            // increase r_pix in steps of 0.1*r_max
            r_pix = 0.05*rGalaxyPlane_pixel_max*root_n;
            xi = r_pix / rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];

            gsl_bspline_eval(xi, _B_incl, _bw_incl);
            gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
            INCL = yi;
            INCL = (INCL*TRparam[0].INCL_MAX_in_degree)*M_PI/180.; // in radian
        }
       
        // calculate fx at _r=0 with PA + INCL given MCMC priors 
        fx_at_end_point = sqrt(pow((-((double)i-XPOS)*sin(PA) + ((double)j-YPOS)*cos(PA)), 2) + pow(((((double)i-XPOS)*cos(PA) + ((double)j-YPOS)*sin(PA))/cos(INCL)), 2)) - r_pix;

        if(fx_at_starting_point*fx_at_end_point < 0)
        {
            x_hi_init = r_pix;
            break;
        }
    }

    rGalaxyPlane_pixel = gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver(i, j, TRparam, 0, x_hi_init, 1E5, "brent");
    //rGalaxyPlane_pixel = gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver(i, j, TRparam, 0, x_hi_init, 1E5, "bisection");
    //rGalaxyPlane_pixel = gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver(i, j, TRparam, 0, x_hi_init, 1E5, "falsepos");
    TRparam[0].rGalaxyPlane_pixel = rGalaxyPlane_pixel;
    x = (rGalaxyPlane_pixel/rGalaxyPlane_pixel_max)*1.0*M_PI;

    /* Calculate Einasto halo rotation velocity */
    r_pc = rGalaxyPlane_pixel * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // in pc
    _x = r_pc/r_2;

    iGamma =  gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);
    e = 2.71828; // natural number
    //VROT = sqrt((4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/rGalaxyPlane_pixel)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma));
    VROT = sqrt((4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/r_pc)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma));
    TRparam[0].vrot0[0] = VROT;

    // PA
    if(strcmp(TRparam[0].pa_function, "bspline") == 0) // b-spline
    {
        // load the spline efficients generated by MCMC
        for(i0=0; i0<TRparam[0].n_coeffs_bspline_pa; i0++)
        {
            _c_pa->data[i0] = TRparam[0]._p_bs[i0];
        }
            
        xi = rGalaxyPlane_pixel/rGalaxyPlane_pixel_max;
        if(xi < x_dat[0]) xi = x_dat[0];
        if(xi > x_dat[n-1]) xi = x_dat[n-1];

        gsl_bspline_eval(xi, _B_pa, _bw_pa);
        gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
        PA = yi*TRparam[0].PA_MAX_in_degree;
    }

    
    // INCL
    if(strcmp(TRparam[0].incl_function, "bspline") == 0) // b-spline
    {
        // load the spline efficients generated by MCMC
        for(i0=0; i0<TRparam[0].n_coeffs_bspline_incl; i0++)
        {
            _c_incl->data[i0] = TRparam[0]._i_bs[i0];
        }
            
        xi = rGalaxyPlane_pixel/rGalaxyPlane_pixel_max;
        if(xi < x_dat[0]) xi = x_dat[0];
        if(xi > x_dat[n-1]) xi = x_dat[n-1];

        gsl_bspline_eval(xi, _B_incl, _bw_incl);
        gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
        INCL = yi*TRparam[0].INCL_MAX_in_degree; // in degree
    }

    // VRAD
    if(strcmp(TRparam[0].vrad_function, "bspline") == 0 && TRparam[0].vrad_nbreak_bspline > 1 && TRparam[0].vrad_order_bspline >= 0) // B-spline
    {
        for (i0=0; i0<n; i0++)
        {
            _ring_temp = TRparam[0].ring_radius[i0]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].vrad_temp[i0]/TRparam[0].vrad_max; // extrapolate inward
            _e_temp = TRparam[0].vrad_temp_e[i0]/TRparam[0].vrad_max; // extrapolate inward

            x_dat[i0] = _ring_temp;
            y_dat[i0] = _value_temp;
            e_dat[i0] = _e_temp;

            gsl_vector_set (_x_vrad, i0, x_dat[i0]);
            gsl_vector_set (_y_vrad, i0, y_dat[i0]);
            gsl_vector_set (_w_vrad, i0, 1/e_dat[i0]);
        }

        /* use uniform breakpoints on [0, 1] */
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_vrad);

        /* construct the fit matrix _X */
        for (i0=0; i0<n; i0++)
        {
            xi = gsl_vector_get(_x_vrad, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);

            /* fill in row i of _X */
            for (j0=0; j0<TRparam[0].n_coeffs_bspline_vrad; j0++)
            {
                //Bj = gsl_vector_get(_B, j);
                Bj = 0.01;
                gsl_matrix_set(_X_vrad, i0, j0, Bj);
            }
        }


        /* construct the fit matrix _cov */
        for (i0=0; i0<TRparam[0].n_coeffs_bspline_vrad; i0++)
        {
            xi = gsl_vector_get(_x_vrad, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);

            /* fill in row i of X */
            for (j0=0; j0<TRparam[0].n_coeffs_bspline_vrad; j0++)
            {
                // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                //Bj = gsl_vector_get(_B, _j);
                Bj = 0.01;
                gsl_matrix_set(_cov_vrad, i0, j0, Bj);
            }
        }

        // load the spline efficients generated by MCMC
        for(i0=0; i0<TRparam[0].n_coeffs_bspline_vrad; i0++)
        {
            _c_vrad->data[i0] = TRparam[0]._vr_bs[i0];
        }
            
        xi = rGalaxyPlane_pixel/rGalaxyPlane_pixel_max;
        if(xi < x_dat[0]) xi = x_dat[0];
        if(xi > x_dat[n-1]) xi = x_dat[n-1];

        gsl_bspline_eval(xi, _B_vrad, _bw_vrad);
        gsl_multifit_linear_est(_B_vrad, _c_vrad, _cov_vrad, &yi, &yerr);
        VRAD = yi*TRparam[0].vrad_max;
    }
    else
        VRAD = 0; // non-fitting for VRAD


    // set PERIMETERS
/*
    if(strcmp(TRparam[0].incl_function, "poly-sersic") == 0 && TRparam[0].SersicPoly_order == 0) // calculate perimeters based on ellipse fit 
    {
        a = TRparam[0].ellipse_semi_mx_boxfiltered;
        b = a * cos(TRparam[0].ellipse_incl_boxfiltered*M_PI/180.);
        perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
        a = sqrt(pow((TRparam[0].ellipse_xpos_boxfiltered-i), 2) + pow((TRparam[0].ellipse_ypos_boxfiltered-j)/cos(TRparam[0].ellipse_incl_boxfiltered*M_PI/180.), 2));
        b = a * cos(TRparam[0].ellipse_incl_boxfiltered*M_PI/180.);
        perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
    }
    else if(strcmp(TRparam[0].incl_function, "fourier") == 0 && TRparam[0].incl_order_Fourier == 0) // calculate perimeters based on ellipse fit 
    {
        a = TRparam[0].rGalaxyPlane_pixel_max;
        b = a * cos(INCL*M_PI/180.);
        perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
        a = TRparam[0].rGalaxyPlane_pixel;
        b = a * cos(INCL*M_PI/180.);
        perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
    }
    else if(strcmp(TRparam[0].incl_function, "bspline") == 0) // calculate perimeters based on ellipse fit 
    {
        //a = TRparam[0].rGalaxyPlane_pixel_max;
        //b = a * cos(INCL*M_PI/180.);
        //perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
        //a = TRparam[0].rGalaxyPlane_pixel;
        //b = a * cos(INCL*M_PI/180.);
        //perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));

        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, TRparam[0].Nrings);
        gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl, TRparam[0].Nrings);

        a = TRparam[0].rGalaxyPlane_pixel_max;
        b = a * cos(TRparam[0].incl[TRparam[0].Nrings-1]*M_PI/180.);
        perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));

        a = TRparam[0].rGalaxyPlane_pixel;            
        if(a >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
        {
            incl_intp = TRparam[0].incl[TRparam[0].Nrings-1]; // outermost incl
        }
        else if(a <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
        {
            incl_intp = TRparam[0].incl[0]; // outermost incl
        }
        else
        {
            incl_intp = gsl_spline_eval (spline, a, acc);
        }

        b = a * cos(incl_intp*M_PI/180.);
        perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
    else if(strcmp(TRparam[0].incl_function, "poly-sersic") == 0 && TRparam[0].SersicPoly_order != 0) // use the constant incl derived with pa_order=0 & incl_order=0 options
    {
        a = TRparam[0].ellipse_semi_mx_boxfiltered;
        b = a * cos(TRparam[0]._i[0]*TRparam[0].INCL_MAX_in_degree*M_PI/180.); // use the constant incl derived with pa_order=0 & incl_order=0 options
        perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
        a = sqrt(pow((TRparam[0].xposF_EinastoFit-i), 2) + pow((TRparam[0].yposF_EinastoFit-j)/cos(TRparam[0]._i[0]*TRparam[0].INCL_MAX_in_degree*M_PI/180.), 2));
        b = a * cos(TRparam[0]._i[0]*TRparam[0].INCL_MAX_in_degree*M_PI/180.);
        perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
    }
    else if(strcmp(TRparam[0].incl_function, "fourier") == 0 && TRparam[0].incl_order_Fourier != 0) // use the constant incl derived with pa_order=0 & incl_order=0 options
    {
        a = TRparam[0].rGalaxyPlane_pixel_max*5;
        b = a * cos(INCL*M_PI/180.);
        perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
        a = TRparam[0].rGalaxyPlane_pixel;
        b = a * cos(INCL*M_PI/180.);
        perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
    }
    else
    {
        a = TRparam[0].rGalaxyPlane_pixel_max;
        b = a * cos(INCL*M_PI/180.);
        perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
        a = TRparam[0].rGalaxyPlane_pixel;
        b = a * cos(INCL*M_PI/180.);
        perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
    }
    //TRparam[0].perimeter = perimeter/perimeter_max; // this is for perimeter weigthing!
*/


    cos_theta = (((double)j-YPOS)*cos(PA*M_PI/180.) - ((double)i-XPOS)*sin(PA*M_PI/180.))/rGalaxyPlane_pixel;
    sin_theta = (((double)j-YPOS)*sin(PA*M_PI/180.) + ((double)i-XPOS)*cos(PA*M_PI/180.))/(rGalaxyPlane_pixel*cos(INCL*M_PI/180.));

    /* Calculate a model line-of-sight velocity given ring parameters */
    *Vmodel_Einasto = VSYS + VROT*cos_theta*sin(INCL*M_PI/180.) + VRAD*sin_theta*sin(INCL*M_PI/180.);

    free(x_dat);
    free(y_dat);
    free(e_dat);

    if(strcmp(TRparam[0].pa_function, "bspline") == 0) // b-spline
    {
        gsl_bspline_free(_bw_pa);
        gsl_vector_free(_B_pa);
        gsl_vector_free(_x_pa);
        gsl_vector_free(_y_pa);
        gsl_matrix_free(_X_pa);
        gsl_vector_free(_c_pa);
        gsl_vector_free(_w_pa);
        gsl_matrix_free(_cov_pa);
        gsl_multifit_linear_free(_mw_pa);
    }

    if(strcmp(TRparam[0].incl_function, "bspline") == 0) // b-spline
    {
        gsl_bspline_free(_bw_incl);
        gsl_vector_free(_B_incl);
        gsl_vector_free(_x_incl);
        gsl_vector_free(_y_incl);
        gsl_matrix_free(_X_incl);
        gsl_vector_free(_c_incl);
        gsl_vector_free(_w_incl);
        gsl_matrix_free(_cov_incl);
        gsl_multifit_linear_free(_mw_incl);

    } 

    if(strcmp(TRparam[0].vrad_function, "bspline") == 0 && TRparam[0].vrad_nbreak_bspline > 1 && TRparam[0].vrad_order_bspline >= 0) // B-spline
    {
        gsl_bspline_free(_bw_vrad);
        gsl_vector_free(_B_vrad);
        gsl_vector_free(_x_vrad);
        gsl_vector_free(_y_vrad);
        gsl_matrix_free(_X_vrad);
        gsl_vector_free(_c_vrad);
        gsl_vector_free(_w_vrad);
        gsl_matrix_free(_cov_vrad);
        gsl_multifit_linear_free(_mw_vrad);
    } 
    return 0;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void estimate_vrot_error(multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, char *finalfit, char final_fit)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0, j=0;
    int i0=0, j0=0;
    int x, y;
    int Nrings;
    double efr, tol, Ztol, logZero;
    double total_Npoints_allrings=0.;
    double ri=0., ro=0., ring;

    //double *Vrot_err_inaring = malloc(sizeof(double) * 1); // resize later using realloc in the loop below
    double Vrot_err_inaring[50000] = {0}; // resize later using realloc in the loop below
    double Vrad_err_inaring[50000] = {0}; // resize later using realloc in the loop below
    //TRparam[0].Nrings = (int)((TRparam[0].ring_e-TRparam[0].ring_s)/TRparam[0].ring_w)+2;

    // dynamic total_error 1D array
    double *total_error_ring_params = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array
    double *total_error_ring_params_temp = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array

    // xpos
    double *xpos_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *xpos_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    // xpos error
    double *xpos_err = malloc(sizeof(double) * TRparam[0].Nrings);

    // ypos
    double *ypos_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *ypos_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    // ypos error
    double *ypos_err = malloc(sizeof(double) * TRparam[0].Nrings);

    // vsys
    double *vsys_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vsys_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    // vsys error
    double *vsys_err = malloc(sizeof(double) * TRparam[0].Nrings);

    // pa
    double *pa_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *pa_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    // pa error
    double *pa_err = malloc(sizeof(double) * TRparam[0].Nrings);

    // incl
    double *incl_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *incl_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    // incl error
    double *incl_err = malloc(sizeof(double) * TRparam[0].Nrings);

    // vrot
    double *vrot_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrot_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    // vrot error
    double *vrot_err = malloc(sizeof(double) * TRparam[0].Nrings);

    // vrad
    double *vrad_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrad_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    // vrad error
    double *vrad_err = malloc(sizeof(double) * TRparam[0].Nrings);

    // incomplete Gamma function
    double _e = 2.71828; // natural number
    double GG = 4.302 * pow(10, -3); // pc Msol^-1 (km/s)^2

    // error propagation for einasto halo model
    double _r, _x, _n, r_2, rho_2, iGamma, vEIN, vEIN_bracket_inside, pdiGamma_pdn, pdiGamma_pdr_2, pdvEIN_pdn, dn, pdvEIN_pdr_2, dr_2, pdvEIN_pdrho_2, drho_2, dvEIN, dvEIN_cor, dvEIN_uncor;
    double cov_n_r_2, cov_r_2_rho_2, cov_n_rho_2;

    double delV_apprec_asymmetry;
    TRparam[0].final_fit = final_fit;

    for(i=0; i<TRparam[0].Nrings; i++)
    {
        ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
        ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        // VROT
        if(side == 0) // derive vrotF_e based on the einasto fit 
        {
            _n = TRparam[0]._n;
            _n_einasto_global = _n; // global variable
            r_2 = TRparam[0].r_2;
            rho_2 = TRparam[0].rho_2;
            _r = ring * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // in pc
            // calculate lower incomplete gamma function
            _x = _r/r_2;
            _e = 2.71828;
            iGamma =  gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);

            // einasto halo rotation velocity
            vEIN_bracket_inside = (4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/_r)*(pow(_e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma);
            vEIN = sqrt(vEIN_bracket_inside);

            // partial derivatives of incomplete gamma function with respect to _n & r_2
            pdiGamma_pdn = pdiGamma_pdn_gsl_numerical_integral(0, _r/r_2);
            pdiGamma_pdr_2 = pow(_e, -1.0*_r/r_2)*pow(_r/r_2, 3*_n-1)*(-1.0*_r*pow(r_2, -2));

            // _n
            pdvEIN_pdn = (1.0/2.0)*(1.0/sqrt(vEIN_bracket_inside))*4.0*M_PI*GG*(pow(r_2, 3)/_r)*rho_2*pow(_e, 2*_n)*pow(2*_n, -3*_n) \
                        *(iGamma + 2*_n*iGamma - 3*_n*(log(2)+1+log(_n))*iGamma + _n*pdiGamma_pdn);
            dn = TRparam[0]._ne;
        
            // r_2  
            pdvEIN_pdr_2 = (1.0/2.0)*(1.0/sqrt(vEIN_bracket_inside))*4.0*M_PI*GG*_n*(pow(r_2, 2)/_r)*rho_2*pow(_e, 2*_n)*pow(2*_n, -3*_n) \
                        *(3*iGamma + r_2*pdiGamma_pdr_2);
            dr_2 = TRparam[0].r_2e;

            // rho_2
            pdvEIN_pdrho_2 = (1.0/2.0)*(1.0/sqrt(vEIN_bracket_inside))*4.0*M_PI*GG*_n*(pow(r_2, 3)/_r)*pow(_e, 2*_n)*pow(2*_n, -3*_n)*iGamma;
            drho_2 = TRparam[0].rho_2e;

            if(final_fit == 'Y')
            {
                cov_n_r_2 = read_einasto_posteriors_and_calculate_cov(TRparam, "n", "r_2");
                cov_r_2_rho_2 = read_einasto_posteriors_and_calculate_cov(TRparam, "r_2", "rho_2");
                cov_n_rho_2 = read_einasto_posteriors_and_calculate_cov(TRparam, "n", "rho_2");
                
            // total error propagated from the einasto fit
                dvEIN_cor = sqrt( pow(pdvEIN_pdn*dn, 2) + pow(pdvEIN_pdr_2*dr_2, 2) + pow(pdvEIN_pdrho_2*drho_2, 2) \
                         + 2*pdvEIN_pdn*pdvEIN_pdr_2*cov_n_r_2 \
                         + 2*pdvEIN_pdr_2*pdvEIN_pdrho_2*cov_r_2_rho_2 \
                         + 2*pdvEIN_pdn*pdvEIN_pdrho_2*cov_n_rho_2 );

                //dvEIN = sqrt( pow(pdvEIN_pdn*dn, 2) + pow(pdvEIN_pdr_2*dr_2, 2) + pow(pdvEIN_pdrho_2*drho_2, 2) \
                //       + 2*pdvEIN_pdn*pdvEIN_pdr_2*dn*dr_2 \
                //       + 2*pdvEIN_pdr_2*pdvEIN_pdrho_2*dr_2*drho_2 \
                //       + 2*pdvEIN_pdn*pdvEIN_pdrho_2*dn*drho_2 );

                dvEIN_uncor = sqrt( pow(pdvEIN_pdn*dn, 2) + pow(pdvEIN_pdr_2*dr_2, 2) + pow(pdvEIN_pdrho_2*drho_2, 2));

                // uncertainties by geometry
                // uncertainties by dispersions
                // uncertainties by asymmetry
                delV_apprec_asymmetry = (TRparam[0].vrot_app[i] - TRparam[0].vrot_rec[i])/4.0; // see de Blok et al. (2008)


                if(isinf(dvEIN_cor) || isnan(dvEIN_cor) || dvEIN_cor > TRparam[0].vrot_dispersion_error[i]) // add TR fit + Einasto fit errors quadratically
                {
                    //TRparam[0].vrot_e[i] = sqrt(dvEIN_uncor*dvEIN_uncor \
                                            + TRparam[0].vrot_e[i]*TRparam[0].vrot_e[i] \
                                            + delV_apprec_asymmetry*delV_apprec_asymmetry \
                                            + TRparam[0].dispersion_VLOS*TRparam[0].dispersion_VLOS);
                    TRparam[0].vrot_e[i] = sqrt(dvEIN_uncor*dvEIN_uncor \
                                            + delV_apprec_asymmetry*delV_apprec_asymmetry \
                                            + TRparam[0].dispersion_VLOS*TRparam[0].dispersion_VLOS);

                    // uncertainties by geometry
                    TRparam[0].vrot_einasto_error[i] = fabs(dvEIN_uncor);
                }
                else
                {
                    //TRparam[0].vrot_e[i] = sqrt(dvEIN_cor*dvEIN_cor \
                                            + TRparam[0].vrot_e[i]*TRparam[0].vrot_e[i] \
                                            + delV_apprec_asymmetry*delV_apprec_asymmetry \
                                            + TRparam[0].dispersion_VLOS*TRparam[0].dispersion_VLOS);
                    TRparam[0].vrot_e[i] = sqrt(dvEIN_cor*dvEIN_cor \
                                            + delV_apprec_asymmetry*delV_apprec_asymmetry \
                                            + TRparam[0].dispersion_VLOS*TRparam[0].dispersion_VLOS);

                    // uncertainties by geometry
                    TRparam[0].vrot_einasto_error[i] = fabs(dvEIN_cor);
                }

                // uncertainties by asymmetry
                TRparam[0].vrot_asymmetry_error[i] = fabs(delV_apprec_asymmetry);
                // uncertainties by dispersions
                TRparam[0].vrot_dispersion_error[i] = fabs(TRparam[0].dispersion_VLOS);
            }
        }
    }


    free(total_error_ring_params);
    free(total_error_ring_params_temp);

    free(xpos_tr);
    free(xpos_tr_err);
    free(xpos_err);

    free(ypos_tr);
    free(ypos_tr_err);
    free(ypos_err);

    free(vsys_tr);
    free(vsys_tr_err);
    free(vsys_err);

    free(pa_tr);
    free(pa_tr_err);
    free(pa_err);

    free(incl_tr);
    free(incl_tr_err);
    free(incl_err);

    free(vrot_tr);
    free(vrot_tr_err);
    free(vrot_err);

    free(vrad_tr);
    free(vrad_tr_err);
    free(vrad_err);

}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double pdiGamma_pdn_gsl_numerical_integral(double a, double b) // from a to b
{
    double result, error;
    double alpha = 1.0;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    gsl_function FF;
    FF.function = &lower_inc_gamma_variant;
    FF.params = &alpha;

    gsl_integration_qags(&FF, a, b, 0, 1e-7, 1000, w, &result, &error);

    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

    gsl_integration_workspace_free(w);

    return result;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double read_einasto_posteriors_and_calculate_cov(TR_ringParameters *TRparam, char *param1, char *param2)
{
    int i=0, j=0, n_lines;
    int n_params=0;
    int status;
    FILE *ifp1;
    double cov_n_r_2, cov_r_2_rho_2, cov_n_rho_2, cov_param12;

    ifp1 = fopen(TRparam[0].txtfile_multinest_output, "rwb");
    n_lines = countlines(ifp1);

    // dynamic mt_post array
    double **mt_post = (double**)malloc(sizeof(double *) * n_lines); // n_lines x n_params

    // update the number of free parameters to fit : TRparam[0].n_freeParams
    set_nfree_params_einasto_halofit_multinest_student(TRparam[0].xpos_fix, TRparam[0].ypos_fix, TRparam[0].vsys_fix, TRparam[0].pa_fix, TRparam[0].incl_fix, TRparam[0]._n_fix, TRparam[0].r_2_fix, TRparam[0].rho_2_fix, TRparam[0].vrad_fix, TRparam[0].sigma_factor_fix, TRparam);


    for(i=0; i<n_lines; i++)
    {
        mt_post[i] = (double*)malloc(sizeof(double) * (TRparam[0].n_freeParams+2));
        for(j=0; j<TRparam[0].n_freeParams+2; j++)
        {
            mt_post[i][j] = 0;
        }
    }
    // dynamic mt_post temp array
    double *mt_post_t = (double*)malloc(sizeof(double) * n_lines); // n_lines x n_params

    // einasto n
    double *mt_post_n = (double*)malloc(sizeof(double) * n_lines); // n_lines x n_params
    // einasto r_2
    double *mt_post_r_2 = (double*)malloc(sizeof(double) * n_lines); // n_lines x n_params
    // einasto rho_2
    double *mt_post_rho_2 = (double*)malloc(sizeof(double) * n_lines); // n_lines x n_params

    ifp1 = fopen(TRparam[0].txtfile_multinest_output, "rb");
    //  Reading the input DATA fileth
    for(i=0; i<n_lines; i++)
    {
        for(j=0; j<TRparam[0].n_freeParams+2; j++)
        {
            status = fscanf(ifp1, "%lf", &mt_post[i][j]);
            //printf("%e ", mt_post[i][j]);
        }
    }

    n_params=0;
    // xpos
    if(TRparam[0].xpos_fix == 'T')
    {
        for(i=0; i<n_lines; i++)
        {
            mt_post_t[i] = mt_post[i][2+n_params];
        }
        n_params++;
    }

    // ypos
    if(TRparam[0].ypos_fix == 'T')
    {
        for(i=0; i<n_lines; i++)
        {
            mt_post_t[i] = mt_post[i][2+n_params];
        }
        n_params++;
    }

    // vsys
    if(TRparam[0].vsys_fix == 'T')
    {
        for(i=0; i<n_lines; i++)
        {
            mt_post_t[i] = mt_post[i][2+n_params];
        }
        n_params++;
    }

    // pa
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(i=0; i<TRparam[0].n_coeffs_bspline_pa; i++)
            {
                for(j=0; j<n_lines; j++)
                {
                    mt_post_t[j] = mt_post[j][2+n_params+i];
                }
            }
        }
        n_params += TRparam[0].n_coeffs_bspline_pa;
    }

    // incl
    if(TRparam[0].incl_fix == 'T')
    {
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            for(i=0; i<TRparam[0].n_coeffs_bspline_incl; i++)
            {
                for(j=0; j<n_lines; j++)
                {
                    mt_post_t[j] = mt_post[j][2+n_params+i];
                }
            }
        }
        n_params += TRparam[0].n_coeffs_bspline_incl;
    }

    // Einasto _n
    if(TRparam[0]._n_fix == 'T')
    {
        for(i=0; i<n_lines; i++)
        {
            mt_post_n[i] = mt_post[i][2+n_params];
        }
        n_params++;
    }

    // Einasto r_2
    if(TRparam[0].r_2_fix == 'T')
    {
        for(i=0; i<n_lines; i++)
        {
            mt_post_r_2[i] = mt_post[i][2+n_params];
            //printf("%f\n", mt_post_r_2[i]);
        }

        n_params++;
    }

    // Einasto rho_2
    if(TRparam[0].rho_2_fix == 'T')
    {
        for(i=0; i<n_lines; i++)
        {
            mt_post_rho_2[i] = mt_post[i][2+n_params];
            //printf("%f\n", mt_post_rho_2[i]);
        }
        n_params++;
    }

    // VRAD
    if(TRparam[0].vrad_fix == 'T')
    {
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
        {
            for(i=0; i<TRparam[0].n_coeffs_bspline_vrad; i++)
            {
                for(j=0; j<n_lines; j++)
                {
                    mt_post_t[j] = mt_post[j][2+n_params+i];
                }
            }
        }
        n_params += TRparam[0].n_coeffs_bspline_vrad;
    }

    // sigma_factor : only if sigma_factor is fitted
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)
    {
        for(i=0; i<n_lines; i++)
        {
            mt_post_t[i] = mt_post[i][2+n_params];
        }
    }

    if(strcmp(param1, "n") == 0 && strcmp(param2, "r_2") == 0)
    {
        cov_n_r_2 = gsl_stats_covariance_m(mt_post_n, 1, mt_post_r_2, 1, n_lines, TRparam[0]._n, TRparam[0].r_2);
        //cov_n_r_2 = gsl_stats_covariance(mt_post_n, 1, mt_post_r_2, 1, n_lines);
        cov_param12 = cov_n_r_2;
    }
    else if(strcmp(param1, "n") == 0 && strcmp(param2, "rho_2") == 0)
    {
        cov_n_rho_2 = gsl_stats_covariance_m(mt_post_n, 1, mt_post_rho_2, 1, n_lines, TRparam[0]._n, TRparam[0].rho_2);
        //cov_n_rho_2 = gsl_stats_covariance(mt_post_n, 1, mt_post_rho_2, 1, n_lines);
        cov_param12 = cov_n_rho_2;
    }
    else if(strcmp(param1, "r_2") == 0 && strcmp(param2, "rho_2") == 0)
    {
        cov_r_2_rho_2 = gsl_stats_covariance_m(mt_post_r_2, 1, mt_post_rho_2, 1, n_lines, TRparam[0].r_2, TRparam[0].rho_2);
        //cov_r_2_rho_2 = gsl_stats_covariance(mt_post_r_2, 1, mt_post_rho_2, 1, n_lines);
        cov_param12 = cov_r_2_rho_2;
    }

    free(mt_post);
    free(mt_post_t);
    free(mt_post_n);
    free(mt_post_r_2);
    free(mt_post_rho_2);
    fclose(ifp1);
    return cov_param12;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double lower_inc_gamma_variant(double x, void * params)
{
    double alpha = *(double *) params;

    // return the function inside the bracket of the variant lower incomplete gamma function
    // : alpha is a scale parameter used in the GSL?
    // _n_einasto_global is a global parameter being set with the derived _n of einasto fit 

    return alpha*3*log(x)*exp(-x)*pow(x, 3.0*_n_einasto_global-1);

}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// --- End of line
