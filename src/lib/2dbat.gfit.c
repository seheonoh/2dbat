#include "2dbat.gfit.h"

// 2DBAT user defined functions
// Gaussian fit related


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Gfit_multinest(multinest_paramters *multinest_param, TR_ringParameters *TRparam)
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
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;
    strcpy(root, "gfit.");

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    //ndims = 3*TRparam[0].n_gauss+1;
    //nPar = 3*TRparam[0].n_gauss+1;
    //nClsPar = 3*TRparam[0].n_gauss+1;
    ndims = 4; 
    nPar = 4;
    nClsPar = 4;

    int pWrap[ndims];
    for(i=0; i<ndims; i++)
    {
        pWrap[i] = multinest_param[0].pWrap[i];
    }

    /* Calling multinest */
    run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_gfit, dumper_Gfits, TRparam);

}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_gfit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0;
    int n_gauss;

    double chi2 = 0.0;
    double sigma_param = 1;
    double logsum_errors = 0.;
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double npoints = 0.0;
    double gauss_model=0;

    n_gauss = 1;

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. priors for coefficients: x1 ~ x2

    for(i=0; i<n_gauss; i++)
    {
        if(i==0)
        {
            Cube[i] = TRparam[0].g01 + Cube[i]*(TRparam[0].g02-TRparam[0].g01); // coefficient g0
            Cube[1+3*i] = TRparam[0].gA1 + Cube[1+3*i]*(TRparam[0].gA2-TRparam[0].gA1); // coefficients A
            Cube[2+3*i] = TRparam[0].gS1 + Cube[2+3*i]*(TRparam[0].gS2-TRparam[0].gS1); // coefficients Sigma
            Cube[3+3*i] = TRparam[0].gX1 + Cube[3+3*i]*(TRparam[0].gX2-TRparam[0].gX1); // coefficients X
        }
        else
        {
            Cube[1+3*i] = TRparam[0].gA1 + Cube[1+3*i]*(TRparam[0].gA2-TRparam[0].gA1); // coefficients A
            Cube[2+3*i] = TRparam[0].gS1 + Cube[2+3*i]*(TRparam[0].gS2-TRparam[0].gS1); // coefficients Sigma
            Cube[3+3*i] = TRparam[0].gX1 + Cube[3+3*i]*(TRparam[0].gX2-TRparam[0].gX1); // coefficients X
        }
    }

    // 8. Total number of n_hist 
    npoints = (double)TRparam[0].n_hist_post;
    GLLhood0 = -(npoints/2.0)*log(2.0*M_PI); // for individual error

    chi2 = 0.;
    logsum_errors = 0.;
    for(i=0; i<TRparam[0].n_hist_post; i++)
    {
        TRparam[0].hist_ye[i] = 1;
        gauss_model = gauss_function(Cube, n_gauss, TRparam[0].hist_x[i]);
        logsum_errors += log(TRparam[0].hist_ye[i]);
        chi2 += pow(((TRparam[0].hist_y[i] - gauss_model)/TRparam[0].hist_ye[i]), 2);
//printf("%f %f\n", TRparam[0].hist_y[i], gauss_model);
    }
    slhood = GLLhood0 - logsum_errors - chi2/2.0;
    *lnew = slhood;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dumper_Gfits(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam)
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

/*
    for(i=0; i<TRparam[0].n_gauss; i++)
    {
        //save the best fit for Gaussian coefficients
        if(i == 0) // for a single gaussian fit
        {
            TRparam[0].g_param[i] = paramConstr[0][*nPar*2+i];
            TRparam[0].g_param[i+1] = paramConstr[0][*nPar*2+i+1];
            TRparam[0].g_param[i+2] = paramConstr[0][*nPar*2+i+2];
            TRparam[0].g_param[i+3] = paramConstr[0][*nPar*2+i+3];
        }
        else
        {
            TRparam[0].g_param[3*i+1] = paramConstr[0][*nPar*2+3*i+1];
            TRparam[0].g_param[3*i+2] = paramConstr[0][*nPar*2+3*i+2];
            TRparam[0].g_param[3*i+3] = paramConstr[0][*nPar*2+3*i+3];
        }
    }

    // sigma
    for(i=0; i<TRparam[0].n_gauss; i++)
    {
        //save the best fit for Gaussian coefficients
        if(TRparam[0].n_gauss == 1) // for a single gaussian fit
        {
            TRparam[0].g_param[i+4] = paramConstr[0][*nPar*2+i+4];
            TRparam[0].g_param[i+5] = paramConstr[0][*nPar*2+i+5];
            TRparam[0].g_param[i+6] = paramConstr[0][*nPar*2+i+6];
            TRparam[0].g_param[i+7] = paramConstr[0][*nPar*2+i+7];
        }
        else if(TRparam[0].n_gauss == 2) // for double gaussian
        {
            TRparam[0].g_param[i+7] = paramConstr[0][*nPar*2+i+4];
            TRparam[0].g_param[i+8] = paramConstr[0][*nPar*2+i+5];
            TRparam[0].g_param[i+9] = paramConstr[0][*nPar*2+i+6];
            TRparam[0].g_param[i+7] = paramConstr[0][*nPar*2+i+7];
        }
    }
*/

    //save the best fit for Gaussian coefficients
    if(TRparam[0].n_gauss == 1) // for a single gaussian fit
    {
        TRparam[0].g_param[0] = paramConstr[0][*nPar*2+0]; // background
        TRparam[0].g_param[1] = paramConstr[0][*nPar*2+1]; // amplitude
        TRparam[0].g_param[2] = paramConstr[0][*nPar*2+2]; // sigma
        TRparam[0].g_param[3] = paramConstr[0][*nPar*2+3]; // centre

        // standard deviation of the parameters
        TRparam[0].g_param[4] = paramConstr[0][*nPar*1+0];
        TRparam[0].g_param[5] = paramConstr[0][*nPar*1+1];
        TRparam[0].g_param[6] = paramConstr[0][*nPar*1+2];
        TRparam[0].g_param[7] = paramConstr[0][*nPar*1+3];
    }
    else if(TRparam[0].n_gauss == 2) // for double gaussian fit
    {
        TRparam[0].g_param[0] = paramConstr[0][*nPar*2+0]; // primary background
        TRparam[0].g_param[1] = paramConstr[0][*nPar*2+1]; // primary amplitude
        TRparam[0].g_param[2] = paramConstr[0][*nPar*2+2]; // primary sigma
        TRparam[0].g_param[3] = paramConstr[0][*nPar*2+3]; // primary centre
        TRparam[0].g_param[4] = paramConstr[0][*nPar*2+4]; // secondary amplitude
        TRparam[0].g_param[5] = paramConstr[0][*nPar*2+5]; // secondary sigma
        TRparam[0].g_param[6] = paramConstr[0][*nPar*2+6]; // secondary centre

        TRparam[0].g_param[7] = paramConstr[0][*nPar*1+0]; // primary background
        TRparam[0].g_param[8] = paramConstr[0][*nPar*1+1]; // primary amplitude
        TRparam[0].g_param[9] = paramConstr[0][*nPar*1+2]; // primary sigma
        TRparam[0].g_param[10] = paramConstr[0][*nPar*1+3]; // primary centre
        TRparam[0].g_param[11] = paramConstr[0][*nPar*1+4]; // secondary amplitude
        TRparam[0].g_param[12] = paramConstr[0][*nPar*1+5]; // secondary sigma
        TRparam[0].g_param[13] = paramConstr[0][*nPar*1+6]; // secondary centre
    }

    TRparam[0].maxLogLikeF = *maxLogLike;
    //TRparam[0].logZF = *logZ;
    //TRparam[0].logZerrF = *logZerr;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gauss_function(double *Cube, int n_gauss, double x)
{
    int i=0;
    double gauss_model = 0.;

    // G0 = Cube[0] : background
    // A0 = Cube[1+3*n_gauss], sigma = Cube[2+3*n_gauss], x0 = Cube[3+3*n_gauss];

    gauss_model += Cube[0];
    for(i=0; i<n_gauss; i++)
    {
        gauss_model += (Cube[1+3*i]/(sqrt(2.0*M_PI)*Cube[2+3*i])) * exp(-0.5*pow((x-Cube[3+3*i])/Cube[2+3*i], 2));
    }

    return gauss_model;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int gsl_nonlinearfit_gauss(TR_ringParameters *TRparam, int nn, double g0_init, double gA_init, double gX_init, double gS_init, double *g0, double *gA, double *gX, double *gS, int *result)
{
    int i, j, iter = 0, n_iter=0;
    int n = nn;
    int status;
    int available_N=0;
    int poly_order;
    const gsl_multifit_fdfsolver_type *Tol;
    gsl_multifit_fdfsolver *s;
    double ring_s, ring_w, rGalaxyPlane_pixel_max;

    int p = 4;
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double ring_radius[n], y[n], sigma[n];
    gsl_multifit_function_fdf f;
    double x_init[4] = {g0_init, gA_init, gX_init, gS_init};
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    struct fourier_inputdata d = {poly_order, n, ring_radius, y, sigma, ring_s, ring_w, rGalaxyPlane_pixel_max};

    f.f = &gauss_f;
    f.df = &gauss_df;
    f.fdf = &gauss_fdf;

    j=0;
    available_N=0;
    for(i=0; i<n; i++)
    {
        ring_radius[j] = TRparam[0].hist_x[i];
        if(TRparam[0].hist_y[i] == 1E9) // if garbage then skip
        {
            continue;
        }
        else
        {
            y[j] = TRparam[0].hist_y[i];
            sigma[j] = 0.1; // useless
            j++; // increase j index
            available_N++;
        }
    }

    n = available_N; // update n;
    d.n = n;
    f.n = n;
    f.p = p;

    f.params = &d;
    Tol = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (Tol, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

    while(1)
    {
        do
        {
            iter++;
            status = gsl_multifit_fdfsolver_iterate (s);
            status = gsl_multifit_test_delta (s->dx, s->x, 1e-10, 1e-10);
        }
        while (status == GSL_CONTINUE && iter < 1E5);

        if(status == GSL_SUCCESS)
        {
            *result = 1;
            break;
        }
        else
        {
            *result = 0;
            break;
        }
    }

    if(*result == 0)
    {
        *g0 = 0;
        *gA = 0;
        *gX = 0;
        *gS = 0;
        gsl_multifit_fdfsolver_free (s);
        gsl_matrix_free (covar);
        return 0;
    }
    
    gsl_multifit_covar (s->J, 0.0, covar);

    #define FIT(i) gsl_vector_get(s->x, i)
    #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    {
        double chi = gsl_blas_dnrm2(s->f);
        double dof = n - p;
        double c = GSL_MAX_DBL(1, chi / sqrt(dof));
            
        if(isnan(ERR(0)) || isinf(ERR(0)) || isinf(c) || isnan(c) || c*ERR(0) > 1E3 || c*ERR(1) > 1E3 ||
           isnan(ERR(1)) || isinf(ERR(1)) || c*ERR(2) > 1E3 || c*ERR(3) > 1E3 ||
           isnan(ERR(2)) || isinf(ERR(2)) || 
           isnan(ERR(3)) || isinf(ERR(3)))
        {
            printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
            printf ("_g0      = %.5f +/- nan\n", FIT(0));
            printf ("_gA    = %.5f +/- nan\n", FIT(1));
            printf ("_gX    = %.5f +/- nan\n", FIT(2));
            printf ("_gS    = %.5f +/- nan\n", FIT(3));

            // failed to fit
            *g0 = 0;
            *gA = 0;
            *gX = 0;
            *gS = 0;
        }
        else
        {
            if(status == GSL_SUCCESS)
            {
                printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
                printf ("_g0    = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
                printf ("_gA    = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
                printf ("_gX    = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
                printf ("_gS    = %.5f +/- %.5f\n", FIT(3), c*ERR(3));

                *g0 = FIT(0);
                *gA = FIT(1);
                *gX = FIT(2);
                *gS = FIT(3);
            }
            else
            {
                *g0 = 0;
                *gA = 0;
                *gX = 0;
                *gS = 0;
            }
        }
    }
    printf ("status = %s\n", gsl_strerror (status));

    //for(i=0; i<n; i++)
    //{
    //    printf("%f %f %f\n", TRparam[0].hist_x[i], TRparam[0].hist_y[i], FIT(0) + (FIT(1)/(sqrt(2.0*M_PI)*FIT(3))) * exp(-0.5*pow((TRparam[0].hist_x[i]-FIT(2))/FIT(3), 2)));
    //}

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    return 0;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int gauss_f(const gsl_vector * x, void *fourier_inputdata, gsl_vector * f)
{
    int i;
    int n = ((struct fourier_inputdata *)fourier_inputdata)->n;
    double *hist_x = ((struct fourier_inputdata *)fourier_inputdata)->ring_radius;
    double *hist_y = ((struct fourier_inputdata *)fourier_inputdata)->y;
    double *sigma = ((struct fourier_inputdata *) fourier_inputdata)->sigma;

    double _g0 = gsl_vector_get (x, 0);
    double _gA = gsl_vector_get (x, 1);
    double _gX = gsl_vector_get (x, 2);
    double _gS = gsl_vector_get (x, 3);

    for (i = 0; i < n; i++)
    {
        double Yi = _g0 + (_gA/(sqrt(2.0*M_PI)*_gS)) * exp(-0.5*pow((hist_x[i]-_gX)/_gS, 2));
        gsl_vector_set (f, i, (Yi - hist_y[i])/sigma[i]);
    }
    return GSL_SUCCESS;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int gauss_df(const gsl_vector *x, void *fourier_inputdata, gsl_matrix *J)
{
    int i;
    int n = ((struct fourier_inputdata *)fourier_inputdata)->n;
    double *hist_x = ((struct fourier_inputdata *)fourier_inputdata)->ring_radius;
    double *sigma = ((struct fourier_inputdata *) fourier_inputdata)->sigma;

    double _g0 = gsl_vector_get (x, 0);
    double _gA = gsl_vector_get (x, 1);
    double _gX = gsl_vector_get (x, 2);
    double _gS = gsl_vector_get (x, 3);

    for (i = 0; i < n; i++)
    {
        /* Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = (Yi - yi)/sigma[i],      */
        /*       Yi = A * exp(-lambda * i) + b  */
        /* and the xj are the parameters (A,lambda,b) */
        double s = sigma[i];
        double dgauss_d0, dgauss_dA, dgauss_dX, dgauss_dS;

        dgauss_d0 = 1;
        dgauss_dA = (1/(sqrt(2.0*M_PI)*_gS)) * exp(-0.5*pow((hist_x[i]-_gX)/_gS, 2));
        dgauss_dX = (_gA/(sqrt(2.0*M_PI)*_gS)) * exp(-0.5*pow((hist_x[i]-_gX)/_gS, 2))*((hist_x[i]-_gX)/_gS);
        dgauss_dS = -1.0*(_gA/(sqrt(2.0*M_PI)*pow(_gS, 2))) * exp(-0.5*pow((hist_x[i]-_gX)/_gS, 2)) + (_gA/(sqrt(2.0*M_PI)*_gS)) * exp(-0.5*pow((hist_x[i]-_gX)/_gS, 2))*pow((hist_x[i]-_gX), 2)*pow(_gS, -3);

        gsl_matrix_set (J, i, 0, dgauss_d0/s);
        gsl_matrix_set (J, i, 1, dgauss_dA/s);
        gsl_matrix_set (J, i, 2, dgauss_dX/s);
        gsl_matrix_set (J, i, 3, dgauss_dS/s);
    }
    return GSL_SUCCESS;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int gauss_fdf(const gsl_vector * x, void *fourier_inputdata, gsl_vector * f, gsl_matrix * J)
{
    gauss_f (x, fourier_inputdata, f);
    gauss_df (x, fourier_inputdata, J);

    return GSL_SUCCESS;
}


// --- End of line



