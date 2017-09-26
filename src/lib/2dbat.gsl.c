
#include "2dbat.gsl.h"

// 2DBAT user defined functions
// GSL related

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void robust_mean_std(double *input, int n, double *robust_mean, double *robust_std)
{
    int i, n_filtered=0;
    int n_loop=0;
    double robust_mean_t;
    double robust_std_t;
    double filter_l, filter_u;

    // 1. Derive individual errors : robust mode values based on histograms
    if(n == 0)
    {
        //printf("The number of input array is %d\n", n);
        //printf("Put 1E9 for now\n");
        robust_mean_t = 1E9;
        robust_std_t = 1E9;
        return;
    }
    else if(n == 1)
    {
        //printf("The number of input array is %d\n", n);
        robust_mean_t = input[0];
        robust_std_t = 1E9;
        return;
    }

    gsl_sort(input, 1, n);
    robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
    robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);

    // first pass with 3 sigma level 
    filter_l = robust_mean_t - 3.0*robust_std_t;
    filter_u = robust_mean_t + 3.0*robust_std_t;
    if(filter_l >= filter_u)
    {
        filter_l = filter_l - 1E9;
        filter_u = filter_l + 1E9;  
    }

    for(i=0; i<n; i++)
    {
        if(input[i] < filter_l || input[i] > filter_u)
        {
            input[i] = robust_mean_t;
        }
    }

    gsl_sort(input, 1, n);
    robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
    robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);
    n_loop = 0;
    while(1)
    {
        n_loop++;
        n_filtered=0;
        filter_l = robust_mean_t - 3*robust_std_t;
        filter_u = robust_mean_t + 3*robust_std_t;
        for(i=0; i<n; i++)
        {
            if(input[i] < filter_l || input[i] > filter_u)
            {
                n_filtered++;
                input[i] = robust_mean_t;
            }
        }
        if(n_filtered == 0)
            break;

        if(n_loop > 1E3)
            break;

        gsl_sort(input, 1, n);
        robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
        robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);
    }

    if(robust_std_t == 0)
        robust_std_t = 1E9;

    *robust_mean = robust_mean_t;
    *robust_std = robust_std_t;
    //printf("mean_input: %f std_input: %f\n", robust_mean_t, robust_std_t);
    return;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void robust_mean_std_e(double *input, int n, double *robust_mean, double *robust_std)
{
    int i, n_filtered=0;
    int n_loop=0;
    double robust_mean_t;
    double robust_std_t;
    double filter_l, filter_u;

    // 1. Derive individual errors : robust mode values based on histograms
    if(n == 0)
    {
        //printf("The number of input array is %d\n", n);
        //printf("Put 1E9 for now\n");
        robust_mean_t = 1E9;
        robust_std_t = 1E9;
        return;
    }
    else if(n == 1)
    {
        //printf("The number of input array is %d\n", n);
        robust_mean_t = input[0];
        robust_std_t = 1E9;
        return;
    }

    gsl_sort(input, 1, n);
    robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
    //robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);
    robust_std_t = fabs(gsl_stats_quantile_from_sorted_data(input, 1, n, 0.75) - gsl_stats_quantile_from_sorted_data(input, 1, n, 0.25))/1.0;

    // first pass with 1 sigma level 
    //filter_l = robust_mean_t - 1.0*robust_std_t;
    filter_l = 0; // as error is always positive & its lower limit is zero
    filter_u = robust_mean_t + 2.0*robust_std_t;
    if(filter_l >= filter_u)
    {
        filter_l = filter_l - 1E9;
        filter_u = filter_l + 1E9;  
    }

    for(i=0; i<n; i++)
    {
        if(input[i] < filter_l || input[i] > filter_u)
        {
            input[i] = robust_mean_t;
        }
    }

    gsl_sort(input, 1, n);
    robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
    //robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);
    robust_std_t = fabs(gsl_stats_quantile_from_sorted_data(input, 1, n, 0.75) - gsl_stats_quantile_from_sorted_data(input, 1, n, 0.25))/1.0;

    n_loop = 0;
    while(1)
    {
        n_loop++;
        n_filtered=0;
        //filter_l = robust_mean_t - 1*robust_std_t;
        filter_l = 0;
        filter_u = robust_mean_t + 2*robust_std_t;
        for(i=0; i<n; i++)
        {
            if(input[i] < filter_l || input[i] > filter_u)
            {
                n_filtered++;
                input[i] = robust_mean_t;
            }
        }
        if(n_filtered == 0)
            break;

        if(n_loop > 1E3)
            break;

        gsl_sort(input, 1, n);
        robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
        //robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);
        robust_std_t = fabs(gsl_stats_quantile_from_sorted_data(input, 1, n, 0.75) - gsl_stats_quantile_from_sorted_data(input, 1, n, 0.25))/1.0;
    }

    if(robust_std_t == 0)
        robust_std_t = 1E9;

    *robust_mean = robust_mean_t;
    *robust_std = robust_std_t;
    //printf("mean_input: %f std_input: %f\n", robust_mean_t, robust_std_t);
    return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void robust_mean_std_histogram(double *input, int n, double *robust_mean, double *robust_std)
{
    gsl_histogram *h_input;
    double lower_bound_input, upper_bound_input;
    double FD_h, IQR;
    double ll_bin, ul_bin;
    int FD_nbins, h_loop;
    int i=0;

    // 1. Derive individual errors : robust mode values based on histograms
    if(n == 0)
    {
        printf("The number of input array is %d\n", n);
        printf("Put 1E9 for now\n");
        *robust_mean = 1E9;
        *robust_std = 1E9;
        return;
    }
    else if(n == 1)
    {
        printf("The number of input array is %d\n", n);
        *robust_mean = input[0];
        *robust_std = 1E9;
        return;
    }

    gsl_sort(input, 1, n);
    lower_bound_input = input[0];
    upper_bound_input = input[n-1];

    if(lower_bound_input >= upper_bound_input)
    {
        lower_bound_input = lower_bound_input - 2*(fabs(lower_bound_input)+1);
        upper_bound_input = lower_bound_input + 2*(fabs(lower_bound_input)+1);
    }
    // find an optimal n_bins following Freedman-Diaconis rule
    //IQR = gsl_stats_quantile_from_sorted_data(input, 1, n, 0.75)  - gsl_stats_quantile_from_sorted_data(input, 1, n, 0.25);
    //FD_h = 2.0*fabs(IQR)*pow(n, -1/3.);
    //FD_nbins = (int)((upper_bound_input-lower_bound_input)/FD_h);

    // optimal # of histogram bin : Shimazaki and Shinomoto. Neural Comput, 2007, 19(6), 1503-1527
    FD_nbins = sshist(input, n);

    if(FD_nbins > 1E3)
    {
        h_loop=0;
        while(1)
        {
            h_loop++;
            upper_bound_input = gsl_stats_quantile_from_sorted_data(input, 1, n, 0.75-0.01*h_loop);
            IQR = upper_bound_input  - gsl_stats_quantile_from_sorted_data(input, 1, n, 0.25);
            FD_h = 2.0*fabs(IQR)*pow(n, -1/3.);
            FD_nbins = (int)((upper_bound_input-lower_bound_input)/FD_h);
            if(FD_nbins < 1E3)
            {
                break;
            }
            if(h_loop > 50)
            {
                printf("Check histogram FD_nbins in find_optimal_sigma_factor_for_likelihood()\n");
            }
        }
    }

    if(FD_nbins < 3) FD_nbins = 3;
    h_input = gsl_histogram_alloc(FD_nbins);
    gsl_histogram_set_ranges_uniform(h_input, lower_bound_input, upper_bound_input);

    // assign histograms
    for (i=0; i<n; i++)
    {
        gsl_histogram_increment(h_input, input[i]);
    }

    // 2. Find mode values of parameter errors's histograms
    // mean & std of ring parameters histogram 
    gsl_histogram_get_range(h_input, gsl_histogram_max_bin(h_input), &ll_bin, &ul_bin);
    *robust_mean = (ll_bin+ul_bin)/2.0;
    *robust_std = gsl_histogram_sigma(h_input);
    if(*robust_std == 0)
        *robust_std = gsl_stats_sd_m(input, 1, n, *robust_mean);

    //printf("mean_input: %f std_input: %f\n", *robust_mean, *robust_std);
    gsl_histogram_free(h_input);

    return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void robust_mean_std_histogram_ac(double *input, double *input_err, int n, double *robust_mean_ac, double *robust_std_ac)
{
    gsl_histogram *h_input_ac;
    double lower_bound_input, upper_bound_input;
    double FD_h, IQR;
    double ll_bin, ul_bin;
    int FD_nbins, h_loop;
    int i=0;

    // 1. Derive individual errors : robust mode values based on histograms
    if(n == 0)
    {
        printf("The number of input array is %d\n", n);
        printf("Put 1E9 for now\n");
        *robust_mean_ac = 1E9;
        *robust_std_ac = 1E9;
        return;
    }
    else if(n == 1)
    {
        printf("The number of input array is %d\n", n);
        *robust_mean_ac = input[0];
        *robust_std_ac = 1E9;
        return;
    }

    gsl_sort(input, 1, n);
    lower_bound_input = input[0];
    upper_bound_input = input[n-1];

    if(lower_bound_input >= upper_bound_input)
    {
        lower_bound_input = lower_bound_input - 2*(fabs(lower_bound_input)+1);
        upper_bound_input = lower_bound_input + 2*(fabs(lower_bound_input)+1);
    }
    // find an optimal n_bins following Freedman-Diaconis rule
    //IQR = gsl_stats_quantile_from_sorted_data(input, 1, n, 0.75)  - gsl_stats_quantile_from_sorted_data(input, 1, n, 0.25);
    //FD_h = 2.0*fabs(IQR)*pow(n, -1/3.);
    //FD_nbins = (int)((upper_bound_input-lower_bound_input)/FD_h);

    // optimal # of histogram bin : Shimazaki and Shinomoto. Neural Comput, 2007, 19(6), 1503-1527
    FD_nbins = sshist(input, n);

    if(FD_nbins > 1E3)
    {
        h_loop=0;
        while(1)
        {
            h_loop++;
            upper_bound_input = gsl_stats_quantile_from_sorted_data(input, 1, n, 0.75-0.01*h_loop);
            IQR = upper_bound_input  - gsl_stats_quantile_from_sorted_data(input, 1, n, 0.25);
            FD_h = 2.0*fabs(IQR)*pow(n, -1/3.);
            FD_nbins = (int)((upper_bound_input-lower_bound_input)/FD_h);
            if(FD_nbins < 1E3)
            {
                break;
            }
            if(h_loop > 50)
            {
                printf("Check histogram FD_nbins in find_optimal_sigma_factor_for_likelihood()\n");
            }
        }
    }

    if(FD_nbins < 3) FD_nbins = 3;
    h_input_ac = gsl_histogram_alloc(FD_nbins);
    gsl_histogram_set_ranges_uniform(h_input_ac, lower_bound_input, upper_bound_input);

    // assign histograms
    for (i=0; i<n; i++)
    {
        gsl_histogram_accumulate(h_input_ac, input[i], 1/input_err[i]);
    }

    // 2. Find mode values of parameter errors's histograms
    // mean & std of ring parameters histogram 
    gsl_histogram_get_range(h_input_ac, gsl_histogram_max_bin(h_input_ac), &ll_bin, &ul_bin);
    *robust_mean_ac = (ll_bin+ul_bin)/2.0;
    *robust_mean_ac = gsl_histogram_mean(h_input_ac);
    *robust_std_ac = gsl_histogram_sigma(h_input_ac);
    if(*robust_std_ac == 0)
        *robust_std_ac = gsl_stats_sd_m(input, 1, n, *robust_mean_ac);

    //printf("mean_input: %f std_input: %f\n", *robust_mean_ac, *robust_std_ac);
    gsl_histogram_free(h_input_ac);

    return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GSL nonlinear solver: rGalaxyplane with PA INCL bspline approximation
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver(int i, int j, TR_ringParameters *TRparam, double x_lo, double x_hi, int max_iter, char *solver)
{
    int status, m=0, i0=0, j0=0;
    int iter = 0;
    const gsl_root_fsolver_type *Tgsl;
    gsl_root_fsolver *s;
    double r = 0, r_expected = 10;
    gsl_function Fgsl;

    struct rGalaxyPlane_pixel_TR_nonlinearEquation_params params;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // save the current ring parameters from MCMC to the gsl struct param

    if(strcmp(TRparam[0].pa_function, "bspline") == 0)
    {
        params.n_coeffs_bspline_pa = TRparam[0].n_coeffs_bspline_pa;
        params.pa_nbreak_bspline = TRparam[0].pa_nbreak_bspline;
    }

    if(strcmp(TRparam[0].incl_function, "bspline") == 0)
    {
        params.n_coeffs_bspline_incl = TRparam[0].n_coeffs_bspline_incl;
        params.incl_nbreak_bspline = TRparam[0].incl_nbreak_bspline;
    }

    strcpy(params.pa_function, TRparam[0].pa_function);
    strcpy(params.incl_function, TRparam[0].incl_function);
    params.N_reliable_rings = TRparam[0].N_reliable_rings;
    params.rmax = TRparam[0].rGalaxyPlane_pixel_max;
    params.rGalaxyPlane_pixel_max = TRparam[0].rGalaxyPlane_pixel_max;
    params.PA_MAX_in_degree = TRparam[0].PA_MAX_in_degree;
    params.INCL_MAX_in_degree = TRparam[0].INCL_MAX_in_degree;
  
    for (i0=0; i0<params.N_reliable_rings; i0++)
    {
        params.ring_radius[i0] = TRparam[0].ring_radius[i0]; 
        params.pa_temp[i0] = TRparam[0].pa_temp[i0]; 
        params.pa_temp_e[i0] = TRparam[0].pa_temp_e[i0]; 
        params.incl_temp[i0] = TRparam[0].incl_temp[i0]; 
        params.incl_temp_e[i0] = TRparam[0].incl_temp_e[i0]; 
    }

    // 1. sky position (i, j)
    params.i = i;
    params.j = j;
    // 2. XPOS YPOS
    params.xpos = TRparam[0].xposF;
    params.ypos = TRparam[0].yposF;
    // 3. PA
    if(strcmp(TRparam[0].pa_function, "bspline") == 0)
    {
        for(m=0; m<TRparam[0].n_coeffs_bspline_pa; m++)
        {
            params._p_bs[m] = TRparam[0]._p_bs[m];
        } 
    }

    // 4. INCL
    if(strcmp(TRparam[0].incl_function, "bspline") == 0)
    {
        for(m=0; m<TRparam[0].n_coeffs_bspline_incl; m++)
        {
            params._i_bs[m] = TRparam[0]._i_bs[m];
        } 
    }

    // 5. Rmax
    params.rmax = TRparam[0].rGalaxyPlane_pixel_max;

    Fgsl.function = &gsl_rGalaxyPlane_pixel_TR_nonlinearEquation;
    Fgsl.params = &params;

    if((strcmp(solver, "brent") == 0))
        Tgsl = gsl_root_fsolver_brent;
    
    if((strcmp(solver, "bisection") == 0))
        Tgsl = gsl_root_fsolver_bisection;

    if((strcmp(solver, "falsepos") == 0))
        Tgsl = gsl_root_fsolver_falsepos;

    s = gsl_root_fsolver_alloc (Tgsl);
    gsl_root_fsolver_set (s, &Fgsl, x_lo, x_hi);

//    printf ("using %s method\n", gsl_root_fsolver_name (s));
//    printf ("%5s [%9s, %9s] %9s %10s %9s\n",
//            "iter", "lower", "upper", "root",
//            "err", "err(est)");

    do
      {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 1E-10);

        if (status == GSL_SUCCESS) ;
        //  printf ("Converged:\n");

//        printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//                iter, x_lo, x_hi,
//                r, r - r_expected,
//                x_hi - x_lo);
      }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
//    return status;
    return r;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver_given_paincl(int i, int j, double _xpos, double _ypos, double _pa, double _incl, TR_ringParameters *TRparam, double x_lo, double x_hi, int max_iter, char *solver)
{
    int status, m=0, i0=0, j0=0;
    int iter = 0;
    const gsl_root_fsolver_type *Tgsl;
    gsl_root_fsolver *s;
    double r = 0;
    gsl_function Fgsl;

    struct rGalaxyPlane_pixel_TR_nonlinearEquation_params params;

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // save the current ring parameters from MCMC to the gsl struct param

    strcpy(params.pa_function, TRparam[0].pa_function);
    strcpy(params.incl_function, TRparam[0].incl_function);
    params.N_reliable_rings = TRparam[0].N_reliable_rings;
    params.rmax = TRparam[0].rGalaxyPlane_pixel_max;
    params.rGalaxyPlane_pixel_max = TRparam[0].rGalaxyPlane_pixel_max;
    params.PA_MAX_in_degree = TRparam[0].PA_MAX_in_degree;
    params.INCL_MAX_in_degree = TRparam[0].INCL_MAX_in_degree;
  
    for (i0=0; i0<params.N_reliable_rings; i0++)
    {
        params.ring_radius[i0] = TRparam[0].ring_radius[i0]; 
        params.pa_temp[i0] = TRparam[0].pa_temp[i0]; 
        params.pa_temp_e[i0] = TRparam[0].pa_temp_e[i0]; 
        params.incl_temp[i0] = TRparam[0].incl_temp[i0]; 
        params.incl_temp_e[i0] = TRparam[0].incl_temp_e[i0]; 
    }

    // 1. sky position (i, j)
    params.i = i;
    params.j = j;
    // 2. XPOS YPOS
    params.xpos = _xpos;
    params.ypos = _ypos;

    // given pa + incl
    params.pa = _pa;
    params.incl = _incl;

    // 5. Rmax
    params.rmax = TRparam[0].rGalaxyPlane_pixel_max;

    Fgsl.function = &gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_given_paincl;
    Fgsl.params = &params;

    if((strcmp(solver, "brent") == 0))
        Tgsl = gsl_root_fsolver_brent;
    
    if((strcmp(solver, "bisection") == 0))
        Tgsl = gsl_root_fsolver_bisection;

    if((strcmp(solver, "falsepos") == 0))
        Tgsl = gsl_root_fsolver_falsepos;

    s = gsl_root_fsolver_alloc (Tgsl);
    gsl_root_fsolver_set (s, &Fgsl, x_lo, x_hi);

//    printf ("using %s method\n", gsl_root_fsolver_name (s));
//    printf ("%5s [%9s, %9s] %9s %10s %9s\n",
//            "iter", "lower", "upper", "root",
//            "err", "err(est)");

    do
      {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 1E-5);

        if (status == GSL_SUCCESS) ;
        //  printf ("Converged:\n");

//        printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//                iter, x_lo, x_hi,
//                r, r - r_expected,
//                x_hi - x_lo);
      }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
//    return status;
    return r;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation(double x, void *params)
{
    int m=0, n=0;
    int i0=0, j0=0;
    int bspline_order;
    struct rGalaxyPlane_pixel_TR_nonlinearEquation_params *p
      = (struct rGalaxyPlane_pixel_TR_nonlinearEquation_params *) params;

    int PA_fourier_order, pa_order, PA_bspline_order, incl_order, INCL_fourier_order, INCL_bspline_order, SersicPoly_order;
    double _x0, fx;
    double PA=0, INCL=0;

    int ncoeffs;
    int nbreak;
    double xi, yi, yerr, Bj;
    double chisq, Rsq, dof, tss;
    double _ring_temp, _value_temp, _e_temp;

    gsl_bspline_workspace *_bw_pa, *_bw_incl;
    gsl_vector *_B_pa, *_B_incl;
    gsl_vector *_c_pa, *_c_incl, *_w_pa, *_w_incl;
    gsl_vector *_x_pa, *_x_incl, *_y_pa, *_y_incl;
    gsl_matrix *_X_pa, *_X_incl, *_cov_pa, *_cov_incl;
    gsl_multifit_linear_workspace *_mw_pa, *_mw_incl;

    n = p->N_reliable_rings;
    // dynamic 1D array
    double *x_dat = malloc(sizeof(double) * n);
    double *y_dat = malloc(sizeof(double) * n);
    double *e_dat = malloc(sizeof(double) * n);

    // B-spline variables
    /* allocate a cubic bspline workspace (k = 4) */

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


    if(strcmp(p->pa_function, "fourier") == 0)
        PA_fourier_order = p->PA_fourier_order;
    else if(strcmp(p->pa_function, "poly") == 0)
        pa_order = p->pa_order;
    else if(strcmp(p->pa_function, "bspline") == 0)
        PA_bspline_order = p->n_coeffs_bspline_pa;

    if(strcmp(p->incl_function, "fourier") == 0)
        INCL_fourier_order = p->INCL_fourier_order;
    else if(strcmp(p->incl_function, "poly-sersic") == 0)
        SersicPoly_order = p->SersicPoly_order;
    else if(strcmp(p->incl_function, "poly") == 0)
        incl_order = p->incl_order;
    else if(strcmp(p->incl_function, "bspline") == 0)
        INCL_bspline_order = p->n_coeffs_bspline_incl;

    int i = p->i;
    double xpos = p->xpos;
    int j = p->j;
    double ypos = p->ypos;

    if(strcmp(p->pa_function, "fourier") == 0)
    {
        _x0 = (x/(p->rmax))*M_PI; // normalised one: r' = (r_pix / r_max) x PI
        PA = p->_p0;
        for(m=0; m<PA_fourier_order; m++)
        {
            PA += p->_p_a[m]*cos((m+1)*p->_p_w*_x0) + p->_p_b[m]*sin((m+1)*p->_p_w*_x0);
        }
        PA = (PA*p->PA_MAX_in_degree)*M_PI/180.; // in radian
    }
    else if(strcmp(p->pa_function, "poly") == 0)
    {
        _x0 = x/p->rmax; // normalised one: r' = (r_pix / r_max)
        PA = 0;
        for(m=0; m<pa_order+1; m++)
        {
            PA += p->_p[m]*pow(_x0, m);
        }
        PA = (PA*p->PA_MAX_in_degree)*M_PI/180.; // in radian
    }
    else if(strcmp(p->pa_function, "bspline") == 0)
    {
        for (i0=0; i0<n; i0++)
        {
            _ring_temp = p->ring_radius[i0]/p->rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = p->pa_temp[i0]/p->PA_MAX_in_degree; // extrapolate inward
            _e_temp = p->pa_temp_e[i0]/p->PA_MAX_in_degree; // extrapolate inward

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
            for (j0=0; j0<p->n_coeffs_bspline_pa; j0++)
            {
                Bj = gsl_vector_get(_B_pa, j0);
                gsl_matrix_set(_X_pa, i0, j0, Bj);
            }
        }

        /* construct the fit matrix _cov */
        for (i0=0; i0<p->n_coeffs_bspline_pa; i0++)
        {
            xi = gsl_vector_get(_x_pa, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_pa, _bw_pa);

            /* fill in row i of X */
            for (j0=0; j0<p->n_coeffs_bspline_pa; j0++)
            {
                // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                //Bj = gsl_vector_get(_B, _j);
                Bj = 0.01;
                gsl_matrix_set(_cov_pa, i0, j0, Bj);
            }
        }


        // load the spline efficients generated by MCMC
        for(i0=0; i0<p->n_coeffs_bspline_pa; i0++)
        {
            _c_pa->data[i0] = p->_p_bs[i0];
        }

        xi = x/p->rmax; // normalised one: r' = (r_pix / r_max)
        if(xi < x_dat[0]) xi = x_dat[0];
        if(xi > x_dat[n-1]) xi = x_dat[n-1];

        gsl_bspline_eval(xi, _B_pa, _bw_pa);
        gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
        PA = yi;
        PA = (PA*p->PA_MAX_in_degree)*M_PI/180.; // in radian
    }

    if(strcmp(p->incl_function, "bspline") == 0)
    {
        for (i0=0; i0<n; i0++)
        {
            _ring_temp = p->ring_radius[i0]/p->rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = p->incl_temp[i0]/p->INCL_MAX_in_degree; // extrapolate inward
            _e_temp = p->incl_temp_e[i0]/p->INCL_MAX_in_degree; // extrapolate inward

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
            for (j0=0; j0<p->n_coeffs_bspline_incl; j0++)
            {
                Bj = gsl_vector_get(_B_incl, j0);
                gsl_matrix_set(_X_incl, i0, j0, Bj);
            }
        }

        /* construct the fit matrix _cov */
        for (i0=0; i0<p->n_coeffs_bspline_incl; i0++)
        {
            xi = gsl_vector_get(_x_incl, i0);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, _B_incl, _bw_incl);

            /* fill in row i of X */
            for (j0=0; j0<p->n_coeffs_bspline_incl; j0++)
            {
                // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                //Bj = gsl_vector_get(_B, _j);
                Bj = 0.01;
                gsl_matrix_set(_cov_incl, i0, j0, Bj);
            }
        }


        // load the spline efficients generated by MCMC
        for(i0=0; i0<p->n_coeffs_bspline_incl; i0++)
        {
            _c_incl->data[i0] = p->_i_bs[i0];
        }

        xi = x/p->rmax; // normalised one: r' = (r_pix / r_max)
        if(xi < x_dat[0]) xi = x_dat[0];
        if(xi > x_dat[n-1]) xi = x_dat[n-1];

        gsl_bspline_eval(xi, _B_incl, _bw_incl);
        gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
        INCL = yi;
        INCL = (INCL*p->INCL_MAX_in_degree)*M_PI/180.; // in radian
    }

    fx = sqrt(pow((-((double)i-xpos)*sin(PA) + ((double)j-ypos)*cos(PA)), 2) + pow(((((double)i-xpos)*cos(PA) + ((double)j-ypos)*sin(PA))/cos(INCL)), 2)) - x;

    free(x_dat);
    free(y_dat);
    free(e_dat);

    gsl_bspline_free(_bw_pa);
    gsl_vector_free(_B_pa);
    gsl_vector_free(_x_pa);
    gsl_vector_free(_y_pa);
    gsl_matrix_free(_X_pa);
    gsl_vector_free(_c_pa);
    gsl_vector_free(_w_pa);
    gsl_matrix_free(_cov_pa);
    gsl_multifit_linear_free(_mw_pa);

    gsl_bspline_free(_bw_incl);
    gsl_vector_free(_B_incl);
    gsl_vector_free(_x_incl);
    gsl_vector_free(_y_incl);
    gsl_matrix_free(_X_incl);
    gsl_vector_free(_c_incl);
    gsl_vector_free(_w_incl);
    gsl_matrix_free(_cov_incl);
    gsl_multifit_linear_free(_mw_incl);

    //free(p);

    return fx;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_given_paincl(double x, void *params)
{
    int m=0, n=0;
    int i0=0, j0=0;
    int bspline_order;
    struct rGalaxyPlane_pixel_TR_nonlinearEquation_params *p
      = (struct rGalaxyPlane_pixel_TR_nonlinearEquation_params *) params;

//    struct rGalaxyPlane_pixel_TR_nonlinearEquation_params *p = malloc(sizeof(struct rGalaxyPlane_pixel_TR_nonlinearEquation_params));
//   p = params;

    double _x0, fx;
    double PA=0, INCL=0;

    int i = p->i;
    double xpos = p->xpos;
    int j = p->j;
    double ypos = p->ypos;

    PA = (p->pa)*M_PI/180.; // in radian
    INCL = (p->incl)*M_PI/180.; // in radian

    fx = sqrt(pow((-((double)i-xpos)*sin(PA) + ((double)j-ypos)*cos(PA)), 2) + pow(((((double)i-xpos)*cos(PA) + ((double)j-ypos)*sin(PA))/cos(INCL)), 2)) - x;

    return fx;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int gsl_nonlinearfit_Bspline_filtering(char *param, TR_ringParameters *TRparam)
{
    int i, j, n, status;
    int bspline_order;
    int N_reliable_rings=0;
    double ri, ro, ring_s, ring_w, rGalaxyPlane_pixel_max;
    float rbm, std, rbm_e, std_e;
    double xi, yi, yerr;
    double y_inner, y_outer;
    double _ring_temp, _value_temp, _e_temp;
    double scale_pa_e_to_total_error;
    double scale_incl_e_to_total_error;
    double scale_vrot_e_to_total_error;
    double gsl_mean_pa, gsl_std_pa, gsl_max_pa, gsl_min_pa;
    double gsl_mean_incl, gsl_std_incl, gsl_max_incl, gsl_min_incl;
    double gsl_mean_vrot, gsl_std_vrot, gsl_max_vrot, gsl_min_vrot;

    int ncoeffs, nbreak;
    gsl_bspline_workspace *bw;
    gsl_vector *B;
    double chisq, Rsq, dof, tss;
    gsl_vector *c, *w;
    gsl_vector *x, *y;
    gsl_matrix *X, *cov;
    gsl_multifit_linear_workspace *mw;

    const gsl_rng_type * type;
    gsl_rng * r;
    gsl_rng_env_setup();
    type = gsl_rng_default;
    r = gsl_rng_alloc (type);

    // gsl vector for ring params : this is for calculating max/min of errors. See below.
    gsl_vector *gsl_xpos_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_ypos_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vsys_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_pa_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_incl_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vrot_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_total_error = gsl_vector_alloc(TRparam[0].Nrings);

    // dynamic total_error 1D array
    double *total_error_ring_params = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array

    // dynamic pa_temp 3D array
    double **pa_temp = (double**)malloc(sizeof(double *) * TRparam[0].Nrings); // 3D array: 999 x 3
    for(i=0; i<TRparam[0].Nrings; i++)
    {
        pa_temp[i] = (double*)malloc(sizeof(double) * 3);
        pa_temp[i][0] = 1E3;
        pa_temp[i][1] = 1E3;
        pa_temp[i][2] = 1E3;
    }
    // dynamic incl_temp 3D array
    double **incl_temp = (double**)malloc(sizeof(double *) * TRparam[0].Nrings); // 3D array: 999 x 3
    for(i=0; i<TRparam[0].Nrings; i++)
    {
        incl_temp[i] = (double*)malloc(sizeof(double) * 3);
        incl_temp[i][0] = 1E3;
        incl_temp[i][1] = 1E3;
        incl_temp[i][2] = 1E3;
    }
    // dynamic vrot_temp 3D array
    double **vrot_temp = (double**)malloc(sizeof(double *) * TRparam[0].Nrings); // 3D array:  x 3
    for(i=0; i<TRparam[0].Nrings; i++)
    {
        vrot_temp[i] = (double*)malloc(sizeof(double) * 3);
        vrot_temp[i][0] = 1E3;
        vrot_temp[i][1] = 1E3;
        vrot_temp[i][2] = 1E3;
    }

    double *pa_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double *incl_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrot_err = malloc(sizeof(double) * TRparam[0].Nrings);

    ring_s = TRparam[0].ring_w;
    ring_w = TRparam[0].ring_w;
    rGalaxyPlane_pixel_max = TRparam[0].rGalaxyPlane_pixel_max;
    N_reliable_rings = 0;

    if(strcmp(param, "PA") == 0)
    {
        trfit_multinest_ellipsefit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "vrot", 'T', "vrad", 'F', "sigmafactor", 'T', multinest_param, TRparam, 0, 0, "finalfit", 'N'); // both sides

        for(j=0; j<TRparam[0].Nrings; j++)
        {
            ri = TRparam[0].ring_s + j*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
            ro = TRparam[0].ring_s + j*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
            if (ri < 0.0) ri = 0.0;

            if(TRparam[0].pa[j] > 0 &&  TRparam[0].pa[j] < TRparam[0].PA_MAX_in_degree)
            {
                pa_temp[j][0] = (ri+ro)/2.0;
                pa_temp[j][1] = TRparam[0].pa[j];
                pa_temp[j][2] = TRparam[0].pa_e[j];
                pa_err[j] = TRparam[0].pa_e[j];
            }
            else
            {
                pa_temp[j][0] = (ri+ro)/2.0;
                pa_temp[j][1] = TRparam[0].pa[j];
                pa_temp[j][2] = TRparam[0].PA_MAX_in_degree; // large error
                pa_err[j] = TRparam[0].PA_MAX_in_degree;
            }

            if(TRparam[0].incl[j] > 0 &&  TRparam[0].incl[j] < TRparam[0].INCL_MAX_in_degree)
            {
                incl_temp[j][0] = (ri+ro)/2.0;
                incl_temp[j][1] = TRparam[0].incl[j];
                incl_temp[j][2] = TRparam[0].incl_e[j];
                incl_err[j] = TRparam[0].incl_e[j];
            }
            else
            {
                incl_temp[j][0] = (ri+ro)/2.0;
                incl_temp[j][1] = TRparam[0].incl[j];
                incl_temp[j][2] = TRparam[0].INCL_MAX_in_degree; // large error
                incl_err[j] = TRparam[0].INCL_MAX_in_degree;
            }

            if(TRparam[0].vrot[j] > 0) 
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = TRparam[0].vrot_e[j];
                vrot_err[j] = TRparam[0].vrot_e[j];
            }
            else
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = 1E3;
                vrot_err[j] = 1E3;
            }

            // copy ring params to gsl vectors for finding min/max of errors 
            gsl_vector_set(gsl_xpos_e, j, TRparam[0].xpos_e[j]);
            gsl_vector_set(gsl_ypos_e, j, TRparam[0].ypos_e[j]);
            gsl_vector_set(gsl_vsys_e, j, TRparam[0].vsys_e[j]);
            gsl_vector_set(gsl_pa_e, j, pa_temp[j][2]);
            gsl_vector_set(gsl_incl_e, j, TRparam[0].incl_e[j]);
            gsl_vector_set(gsl_vrot_e, j, TRparam[0].vrot_e[j]);

        }
        gsl_mean_pa = gsl_stats_mean(pa_err, 1, TRparam[0].Nrings);
        gsl_std_pa = sqrt(gsl_stats_variance(pa_err, 1, TRparam[0].Nrings));
        gsl_min_pa = gsl_stats_min(pa_err, 1, TRparam[0].Nrings);
        gsl_max_pa = gsl_stats_max(pa_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_pa) || isinf(gsl_min_pa) || gsl_min_pa <= 0)
            gsl_min_pa = 1E3; // large value

        gsl_mean_incl = gsl_stats_mean(incl_err, 1, TRparam[0].Nrings);
        gsl_std_incl = sqrt(gsl_stats_variance(incl_err, 1, TRparam[0].Nrings));
        gsl_min_incl = gsl_stats_min(incl_err, 1, TRparam[0].Nrings);
        gsl_max_incl = gsl_stats_max(incl_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_incl) || isinf(gsl_min_incl) || gsl_min_incl <= 0)
            gsl_min_incl = 1E3; // large value

        gsl_mean_vrot = gsl_stats_mean(vrot_err, 1, TRparam[0].Nrings);
        gsl_std_vrot = sqrt(gsl_stats_variance(vrot_err, 1, TRparam[0].Nrings));
        gsl_min_vrot = gsl_stats_min(vrot_err, 1, TRparam[0].Nrings);
        gsl_max_vrot = gsl_stats_max(vrot_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrot) || isinf(gsl_min_vrot) || gsl_min_vrot <= 0)
            gsl_min_vrot = 1E3; // small value

        for(j=0; j<TRparam[0].Nrings; j++)
        {
            if(TRparam[0].pa_e[j] > 0)
            {
                pa_err[j] = TRparam[0].pa_e[j]/gsl_min_pa;
            }
            else
            {
                pa_err[j] = 1E3; // large value
            }

            if(TRparam[0].incl_e[j] > 0)
            {
                incl_err[j] = TRparam[0].incl_e[j]/gsl_min_incl;
            }
            else
            {
                incl_err[j] = 1E3; // large value
            }

            if(TRparam[0].vrot_e[j] > 0) 
            {
                vrot_err[j] = TRparam[0].vrot_e[j]/gsl_min_vrot;
            }
            else
            {
                vrot_err[j] = 1E3;
            }

            // total errors (PA_e+INCL_e+VROT_e): normalised ones to their minimums
            total_error_ring_params[j] = sqrt(pow(pa_err[j], 2) + \
                                                pow(incl_err[j], 2) + \
                                                pow(vrot_err[j], 2));

            if(isinf(total_error_ring_params[j]) || isnan(total_error_ring_params[j]) || total_error_ring_params[j] <= 0)
                total_error_ring_params[i] = 1E3; // large value
        }
    }
    else if(strcmp(param, "INCL") == 0)
    {
        trfit_multinest_ellipsefit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "vrot", 'T', "vrad", 'T', "sigmafactor", 'T', multinest_param, TRparam, 0, 0, "finalfit", 'N'); // both sides

        for(j=0; j<TRparam[0].Nrings; j++)
        {
            ri = TRparam[0].ring_s + j*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
            ro = TRparam[0].ring_s + j*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
            if (ri < 0.0) ri = 0.0;

            if(TRparam[0].pa[j] > 0 &&  TRparam[0].pa[j] < TRparam[0].PA_MAX_in_degree)
            {
                pa_temp[j][0] = (ri+ro)/2.0;
                pa_temp[j][1] = TRparam[0].pa[j];
                pa_temp[j][2] = TRparam[0].pa_e[j];
                pa_err[j] = TRparam[0].pa_e[j];
            }
            else
            {
                pa_temp[j][0] = (ri+ro)/2.0;
                pa_temp[j][1] = TRparam[0].pa[j];
                pa_temp[j][2] = TRparam[0].PA_MAX_in_degree; // large error
                pa_err[j] = TRparam[0].PA_MAX_in_degree;
            }

            if(TRparam[0].incl[j] > 0 &&  TRparam[0].incl[j] < TRparam[0].INCL_MAX_in_degree)
            {
                incl_temp[j][0] = (ri+ro)/2.0;
                incl_temp[j][1] = TRparam[0].incl[j];
                incl_temp[j][2] = TRparam[0].incl_e[j];
                incl_err[j] = TRparam[0].incl_e[j];
            }
            else
            {
                incl_temp[j][0] = (ri+ro)/2.0;
                incl_temp[j][1] = TRparam[0].incl[j];
                incl_temp[j][2] = TRparam[0].INCL_MAX_in_degree; // large error
                incl_err[j] = TRparam[0].INCL_MAX_in_degree;
            }

            if(TRparam[0].vrot[j] > 0) 
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = TRparam[0].vrot_e[j];
                vrot_err[j] = TRparam[0].vrot_e[j];
            }
            else
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = 1E3;
                vrot_err[j] = 1E3;
            }

            // copy ring params to gsl vectors for finding min/max of errors 
            gsl_vector_set(gsl_xpos_e, j, TRparam[0].xpos_e[j]);
            gsl_vector_set(gsl_ypos_e, j, TRparam[0].ypos_e[j]);
            gsl_vector_set(gsl_vsys_e, j, TRparam[0].vsys_e[j]);
            gsl_vector_set(gsl_pa_e, j, pa_temp[j][2]);
            gsl_vector_set(gsl_incl_e, j, TRparam[0].incl_e[j]);
            gsl_vector_set(gsl_vrot_e, j, TRparam[0].vrot_e[j]);

        }
        gsl_mean_pa = gsl_stats_mean(pa_err, 1, TRparam[0].Nrings);
        gsl_std_pa = sqrt(gsl_stats_variance(pa_err, 1, TRparam[0].Nrings));
        gsl_min_pa = gsl_stats_min(pa_err, 1, TRparam[0].Nrings);
        gsl_max_pa = gsl_stats_max(pa_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_pa) || isinf(gsl_min_pa) || gsl_min_pa <= 0)
            gsl_min_pa = 1E3; // large value

        gsl_mean_incl = gsl_stats_mean(incl_err, 1, TRparam[0].Nrings);
        gsl_std_incl = sqrt(gsl_stats_variance(incl_err, 1, TRparam[0].Nrings));
        gsl_min_incl = gsl_stats_min(incl_err, 1, TRparam[0].Nrings);
        gsl_max_incl = gsl_stats_max(incl_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_incl) || isinf(gsl_min_incl) || gsl_min_incl <= 0)
            gsl_min_incl = 1E3; // large value

        gsl_mean_vrot = gsl_stats_mean(vrot_err, 1, TRparam[0].Nrings);
        gsl_std_vrot = sqrt(gsl_stats_variance(vrot_err, 1, TRparam[0].Nrings));
        gsl_min_vrot = gsl_stats_min(vrot_err, 1, TRparam[0].Nrings);
        gsl_max_vrot = gsl_stats_max(vrot_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrot) || isinf(gsl_min_vrot) || gsl_min_vrot <= 0)
            gsl_min_vrot = 1E3; // small value

        for(j=0; j<TRparam[0].Nrings; j++)
        {
            if(TRparam[0].pa_e[j] > 0)
            {
                pa_err[j] = TRparam[0].pa_e[j]/gsl_min_pa;
            }
            else
            {
                pa_err[j] = 1E3; // large value
            }

            if(TRparam[0].incl_e[j] > 0)
            {
                incl_err[j] = TRparam[0].incl_e[j]/gsl_min_incl;
            }
            else
            {
                incl_err[j] = 1E3; // large value
            }

            if(TRparam[0].vrot_e[j] > 0) 
            {
                vrot_err[j] = TRparam[0].vrot_e[j]/gsl_min_vrot;
            }
            else
            {
                vrot_err[j] = 1E3;
            }

            // total errors (PA_e+INCL_e+VROT_e): normalised ones to their minimums
            total_error_ring_params[j] = sqrt(pow(pa_err[j], 2) + \
                                                pow(incl_err[j], 2) + \
                                                pow(vrot_err[j], 2));

            if(isinf(total_error_ring_params[j]) || isnan(total_error_ring_params[j]) || total_error_ring_params[j] <= 0)
                total_error_ring_params[i] = 1E3; // large value
        }
    }
    else if(strcmp(param, "VROT-ellipsefit_rings") == 0)
    {
        trfit_multinest_ellipsefit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", 'F', "sigmafactor", 'T', multinest_param, TRparam, 0, 0, "finalfit", 'N'); // both sides
        for(j=0; j<TRparam[0].Nrings; j++)
        {
            ri = TRparam[0].ring_s + j*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
            ro = TRparam[0].ring_s + j*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
            if (ri < 0.0) ri = 0.0;

            if(TRparam[0].vrot[j] > 0) 
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = TRparam[0].vrot_e[j];
                vrot_err[j] = TRparam[0].vrot_e[j];
            }
            else
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = 1E3;
                vrot_err[j] = 1E3;
            }

            // copy ring params to gsl vectors for finding min/max of errors 
            gsl_vector_set(gsl_vrot_e, j, TRparam[0].vrot_e[j]);
        }

        gsl_mean_vrot = gsl_stats_mean(vrot_err, 1, TRparam[0].Nrings);
        gsl_std_vrot = sqrt(gsl_stats_variance(vrot_err, 1, TRparam[0].Nrings));
        gsl_min_vrot = gsl_stats_min(vrot_err, 1, TRparam[0].Nrings);
        gsl_max_vrot = gsl_stats_max(vrot_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrot) || isinf(gsl_min_vrot) || gsl_min_vrot <= 0)
            gsl_min_vrot = 1E3; // small value

        for(j=0; j<TRparam[0].Nrings; j++)
        {
            if(TRparam[0].vrot_e[j] > 0) 
            {
                vrot_err[j] = TRparam[0].vrot_e[j]/gsl_min_vrot;
            }
            else
            {
                vrot_err[j] = 1E3;
            }

            // total errors (PA_e+INCL_e+VROT_e): normalised ones to their minimums
            total_error_ring_params[j] = vrot_err[j];

            if(isinf(total_error_ring_params[j]) || isnan(total_error_ring_params[j]) || total_error_ring_params[j] <= 0)
                total_error_ring_params[i] = 1E3; // large value
        }
    }
    else if(strcmp(param, "VROT-einasto_halofit_rings") == 0)
    {
        trfit_multinest_trfit_rings_normal("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", 'T', "sigmafactor", 'T', multinest_param, TRparam, 999, "finalfit", 'N');
        //trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", 'T', "sigmafactor", 'T', multinest_param, TRparam, 999, "finalfit", 'N');

        for(j=0; j<TRparam[0].Nrings; j++)
        {
            ri = TRparam[0].ring_s + j*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
            ro = TRparam[0].ring_s + j*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
            if (ri < 0.0) ri = 0.0;

            if(TRparam[0].vrot[j] > 0) 
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = TRparam[0].vrot_e[j];
                vrot_err[j] = TRparam[0].vrot_e[j];
            }
            else
            {
                vrot_temp[j][0] = (ri+ro)/2.0;
                vrot_temp[j][1] = TRparam[0].vrot[j];
                vrot_temp[j][2] = 1E3;
                vrot_err[j] = 1E3;
            }

            // copy ring params to gsl vectors for finding min/max of errors 
            gsl_vector_set(gsl_vrot_e, j, TRparam[0].vrot_e[j]);
        }

        gsl_mean_vrot = gsl_stats_mean(vrot_err, 1, TRparam[0].Nrings);
        gsl_std_vrot = sqrt(gsl_stats_variance(vrot_err, 1, TRparam[0].Nrings));
        gsl_min_vrot = gsl_stats_min(vrot_err, 1, TRparam[0].Nrings);
        gsl_max_vrot = gsl_stats_max(vrot_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrot) || isinf(gsl_min_vrot) || gsl_min_vrot <= 0)
            gsl_min_vrot = 1E3; // small value

        for(j=0; j<TRparam[0].Nrings; j++)
        {
            if(TRparam[0].vrot_e[j] > 0) 
            {
                vrot_err[j] = TRparam[0].vrot_e[j]/gsl_min_vrot;
            }
            else
            {
                vrot_err[j] = 1E3;
            }

            // total errors (PA_e+INCL_e+VROT_e): normalised ones to their minimums
            total_error_ring_params[j] = vrot_err[j];

            if(isinf(total_error_ring_params[j]) || isnan(total_error_ring_params[j]) || total_error_ring_params[j] <= 0)
                total_error_ring_params[i] = 1E3; // large value
        }
    }
    else
    {
        printf("Invailid ring parameters, %s: Check gsl_nonlinearfit_polynomial()\n", param);
        exit(1);
    }

    /*
    qsort(pa_temp, 999, sizeof(pa_temp[0]), array1D_comp);
    qsort(incl_temp, 999, sizeof(incl_temp[0]), array1D_comp);
    qsort(vrot_temp, 999, sizeof(vrot_temp[0]), array1D_comp);
    */
    TRparam[0].N_reliable_rings = TRparam[0].Nrings;
    n = TRparam[0].Nrings;

    for(i=0; i<n; i++)
    {
        if(strcmp(param, "PA") == 0)
        {
            TRparam[0].ring_radius[i] = pa_temp[i][0];
            TRparam[0].pa_temp[i] = pa_temp[i][1];
            TRparam[0].pa_temp_e[i] = pa_temp[i][2];
        }
        else if(strcmp(param, "INCL") == 0)
        {
            TRparam[0].ring_radius[i] = incl_temp[i][0];
            TRparam[0].incl_temp[i] = incl_temp[i][1];
            TRparam[0].incl_temp_e[i] = incl_temp[i][2];
        }
        else if(strcmp(param, "VROT-ellipsefit_rings") == 0 || strcmp(param, "VROT-einasto_halofit_rings") == 0)
        {
            TRparam[0].ring_radius[i] = vrot_temp[i][0];
            TRparam[0].vrot_temp[i] = vrot_temp[i][1];
            TRparam[0].vrot_temp_e[i] = vrot_temp[i][2];
        }
    }

    // dynamic 1D array
    double *x_dat = malloc(sizeof(double) * n);
    double *y_dat = malloc(sizeof(double) * n);
    double *e_dat = malloc(sizeof(double) * n);
    double *sigma = malloc(sizeof(double) * n);
    float *residual = malloc(sizeof(float) * n);
    float *error = malloc(sizeof(float) * n);
    double *Y_spline = malloc(sizeof(double) * n);

    // B-spline variables
    /* allocate a cubic bspline workspace (k = 4) */
    if(strcmp(param, "PA") == 0)
    {
        bspline_order = TRparam[0].pa_order_bspline;
        nbreak = TRparam[0].pa_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_pa = ncoeffs;
    }
    else if(strcmp(param, "INCL") == 0)
    {
        bspline_order = TRparam[0].incl_order_bspline;
        nbreak = TRparam[0].incl_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_incl = ncoeffs;
    }
    else if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
    {
        bspline_order = 3;
        nbreak = 2;
        ncoeffs = 4;
    }

    bw = gsl_bspline_alloc(bspline_order+1, nbreak);
    B = gsl_vector_alloc(ncoeffs);

    x = gsl_vector_alloc(n);
    y = gsl_vector_alloc(n);
    X = gsl_matrix_alloc(n, ncoeffs);
    c = gsl_vector_alloc(ncoeffs);
    w = gsl_vector_alloc(n);
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    mw = gsl_multifit_linear_alloc(n, ncoeffs);

    for (i=0; i<n; i++)
    {
        if(strcmp(param, "PA") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].pa[i]/TRparam[0].PA_MAX_in_degree; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = total_error_ring_params[i]; // as derived from above
        }
        else if(strcmp(param, "INCL") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].incl[i]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward
            _e_temp = TRparam[0].incl_e[i]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = total_error_ring_params[i]; // as derived from above
        }
        else if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].vrot[i]; // normalise with mean as the maximum can be very large!

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = total_error_ring_params[i]; // as derived from above
        }
        else
        {
            printf("Approximation for %s is not available. Only PA or INCL is available\n", param);
            exit(1);
        }
    }

    for (i=0; i<n; i++)
    {
          gsl_vector_set (x, i, x_dat[i]);
          gsl_vector_set (y, i, y_dat[i]);
          gsl_vector_set (w, i, 1/e_dat[i]);
    }

    /* use uniform breakpoints on [0, 15] */
    gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], bw);

    /* construct the fit matrix X */
    for (i=0; i<n; i++)
    {
      double xi = gsl_vector_get(x, i);
      if(xi < x_dat[0]) xi = x_dat[0];
      if(xi > x_dat[n-1]) xi = x_dat[n-1];
      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);

      /* fill in row i of X */
      for (j=0; j<ncoeffs; j++)
      {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
      }
    }

    /* do the fit */
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

    dof = n - ncoeffs;
    tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
    Rsq = 1.0 - chisq / tss;

    fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", chisq / dof, Rsq);

/*
    // 1. extrapolation to the inner from the first spline derived 
    xi = x_dat[(int)(n*0.1)];
    gsl_bspline_eval(xi, B, bw);
    gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
    y_inner = yi;
    // 2. extrapolation to the outer from the last spline derived 
    xi = x_dat[n-1];
    gsl_bspline_eval(xi, B, bw);
    gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
    y_outer = yi;
*/
  
    /* output the smoothed curve */
    {
        if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
        {
            for (i=0; i<n; i++)
            {
                if(i < (int)n*0.0) // it was 0.1 for the old version.
                {
                    Y_spline[i] = 0.5;
                }
                else
                {
                    xi = x_dat[i];
                    gsl_bspline_eval(xi, B, bw);
                    gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
                    Y_spline[i] = yi;
                }
            }
        }
        else if(strcmp(param, "PA") == 0 || strcmp(param, "INCL") == 0)
        {
            for (i=0; i<n; i++)
            {
                xi = x_dat[i];
                gsl_bspline_eval(xi, B, bw);
                gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
                Y_spline[i] = yi;
            //    printf("inclnew: %f %f %f\n", x_dat[i], y_dat[i], Y_spline[i]);
            }
        }

    }

    for(i=0; i<n; i++)
    {
        residual[i] = 999;
        error[i] = 999;
    }

    for (i=0; i<n; i++)
    {
        if(strcmp(param, "PA") == 0)
        {
            sigma[i] = e_dat[i]; // useless
            residual[i] = y_dat[i] - Y_spline[i];
            error[i] = TRparam[0].pa_e[i];
        }
        else if(strcmp(param, "INCL") == 0)
        {
            sigma[i] = e_dat[i]; // useless
            residual[i] = y_dat[i] - Y_spline[i];
            error[i] = TRparam[0].incl_e[i];
        }
        else if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
        {
            sigma[i] = e_dat[i]; // useless
            residual[i] = y_dat[i] - Y_spline[i];
            error[i] = TRparam[0].vrot_e[i];
        }
        else
        {
            printf("Approximation for %s is not available. Only PA or INCL is available\n", param);
            exit(1);
        }
    }
    get_rbm(residual, n, -1E3, 1E3, &rbm, &std);
    get_rbm(error, n, -1E3, 1E3, &rbm_e, &std_e);


    TRparam[0].N_reliable_rings = 0;
    if(strcmp(param, "PA") == 0)
    {
        for(i=0; i<n; i++) // for the intrinsic TR rings ignoring the edge points
        {
            if(residual[i] < rbm + 5*std && error[i] < rbm_e+5*std_e)
            {
                TRparam[0].ring_radius[TRparam[0].N_reliable_rings] = x_dat[i]*TRparam[0].rGalaxyPlane_pixel_max; // restore the galaxy radius in pix
                TRparam[0].pa_temp[TRparam[0].N_reliable_rings] = y_dat[i]*TRparam[0].PA_MAX_in_degree;
                if(!isinf(sigma[i]) && !isnan(sigma[i]))
                    TRparam[0].pa_temp_e[TRparam[0].N_reliable_rings] = sigma[i]*TRparam[0].PA_MAX_in_degree;
                else
                    TRparam[0].pa_temp_e[TRparam[0].N_reliable_rings] = 9999; // useless

                TRparam[0].N_reliable_rings++;
            }
            else
            {
                TRparam[0].ring_radius[TRparam[0].N_reliable_rings] = x_dat[i]*TRparam[0].rGalaxyPlane_pixel_max;
               // TRparam[0].pa_temp[TRparam[0].N_reliable_rings] = Y_spline[i]*TRparam[0].PA_MAX_in_degree;
                TRparam[0].pa_temp[TRparam[0].N_reliable_rings] = (Y_spline[i] + gsl_ran_gaussian(r, std*0.5))*TRparam[0].PA_MAX_in_degree;

                if(!isinf(sigma[i]) && !isnan(sigma[i]))
                    TRparam[0].pa_temp_e[TRparam[0].N_reliable_rings] = sigma[i]*TRparam[0].PA_MAX_in_degree;
                else
                    TRparam[0].pa_temp_e[TRparam[0].N_reliable_rings] = 9999; // useless

                TRparam[0].N_reliable_rings++;
            }
        }
    }
    else if(strcmp(param, "INCL") == 0)
    {
        for(i=0; i<n; i++) // for the intrinsic TR rings ignoring the edge points
        {
            if(residual[i] < rbm + 5*std && error[i] < rbm_e+3*std_e)
            {
                TRparam[0].ring_radius[TRparam[0].N_reliable_rings] = x_dat[i]*TRparam[0].rGalaxyPlane_pixel_max;
                TRparam[0].incl_temp[TRparam[0].N_reliable_rings] = y_dat[i]*TRparam[0].INCL_MAX_in_degree;
                TRparam[0].incl_temp_e[TRparam[0].N_reliable_rings] = sigma[i]*TRparam[0].INCL_MAX_in_degree;

                if(!isinf(sigma[i]) && !isnan(sigma[i]))
                    TRparam[0].incl_temp_e[TRparam[0].N_reliable_rings] = sigma[i]*TRparam[0].INCL_MAX_in_degree;
                else
                    TRparam[0].incl_temp_e[TRparam[0].N_reliable_rings] = 9999; // useless

                TRparam[0].N_reliable_rings++;
            }
            else
            {
                TRparam[0].ring_radius[TRparam[0].N_reliable_rings] = x_dat[i]*TRparam[0].rGalaxyPlane_pixel_max;
                TRparam[0].incl_temp[TRparam[0].N_reliable_rings] = (Y_spline[i] + gsl_ran_gaussian(r, std*0.5))*TRparam[0].INCL_MAX_in_degree;

                if(!isinf(sigma[i]) && !isnan(sigma[i]))
                    TRparam[0].incl_temp_e[TRparam[0].N_reliable_rings] = sigma[i]*TRparam[0].INCL_MAX_in_degree;
                else
                    TRparam[0].incl_temp_e[TRparam[0].N_reliable_rings] = 9999; // useless

                TRparam[0].N_reliable_rings++;
            }
        }
    }
    else if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
    {
        for(i=0; i<n; i++) // for the intrinsic TR rings ignoring the edge points
        {
            if(residual[i] < rbm + 5*std && error[i] < rbm_e+3*std_e)
            {
                TRparam[0].ring_radius[TRparam[0].N_reliable_rings] = x_dat[i]*TRparam[0].rGalaxyPlane_pixel_max;
                TRparam[0].vrot_temp[TRparam[0].N_reliable_rings] = y_dat[i];
                TRparam[0].vrot_temp_e[TRparam[0].N_reliable_rings] = sigma[i];

                if(!isinf(sigma[i]) && !isnan(sigma[i]))
                    TRparam[0].vrot_temp_e[TRparam[0].N_reliable_rings] = sigma[i];
                else
                    TRparam[0].vrot_temp_e[TRparam[0].N_reliable_rings] = 9999; // useless

                TRparam[0].N_reliable_rings++;
            }
            else
            {
                TRparam[0].ring_radius[TRparam[0].N_reliable_rings] = x_dat[i]*TRparam[0].rGalaxyPlane_pixel_max;
                TRparam[0].vrot_temp[TRparam[0].N_reliable_rings] = 1E9; // garbage

                if(!isinf(sigma[i]) && !isnan(sigma[i]))
                    TRparam[0].vrot_temp_e[TRparam[0].N_reliable_rings] = sigma[i];
                else
                    TRparam[0].vrot_temp_e[TRparam[0].N_reliable_rings] = 9999; // useless

                TRparam[0].N_reliable_rings++;
            }
        }
    }

    free(pa_temp);
    free(incl_temp);
    free(vrot_temp);
    free(total_error_ring_params);
    free(pa_err);
    free(incl_err);
    free(vrot_err);

    free(x_dat);
    free(y_dat);
    free(e_dat);
    free(sigma);
    free(residual);
    free(error);
    free(Y_spline);

    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(w);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);

    gsl_vector_free(gsl_xpos_e);
    gsl_vector_free(gsl_ypos_e);
    gsl_vector_free(gsl_vsys_e);
    gsl_vector_free(gsl_pa_e);
    gsl_vector_free(gsl_incl_e);
    gsl_vector_free(gsl_vrot_e);

    return 0;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int bsplinefit_set_unipriors(char *param, TR_ringParameters *TRparam)
{
    int i, j, n;
    int bspline_order;
    double _ring_temp, _value_temp, _e_temp, rGalaxyPlane_pixel_max;
    double xi, yi, yerr;
    double y_inner, y_outer;

    int ncoeffs, nbreak, n_loop;
    gsl_bspline_workspace *bw;
    gsl_vector *B;
    double chisq, Rsq, dof, tss;
    gsl_vector *c, *c_temp, *w;
    gsl_vector *x, *y;
    gsl_matrix *X, *cov;
    gsl_multifit_linear_workspace *mw;

    rGalaxyPlane_pixel_max = TRparam[0].rGalaxyPlane_pixel_max;

    //n = TRparam[0].N_reliable_rings;
    n = TRparam[0].Nrings;

    // dynamic 1D array
    double *x_dat = malloc(sizeof(double) * n);
    double *y_dat = malloc(sizeof(double) * n);
    double *y_dat_t = malloc(sizeof(double) * n);
    double *e_dat = malloc(sizeof(double) * n);
    double *e_dat_t = malloc(sizeof(double) * n);
    double *sigma = malloc(sizeof(double) * n);
    float *error = malloc(sizeof(float) * n);
    double *Y_spline = malloc(sizeof(double) * n);
    double *Y_spline_err = malloc(sizeof(double) * n);
    double *Y_spline_err_t = malloc(sizeof(double) * n);
    double hist_mean_Y_spline_err, hist_std_Y_spline_err;

    double *xpos = malloc(sizeof(double) * TRparam[0].Nrings);
    double *xpos_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_xpos_err, hist_std_xpos_err;
    double mean_xpos, std_xpos, mean_xpos_e, std_xpos_e;

    double *ypos = malloc(sizeof(double) * TRparam[0].Nrings);
    double *ypos_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_ypos_err, hist_std_ypos_err;
    double mean_ypos, std_ypos, mean_ypos_e, std_ypos_e;

    double *vsys = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vsys_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vsys_err, hist_std_vsys_err;
    double mean_vsys, std_vsys, mean_vsys_e, std_vsys_e;

    double *pa = malloc(sizeof(double) * TRparam[0].Nrings);
    double *pa_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_pa_err, hist_std_pa_err;
    double mean_pa, std_pa, mean_pa_e, std_pa_e;

    double *incl = malloc(sizeof(double) * TRparam[0].Nrings);
    double *incl_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_incl_err, hist_std_incl_err;
    double mean_incl, std_incl, mean_incl_e, std_incl_e;

    double *vrot = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrot_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vrot_err, hist_std_vrot_err;
    double hist_mean_e_dat, hist_std_e_dat;

    double *vrad = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrad_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vrad_err, hist_std_vrad_err;

    double *total_error_norm = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array
    double mean_vrot, std_vrot, mean_vrot_e, std_vrot_e;
    double mean_vrad, std_vrad, mean_vrad_e, std_vrad_e;

    double robust_mean_y, robust_std_y;

    int filter_flag = 0, outliered=0;

    // gsl statistics
    gsl_vector *gsl_xpos = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_xpos_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_ypos = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_ypos_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vsys = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vsys_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_pa = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_pa_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_incl = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_incl_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vrot = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vrot_e = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vrad = gsl_vector_alloc(TRparam[0].Nrings);
    gsl_vector *gsl_vrad_e = gsl_vector_alloc(TRparam[0].Nrings);
    

    double gsl_max_vrot, gsl_min_vrot;
    double gsl_max_vrot_e, gsl_min_vrot_e;
    double gsl_max_vrad, gsl_min_vrad;
    double gsl_max_vrad_e, gsl_min_vrad_e;

    double xpos_L, xpos_U, xpos_eL, xpos_eU;
    double ypos_L, ypos_U, ypos_eL, ypos_eU;
    double vsys_L, vsys_U, vsys_eL, vsys_eU;
    double pa_L, pa_U, pa_eL, pa_eU;
    double incl_L, incl_U, incl_eL, incl_eU;
    double vrot_L, vrot_U, vrot_eL, vrot_eU;
    double vrad_L, vrad_U, vrad_eL, vrad_eU;

    double x0, y0, e0, x1, y1, e1, x2, y2, e2;



    // B-spline variables
    if(strcmp(param, "PA") == 0)
    {
        bspline_order = TRparam[0].pa_order_bspline;
        nbreak = TRparam[0].pa_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_pa = ncoeffs;
    }
    else if(strcmp(param, "INCL") == 0)
    {
        bspline_order = TRparam[0].incl_order_bspline;
        nbreak = TRparam[0].incl_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_incl = ncoeffs;
    }
    else if(strcmp(param, "VRAD") == 0)
    {
        bspline_order = TRparam[0].vrad_order_bspline;
        nbreak = TRparam[0].vrad_nbreak_bspline;
        ncoeffs = nbreak - 1 + bspline_order;
        TRparam[0].n_coeffs_bspline_vrad = ncoeffs;
    }
    else if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
    {
        // 1 break cubic spline fitting for vrot
        // : cubic
        //bspline_order = 3;
        //nbreak = 2;
        //ncoeffs = 4;
        // : quadratic
        bspline_order = 2;
        nbreak = 2;
        ncoeffs = 3;
    }
    else if(strcmp(param, "XPOS") == 0 || strcmp(param, "YPOS") == 0 || strcmp(param, "VSYS") == 0)
    {
        // put default values : not used
        // 1 break cubic spline fitting for vrot
        bspline_order = 3;
        nbreak = 2;
        ncoeffs = 4;
    }

    bw = gsl_bspline_alloc(bspline_order+1, nbreak);
    B = gsl_vector_alloc(ncoeffs);

    x = gsl_vector_alloc(n);
    y = gsl_vector_alloc(n);
    X = gsl_matrix_alloc(n, ncoeffs);
    c = gsl_vector_alloc(ncoeffs);
    c_temp = gsl_vector_alloc(ncoeffs);
    w = gsl_vector_alloc(n);
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    mw = gsl_multifit_linear_alloc(n, ncoeffs);

    if(strcmp(param, "XPOS") == 0 || strcmp(param, "YPOS") == 0 || strcmp(param, "VSYS") == 0)
    {
        for (i=0; i<n; i++)
        {
            // XPOS
            if(!isinf(TRparam[0].xpos[i]) && !isnan(TRparam[0].xpos[i]) &&
               !isinf(TRparam[0].xpos_e[i]) && !isnan(TRparam[0].xpos_e[i]) && TRparam[0].xpos_e[i] > 0)
            {
                xpos[i] = TRparam[0].xpos[i];
                xpos_err[i] = TRparam[0].xpos_e[i];
            }
            else
            {
                xpos[i] = 1E3; // large value
                xpos_err[i] = 1E3; // large value
            }

            // YPOS
            if(!isinf(TRparam[0].ypos[i]) && !isnan(TRparam[0].ypos[i]) &&
               !isinf(TRparam[0].ypos_e[i]) && !isnan(TRparam[0].ypos_e[i]) && TRparam[0].ypos_e[i] > 0)
            {
                ypos[i] = TRparam[0].ypos[i];
                ypos_err[i] = TRparam[0].ypos_e[i];
            }
            else
            {
                ypos[i] = 1E3; // large value
                ypos_err[i] = 1E3; // large value
            }

            // VSYS
            if(!isinf(TRparam[0].vsys[i]) && !isnan(TRparam[0].vsys[i]) &&
               !isinf(TRparam[0].vsys_e[i]) && !isnan(TRparam[0].vsys_e[i]) && TRparam[0].vsys_e[i] > 0)
            {
                vsys[i] = TRparam[0].vsys[i];
                vsys_err[i] = TRparam[0].vsys_e[i];
            }
            else
            {
                vsys[i] = 1E3; // large value
                vsys_err[i] = 1E3; // large value
            }

            // VRAD
            if(!isinf(TRparam[0].vrad[i]) && !isnan(TRparam[0].vrad[i]) &&
               !isinf(TRparam[0].vrad_e[i]) && !isnan(TRparam[0].vrad_e[i]) && TRparam[0].vrad_e[i] > 0)
            {
                vrad[i] = TRparam[0].vrad[i];
                vrad_err[i] = TRparam[0].vrad_e[i];
            }
            else
            {
                vrad[i] = 1E3; // large value
                vrad_err[i] = 1E3; // large value
            }
        
            // VROT
            if(!isinf(TRparam[0].vrot[i]) && !isnan(TRparam[0].vrot[i]) &&
               !isinf(TRparam[0].vrot_e[i]) && !isnan(TRparam[0].vrot_e[i]) && TRparam[0].vrot_e[i] > 0)
            {
                vrot[i] = TRparam[0].vrot[i];
                vrot_err[i] = TRparam[0].vrot_e[i];
            }
            else
            {
                vrot[i] = 1E3; // large value
                vrot_err[i] = 1E3; // large value
            }

            // copy ring params to gsl vectors for finding min/max of errors
            gsl_vector_set(gsl_xpos, i, xpos[i]);
            gsl_vector_set(gsl_xpos_e, i, xpos_err[i]);

            gsl_vector_set(gsl_ypos, i, ypos[i]);
            gsl_vector_set(gsl_ypos_e, i, ypos_err[i]);

            gsl_vector_set(gsl_vsys, i, vsys[i]);
            gsl_vector_set(gsl_vsys_e, i, vsys_err[i]);

            gsl_vector_set(gsl_vrot, i, vrot[i]);
            gsl_vector_set(gsl_vrot_e, i, vrot_err[i]);

            gsl_vector_set(gsl_vrad, i, vrad[i]);
            gsl_vector_set(gsl_vrad_e, i, vrad_err[i]);
        }

        // Calculate mean & std 
        // XPOS
        robust_mean_std(xpos, n, &mean_xpos, &std_xpos);
        robust_mean_std_e(xpos_err, n, &mean_xpos_e, &std_xpos_e);
        // YPOS
        robust_mean_std(ypos, n, &mean_ypos, &std_ypos);
        robust_mean_std_e(ypos_err, n, &mean_ypos_e, &std_ypos_e);
        // VSYS
        robust_mean_std(vsys, n, &mean_vsys, &std_vsys);
        robust_mean_std_e(vsys_err, n, &mean_vsys_e, &std_vsys_e);
        // VROT
        robust_mean_std(vrot, n, &mean_vrot, &std_vrot);
        robust_mean_std_e(vrot_err, n, &mean_vrot_e, &std_vrot_e);
        // VRAD
        robust_mean_std(vrad, n, &mean_vrad, &std_vrad);
        robust_mean_std_e(vrad_err, n, &mean_vrad_e, &std_vrad_e);

    }
    else if(strcmp(param, "PA") == 0)
    {
        for (i=0; i<n; i++)
        {
            // PA
            if(!isinf(TRparam[0].pa[i]) && !isnan(TRparam[0].pa[i]) &&
               !isinf(TRparam[0].pa_e[i]) && !isnan(TRparam[0].pa_e[i]) && TRparam[0].pa_e[i] > 0)
            {
                pa[i] = TRparam[0].pa[i];
                pa_err[i] = TRparam[0].pa_e[i];
            }
            else
            {
                pa[i] = 1E3; // large value
                pa_err[i] = 1E3; // large value
            }
        
            // INCL
            if(!isinf(TRparam[0].incl[i]) && !isnan(TRparam[0].incl[i]) &&
               !isinf(TRparam[0].incl_e[i]) && !isnan(TRparam[0].incl_e[i]) && TRparam[0].incl_e[i] > 0)
            {
                incl[i] = TRparam[0].incl[i];
                incl_err[i] = TRparam[0].incl_e[i];
            }
            else
            {
                incl[i] = 1E3; // large value
                incl_err[i] = 1E3; // large value
            }
        
            // VROT
            if(!isinf(TRparam[0].vrot[i]) && !isnan(TRparam[0].vrot[i]) &&
               !isinf(TRparam[0].vrot_e[i]) && !isnan(TRparam[0].vrot_e[i]) && TRparam[0].vrot_e[i] > 0)
            {
                vrot[i] = TRparam[0].vrot[i];
                vrot_err[i] = TRparam[0].vrot_e[i];
            }
            else
            {
                vrot[i] = 1E3; // large value
                vrot_err[i] = 1E3; // large value
            }

            // VRAD
            if(!isinf(TRparam[0].vrad[i]) && !isnan(TRparam[0].vrad[i]) &&
               !isinf(TRparam[0].vrad_e[i]) && !isnan(TRparam[0].vrad_e[i]) && TRparam[0].vrad_e[i] > 0)
            {
                vrad[i] = TRparam[0].vrad[i];
                vrad_err[i] = TRparam[0].vrad_e[i];
            }
            else
            {
                vrad[i] = 1E3; // large value
                vrad_err[i] = 1E3; // large value
            }

            // copy ring params to gsl vectors for finding min/max of errors
            gsl_vector_set(gsl_pa, i, pa[i]);
            gsl_vector_set(gsl_pa_e, i, pa_err[i]);

            gsl_vector_set(gsl_incl, i, incl[i]);
            gsl_vector_set(gsl_incl_e, i, incl_err[i]);

            gsl_vector_set(gsl_vrot, i, vrot[i]);
            gsl_vector_set(gsl_vrot_e, i, vrot_err[i]);

            gsl_vector_set(gsl_vrad, i, vrad[i]);
            gsl_vector_set(gsl_vrad_e, i, vrad_err[i]);
        }

        // Calculate mean & std of pa_e[j] & incl_e[j] histograms to match with total error below
        // PA
        robust_mean_std(pa, n, &mean_pa, &std_pa);
        robust_mean_std_e(pa_err, n, &mean_pa_e, &std_pa_e);
        // INCL
        robust_mean_std(incl, n, &mean_incl, &std_incl);
        robust_mean_std_e(incl_err, n, &mean_incl_e, &std_incl_e);
        // VROT
        robust_mean_std(vrot, n, &mean_vrot, &std_vrot);
        robust_mean_std_e(vrot_err, n, &mean_vrot_e, &std_vrot_e);
        // VRAD
        robust_mean_std(vrad, n, &mean_vrad, &std_vrad);
        robust_mean_std_e(vrad_err, n, &mean_vrad_e, &std_vrad_e);
    }
    else if(strcmp(param, "INCL") == 0)
    {
        // if pa_function is bspilne, TRfits was already made with PA+INCL free
        if(strcmp(TRparam[0].pa_function, "bspline") != 0)
        {
            //trfit_multinest_trfit_rings_normal("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "vrot", 'T', "sigmafactor", 'T', multinest_param, TRparam, 0, "finalfit", 'N'); // both sidesK
        }
        for (i=0; i<n; i++)
        {
            // PA
            if(!isinf(TRparam[0].pa[i]) && !isnan(TRparam[0].pa[i]) &&
               !isinf(TRparam[0].pa_e[i]) && !isnan(TRparam[0].pa_e[i]) && TRparam[0].pa_e[i] > 0)
            {
                pa[i] = TRparam[0].pa[i];
                pa_err[i] = TRparam[0].pa_e[i];
            }
            else
            {
                pa[i] = 1E3; // large value
                pa_err[i] = 1E3; // large value
            }
        
            // INCL
            if(!isinf(TRparam[0].incl[i]) && !isnan(TRparam[0].incl[i]) &&
               !isinf(TRparam[0].incl_e[i]) && !isnan(TRparam[0].incl_e[i]) && TRparam[0].incl_e[i] > 0)
            {
                incl[i] = TRparam[0].incl[i];
                incl_err[i] = TRparam[0].incl_e[i];
            }
            else
            {
                incl[i] = 1E3; // large value
                incl_err[i] = 1E3; // large value
            }
        
            // VROT
            if(!isinf(TRparam[0].vrot[i]) && !isnan(TRparam[0].vrot[i]) &&
               !isinf(TRparam[0].vrot_e[i]) && !isnan(TRparam[0].vrot_e[i]) && TRparam[0].vrot_e[i] > 0)
            {
                vrot[i] = TRparam[0].vrot[i];
                vrot_err[i] = TRparam[0].vrot_e[i];
            }
            else
            {
                vrot[i] = 1E3; // large value
                vrot_err[i] = 1E3; // large value
            }

            // VRAD
            if(!isinf(TRparam[0].vrad[i]) && !isnan(TRparam[0].vrad[i]) &&
               !isinf(TRparam[0].vrad_e[i]) && !isnan(TRparam[0].vrad_e[i]) && TRparam[0].vrad_e[i] > 0)
            {
                vrad[i] = TRparam[0].vrad[i];
                vrad_err[i] = TRparam[0].vrad_e[i];
            }
            else
            {
                vrad[i] = 1E3; // large value
                vrad_err[i] = 1E3; // large value
            }

            // copy ring params to gsl vectors for finding min/max of errors
            gsl_vector_set(gsl_pa, i, pa[i]);
            gsl_vector_set(gsl_pa_e, i, pa_err[i]);

            gsl_vector_set(gsl_incl, i, incl[i]);
            gsl_vector_set(gsl_incl_e, i, incl_err[i]);

            gsl_vector_set(gsl_vrot, i, vrot[i]);
            gsl_vector_set(gsl_vrot_e, i, vrot_err[i]);

            gsl_vector_set(gsl_vrad, i, vrad[i]);
            gsl_vector_set(gsl_vrad_e, i, vrad_err[i]);
        }

        // Calculate mean & std of pa_e[j] & incl_e[j] histograms to match with total error below
        // PA
        robust_mean_std(pa, n, &mean_pa, &std_pa);
        robust_mean_std_e(pa_err, n, &mean_pa_e, &std_pa_e);
        // INCL
        robust_mean_std(incl, n, &mean_incl, &std_incl);
        robust_mean_std_e(incl_err, n, &mean_incl_e, &std_incl_e);
        // VROT
        robust_mean_std(vrot, n, &mean_vrot, &std_vrot);
        robust_mean_std_e(vrot_err, n, &mean_vrot_e, &std_vrot_e);
        // VRAD
        robust_mean_std(vrad, n, &mean_vrad, &std_vrad);
        robust_mean_std_e(vrad_err, n, &mean_vrad_e, &std_vrad_e);
    }
    else if(strcmp(param, "VRAD") == 0)
    {
        for (i=0; i<n; i++)
        {
            // VRAD
            if(!isinf(TRparam[0].vrad[i]) && !isnan(TRparam[0].vrad[i]) &&
               !isinf(TRparam[0].vrad_e[i]) && !isnan(TRparam[0].vrad_e[i]) && TRparam[0].vrad_e[i] > 0)
            {
                vrad[i] = TRparam[0].vrad[i];
                vrad_err[i] = TRparam[0].vrad_e[i];
            }
            else
            {
                vrad[i] = 1E3; // large value
                vrad_err[i] = 1E3; // large value
            }
        
            // VROT
            if(!isinf(TRparam[0].vrot[i]) && !isnan(TRparam[0].vrot[i]) &&
               !isinf(TRparam[0].vrot_e[i]) && !isnan(TRparam[0].vrot_e[i]) && TRparam[0].vrot_e[i] > 0)
            {
                vrot[i] = TRparam[0].vrot[i];
                vrot_err[i] = TRparam[0].vrot_e[i];
            }
            else
            {
                vrot[i] = 1E3; // large value
                vrot_err[i] = 1E3; // large value
            }

            // copy ring params to gsl vectors for finding min/max of errors
            gsl_vector_set(gsl_vrad, i, vrad[i]);
            gsl_vector_set(gsl_vrad_e, i, vrad_err[i]);

            gsl_vector_set(gsl_vrot, i, vrot[i]);
            gsl_vector_set(gsl_vrot_e, i, vrot_err[i]);
        }

        // Calculate mean & std of vrad_e[j] & vrot_e[j] histograms to match with total error below
        // VRAD
        robust_mean_std(vrad, n, &mean_vrad, &std_vrad);
        robust_mean_std_e(vrad_err, n, &mean_vrad_e, &std_vrad_e);
        gsl_max_vrad = fabs(gsl_stats_max(vrad, 1, n));
        gsl_min_vrad = fabs(gsl_stats_min(vrad, 1, n));
        if(gsl_max_vrad < gsl_min_vrad)
            gsl_max_vrad = gsl_min_vrad;

        TRparam[0].vrad_max = gsl_max_vrad;

        // VROT
        robust_mean_std(vrot, n, &mean_vrot, &std_vrot);
        robust_mean_std_e(vrot_err, n, &mean_vrot_e, &std_vrot_e);
    }
    else if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
    {
        for (i=0; i<n; i++)
        {
            // VROT
            if(!isinf(TRparam[0].vrot[i]) && !isnan(TRparam[0].vrot[i]) &&
               !isinf(TRparam[0].vrot_e[i]) && !isnan(TRparam[0].vrot_e[i]) && TRparam[0].vrot_e[i] > 0)
            {
                vrot[i] = TRparam[0].vrot[i];
                vrot_err[i] = TRparam[0].vrot_e[i];
            }
            else
            {
                vrot[i] = 1E3; // large value
                vrot_err[i] = 1E3; // large value
            }

            // VRAD
            if(!isinf(TRparam[0].vrad[i]) && !isnan(TRparam[0].vrad[i]) &&
               !isinf(TRparam[0].vrad_e[i]) && !isnan(TRparam[0].vrad_e[i]) && TRparam[0].vrad_e[i] > 0)
            {
                vrad[i] = TRparam[0].vrad[i];
                vrad_err[i] = TRparam[0].vrad_e[i];
            }
            else
            {
                vrad[i] = 1E3; // large value
                vrad_err[i] = 1E3; // large value
            }

            // copy ring params to gsl vectors for finding min/max of errors
            gsl_vector_set(gsl_vrot_e, i, vrot_err[i]);
            gsl_vector_set(gsl_vrad_e, i, vrad_err[i]);
        }
 
        // VROT
        robust_mean_std(vrot, n, &mean_vrot, &std_vrot);
        robust_mean_std_e(vrot_err, n, &mean_vrot_e, &std_vrot_e);

        // VRAD
        robust_mean_std(vrad, n, &mean_vrad, &std_vrad);
        robust_mean_std_e(vrad_err, n, &mean_vrad_e, &std_vrad_e);
    }
    else
    {
        printf("B-spline approximation for %s is not available. Only PA, INCL or VRAD is available\n", param);
        exit(1);
    }

    for (i=0; i<n; i++)
    {
        if(strcmp(param, "XPOS") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].xpos[i]; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            y_dat_t[i] = _value_temp;
            e_dat[i] = TRparam[0].xpos_e[i];
        }
        else if(strcmp(param, "YPOS") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].ypos[i]; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            y_dat_t[i] = _value_temp;
            e_dat[i] = TRparam[0].ypos_e[i]; // extrapolate inward
        }
        else if(strcmp(param, "VSYS") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].vsys[i]; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            y_dat_t[i] = _value_temp;
            e_dat[i] = TRparam[0].vsys_e[i]; // extrapolate inward
        }
        else if(strcmp(param, "PA") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].pa[i]/TRparam[0].PA_MAX_in_degree; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            y_dat_t[i] = _value_temp;
            e_dat[i] = TRparam[0].pa_e[i]/TRparam[0].PA_MAX_in_degree; // extrapolate inward
        }
        else if(strcmp(param, "INCL") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].incl[i]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            y_dat_t[i] = _value_temp;
            e_dat[i] = TRparam[0].incl_e[i]/TRparam[0].INCL_MAX_in_degree; // extrapolate inward
        }
        else if(strcmp(param, "VRAD") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].vrad[i]/TRparam[0].vrad_max; // normalised to 1

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            y_dat_t[i] = _value_temp;
            e_dat[i] = TRparam[0].vrad_e[i]/TRparam[0].vrad_max; // normalised to 1
        }
        else if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].vrot[i]; // extrapolate inward

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            y_dat_t[i] = _value_temp;
            e_dat[i] = TRparam[0].vrot_e[i]; // extrapolate inward
        }
        else
        {
            printf("B-spline approximation for %s is not available. Only PA, INCL or VRAD is available\n", param);
            exit(1);
        }
    }
    // robust mean & std of the parameter
    robust_mean_std(y_dat_t, n, &robust_mean_y, &robust_std_y);


    // params & their error limits
    xpos_L = mean_xpos - 9*std_xpos;
    xpos_U = mean_xpos + 9*std_xpos;
    xpos_eL = 0;
    xpos_eU = mean_xpos_e + 3*std_xpos_e;

    ypos_L = mean_ypos - 9*std_ypos;
    ypos_U = mean_ypos + 9*std_ypos;
    ypos_eL = 0;
    ypos_eU = mean_ypos_e + 3*std_ypos_e;

    vsys_L = mean_vsys - 9*std_vsys;
    vsys_U = mean_vsys + 9*std_vsys;
    vsys_eL = 0;
    vsys_eU = mean_vsys_e + 3*std_vsys_e;

    pa_L = mean_pa - 9*std_pa;
    if(pa_L < 0) pa_L = 0;
    pa_U = mean_pa + 9*std_pa;
    if(pa_U > TRparam[0].PA_MAX_in_degree) pa_U = TRparam[0].PA_MAX_in_degree;

    pa_eL = 0;
    pa_eU = mean_pa_e + 3*std_pa_e;

    incl_L = mean_incl - 9*std_incl;
    if(incl_L < 0) incl_L = 0;
    incl_U = mean_incl + 9*std_incl;
    if(incl_U > TRparam[0].INCL_MAX_in_degree) incl_U = TRparam[0].INCL_MAX_in_degree;
    incl_eL = 0;
    incl_eU = mean_incl_e + 3*std_incl_e;

    vrot_L = 0;
    vrot_U = mean_vrot + 9*std_vrot;
    vrot_eL = 0;
    vrot_eU = mean_vrot_e + 9*std_vrot_e;

    vrad_L = mean_vrad - 9*std_vrad;
    vrad_U = mean_vrad + 9*std_vrad;
    vrad_eL = 0;
    vrad_eU = mean_vrad_e + 3*std_vrad_e;

    // Filtering out outliers based on the total error + the scatter of parameters
    n_loop = 0;
    outliered = 0;
    for (i=0; i<n; i++)
    {
        outliered = ringOutliered(TRparam, i, xpos_L, xpos_U, xpos_eL, xpos_eU,
                          ypos_L, ypos_U, ypos_eL, ypos_eU,
                          vsys_L, vsys_U, vsys_eL, vsys_eU,
                          pa_L,   pa_U,   pa_eL,   pa_eU,
                          incl_L, incl_U, incl_eL, incl_eU,
                          vrot_L, vrot_U, vrot_eL, vrot_eU,
                          vrad_L, vrad_U, vrad_eL, vrad_eU);

        if(outliered == 0)
        {
            gsl_vector_set(x, i, x_dat[i]);
            gsl_vector_set(y, i, y_dat[i]);
            gsl_vector_set(w, i, 1/e_dat[i]);
            filter_flag = 0;
        }
        else
        {
            if(i==0) // first ring 
            {
                // find the next available point
                n_loop = 0;
                while(1)
                {
                    n_loop++;
                    if(i+n_loop > n-1)
                    {
                        printf("n_loop exceeds the number of rings. Check bsplinefit_set_unipriors() for %s\n", param);
                        exit(1);
                    }
                    outliered = ringOutliered(TRparam, i+n_loop, xpos_L, xpos_U, xpos_eL, xpos_eU,
                          ypos_L, ypos_U, ypos_eL, ypos_eU,
                          vsys_L, vsys_U, vsys_eL, vsys_eU,
                          pa_L,   pa_U,   pa_eL,   pa_eU,
                          incl_L, incl_U, incl_eL, incl_eU,
                          vrot_L, vrot_U, vrot_eL, vrot_eU,
                          vrad_L, vrad_U, vrad_eL, vrad_eU);

                    if(outliered == 0)
                    {
                        y_dat[i] = y_dat[i+n_loop];
                        if(filter_flag == 0)
                            e_dat[i] = 2.0*e_dat[i+n_loop]; // 2 times the error to compensate for uncertainties
                        else if(filter_flag == 1)
                            e_dat[i] = e_dat[i+n_loop];

                        gsl_vector_set (x, i, x_dat[i]);
                        gsl_vector_set (y, i, y_dat[i]);
                        gsl_vector_set (w, i, 1/e_dat[i]);
                        
                        // found the first next available point
                        x1 = x_dat[i+n_loop];
                        y1 = y_dat[i+n_loop];
                        e1 = e_dat[i+n_loop];
            
                        if(strcmp(param, "XPOS") == 0 || strcmp(param, "YPOS") == 0 || strcmp(param, "VSYS") == 0 || strcmp(param, "PA") == 0 || strcmp(param, "INCL") == 0)
                            break;
                    }
                    else
                    {
                        while(1)
                        {
                            n_loop++;
                            if(i+n_loop > n-1)
                            {
                                printf("n_loop exceeds the number of rings. Check bsplinefit_set_unipriors() for %s\n", param);
                                exit(1);
                            }
                            outliered = ringOutliered(TRparam, i+n_loop, xpos_L, xpos_U, xpos_eL, xpos_eU,
                                  ypos_L, ypos_U, ypos_eL, ypos_eU,
                                  vsys_L, vsys_U, vsys_eL, vsys_eU,
                                  pa_L,   pa_U,   pa_eL,   pa_eU,
                                  incl_L, incl_U, incl_eL, incl_eU,
                                  vrot_L, vrot_U, vrot_eL, vrot_eU,
                                  vrad_L, vrad_U, vrad_eL, vrad_eU);

                            if(outliered == 0)
                            {
                                y_dat[i] = y_dat[i+n_loop];
                                if(filter_flag == 0)
                                    e_dat[i] = 2.0*e_dat[i+n_loop]; // 2 times the error to compensate for uncertainties
                                else if(filter_flag == 1)
                                    e_dat[i] = e_dat[i+n_loop];

                                gsl_vector_set (x, i, x_dat[i]);
                                gsl_vector_set (y, i, y_dat[i]);
                                gsl_vector_set (w, i, 1/e_dat[i]);
                                
                                // found the first next available point
                                x1 = x_dat[i+n_loop];
                                y1 = y_dat[i+n_loop];
                                e1 = e_dat[i+n_loop];
                                break;
                            }
                        }
                        if(strcmp(param, "XPOS") == 0 || strcmp(param, "YPOS") == 0 || strcmp(param, "VSYS") == 0 || strcmp(param, "PA") == 0 || strcmp(param, "INCL") == 0)
                        break;
                    }

                    // for VROT, find the second next available point
                    if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0 || strcmp(param, "VRAD") == 0)
                    {
                        while(1)
                        {
                            n_loop++;
                            if(i+n_loop > n-1)
                            {
                                // put the first one
                                y_dat[i] = y1;
                                e_dat[i] = e1;  
                                gsl_vector_set (x, i, x_dat[i]);
                                gsl_vector_set (y, i, y_dat[i]);
                                gsl_vector_set (w, i, 1/e_dat[i]);
                                printf("n_loop1 exceeds the number of rings. Check bsplinefit_set_unipriors() for %s\n", param);
                                break;
                            }
                            outliered = ringOutliered(TRparam, i+n_loop, xpos_L, xpos_U, xpos_eL, xpos_eU,
                              ypos_L, ypos_U, ypos_eL, ypos_eU,
                              vsys_L, vsys_U, vsys_eL, vsys_eU,
                              pa_L,   pa_U,   pa_eL,   pa_eU,
                              incl_L, incl_U, incl_eL, incl_eU,
                              vrot_L, vrot_U, vrot_eL, vrot_eU,
                              vrad_L, vrad_U, vrad_eL, vrad_eU);

                            if(outliered == 0)
                            {
                                // found the first next available point
                                x2 = x_dat[i+n_loop];
                                y2 = y_dat[i+n_loop];
                                e2 = e_dat[i+n_loop];
                                x0 = x_dat[0];
                                y0 = (y1-y2)/(x1-x2)*(x0-x1) + y1; // linear extrapolation to the inner point!
                                e0 = (e1+e2)/2.0;

                                y_dat[i] = y0;
                                if(filter_flag == 0)
                                    e_dat[i] = 2.0*e0; // 3 times the error to compensate for uncertainties
                                if(filter_flag == 1)
                                    e_dat[i] = e0; // 3 times the error to compensate for uncertainties
                                
                                gsl_vector_set (x, i, x_dat[i]);
                                gsl_vector_set (y, i, y_dat[i]);
                                gsl_vector_set (w, i, 1/e_dat[i]);
                                break;
                            }
                        }
                        break; // exit the main loop
                    }
                }
            }
            else if(i == n-1) // last ring : the very previous ring already filtered out
            {
                y_dat[i] = y_dat[i-1];
                if(filter_flag == 0)
                    e_dat[i] = 2.0*e_dat[i-1]; // 3 times the error to compensate for uncertainties
                else if(filter_flag == 1)
                    e_dat[i] = e_dat[i-1]; // 3 times the error to compensate for uncertainties

                gsl_vector_set (x, i, x_dat[i]);
                gsl_vector_set (y, i, y_dat[i]);
                gsl_vector_set (w, i, 1/e_dat[i]);
            }
            else
            {
                // find the next available point
                n_loop = 0;
                while(1)
                {
                    n_loop++;
                    if(i+n_loop > n-1)
                    {
                        y_dat[i] = y_dat[i-1]; // put the very previous ring 
                        e_dat[i] = e_dat[i-1]; // put the very previous ring
                        gsl_vector_set (x, i, x_dat[i]);
                        gsl_vector_set (y, i, y_dat[i]);
                        gsl_vector_set (w, i, 1/e_dat[i]);
                        printf("Check bsplinefit_set_unipriors() for %s\n", param);
                        break;
                    }
                    outliered = ringOutliered(TRparam, i+n_loop, xpos_L, xpos_U, xpos_eL, xpos_eU,
                        ypos_L, ypos_U, ypos_eL, ypos_eU,
                        vsys_L, vsys_U, vsys_eL, vsys_eU,
                        pa_L,   pa_U,   pa_eL,   pa_eU,
                        incl_L, incl_U, incl_eL, incl_eU,
                        vrot_L, vrot_U, vrot_eL, vrot_eU,
                        vrad_L, vrad_U, vrad_eL, vrad_eU);

                    if(outliered == 0)
                    {
                        // i-1 point already filtered out
                        y_dat[i] = (y_dat[i-1] + y_dat[i+n_loop])/2.0; // linear interpolation between the neighbours
                        if(filter_flag == 0)
                            e_dat[i] = 2.0*(e_dat[i-1] + e_dat[i+n_loop])/2.0;// linear interpolation between the neighbours; 5 times the error to compensate for uncertainties
                        else if(filter_flag == 1)
                            e_dat[i] = (e_dat[i-1] + e_dat[i+n_loop])/2.0;// linear interpolation between the neighbours; 5 times the error to compensate for uncertainties
                   
                        gsl_vector_set (x, i, x_dat[i]);
                        gsl_vector_set (y, i, y_dat[i]);
                        gsl_vector_set (w, i, 1/e_dat[i]);
                        break;
                    }
                }
            }
            filter_flag = 1;
        }
    }

    gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], bw);
    /* construct the fit matrix X */
    for (i=0; i<n; i++)
    {
      double xi = gsl_vector_get(x, i);
      e_dat_t[i] = e_dat[i]; // copy for scaling. see below

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);

      /* fill in row i of X */
      for (j=0; j<ncoeffs; j++)
      {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
      }
    }

    if(strcmp(param, "PA") == 0 || strcmp(param, "INCL") == 0 || strcmp(param, "VRAD") == 0)
    {
        /* do the fit */
        gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

        dof = n - ncoeffs;
        tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
        Rsq = 1.0 - chisq / tss;

        //fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", chisq / dof, Rsq);

        robust_mean_std_e(e_dat_t, n, &hist_mean_e_dat, &hist_std_e_dat);
        /* output the smoothed curve */
        {
            for (i=0; i<n; i++)
            {
                xi = x_dat[i];
                gsl_bspline_eval(xi, B, bw);
                gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
                Y_spline[i] = yi;
                Y_spline_err[i] = yerr;
                Y_spline_err_t[i] = yerr;
            }

            // compute rms of errors for setting up the boundaries of PA/INCL/VRAD below
            robust_mean_std(Y_spline_err_t, n, &hist_mean_Y_spline_err, &hist_std_Y_spline_err);
        }

        for(i=0; i<c->size; i++)
        {
            c_temp->data[i] = c->data[i]; // copy the fitted c temp
            if(strcmp(param, "PA") == 0)
            {
                TRparam[0].n_coeffs_bspline_pa = c->size;
                TRparam[0]._p_bs[i] = c->data[i];
                TRparam[0]._p_bs_tr[i] = c->data[i]; // this is used for checking boundary limit later in check_boundary_limit_pa()
                TRparam[0].bspline1pa[i] = c->data[i] - 5.0*hist_mean_Y_spline_err;
                TRparam[0].bspline2pa[i] = c->data[i] + 5.0*hist_mean_Y_spline_err;

                if(TRparam[0].bspline1pa[i] <=0) TRparam[0].bspline1pa[i] = 0.01;
                if(TRparam[0].bspline2pa[i] >=1) TRparam[0].bspline2pa[i] = 0.99;

                if(TRparam[0].bspline1pa[i] >= TRparam[0].bspline2pa[i])
                {
                    TRparam[0].bspline1pa[i] = c->data[i] - 0.1;
                    TRparam[0].bspline2pa[i] = c->data[i] + 0.1;
                    if(TRparam[0].bspline1pa[i] <=0) TRparam[0].bspline1pa[i] = 0.01;
                    if(TRparam[0].bspline2pa[i] >=1) TRparam[0].bspline2pa[i] = 0.99;
                }
                TRparam[0]._p1_tr[i] = TRparam[0].bspline1pa[i];
                TRparam[0]._p2_tr[i] = TRparam[0].bspline2pa[i];

                TRparam[0]._p_bs_t[i] = c->data[i]; // the first bspline coeffs based on the TR fit: this will be updated later in read_mtPost_set_mtPriors1()
                TRparam[0]._p_bs_e_t[i] = 1.0*hist_mean_Y_spline_err; // the first sigma based on the TR fit : this will be kept in read_mtPost_set_mtPriors1()

            }
            else if(strcmp(param, "INCL") == 0)
            {
                TRparam[0].n_coeffs_bspline_incl = c->size;
                TRparam[0]._i_bs[i] = c->data[i];
                TRparam[0]._i_bs_tr[i] = c->data[i]; // this is used for checking boundary limit later in check_boundary_limit_pa()
                TRparam[0].bspline1incl[i] = c->data[i] - 5.0*hist_mean_Y_spline_err;
                TRparam[0].bspline2incl[i] = c->data[i] + 5.0*hist_mean_Y_spline_err;
                if(TRparam[0].bspline1incl[i] <= 0) TRparam[0].bspline1incl[i] = 0.01;
                if(TRparam[0].bspline2incl[i] >= 1) TRparam[0].bspline2incl[i] = 0.99;

                if(TRparam[0].bspline1incl[i] >= TRparam[0].bspline2incl[i])
                {
                    TRparam[0].bspline1incl[i] = c->data[i] - 0.1;
                    TRparam[0].bspline2incl[i] = c->data[i] + 0.1;
                    if(TRparam[0].bspline1incl[i] <=0) TRparam[0].bspline1incl[i] = 0.01;
                    if(TRparam[0].bspline2incl[i] >=1) TRparam[0].bspline2incl[i] = 0.99;
                }
                TRparam[0]._i1_tr[i] = TRparam[0].bspline1incl[i];
                TRparam[0]._i2_tr[i] = TRparam[0].bspline2incl[i];

                TRparam[0]._i_bs_t[i] = c->data[i]; // the first bspline coeffs based on the TR fit: this will be updated later in read_mtPost_set_mtPriors1()
                TRparam[0]._i_bs_e_t[i] = 1.0*hist_mean_Y_spline_err; // the first sigma based on the TR fit : this will be kept in read_mtPost_set_mtPriors1()
            }
            else if(strcmp(param, "VRAD") == 0)
            {
                TRparam[0].n_coeffs_bspline_vrad = c->size;
                TRparam[0]._vr_bs[i] = c->data[i];
                TRparam[0]._vr_bs_tr[i] = c->data[i]; // this is used for checking boundary limit later in check_boundary_limit_pa()
                TRparam[0].bspline1vrad[i] = c->data[i] - 5*hist_mean_Y_spline_err;
                TRparam[0].bspline2vrad[i] = c->data[i] + 5*hist_mean_Y_spline_err;

                if(TRparam[0].bspline1vrad[i] >= TRparam[0].bspline2vrad[i])
                {
                    TRparam[0].bspline1vrad[i] = c->data[i] - 0.1;
                    TRparam[0].bspline2vrad[i] = c->data[i] + 0.1;
                    //if(TRparam[0].bspline1vrad[i] <=0) TRparam[0].bspline1vrad[i] = 1E-2;
                    //if(TRparam[0].bspline2vrad[i] >=1) TRparam[0].bspline2vrad[i] = 0.99;
                }
            }
        }
    }

    if(strcmp(param, "VROT-einasto_halofit_rings") == 0 || strcmp(param, "VROT-ellipsefit_rings") == 0)
    {
        gsl_max_vrot = gsl_stats_max(y_dat, 1, n);
        gsl_max_vrot_e = gsl_stats_max(e_dat, 1, n);

        for (i=0; i<n; i++)
        {
            //TRparam[0].ring_radius[i] = x_dat[i]*TRparam[0].rGalaxyPlane_pixel_max;
            //TRparam[].vrot_temp[i] = y_dat[i];
            //TRparam[0].vrot_temp_e[i] = e_dat[i];
            //TRparam[0].vrot_temp[i] = y_dat[i] / gsl_max_vrot;
            //TRparam[0].vrot_temp_e[i] = e_dat[i] / gsl_max_vrot_e / 100; // 1% of the maximum vrot
        }
    }
    else if(strcmp(param, "XPOS") == 0)
    {
        robust_mean_std(y_dat_t, n, &mean_xpos, &std_xpos);
        if(10*std_xpos < TRparam[0].ring_w)
            std_xpos = TRparam[0].ring_w / 10;

        //printf("xpos-mean:%f xpos-std:%f", mean_xpos, std_xpos);
        TRparam[0].xpos1 = mean_xpos - 5*std_xpos;
        TRparam[0].xpos2 = mean_xpos + 5*std_xpos;
    }
    else if(strcmp(param, "YPOS") == 0)
    {
        robust_mean_std(y_dat_t, n, &mean_ypos, &std_ypos);
        if(10*std_ypos < TRparam[0].ring_w)
            std_ypos = TRparam[0].ring_w / 10;

        //printf("ypos-mean:%f ypos-std:%f", mean_ypos, std_ypos);
        TRparam[0].ypos1 = mean_ypos - 5*std_ypos;
        TRparam[0].ypos2 = mean_ypos + 5*std_ypos;
    }
    else if(strcmp(param, "VSYS") == 0)
    {
        robust_mean_std(y_dat_t, n, &mean_vsys, &std_vsys);
        //printf("vsys-mean:%f vsys-std:%f", mean_vsys, std_vsys);
        TRparam[0].vsys1 = mean_vsys - 5*std_vsys;
        TRparam[0].vsys2 = mean_vsys + 5*std_vsys;

        // put a wider vsys prior to include R-mean derived from Vlos histogram
        if(TRparam[0].LOS_vel_hist_rbm <= TRparam[0].vsys1)
        {
            TRparam[0].vsys1 =  mean_vsys - 3*fabs(TRparam[0].LOS_vel_hist_rbm - TRparam[0].vsys2);
        }
        else if(TRparam[0].LOS_vel_hist_rbm >= TRparam[0].vsys2)
        {
            TRparam[0].vsys2 =  mean_vsys + 3*fabs(TRparam[0].LOS_vel_hist_rbm - TRparam[0].vsys1);
        }
    }

    free(x_dat);
    free(y_dat);
    free(y_dat_t);
    free(e_dat);
    free(e_dat_t);
    free(sigma);
    free(error);
    free(Y_spline);
    free(Y_spline_err);
    free(Y_spline_err_t);

    free(xpos);
    free(ypos);
    free(vsys);
    free(pa);
    free(incl);
    free(vrot);

    free(xpos_err);
    free(ypos_err);
    free(vsys_err);
    free(pa_err);
    free(incl_err);
    free(vrot_err);

    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(w);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);

    return 0;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GIPSY median
double getmedian_g(float *calcarray, int len)
/*------------------------------------------------------------------*/
/* Determine the median of the 'calcarray' by sorting the array.    */
/* (this is not the fastest possible way).                          */
/*------------------------------------------------------------------*/
{
   int n;

   if (len == 0) return(1E90);
   if (len == 1) return(calcarray[0]);

   sortra_c( calcarray, &len );
   n = len / 2;
   if (len%2) return( calcarray[n] );                           /* Odd length */
   else       return( 0.5*(calcarray[n-1]+calcarray[n]) );            /* Even */
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void sortra_c(float *x, int *n)
{
   int ibnd;

   if ((*n) < 2) return;
   ibnd = (*n) - 1;
   do {
      int ixch = -1;
      int j;
      for (j = 0; j < ibnd; j++) {
         if (x[j] > x[j+1]) {
            double xtemp = x[j];
            x[j] = x[j+1];
            x[j+1] = xtemp;
            ixch = j;
         }
      }
      ibnd = ixch;
   } while (ibnd != -1);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void get_rbm(float *array, int arraySize, float filter_lower, float filter_upper, float *rbm_out, float *std_out)
{
    int i=0, j=0, LOOP;
    int N_tot, N_sel, N, iter=0;

    float nelements=0, filter_l=0, filter_u=0;
    float sum_of_value=0, sum_of_square=0, square_of_sum=0;
    float variation=0, degree_of_freedom=0;
    float std_dev=0, pre_std_dev=0;
    float mean=0, pre_mean=0;

    N_tot = arraySize;
    nelements = N_tot;

    while(1)
    {
        iter++;
        if(iter>5000)
        {
            *rbm_out = mean;
            *std_out = std_dev;
            break;
        }
        for(i=0; i<N_tot; i++)
        {
            if(array[i] > filter_upper || array[i] < filter_lower)
            {
                 nelements -= 1;
                 continue;
            }
            sum_of_value += array[i];
            sum_of_square += array[i]*array[i];
        }

        degree_of_freedom = nelements - 1;
        mean = sum_of_value / nelements;
        square_of_sum = sum_of_value*sum_of_value;
        variation = (sum_of_square - (square_of_sum/nelements))/degree_of_freedom;

        if(variation < 0)
            variation = -1*variation;

        std_dev = sqrt(variation);
        if(std_dev == 0)
        {
            *rbm_out = mean;
            *std_out = std_dev;
            break;
        }

        if((fabs(mean - pre_mean)/fabs(pre_mean) < 0.2) && (fabs(std_dev - pre_std_dev)/fabs(pre_std_dev) < 0.2))
        {
            *rbm_out = mean;
            *std_out = std_dev;
            break;
        }
        else
        {
            pre_mean = mean;
            pre_std_dev = std_dev;
            filter_l = mean - 3*std_dev;
            filter_u = mean + 3*std_dev;
            sum_of_value = 0;
            sum_of_square = 0;
            nelements = N_tot;
        }
    }
    return;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void get_rbm_double(double *array, int arraySize, double filter_lower, double filter_upper, double *rbm_out, double *std_out)
{
    int i=0, j=0, LOOP;
    int N_tot, N_sel, N, iter=0;

    double nelements=0, filter_l=0, filter_u=0;
    double sum_of_value=0, sum_of_square=0, square_of_sum=0;
    double variation=0, degree_of_freedom=0;
    double std_dev=0, pre_std_dev=0;
    double mean=0, pre_mean=0;

    N_tot = arraySize;
    nelements = N_tot;

    while(1)
    {
        iter++;
        if(iter>5000)
        {
            *rbm_out = mean;
            *std_out = std_dev;
            break;
        }
        for(i=0; i<N_tot; i++)
        {
            if(array[i] > filter_upper || array[i] < filter_lower)
            {
                 nelements -= 1;
                 continue;
            }
            sum_of_value += array[i];
            sum_of_square += array[i]*array[i];
        }

        degree_of_freedom = nelements - 1;
        mean = sum_of_value / nelements;
        square_of_sum = sum_of_value*sum_of_value;
        variation = (sum_of_square - (square_of_sum/nelements))/degree_of_freedom;

        if(variation < 0)
            variation = -1*variation;

        std_dev = sqrt(variation);
        if(std_dev == 0)
        {
            *rbm_out = mean;
            *std_out = std_dev;
            break;
        }

        if((fabs(mean - pre_mean)/fabs(pre_mean) < 0.2) && (fabs(std_dev - pre_std_dev)/fabs(pre_std_dev) < 0.2))
        {
            *rbm_out = mean;
            *std_out = std_dev;
            break;
        }
        else
        {
            pre_mean = mean;
            pre_std_dev = std_dev;
            filter_l = mean - 3*std_dev;
            filter_u = mean + 3*std_dev;
            sum_of_value = 0;
            sum_of_square = 0;
            nelements = N_tot;
        }
    }
    return;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* calculate a median value from an 1D array */
double getMedian(double* array, size_t arraySize)
{
    size_t center = arraySize / 2; 

    qsort((void *)array, arraySize, sizeof(array[0]), comparisonFunctionDouble);
    if (arraySize % 2 == 1) {
        return array[center];
    } else {
    return (array[center - 1] + array[center]) / 2.0;
    }
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* calculate mean and standard deviation of 2D array: <value> <weight> */
void get_wMean_STD(double* array, double* array_w, int arraySize, double *mean_array, double *std_array)
{
    int i=0;
    double deviation, sigma_wx, sigma_w, sumsqr, weighted_mean, variance, stddeviation;
    sumsqr=0;
    sigma_wx=0;
    sigma_w=0;

    for (i=0; i<arraySize ; i++)
    {
        if(!isinf(array[i]) && !isnan(array[i]) && !isinf(array_w[i]) && !isnan(array_w[i]))
        {
            sigma_wx += array[i]*array_w[i];
            sigma_w += array_w[i];
        }
    }
    weighted_mean = sigma_wx/sigma_w;

    for (i=0; i<arraySize; i++)
    {
        if(!isinf(array[i]) && !isnan(array[i]))
        {
            deviation = array[i] - weighted_mean;
            sumsqr += deviation*deviation;
        }
    }
    variance = sumsqr/(double)arraySize;
    stddeviation = sqrt(variance) ;

    *mean_array = weighted_mean;
    *std_array = stddeviation;
}


//******************************************************************************
// Brent nonlinear equation solver
//******************************************************************************
double zero_rc_one(double a, double b, double machep, double t, double f(double XPOS, double YPOS, double i, double j, double *_p, double *_i, double k, double alpha, double n, double rGalaxyPlane_pixel, double rGalaxyPlane_pixel_max, TR_ringParameters *TRparam), double XPOS, double YPOS, double i, double j, double *_p, double *_i, double k, double alpha, double n, double rGalaxyPlane_pixel_max, TR_ringParameters *TRparam)
/*
  Parameters:

    Input, double A, B, the two endpoints of the change of sign
    interval.

    Input, double MACHEP, an estimate for the relative machine
    precision.

    Input, double T, a positive error tolerance.

    Input, double F ( double x ), the name of a user-supplied
    function which evaluates the function whose zero is being sought.

    Input, char *TITLE, a title for the problem.
*/
{
    double arg;
    int status;
    double value;

//  printf ("solving r_galaxy_plane using a nonlinear equation solver\n");
//  printf ("\n");
//  printf ("    STATUS      X               F(X)\n");
//  printf ("\n");

    status=0;
    for( ; ; )
    {
        zero_rc(a, b, t, &arg, &status, value);

        if(status < 0)
        {
            //printf ("  ZERO_RC returned an error flag!\n");
        //  return 9999999999.99999; // return non-value flag
        //  return 1; // return non-value flag
            return TRparam[0].nax1; // return non-value flag
        //  return 0; // return non-value flag
        //  return 99999.0; // return non-value flag
        //  return sqrt(pow((TRparam[0].ellipse_xpos_boxfiltered-i), 2) + pow((TRparam[0].ellipse_ypos_boxfiltered-j)/cos(TRparam[0].ellipse_incl_boxfiltered*M_PI/180.), 2));
            break;
        }
        value = f(XPOS, YPOS, i, j, _p, _i, k, alpha, n, arg, rGalaxyPlane_pixel_max, TRparam);
        //printf ( "  %8d  %14e  %14e\n", status, arg, value );

        if (status == 0)
        {
            return arg;
            break;
        }
    }
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void zero_rc( double a, double b, double t, double *arg, int *status, double value )

/******************************************************************************/
/*
  Purpose:

    ZERO_RC seeks the root of a function F(X) using reverse communication.

  Discussion:

    The interval [A,B] must be a change of sign interval for F.
    That is, F(A) and F(B) must be of opposite signs.  Then
    assuming that F is continuous implies the existence of at least
    one value C between A and B for which F(C) = 0.

    The location of the zero is determined to within an accuracy
    of 6 * MACHEPS * r8_abs ( C ) + 2 * T.

    The routine is a revised version of the Brent zero finder 
    algorithm, using reverse communication.

    Thanks to Thomas Secretin for pointing out a transcription error in the
    setting of the value of P, 11 February 2013.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 February 2013

  Author:

    John Burkardt

  Reference:

    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.

  Parameters:

    Input, double A, B, the endpoints of the change of sign interval.

    Input, double T, a positive error tolerance.

    Output, double *ARG, the currently considered point.  The user
    does not need to initialize this value.  On return with STATUS positive,
    the user is requested to evaluate the function at ARG, and return
    the value in VALUE.  On return with STATUS zero, ARG is the routine's
    estimate for the function's zero.

    Input/output, int *STATUS, used to communicate between 
    the user and the routine.  The user only sets STATUS to zero on the first 
    call, to indicate that this is a startup call.  The routine returns STATUS
    positive to request that the function be evaluated at ARG, or returns
    STATUS as 0, to indicate that the iteration is complete and that
    ARG is the estimated zero

    Input, double VALUE, the function value at ARG, as requested
    by the routine on the previous call.
*/

{
    static double c;
    static double d;
    static double e;
    static double fa;
    static double fb;
    static double fc;
    double m;
    static double machep;
    double p;
    double q;
    double r;
    double s;
    static double sa;
    static double sb;
    double tol;
    //Input STATUS = 0.
    //Initialize, request F(A).
    if ( *status == 0 )
    {
        machep = r8_epsilon ( );

        sa = a;
        sb = b;
        e = sb - sa;
        d = e;

        *status = 1;
        *arg = a;
        return;
    }
    //Input STATUS = 1.
    //Receive F(A), request F(B).
    else if ( *status == 1 )
    {
        fa = value;
        *status = 2;
        *arg = sb;
        return;
    }
    //Input STATUS = 2
    //Receive F(B).

    else if ( *status == 2 )
    {
        fb = value;
        if ( 0.0 < fa * fb )
        {
            *status = -1;
            return;
        }
        c = sa;
        fc = fa;
    }
    else
    {
        fb = value;
        if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
        {
            c = sa;
            fc = fa;
            e = sb - sa;
            d = e;
        }
    }

/*
  Compute the next point at which a function value is requested.
*/
    if ( r8_abs ( fc ) < r8_abs ( fb ) )
    {
        sa = sb;
        sb = c;
        c = sa;
        fa = fb;
        fb = fc;
        fc = fa;
    }

    tol = 2.0 * machep * r8_abs ( sb ) + t;
    m = 0.5 * ( c - sb );

    if ( r8_abs ( m ) <= tol || fb == 0.0 )
    {
        *status = 0;
        *arg = sb;
        return;
    }

    if ( r8_abs ( e ) < tol || r8_abs ( fa ) <= r8_abs ( fb ) )
    {
        e = m;
        d = e;
    }
    else
    {
        s = fb / fa;

        if ( sa == c )
        {
            p = 2.0 * m * s;
            q = 1.0 - s;
        }
        else
        {
            q = fa / fc;
            r = fb / fc;
            p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
            q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
        }

        if ( 0.0 < p )
        {
            q = - q;
        }
        else
        {
            p = - p;
        }
        s = e;
        e = d;

        if ( 2.0 * p < 3.0 * m * q - r8_abs ( tol * q ) && p < r8_abs ( 0.5 * s * q ) )
        {
          d = p / q;
        }
        else
        {
          e = m;
          d = e;
        }
    }
    sa = sb;
    fa = fb;

    if ( tol < r8_abs ( d ) )
    {
        sb = sb + d;
    }
    else if ( 0.0 < m )
    {
        sb = sb + tol;
    }
    else
    {
        sb = sb - tol;
    }
    *arg = sb;
    *status = *status + 1;

    return;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double r8_abs(double x)
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
    double value;

    if ( 0.0 <= x )
    {
        value = x;
    }
    else
    {
        value = - x;
    }
    return value;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Derive r(i,j) used for computing weight: See the relevant function()
double r_ij_pa_incl_bspline_W(int i, int j, TR_ringParameters *TRparam)
{
    int pa_n=0, pa_nF=0, pa_nS=0;
    int bspline_order;
    int incl_n=0, incl_nF=0, incl_nS=0;
    int vrad_n=0, vrad_nF=0, vrad_nS=0;
    int  n=0, n_freeParams=0;
    int pa_order_, incl_order_, root_n;
    double cos_theta, sin_theta;
    double XPOS, YPOS, VSYS, PA, INCL, VROT, VRAD;
    double _n, r_2, rho_2, _x, e, iGamma;
    double _p[999], _i[999], kapa, alpha, sersic_n;
    //double _p_w=0, _p_a[9999], _p_ae[9999], _p_b[9999], _p_be[9999];
    //double _i_w=0, _i_a[9999], _i_ae[9999], _i_b[9999], _i_be[9999];
    double _p_bs[99];
    double _i_bs[99];
    double _vr_bs[99];
    double pixelScale=0.0;
    double x, rGalaxyPlane_pixel, rGalaxyPlane_arcsec, rGalaxyPlane_pc, perimeter, perimeter_max, a, b;
    double costheta, sintheta; 

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
    double _ring_temp, _value_temp, _e_temp;

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
    XPOS = TRparam[0].xposF;
    // YPOS
    YPOS = TRparam[0].yposF;


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
    return rGalaxyPlane_pixel;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Derive r(i,j) used for computing weight: See the relevant function()
double r_ij_pa_incl_given_W(int i, int j, double _xpos, double _ypos, double _pa, double _incl, TR_ringParameters *TRparam)
{
    int root_n=0;
    double XPOS, YPOS, PA, INCL;
    double x, rGalaxyPlane_pixel;

    // parameters for nonlinear equation solver
    double rGalaxyPlane_pixel_max;
    double fx_at_starting_point=0;
    double fx_at_end_point=0;
    double r_pix;
    double x_hi_init;

    rGalaxyPlane_pixel_max = TRparam[0].rGalaxyPlane_pixel_max;

    // XPOS
    XPOS = _xpos;
    // YPOS
    YPOS = _ypos;

    // I. at starting point PA & INCL:
    PA = _pa*M_PI/180.; // in radian
    INCL = _incl*M_PI/180.; // in radian

    r_pix = 0;
    // calculate fx at _r=0 with PA + INCL given MCMC priors 
    fx_at_starting_point = sqrt(pow((-((double)i-XPOS)*sin(PA) + ((double)j-YPOS)*cos(PA)), 2) + pow(((((double)i-XPOS)*cos(PA) + ((double)j-YPOS)*sin(PA))/cos(INCL)), 2)) - r_pix;


    // II. at end point:
    root_n=0;
    while(1)
    {
        root_n++;
        // increase r_pix in steps of 0.1*r_max
        r_pix = 0.1*rGalaxyPlane_pixel_max*root_n;
        //PA = _pa*M_PI/180.; // in radian
        //INCL = _incl*M_PI/180.; // in radian
       
        // calculate fx at _r=0 with PA + INCL given MCMC priors 
        fx_at_end_point = sqrt(pow((-((double)i-XPOS)*sin(PA) + ((double)j-YPOS)*cos(PA)), 2) + pow(((((double)i-XPOS)*cos(PA) + ((double)j-YPOS)*sin(PA))/cos(INCL)), 2)) - r_pix;

        if(fx_at_starting_point*fx_at_end_point < 0)
        {
            x_hi_init = r_pix;
            break;
        }
    }

    rGalaxyPlane_pixel = gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver_given_paincl(i, j, _xpos, _ypos, _pa, _incl, TRparam, 0, x_hi_init, 1E5, "brent");
    //rGalaxyPlane_pixel = gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver(i, j, TRparam, 0, x_hi_init, 1E5, "bisection");
    //rGalaxyPlane_pixel = gsl_rGalaxyPlane_pixel_TR_nonlinearEquation_solver(i, j, TRparam, 0, x_hi_init, 1E5, "falsepos");

    return rGalaxyPlane_pixel;
}


double r8_epsilon(void)
/******************************************************************************/
/*
    R8_EPSILON returns the R8 round off unit.

  Discussion:
    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but 
      1 = ( 1 + R / 2 )

    Output, double R8_EPSILON, the double precision round-off unit.
*/
{
    double r;

    r = 1.0;
    while ( 1.0 < ( double ) ( 1.0 + r )  )
    {
        r = r / 2.0;
    }

    return (2.0 * r);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// --- End of line




