
#include "2dbat.etc.h"

// 2DBAT user defined functions
// ETC


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int countlines(FILE *ifp1)
{
    char ch;
    int lines=0; 
     
    while(!feof(ifp1))
    {
        ch = fgetc(ifp1);
        if(ch == '\n')
        {
            lines++;
        }
    }
    fclose(ifp1);
    return lines;
} 


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void add_MAP_line_multinest_output(char *multinest_output, TR_ringParameters *TRparam)
{
    int i=0, j=0, n_lines;
    int status;
    int n_params=0;
    int _i, _p, _vr;
    double logZ_dummy=0, maxLogLike_dummy=0;

    FILE *ifp1;

    ifp1 = fopen(multinest_output, "rb");
    n_lines = countlines(ifp1);

    // dynamic mt_post array
    double **mt_post = (double**)malloc(sizeof(double *) * n_lines); // n_lines x n_params
    for(i=0; i<n_lines; i++)
    {
        mt_post[i] = (double*)malloc(sizeof(double) * (TRparam[0].n_freeParams+2));
        for(j=0; j<TRparam[0].n_freeParams+2; j++)
        {
            mt_post[i][j] = 0;
        }
    }

    ifp1 = fopen(multinest_output, "rb");
    //  Reading the input DATA fileth
    for(i=0; i<n_lines; i++)
    {
        for(j=0; j<TRparam[0].n_freeParams+2; j++)
        {
            status = fscanf(ifp1, "%lf", &mt_post[i][j]);
            //printf("%e ", mt_post[i][j]);
        }
    }
    fclose(ifp1);


    //  ADD the first line where the MAP values of the parameters are written into the MULTINEST outputfile : for plotting
    ifp1 = fopen(multinest_output, "wb");
    // loglikelihood & evidence : add 2 null columns
    fprintf(ifp1, "%e\t%e\t", logZ_dummy, maxLogLike_dummy);

    // XPOS
    if(TRparam[0].xpos_fix == 'T')
    {
        fprintf(ifp1, "%e\t", TRparam[0].xposF_EinastoFit);
    }
    // YPOS
    if(TRparam[0].ypos_fix == 'T')
    {
        fprintf(ifp1, "%e\t", TRparam[0].yposF_EinastoFit);
    }
    // VSYS
    if(TRparam[0].vsys_fix == 'T')
    {
        fprintf(ifp1, "%e\t", TRparam[0].vsysF_EinastoFit);
    }
    // PA
    if(TRparam[0].pa_fix == 'T')
    {
        for(_p=0; _p<TRparam[0].n_coeffs_bspline_pa; _p++)
        {
            fprintf(ifp1, "%e\t", TRparam[0]._p_bs[_p]);
        }
    }
    // INCL
    if(TRparam[0].incl_fix == 'T')
    {
        for(_i=0; _i<TRparam[0].n_coeffs_bspline_incl; _i++)
        {
            fprintf(ifp1, "%e\t", TRparam[0]._i_bs[_i]);
        }
    }
    // n
    if(TRparam[0]._n_fix == 'T')
    {
        fprintf(ifp1, "%e\t", TRparam[0]._n);
    }
    // r_2
    if(TRparam[0].r_2_fix == 'T')
    {
        fprintf(ifp1, "%e\t", TRparam[0].r_2);
    }
    // rho_2
    if(TRparam[0].rho_2_fix == 'T')
    {
        fprintf(ifp1, "%e\t", TRparam[0].rho_2);
    }
    // VRAD
    if(TRparam[0].vrad_fix == 'T')
    {
        for(_i=0; _i<TRparam[0].n_coeffs_bspline_vrad; _i++)
        {
            fprintf(ifp1, "%e\t", TRparam[0]._vr_bs[_i]);
        }
    }
    // sigma factor
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)
    {
        fprintf(ifp1, "%e\n", TRparam[0].sigma_factor);
    }

    // e_sigma
    if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0)
    {
        fprintf(ifp1, "%e\n", TRparam[0].e_sigma_fitted);
    }


    // reprint multinest posterior output
    for(i=0; i<n_lines; i++)
    {
        for(j=0; j<TRparam[0].n_freeParams+2; j++)
        {
            fprintf(ifp1, "%e\t", mt_post[i][j]);
        }
        fprintf(ifp1, "\n");
    }

    free(mt_post);
    fclose(ifp1);
    return;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void update_Einasto_params_priors(TR_ringParameters *TRparam, float Ntimes)
{
    // n
    //
    if(fabs((TRparam[0]._n_t-TRparam[0]._n)/TRparam[0]._n) > 0.1 || \
            fabs((TRparam[0]._n-TRparam[0].Ein_n_max)/TRparam[0].Ein_n_max) < 0.2) // not converged
    {
        TRparam[0].Ein_n_min = TRparam[0]._n/Ntimes;
        TRparam[0].Ein_n_max = Ntimes*TRparam[0]._n;
        TRparam[0]._n1_t = TRparam[0].Ein_n_min;
        TRparam[0]._n2_t = TRparam[0].Ein_n_max;
        if(TRparam[0]._n1_t < TRparam[0].Ein_n_min) TRparam[0]._n1_t = TRparam[0].Ein_n_min;
        if(TRparam[0]._n2_t > TRparam[0].Ein_n_max) TRparam[0]._n2_t = TRparam[0].Ein_n_max;
    }

    // r_2
    if(fabs((TRparam[0]._r_2_t-TRparam[0].r_2)/TRparam[0].r_2) > 0.1 || \
            fabs((TRparam[0].r_2-TRparam[0].Ein_r_2_max)/TRparam[0].Ein_r_2_max) < 0.2) // not converged
    {
        TRparam[0].Ein_r_2_min = TRparam[0].r_2/Ntimes;
        TRparam[0].Ein_r_2_max = Ntimes*TRparam[0].r_2;
        TRparam[0].r_21_t = TRparam[0].Ein_r_2_min;
        TRparam[0].r_22_t = TRparam[0].Ein_r_2_max;
        if(TRparam[0].r_21_t < TRparam[0].Ein_r_2_min) TRparam[0].r_21_t = TRparam[0].Ein_r_2_min;
        if(TRparam[0].r_22_t > TRparam[0].Ein_r_2_max) TRparam[0].r_22_t = TRparam[0].Ein_r_2_max;
    }

    // rho_2
    if(fabs((TRparam[0]._rho_2_t-TRparam[0].rho_2)/TRparam[0].rho_2) > 0.1 || \
            fabs((TRparam[0].rho_2-TRparam[0].Ein_rho_2_max)/TRparam[0].Ein_rho_2_max) < 0.2) // not converged
    {
        TRparam[0].Ein_rho_2_min = TRparam[0].rho_2/Ntimes;
        TRparam[0].Ein_rho_2_max = Ntimes*TRparam[0].rho_2;
        TRparam[0].rho_21_t = TRparam[0].Ein_rho_2_min;
        TRparam[0].rho_22_t = TRparam[0].Ein_rho_2_max;
        if(TRparam[0].rho_21_t < TRparam[0].Ein_rho_2_min) TRparam[0].rho_21_t = TRparam[0].Ein_rho_2_min;
        if(TRparam[0].rho_22_t > TRparam[0].Ein_rho_2_max) TRparam[0].rho_22_t = TRparam[0].Ein_rho_2_max;
    }
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void update_multinest_params_dirtyfit(multinest_paramters *multinest_param)
{
    multinest_param[0].is = multinest_param[0]._is_d;                  // run in constant efficiency mode?
    multinest_param[0].ceff = multinest_param[0]._ceff_d;                  // run in constant efficiency mode?
    if(multinest_param[0]._nlive_d < TRparam[0].n_freeParams)
    {
        multinest_param[0]._nlive_d = 2*TRparam[0].n_freeParams;
    }
    multinest_param[0].nlive_einasto_halofit = multinest_param[0]._nlive_einasto_halofit_d;                // number of live points
    multinest_param[0].efr = multinest_param[0]._efr_d;                // set the required efficiency
    multinest_param[0].tol = multinest_param[0]._tol_d;                // tol, defines the stopping criteria
    multinest_param[0].fb = multinest_param[0]._fb_d;                  // need feedback on standard output?
    multinest_param[0].outfile = multinest_param[0]._outfile_d;                // write output files?
    multinest_param[0].maxiter = multinest_param[0]._maxiter_d;                // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
    //TRparam[0].decimX_einasto_halofit = _decimX;
    //TRparam[0].decimY_einasto_halofit = _decimY;
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void update_multinest_params_fullfit(multinest_paramters *multinest_param)
{
    multinest_param[0].is = multinest_param[0]._is;                  // run in constant efficiency mode?
    multinest_param[0].ceff = multinest_param[0]._ceff;                  // run in constant efficiency mode?
    if(multinest_param[0]._nlive < TRparam[0].n_freeParams)
    {
        multinest_param[0]._nlive = 2*TRparam[0].n_freeParams;
    }
    multinest_param[0].nlive_einasto_halofit = multinest_param[0]._nlive_einasto_halofit;                // number of live points
    multinest_param[0].efr = multinest_param[0]._efr;                // set the required efficiency
    multinest_param[0].tol = multinest_param[0]._tol;                // tol, defines the stopping criteria
    multinest_param[0].fb = multinest_param[0]._fb;                  // need feedback on standard output?
    multinest_param[0].outfile = multinest_param[0]._outfile;                // write output files?
    multinest_param[0].maxiter = multinest_param[0]._maxiter;                // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
    //TRparam[0].decimX_einasto_halofit = _decimX;
    //TRparam[0].decimY_einasto_halofit = _decimY;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void update_vrot_prior(TR_ringParameters *TRparam)
{
    int i=0;
    double _n, r_2, rho_2, _r, _x, e, iGamma, GG;

    GG = 4.302 * pow(10, -3); // pc * Mo^-1 * (km/s)^2; Gravitational Constant
    // Update VROT for TR fitting based on the Einasto fit
    _n = TRparam[0]._n;
    r_2 = TRparam[0].r_2; // in pc
    rho_2 = TRparam[0].rho_2;
    //r_2 = pow(10, (-log10(rho_2) - 0.81)/1.61) * pow(10, 3); // in pc

    _r = TRparam[0].ring_radius[TRparam[0].Nrings-1] * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // the outermost ring in pc
    _x = _r/r_2;
    e = 2.71828;
    iGamma =  gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);
    TRparam[0].vrot1 = 0;
    TRparam[0].vrot2 = 3*sqrt((4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/_r)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma));

    //for(i=0; i<TRparam[0].Nrings_intp; i++)
    //{
    //    _r = TRparam[0].ring_intp[i] * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // in pc
    //    _x = _r/r_2;
    //    iGamma =  gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);

        //printf("radius_intp: %f vrot_intp: %f vrot_input: %f vrot_input_e: %f vrot_einasto: %f\n", TRparam[0].ring_intp[i], TRparam[0].vrot_intp[i], TRparam[0].vrot_temp[i], TRparam[0].vrot_e_intp[i], sqrt((4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/_r)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma)));  
    //}
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int write_trfit_results(TR_ringParameters *TRparam, char *bayesFit_outputfile)
{
    FILE *TRfits_output;
    char dirname[1000];
    struct stat st={0};
    int error=0, writefile=1, i=0;
    double ri, ro, ring=0.;
    double vrot_EinastoFit;
    double GG;
    double vrot_error;
    double _n, r_2, rho_2, _x, _r, iGamma, e;

    e = 2.71828;
    GG = 4.302*pow(10, -3);
    sprintf(dirname, "%s/2dbat_output",TRparam[0].wdir);

    //if(stat("wallaby_2D",&st) != 0)
    if(stat(dirname, &st) != 0)
        mkdir(dirname, 0775);

    TRfits_output = fopen(bayesFit_outputfile, "w");
    if(TRfits_output != NULL) //if file opened
        writefile = 1;
    else // if not opened with some reasons..
        error = 1;

    if(error)
    {
        printf("Disc full or no permission\n");
        return 1;
    }

    if(writefile)
    {
        fprintf(TRfits_output,"! 1. Radius (pixel)\n");
        fprintf(TRfits_output,"! 2. VROT_Einasto_fit (km/s)\n");
        fprintf(TRfits_output,"! 3. VROT_both (km/s)\n");
        fprintf(TRfits_output,"! 4. VROT_rec (km/s)\n");
        fprintf(TRfits_output,"! 5. VROT_app (km/s)\n");
        fprintf(TRfits_output,"! 6. VROT_total_error (km/s) = sqrt(VROT_einasto_error**2 + VROT_asym_error**2 + VROT_disp_error**2)\n");
        fprintf(TRfits_output,"! 7. VROT_einasto_error (km/s) : propagated error of VROT_einasto from the errors of n, r_2 rho_2\n");
        fprintf(TRfits_output,"! 8. VROT_asym_error (km/s) : uncertainty by asymmetry = fabs(VROT_app - VROT_rec)/4.0 : see Swaters et al. 1999\n");
        fprintf(TRfits_output,"! 9. VROT_disp_error (km/s) : uncertainty by dispersion : mom2 or her3_sigma used\n");
        fprintf(TRfits_output,"! 10. XPOS (pixel)\n");
        fprintf(TRfits_output,"! 11. XPOS_error (pixel)\n");
        fprintf(TRfits_output,"! 12. YPOS (pixel)\n");
        fprintf(TRfits_output,"! 13. YPOS_error (pixel)\n");
        fprintf(TRfits_output,"! 14. VSYS (km/s)\n");
        fprintf(TRfits_output,"! 15. VSYS_error (km/s)\n");
        fprintf(TRfits_output,"! 16. PA (degree)\n");
        fprintf(TRfits_output,"! 17. PA_error (degree)\n");
        fprintf(TRfits_output,"! 18. INCL (degree)\n");
        fprintf(TRfits_output,"! 19. INCL_error (degree)\n");
        fprintf(TRfits_output,"! 20. _n (Einasto profile) (pixel)\n");
        fprintf(TRfits_output,"! 21. _n_error (Einasto profile) (pixel)\n");
        fprintf(TRfits_output,"! 22. r_2 (pc at 1 Mpc distance)\n");
        fprintf(TRfits_output,"! 23. r_2_error (pc at 1 Mpc distance)\n");
        fprintf(TRfits_output,"! 24. rho_2 (Mo/pc^3 at 1 Mpc distance)\n");
        fprintf(TRfits_output,"! 25. rho_2_error (Mo/pc^3 at 1 Mpc distance)\n");
        fprintf(TRfits_output,"! 26. Pixel_Scale (arcsec/pixel)\n");
        fprintf(TRfits_output,"! 27. 1 sigma of student-t distribution's error\n");
        fprintf(TRfits_output,"! 28. NU value of student-t distribution\n");
        fprintf(TRfits_output,"! 29. BIC of Einasto fit\n");
        fprintf(TRfits_output, "! %9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\t%9s\n", \
                "Radius", "VROT_Ein", \
                "VROT_both", \
                "VROT_rec", \
                "VROT_app", "VROT_te", \
                "VROT_ee", "VROT_ae", "VROT_de", \
                "XPOS", "XPOS_e", \
                "YPOS", "YPOS_e", \
                "VSYS", "VSYS_e", \
                "PA", "PA_e", \
                "INCL", "INCL_e", \
                "_n", "_n_error", \
                "r_2", "r_2_e", \
                "rho_2", "rho_2_e", \
                "pixelScale", \
                "e_sig_stu.", \
                "nu_ttu.", \
                "BIC_EinastoFit");

        for(i=0; i<TRparam[0].Nrings; i++)
        {
            //i++;
            if(i == TRparam[0].Nrings) break;

            ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
            ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
            ring = (ri+ro)/2.0;

            //vrot_error = sqrt((TRparam[0].vrot_app[i]-TRparam[0].vrot_rec[i])*(TRparam[0].vrot_app[i]-TRparam[0].vrot_rec[i])/16.0 + TRparam[0].vrot_e[i]*TRparam[0].vrot_e[i]);
            vrot_error = TRparam[0].vrot_e[i];

            _n = TRparam[0]._n;
            r_2 = TRparam[0].r_2;
            rho_2 = TRparam[0].rho_2;
            _r = ring * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265; // in pc
            _x = _r/r_2;
            iGamma =  gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);

            vrot_EinastoFit = sqrt((4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/_r)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma));
                
            if(isinf(vrot_EinastoFit) || isnan(vrot_EinastoFit)) vrot_EinastoFit = 0.0;
            if(isinf(TRparam[0].vrot[i]) || isnan(TRparam[0].vrot[i])) TRparam[0].vrot[i] = 0.0;
            if(isinf(TRparam[0].vrot_rec[i]) || isnan(TRparam[0].vrot_rec[i])) TRparam[0].vrot_rec[i] = 0.0;
            if(isinf(TRparam[0].vrot_app[i]) || isnan(TRparam[0].vrot_app[i])) TRparam[0].vrot_app[i] = 0.0;
            if(isinf(vrot_error) || isnan(vrot_error)) vrot_error = 0.0;
            if(isinf(TRparam[0].xposF_EinastoFit) || isnan(TRparam[0].xposF_EinastoFit)) TRparam[0].xposF_EinastoFit = 0.0;
            if(isinf(TRparam[0].xposF_EinastoFit_e) || isnan(TRparam[0].xposF_EinastoFit_e)) TRparam[0].xposF_EinastoFit_e = 0.0;
            if(isinf(TRparam[0].yposF_EinastoFit) || isnan(TRparam[0].yposF_EinastoFit)) TRparam[0].yposF_EinastoFit = 0.0;
            if(isinf(TRparam[0].yposF_EinastoFit_e) || isnan(TRparam[0].yposF_EinastoFit_e)) TRparam[0].yposF_EinastoFit_e = 0.0;
            if(isinf(TRparam[0].vsysF_EinastoFit) || isnan(TRparam[0].vsysF_EinastoFit)) TRparam[0].vsysF_EinastoFit = 0.0;
            if(isinf(TRparam[0].vsysF_EinastoFit_e) || isnan(TRparam[0].vsysF_EinastoFit_e)) TRparam[0].vsysF_EinastoFit_e = 0.0;
            if(isinf(TRparam[0].pa[i]) || isnan(TRparam[0].pa[i])) TRparam[0].pa[i] = 0.0;
            if(isinf(TRparam[0].pa_e[i]) || isnan(TRparam[0].pa_e[i])) TRparam[0].pa_e[i] = 0.0;
            if(isinf(TRparam[0].incl[i]) || isnan(TRparam[0].incl[i])) TRparam[0].incl[i] = 0.0;
            if(isinf(TRparam[0].incl_e[i]) || isnan(TRparam[0].incl_e[i])) TRparam[0].incl_e[i] = 0.0;
            if(isinf(TRparam[0]._n) || isnan(TRparam[0]._n)) TRparam[0]._n = 0.0;
            if(isinf(TRparam[0]._ne) || isnan(TRparam[0]._ne)) TRparam[0]._ne = 0.0;
            if(isinf(TRparam[0].r_2) || isnan(TRparam[0].r_2)) TRparam[0].r_2 = 0.0;
            if(isinf(TRparam[0].r_2e) || isnan(TRparam[0].r_2e)) TRparam[0].r_2e = 0.0;
            if(isinf(TRparam[0].rho_2) || isnan(TRparam[0].rho_2)) TRparam[0].rho_2 = 0.0;
            if(isinf(TRparam[0].rho_2e) || isnan(TRparam[0].rho_2e)) TRparam[0].rho_2e = 0.0;
            if(isinf(TRparam[0].pixelScale) || isnan(TRparam[0].pixelScale)) TRparam[0].pixelScale = 0.0;
            if(isinf(TRparam[0].e_sigma_student_TR[i]) || isnan(TRparam[0].e_sigma_student_TR[i])) TRparam[0].e_sigma_student_TR[i] = 0.0;
            if(isinf(TRparam[0].einastofit_BIC) || isnan(TRparam[0].einastofit_BIC)) TRparam[0].einastofit_BIC = 0.0;

            fprintf(TRfits_output,
"%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", \
                    ring, vrot_EinastoFit, \
                    TRparam[0].vrot[i], \
                    TRparam[0].vrot_rec[i], \
                    TRparam[0].vrot_app[i], vrot_error, \
                    TRparam[0].vrot_einasto_error[i], TRparam[0].vrot_asymmetry_error[i], TRparam[0].vrot_dispersion_error[i], \
                    TRparam[0].xpos[i], TRparam[0].xpos_e[i], \
                    TRparam[0].ypos[i], TRparam[0].ypos_e[i], \
                    TRparam[0].vsys[i], TRparam[0].vsys_e[i], \
                    TRparam[0].pa[i], TRparam[0].pa_e[i], \
                    TRparam[0].incl[i], TRparam[0].incl_e[i], \
                    TRparam[0]._n, TRparam[0]._ne, \
                    TRparam[0].r_2, TRparam[0].r_2e, \
                    TRparam[0].rho_2, TRparam[0].rho_2e, \
                    TRparam[0].pixelScale, \
                    TRparam[0].e_sigma_student_TR[i], \
                    TRparam[0]._nu_studenT, \
                    TRparam[0].einastofit_BIC);
        }
        fclose(TRfits_output);
    }
    return 0;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void read_user_input_init_params(TR_ringParameters *TRparam, multinest_paramters *multinest_param, filename_2dbat *fname, int argc, char *argv[])
{
    int i=0;
    FILE *file_exist;

    // -----------------------------------------------------------------------
    // A. 2D maps ----------------------------
    // [1. WORKING DIRECTORY= /wdir]
    sprintf(TRparam[0].wdir, "%s", argv[1]);
    // [2. 2D VELOCITY FIELD= vf.fits]
    sprintf(fname[0].fitsfile_2Dinput_VF, "%s/%s", argv[1], argv[2]);
    // [3. 2D VELOCITY FIELD ERROR= vf_e.fits]
    sprintf(fname[0].fitsfile_2Dinput_VF_e, "%s/%s", argv[1], argv[3]);
    if((file_exist = fopen(fname[0].fitsfile_2Dinput_VF_e, "r")) == NULL)
    {
        TRparam[0].vf_e_user = atof(argv[3]); // vf_e supplied by user in km/s
    }
    // [4. 2D mom0= mom0.fits]
    sprintf(fname[0].fitsfile_mom0, "%s/%s", argv[1], argv[4]);
    // [5. 2D mom2 = mom2.fits]
    sprintf(fname[0].fitsfile_mom2, "%s/%s", argv[1], argv[5]);
    if((file_exist = fopen(fname[0].fitsfile_mom2, "r")) == NULL)
    {
        TRparam[0].vdisp_user = atof(argv[5]); // vdisp supplied by user in km/s
    }
    // -----------------------------------------------------------------------
    // B. GRID option ------------------------
    // [6. Grid_X= 5]
    TRparam[0].decimX_einasto_halofit = atof(argv[6]); // sampling for x
    // [7. Grid_Y= 5]
    TRparam[0].decimY_einasto_halofit = atof(argv[7]); // sampling for y
    // [8. Grid_X_dirty= 5]
    TRparam[0].decimX_einasto_halofit_d = atof(argv[8]); // sampling for x
    // [9. Grid_Y_dirty= 5]
    TRparam[0].decimY_einasto_halofit_d = atof(argv[9]); // sampling for y
    // [10. Grid_X_TRfit= 0]
    TRparam[0].decimX_trfit = atof(argv[10]); // x bin
    // [11. Grid_Y_TRfit= 0]
    TRparam[0].decimY_trfit = atof(argv[11]); // y bin
    // [12. median_box_x= 3]
    TRparam[0].box_x = atoi(argv[12]); // median box_x
    // [13. median_box_y= 3]
    TRparam[0].box_y = atoi(argv[13]); // median box_y
    // [14. use_allPixels= 0]
    TRparam[0].use_allPixels = atoi(argv[14]); // use all pixels?

    // -----------------------------------------------------------------------
    //C. PA / INCL splines ------------------------
    // [15. RING_width= 4]
    TRparam[0].ring_w = atof(argv[15]);
    // [16. PA_function= \"bspline\"]
    strcpy(TRparam[0].pa_function, argv[16]); // "bspline"
    // [17. PA_Bspline_section= 2]
    TRparam[0].pa_nbreak_bspline = atoi(argv[17])+1;
    // [18. PA_Bspline_order= 0]
    TRparam[0].pa_order_bspline = atoi(argv[18]);
    // [19. INCL_function= \"bspline\"]
    strcpy(TRparam[0].incl_function, argv[19]); // "bspline"
    // [20. INCL_Bspline_section= 2]
    TRparam[0].incl_nbreak_bspline = atoi(argv[20])+1;
    // [21. INCL_Bspline_order= 0]
    TRparam[0].incl_order_bspline = atoi(argv[21]);
    // [22. vrad_function= \"bspline\"]
    strcpy(TRparam[0].vrad_function, argv[22]); // "bspline"
    // [23. vrad_Bspline_section= 2]
    TRparam[0].vrad_nbreak_bspline = atoi(argv[23])+1;
    // [24. vrad_Bspline_order= 0]
    TRparam[0].vrad_order_bspline = atoi(argv[24]);

    // -----------------------------------------------------------------------
    //D. Weights -----------------------------
    // [25. NU (student-t)= [1 ~ 30] : [NU=30 for Normal] [NU=1 for removing outliers]
    TRparam[0]._nu_studenT = atof(argv[25]); // radial weighting order
    // [26. Radial_WEIGHT= 0, 1 or 2 : 1/R^WEIGHT]
    TRparam[0].rwpow = atof(argv[26]); // radial weighting order
    // [27. Cosine_WEIGHT= 0, 1 or 2 : |cos(theta)|^WEIGHT]
    TRparam[0].wpow = atof(argv[27]); // cosine weight order as used in rotcur
    // [28. Free_Angle= 10]
    TRparam[0].free_angle = atof(argv[28]); // free angle within which velocities are discarded
    // [29. velocity dispersion value in (km/s)] : if 0 : mom2 will be used
    TRparam[0].e_sigma = atof(argv[29]); // Constant VF_e mode or derive optimal sigma?
    // [30. Sigma_factor] : a scale factor for vlos_e fit? T/F
    TRparam[0].sigma_factor_fix =  argv[30][0];

    // -----------------------------------------------------------------------
    //E. MULTINEST parameters : full run ----------------
    //multinest_param = (multinest_paramters *)malloc(sizeof (multinest_paramters) * 1);
    // [31. is= 0 or 1]
    multinest_param[0].is = atoi(argv[31]);                 // important nested sampling?
    // [32. ceff= 0 or 1]
    multinest_param[0].ceff = atoi(argv[32]);                   // run in constant efficiency mode?
    // [33. nlive= 50]
    multinest_param[0].nlive = atoi(argv[33]);              // number of live points
    multinest_param[0].nlive_einasto_halofit = atoi(argv[33]);              // number of live points
    // [34. efr= 0.8]
    multinest_param[0].efr = atof(argv[34]);                // set the required efficiency
    // [35. tol= 0.3]
    multinest_param[0].tol = atof(argv[35]);                // tol, defines the stopping criteria
    // [36. fb= 0]
    multinest_param[0].fb = atoi(argv[36]);                 // need feedback on standard output?
    // [37. outfile= 0]
    multinest_param[0].outfile = atoi(argv[37]);                // write output files?
    // [38. maxiter= 0]
    multinest_param[0].maxiter = atoi(argv[38]);                // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 

    //F. MULTINEST parameters : copy fullfit params ----------------
    multinest_param[0]._is = multinest_param[0].is;
    multinest_param[0]._ceff = multinest_param[0].ceff;
    multinest_param[0]._nlive_einasto_halofit = multinest_param[0].nlive_einasto_halofit;
    multinest_param[0]._efr = multinest_param[0].efr;
    multinest_param[0]._tol = multinest_param[0].tol;
    multinest_param[0]._fb = multinest_param[0].fb;
    multinest_param[0]._outfile = multinest_param[0].outfile;
    multinest_param[0]._maxiter = multinest_param[0].maxiter;

    // -----------------------------------------------------------------------
    //F. MULTINEST parameters : dirty run ----------------
    // [39. is= 0 or 1]
    multinest_param[0]._is_d = atoi(argv[39]);                 // important nested sampling?
    // [40. ceff= 0 or 1]
    multinest_param[0]._ceff_d = atoi(argv[40]);                   // run in constant efficiency mode?
    // [41. nlive= 50]
    multinest_param[0]._nlive_einasto_halofit_d = atoi(argv[41]);              // number of live points
    multinest_param[0]._nlive_d = atoi(argv[41]);              // number of live points
    // [42. efr= 0.8]
    multinest_param[0]._efr_d = atof(argv[42]);                // set the required efficiency
    // [43. tol= 0.3]
    multinest_param[0]._tol_d = atof(argv[43]);                // tol, defines the stopping criteria
    // [44. fb= 0]
    multinest_param[0]._fb_d = atoi(argv[44]);                 // need feedback on standard output?
    // [45. outfile= 0]
    multinest_param[0]._outfile_d = atoi(argv[45]);                // write output files?
    // [46. maxiter= 0]
    multinest_param[0]._maxiter_d = atoi(argv[46]);                // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
    // [47. nfilter= 0 | 1 | 2 | ...]
    TRparam[0]._nfilter = atoi(argv[47])+1;                // no. of iteration for filtering out outlying pixels which are either below or above 3*std with respect to the first TR model

    TRparam[0].decim_x0 = 0;
    TRparam[0].decim_y0 = 0;

    // default multinest parameters
    multinest_param[0].mmodal = 1;                  // do mode separation?
    multinest_param[0].ndims = 999;                 // dimensionality (no. of free parameters)
    multinest_param[0].nPar = 999;                  // total no. of parameters including free & derived parameters
    multinest_param[0].nClsPar = 999;
    multinest_param[0].updInt = 100;                // after how many iterations feedback is required & the output files should be updated
                            // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
    multinest_param[0].Ztol = -1E90;                // all the modes with logZ < Ztol are ignored
    multinest_param[0].maxModes = 100;              // expected max no. of modes (used only for memory allocation)
    for(i = 0; i < multinest_param[0].ndims; i++) multinest_param[0].pWrap[i] = 0;

    //strcpy(multinest_param[0].root, "%s/2dbat_output/multinest.");        // root for output files
    sprintf(multinest_param[0].root, "%s", argv[1]);        // root for output files
    multinest_param[0].seed = -1;       //          // random no. generator seed, if < 0 then take the seed from system clock
    multinest_param[0].resume = 0;                  // resume from a previous job?
    multinest_param[0].initMPI = 0;             // initialize MPI routines?, relevant only if compiling with MPI
                            // set it to F if you want your main program to handle MPI initialization
    multinest_param[0].logZero = -DBL_MAX;          // points with loglike < logZero will be ignored by MultiNest
                            // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied

    // C-2. File and velocity field names 
    sprintf(fname[0].bayesFit_outputfile, "%s/2dbat_output/%s.trfit.dat", argv[1], argv[2]);
    sprintf(fname[0].bayesFit_outputfile_allfree, "%s/2dbat_output/%s.trfit.allfree.dat", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_2Dinput_VF_boxfiltered, "%s/2dbat_output/%s.med_filtered.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_2Dinput_VF_boxfiltered_sigma, "%s/2dbat_output/%s.vf_e.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_2Dinput_VF_boxfiltered_sigma_ew, "%s/2dbat_output/%s.vf_ew.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_HI_VF_sigma_geo_radial_angle_w, "%s/2dbat_output/%s.geo_radial_angle_w.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_2Dinput_VF_boxfiltered_decim0, "%s/2dbat_output/%s.med_filtered.decim.x0.y0.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_2Dinput_VF_boxfiltered_decim_user, "%s/2dbat_output/%s.med_filtered.decim.x%d.y%d.fits", argv[1], argv[2], (int)TRparam[0].decimX_einasto_halofit, (int)TRparam[0].decimY_einasto_halofit);
    sprintf(fname[0].fitsfile_trfit_model, "%s/2dbat_output/%s.trmodel.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_einasto_modelvf, "%s/2dbat_output/%s.einastomodel.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_res_input_minus_trfit, "%s/2dbat_output/%s.res.input_m_trmodel.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_res_input_minus_einastofit, "%s/2dbat_output/%s.res.input_m_einastomodel.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_res_trfit_minus_einastofit, "%s/2dbat_output/%s.res.trmodel_m_einastomodel.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_einasto_halofit_geo_radial_w, "%s/2dbat_output/%s.einasto_halofit.geo_w.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_trfit_w, "%s/2dbat_output/%s.trfit.w.fits", argv[1], argv[2]);
    sprintf(fname[0].fitsfile_fract_Navail_Nall, "%s/2dbat_output/%s.trfit.navail_nall.fits", argv[1], argv[2]);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void usage_2dbat()
{
    printf("\n+ -------------------------------------------------------------------------------------- +\n");
    printf("+ 2D Bayesian Automated Tilted-ring fitter (2DBAT)                                       +\n");
    printf("+ by SE-HEON OH (KASI/ICRAR) + WALLABY KINEMATICS WORKING GROUP                          +\n");
    printf("+ -------------------------------------------------------------------------------------- +\n");
    printf("+                                                                                        +\n");
    printf("+ Development history                                                                    +\n");
    printf("+ : V.1.0 8/Dec/2016                                                                     +\n");
    printf("+                                                                                        +\n");
    printf("+ Usage                                                                                  +\n");
    printf("+ : mpirun -np [0. N-cores= 8] ./2dbat                                                   +\n");
    printf("+                                                                                        +\n");
    printf("+ A. 2D maps --------------------------------------------------------------------------- +\n");
    printf("+  [1. WORKING DIRECTORY= /wdir]                                                         +\n");
    printf("+  [2. 2D VELOCITY FIELD= vf.fits]                                                       +\n");
    printf("+  [3. 2D VELOCITY FIELD ERROR= vf_e.fits || vf_e in (km/s)]                             +\n");
    printf("+  [4. 2D mom0= mom0.fits]                                                               +\n");
    printf("+  [5. 2D mom2= mom2.fits || vdisp in (km/s)]                                            +\n");
    printf("+                                                                                        +\n");
    printf("+ B. GRID option ----------------------------------------------------------------------- +\n");
    printf("+  [6. Grid_X_EinastoFit= 5]                                                             +\n");
    printf("+  [7. Grid_Y_EinastoFit= 5]                                                             +\n");
    printf("+  [8. Grid_X_EinastoFit_dirty= 5]                                                       +\n");
    printf("+  [9. Grid_Y_EinastoFit_dirty= 5]                                                       +\n");
    printf("+  [10. Grid_X_TRfit= 0]                                                                 +\n");
    printf("+  [11. Grid_Y_TRfit= 0]                                                                 +\n");
    printf("+  [12. median_box_x= 0]                                                                 +\n");
    printf("+  [13. median_box_y= 0]                                                                 +\n");
    printf("+                                                                                        +\n");
    printf("+ C. Tilted rings ---------------------------------------------------------------------- +\n");
    printf("+  [14. RING_width= 4]                                                                   +\n");
    printf("+  [15. PA_function= \"bspline\"]                                                          +\n");
    printf("+  [16. PA_Bspline_section= 1]                                                           +\n");
    printf("+  [17. PA_Bspline_order= 0]                                                             +\n");
    printf("+  [19. INCL_function= \"bspline\"]                                                        +\n");
    printf("+  [20. INCL_Bspline_section= 1]                                                         +\n");
    printf("+  [21. INCL_Bspline_order= 0]                                                           +\n");
    printf("+  [22. vrad_function= \"bspline\"]                                                        +\n");
    printf("+  [23. vrad_Bspline_section= 0]                                                         +\n");
    printf("+  [24. vrad_Bspline_order= 0]                                                           +\n");
    printf("+                                                                                        +\n");
    printf("+ D. Weights --------------------------------------------------------------------------- +\n");
    printf("+  [25. NU (student-t)= [1 ~ 30] : [NU=30 for Normal] [NU=1 for removing outliers]       +\n");
    printf("+  [26. Radial_WEIGHT= 0, 1 or 2 : 1/R^WEIGHT]                                           +\n");
    printf("+  [27. Cosine_WEIGHT= 0, 1 or 2 : |cos(theta)|^WEIGHT]                                  +\n");
    printf("+  [28. Free_Angle= 10]                                                                  +\n");
    printf("+  [29. Velocity_dispersion= 0] if 0, mom2 will be used                                  +\n");
    printf("+  [30. sigma_factor= F] : T/F                                                           +\n");
    printf("+                                                                                        +\n");
    printf("+ E. MULTINEST parameters : full fit --------------------------------------------------- +\n");
    printf("+  [31. is= 0 or 1]                                                                      +\n");
    printf("+  [32. ceff= 0 or 1]                                                                    +\n");
    printf("+  [33. nlive= 50]                                                                       +\n");
    printf("+  [34. efr= 0.8]                                                                        +\n");
    printf("+  [35. tol= 0.3]                                                                        +\n");
    printf("+  [36. fb= 0]                                                                           +\n");
    printf("+  [37. outfile= 0]                                                                      +\n");
    printf("+  [38. maxiter= 0]                                                                      +\n");
    printf("+                                                                                        +\n");
    printf("+ F. MULTINEST parameters : dirty fit -------------------------------------------------- +\n");
    printf("+  [39. is= 0 or 1]                                                                      +\n");
    printf("+  [40. ceff= 0 or 1]                                                                    +\n");
    printf("+  [41. nlive= 50]                                                                       +\n");
    printf("+  [42. efr= 0.8]                                                                        +\n");
    printf("+  [43. tol= 0.3]                                                                        +\n");
    printf("+  [44. fb= 0]                                                                           +\n");
    printf("+  [45. outfile= 0]                                                                      +\n");
    printf("+  [46. maxiter= 0]                                                                      +\n");
    printf("+ G. No. of iterations for filtering out outlying pixels based on the first TR model --- +\n");
    printf("+  [47. nfilter= 0]                                                                      +\n");
    printf("+ -------------------------------------------------------------------------------------- +\n");
    printf("+                                                                                        +\n");
    printf("+ ! EXAMPLE ---------------------------------------------------------------------------- +\n");
    printf(" mpirun -np 5 ./2dbat \\                                                                  +\n");
    printf(" WDIR \\                                                                                  +\n");
    printf(" INPUT_VF.fits \\                                                                         +\n");
    printf(" INPUT_VF_e.fits \\                                                                       +\n");
    printf(" MOM0.fits \\                                                                             +\n");
    printf(" MOM2.fits \\                                                                             +\n");
    printf(" 2 2 2 2 0 0 9 9 0 \\                                                                     +\n");
    printf(" 4 \"bspline\" 1 0 \"bspline\" 1 0 \"bspline\" 0 0 \\                                           +\n");
    printf(" 1 0 0 10 0 F \\                                                                          +\n");
    printf(" 0 0 100 0.8 0.3 0 0 0 \\                                                                 +\n");
    printf(" 0 0 50 0.8 0.3 0 0 0 \\                                                                   +\n");
    printf(" 0                                                                                       +\n");
    printf("                                                                                         +\n");
    printf("+ -------------------------------------------------------------------------------------- +\n\n");
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void print_priors_info(char *mode, TR_ringParameters *TRparam, int n_node, int iter_find_optimal_priors)
{
    int ic=0;
    int i=0;

    double _r1, _r2, mu, sigma, mu_m_5sig, mu_p_5sig;
    double BIC;

    // BIC computation
    BIC = -2.0*TRparam[0].maxLogLikeF + TRparam[0].n_freeParams*log(TRparam[0].Npoints_in_tilted_ring);
    TRparam[0].einastofit_BIC = BIC;

    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("- %s FIT: %f pixels (Nlive: %d) (BIC: %.5E) (studenT_nu: %.2f) (%d THREADS USED)\n", mode, TRparam[0].Npoints_in_tilted_ring, multinest_param[0].nlive_einasto_halofit, BIC, TRparam[0]._nu_studenT, n_node);
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("      Total number of pixels in the area to fit (N: %f)\n", TRparam[0].Npoints_in_tilted_ring);
    printf("      PA   (bspline-knots:%d bspline-order:%d)\n", TRparam[0].pa_nbreak_bspline-1, TRparam[0].pa_order_bspline);
    printf("      INCL (bspline-knots:%d bspline-order:%d)\n", TRparam[0].incl_nbreak_bspline-1, TRparam[0].incl_order_bspline);
    printf("      VRAD (bspline-knots:%d bspline-order:%d)\n", TRparam[0].vrad_nbreak_bspline-1, TRparam[0].vrad_order_bspline);
    printf("\t++ Updated uniform priors of tilted-ring parametres ++\n");

    // XPOS
    _r1 = TRparam[0].xpos1;
    _r2 = TRparam[0].xpos2;

    mu = TRparam[0].xposF_EinastoFit;
    sigma = TRparam[0].ellipse_semi_mx_boxfiltered / 5.;

    mu_m_5sig = mu - 5*sigma;
    mu_p_5sig = mu + 5*sigma;
    printf("\t- XPOS\t\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);

    // YPOS
    _r1 = TRparam[0].ypos1;
    _r2 = TRparam[0].ypos2;

    mu = TRparam[0].yposF_EinastoFit;
    sigma = TRparam[0].ellipse_semi_mx_boxfiltered / 5.;

    mu_m_5sig = mu - 5*sigma;
    mu_p_5sig = mu + 5*sigma;
    printf("\t- YPOS\t\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);

    // VSYS
    _r1 = TRparam[0].vsys1;
    _r2 = TRparam[0].vsys2;

    mu = TRparam[0].vsysF_EinastoFit;
    sigma = TRparam[0].LOS_hist_Gfit_sigma/2.0;

    mu_m_5sig = mu - 5*sigma;
    mu_p_5sig = mu + 5*sigma;
    printf("\t- VSYS\t\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);

    // PA
    for(ic=0; ic<TRparam[0].n_coeffs_bspline_pa; ic++)
    {
        _r1 = TRparam[0].bspline1pa[ic];
        _r2 = TRparam[0].bspline2pa[ic];

        mu = TRparam[0]._p_bs_t[ic];
        sigma = TRparam[0]._p_bs_e_t[ic];

        mu_m_5sig = mu - 5*sigma;
        mu_p_5sig = mu + 5*sigma;

        if(mu_m_5sig < 0) mu_m_5sig = 0;
        if(mu_p_5sig > 1) mu_p_5sig = 1;
        printf("\t- PA-bsc-%d\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", ic, mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);
    }

    // INCL
    for(ic=0; ic<TRparam[0].n_coeffs_bspline_incl; ic++)
    {
        _r1 = TRparam[0].bspline1incl[ic];
        _r2 = TRparam[0].bspline2incl[ic];

        mu = TRparam[0]._i_bs_t[ic];
        sigma = TRparam[0]._i_bs_e_t[ic];

        mu_m_5sig = mu - 5*sigma;
        mu_p_5sig = mu + 5*sigma;

        if(mu_m_5sig < 0) mu_m_5sig = 0;
        if(mu_p_5sig > 1) mu_p_5sig = 1;
        printf("\t- INCL-bsc-%d\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", ic, mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);
    }

    // VRAD
    for(ic=0; ic<TRparam[0].n_coeffs_bspline_vrad; ic++)
    {
        _r1 = TRparam[0].bspline1vrad[ic];
        _r2 = TRparam[0].bspline2vrad[ic];

        mu = TRparam[0]._vr_bs_t[ic];
        sigma = TRparam[0]._vr_bs_e_t[ic];

        mu_m_5sig = mu - 5*sigma;
        mu_p_5sig = mu + 5*sigma;

        if(mu_m_5sig < 0) mu_m_5sig = 0;
        if(mu_p_5sig > 1) mu_p_5sig = 1;
        printf("\t- VRAD-bsc-%d\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", ic, mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);
    }

    // n
    _r1 = TRparam[0]._n1;
    _r2 = TRparam[0]._n2;

    mu = TRparam[0]._n_t;
    sigma = TRparam[0]._ne_t;

    mu_m_5sig = mu - 5*sigma;
    mu_p_5sig = mu + 5*sigma;
    if(mu_m_5sig < 0) mu_m_5sig = 0;
    printf("\t- _n\t\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);

    // r2
    _r1 = TRparam[0].r_21;
    _r2 = TRparam[0].r_22;

    mu = TRparam[0]._r_2_t;
    sigma = TRparam[0]._r_2e_t;

    mu_m_5sig = mu - 5*sigma;
    mu_p_5sig = mu + 5*sigma;
    if(mu_m_5sig < 0) mu_m_5sig = 0;
    printf("\t- r_2\t\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);

    // rho2
    _r1 = TRparam[0].rho_21;
    _r2 = TRparam[0].rho_22;

    mu = TRparam[0]._rho_2_t;
    sigma = TRparam[0]._rho_2e_t;

    mu_m_5sig = mu - 5*sigma;
    mu_p_5sig = mu + 5*sigma;
    if(mu_m_5sig < 0) mu_m_5sig = 0;
    printf("\t- rho_2\t\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);
    // SIGMA
    _r1 = TRparam[0].sigma_factor1;
    _r2 = TRparam[0].sigma_factor2;

    mu = TRparam[0].sigma_factor_mode;
    sigma = TRparam[0].sigma_factor_std;

    mu_m_5sig = mu - 5*sigma;
    mu_p_5sig = mu + 5*sigma;
    if(mu_m_5sig < 0) mu_m_5sig = 0;
    printf("\t- SIGMAf\t: (mu: %6.5E, sigma: %6.5E) : (g1-g2: %6.5E - %6.5E) : (u1-u2: %6.5E - %6.5E)\n", mu, sigma, mu_m_5sig, mu_p_5sig, _r1, _r2);

    // VROT
    printf("\t- VROT\t\t: (u1-u2: %6.5E - %6.5E)\n", TRparam[0].vrot1, TRparam[0].vrot2);
    printf("\n\n");
}

// --- End of line

