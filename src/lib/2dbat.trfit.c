#include "2dbat.trfit.h"

// 2DBAT user defined functions
// Tilted-ring fit related

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void set_Nrings(TR_ringParameters *TRparam, double *Nall_ring, double *Navail_ring)
{
    int i=0;
    int decimX, decimY;
    int x0, y0, x1, y1, box_x, box_y;
    double ring, _ring_w, _xpos, _ypos, _pa, _incl, ri, ro;
    double N_fract;

    FILE *ftemp;

    reset_HI_VF_fract_navail_nall(TRparam);
    _ring_w = TRparam[0].ring_w;
    _xpos = TRparam[0].xposF;
    _ypos = TRparam[0].yposF;
    _pa = TRparam[0].paF;
    _incl = TRparam[0].inclF;

    // --- a-4. Set tilted-ring radii ---
    TRparam[0].ring_s = TRparam[0].ring_w;
    TRparam[0].ring_e = TRparam[0].ellipse_semi_mx_boxfiltered;
    TRparam[0].Nrings_to_semi_mx = (int)(1.0*(TRparam[0].ellipse_semi_mx_boxfiltered-TRparam[0].ring_s)/TRparam[0].ring_w);
    TRparam[0].Nrings = TRparam[0].Nrings_to_semi_mx;

    for(i=0; i<2*TRparam[0].Nrings_to_semi_mx; i++)
    {
        ri = TRparam[0].ring_s + i*_ring_w - 0.5*_ring_w;
        ro = TRparam[0].ring_s + i*_ring_w + 0.5*_ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        find_Navail_Nall_pixels(_xpos, _ypos, _pa, _incl, ri, ro, TRparam, 0, Nall_ring, Navail_ring);
        if(i >= TRparam[0].Nrings_to_semi_mx-1 && (*Navail_ring)/(*Nall_ring) < 0.1)
        {
            TRparam[0].Nrings = i;
            break;
        }
    }

    ftemp = fopen("paincl.order.txt", "w");
    fprintf(ftemp, "%d %d %d %d\n", TRparam[0].pa_nbreak_bspline-1, TRparam[0].pa_order_bspline, \
                                    TRparam[0].incl_nbreak_bspline-1, TRparam[0].incl_order_bspline);
    fclose(ftemp);

    TRparam[0].sigma_factor1 = 0;
    TRparam[0].sigma_factor2 = 10;

    return;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// print out TR fit results on screen
void print_trfit_results(TR_ringParameters *TRparam)
{
    int i =0;
    double ri, ro;

    printf("!+++ TILTED-RING FIT: (XPOS:%c YPOS:%c VSYS:%c PA:%c INCL:%c VROT:%c VRAD:%c sigmaF:%c)\n", 
            TRparam[0].xpos_fix, TRparam[0].ypos_fix, TRparam[0].vsys_fix,
            TRparam[0].pa_fix, TRparam[0].incl_fix, TRparam[0].vrot_fix, TRparam[0].vrad_fix, TRparam[0].sigma_factor_fix);
    printf("\n!RADIUS\tXPOS\tXPOS_e\tYPOS\tYPOS_e\tVSYS\tVSYS_e\tPA\tPA_e\tINCL\tINCL_e\tVROT\tVROT_et\tVRAD\tVRAD_e\tVRAD_bs\tVRAD_bs_e\tNpix\tNpix\te_sigma_studenT\n");
    printf("!(pix)\t(pix)\t(pix)\t(pix)\t(pix)\t(km/s)\t(km/s)\t(deg)\t(deg)\t(deg)\t(deg)\t(km/s)\t(km/s)\t(km/s)\t(km/s)\t(km/s)\t(km/s)\t\t(%dx%d)\t(0x0)\t\n", (int)TRparam[0].decimX_trfit, (int)TRparam[0].decimX_trfit);

    for(i=0; i<(int)TRparam[0].Nrings; i++)
    {
        ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
        ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
        printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\t%.2f\t%.2f\t%.2f\n", (ri+ro)/2.0, TRparam[0].xpos[i], TRparam[0].xpos_e[i], TRparam[0].ypos[i], TRparam[0].ypos_e[i], TRparam[0].vsys[i], TRparam[0].vsys_e[i], TRparam[0].pa[i], TRparam[0].pa_e[i], TRparam[0].incl[i], TRparam[0].incl_e[i], TRparam[0].vrot[i], TRparam[0].vrot_e[i], TRparam[0].vrad[i], TRparam[0].vrad_e[i], TRparam[0].vrad_einastofit_bs[i], TRparam[0].vrad_einastofit_bs_e[i], TRparam[0].npoints_inaring[i], TRparam[0].npoints_inaring_decim0[i], TRparam[0].e_sigma_student_TR[i]);
    }
    printf("\n\n");
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int ringOutliered(TR_ringParameters *TRparam, int ringN, double xposL, double xposU, double xpos_eL, double xpos_eU,
                           double yposL, double yposU, double ypos_eL, double ypos_eU,
                           double vsysL, double vsysU, double vsys_eL, double vsys_eU,
                           double paL, double paU, double pa_eL, double pa_eU,
                           double inclL, double inclU, double incl_eL, double incl_eU,
                           double vrotL, double vrotU, double vrot_eL, double vrot_eU,
                           double vradL, double vradU, double vrad_eL, double vrad_eU)
{
    int outlier=0;
    if(TRparam[0].xpos_fix == 'T')
    {
        if(TRparam[0].xpos[ringN] < xposL || TRparam[0].xpos[ringN] > xposU || TRparam[0].xpos_e[ringN] < xpos_eL || TRparam[0].xpos_e[ringN] > xpos_eU)
            outlier++;
    }
    if(TRparam[0].ypos_fix == 'T')
    {
        if(TRparam[0].ypos[ringN] < yposL || TRparam[0].ypos[ringN] > yposU || TRparam[0].ypos_e[ringN] < ypos_eL || TRparam[0].ypos_e[ringN] > ypos_eU)
            outlier++;
    }
    if(TRparam[0].vsys_fix == 'T')
    {
        if(TRparam[0].vsys[ringN] < vsysL || TRparam[0].vsys[ringN] > vsysU || TRparam[0].vsys_e[ringN] < vsys_eL || TRparam[0].vsys_e[ringN] > vsys_eU)
            outlier++;
    }
    if(TRparam[0].pa_fix == 'T')
    {
        if(TRparam[0].pa[ringN] < paL || TRparam[0].pa[ringN] > paU || TRparam[0].pa_e[ringN] < pa_eL || TRparam[0].pa_e[ringN] > pa_eU)
            outlier++;
    }
    if(TRparam[0].incl_fix == 'T')
    {
        if(TRparam[0].incl[ringN] < inclL || TRparam[0].incl[ringN] > inclU || TRparam[0].incl_e[ringN] < incl_eL || TRparam[0].incl_e[ringN] > incl_eU)
            outlier++;
    }
    if(TRparam[0].vrot_fix == 'T')
    {
        if(TRparam[0].vrot[ringN] < vrotL || TRparam[0].vrot[ringN] > vrotU || TRparam[0].vrot_e[ringN] < vrot_eL || TRparam[0].vrot_e[ringN] > vrot_eU)
            outlier++;
    }
    if(TRparam[0].vrad_fix == 'T')
    {
        if(TRparam[0].vrad[ringN] < vradL || TRparam[0].vrad[ringN] > vradU || TRparam[0].vrad_e[ringN] < vrad_eL || TRparam[0].vrad_e[ringN] > vrad_eU)
            outlier++;
    }

    return outlier;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_trfit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0, i0=0, j0=0;
    int n_freeParams = 0;

    double chi2 = 0.0;
    double v_sigma_weighted, sigma_factor;
    double logsum_errors = 0.;  
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double NobsPoints = 0.0, NobsPoints_available=0;
    double TRmodels;
    float rbm_SN, std_SN, rbm_sigma, std_sigma;
    float HI_VF_boxFiltered_sigma_temp;
    float HI_VF_boxFiltered_SN_temp;

    double _r, _mu, _sigma;

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. XPOS prior : xpo1 ~ xpo2
    if(TRparam[0].xpos_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].xpos1 + Cube[n_freeParams]*(TRparam[0].xpos2-TRparam[0].xpos1);
        n_freeParams++;
    }
    // 2. YPOS prior : ypo1 ~ ypo2
    if(TRparam[0].ypos_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].ypos1 + Cube[n_freeParams]*(TRparam[0].ypos2-TRparam[0].ypos1);
        n_freeParams++;
    }
    // 3. VSYS prior : vsys1 ~ vsys2
    if(TRparam[0].vsys_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].vsys1 + Cube[n_freeParams]*(TRparam[0].vsys2-TRparam[0].vsys1);
        n_freeParams++;
    }
    // 4. PA prior : pa1 ~ pa2: SEHEON
    if(TRparam[0].pa_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].pa1_for_TRfit*TRparam[0].PA_MAX_in_degree + Cube[n_freeParams]*(TRparam[0].pa2_for_TRfit*TRparam[0].PA_MAX_in_degree-TRparam[0].pa1_for_TRfit*TRparam[0].PA_MAX_in_degree);
        n_freeParams++;
    }
    // 5. INCL prior : incl1 ~ incl2
    if(TRparam[0].incl_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].incl1_for_TRfit*TRparam[0].INCL_MAX_in_degree + Cube[n_freeParams]*(TRparam[0].incl2_for_TRfit*TRparam[0].INCL_MAX_in_degree-TRparam[0].incl1_for_TRfit*TRparam[0].INCL_MAX_in_degree);
        _r = Cube[n_freeParams];
        _mu = TRparam[0].ellipse_incl_boxfiltered;
        _sigma = TRparam[0].incl_EllipseFit_e;
        //Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
        n_freeParams++;
    }
    // 6. VROT prior : vrot1 ~ vrot2
    if(TRparam[0].vrot_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].vrot1 + Cube[n_freeParams]*(TRparam[0].vrot2-TRparam[0].vrot1);
        n_freeParams++;
    }
    // 7. VRAD prior : vrot1 ~ vrot2
    if(TRparam[0].vrad_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].vrad1 + Cube[n_freeParams]*(TRparam[0].vrad2-TRparam[0].vrad1);
        n_freeParams++;
    }
    // 8. sigma_factor prior : sigma_factor1 ~ sigma_factor2
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].final_fit == 'Y')
    {
        Cube[n_freeParams] = TRparam[0].sigma_factor1 + Cube[n_freeParams]*(TRparam[0].sigma_factor2-TRparam[0].sigma_factor1);
        sigma_factor = Cube[n_freeParams];
        n_freeParams++;
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && TRparam[0].final_fit == 'Y')
    {
        Cube[n_freeParams] = TRparam[0].sigma_factor1 + Cube[n_freeParams]*(TRparam[0].sigma_factor2-TRparam[0].sigma_factor1);
        sigma_factor = Cube[n_freeParams];
        n_freeParams++;
    }
    else
        sigma_factor = 1;

    // 7. Total number of points to fit in a given ring
    NobsPoints = (double)TRparam[0].Npoints_in_tilted_ring;
    chi2 = 0.;
    GLLhood0 = 0.;
    logsum_errors = 0.;
    NobsPoints_available = 0;

    for(i=0; i<TRparam[0].Npoints_in_tilted_ring; i++)
    {
        i0 = TRparam[0].tilted_ring[i][0];
        j0 = TRparam[0].tilted_ring[i][1];

        TiltedRingModel(Cube, i0, j0, TRparam, &TRmodels); // Update TRmodels

        if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma > 0 && HI_VF_weight_TRfit[0].data[j0][i0] > 0) // constant vlos_e mode
        { 
            if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0 \
            && !isinf(HI_VF_weight_TRfit[0].data[j0][i0]) \
            && !isnan(HI_VF_weight_TRfit[0].data[j0][i0]))
            {
                // geometry (perimeter + cos(theta)) weighted constant vlos_e
                // HI_VF_weight_TRfit[0].data[j0][i0] is weighted by perimeter + cos(theta) and the normalised (0 ~ 1) 
                v_sigma_weighted = TRparam[0].e_sigma * TRparam[0].scale_factor_const_vlose_w / HI_VF_weight_TRfit[0].data[j0][i0];
                // put sigma_factor which is derived from the Bayesian fit
                logsum_errors += log(v_sigma_weighted);
                chi2 += pow(((HI_VF_boxFiltered[0].data[j0][i0] - TRmodels)/v_sigma_weighted), 2);
                NobsPoints_available += 1;
            }
        }
        else if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0 && HI_VF_weight_TRfit[0].data[j0][i0] > 0) // fitted vlos_e mode : full fit
        {
            if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
            && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0 \
            && !isinf(HI_VF_mom0[0].data[j0][i0]) \
            && !isnan(HI_VF_mom0[0].data[j0][i0]) \
            && HI_VF_mom0[0].data[j0][i0] > 0.0) // no blank
            {
                v_sigma_weighted = HI_VF_boxFiltered_sigma[0].data[j0][i0] * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[j0][i0] / HI_VF_mom0[0].data[j0][i0];

                // put sigma_factor which is derived from the Bayesian fit
                v_sigma_weighted = sigma_factor*v_sigma_weighted;
                logsum_errors += log(v_sigma_weighted);
                chi2 += pow(((HI_VF_boxFiltered[0].data[j0][i0] - TRmodels)/v_sigma_weighted), 2);
                NobsPoints_available += 1;
            }
        }
        else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && HI_VF_weight_TRfit[0].data[j0][i0] > 0) // fitted vlos_e mode : partial fit
        {
            if(TRparam[0].final_fit == 'Y')
            {
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0 \
                && !isinf(HI_VF_mom0[0].data[j0][i0]) \
                && !isnan(HI_VF_mom0[0].data[j0][i0]) \
                && HI_VF_mom0[0].data[j0][i0] > 0.0) // no blank
                {
                    v_sigma_weighted = sigma_factor*HI_VF_boxFiltered_sigma[0].data[j0][i0] * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[j0][i0] / HI_VF_mom0[0].data[j0][i0];

                    // put sigma_factor which is derived from the Bayesian fit
                    //v_sigma_weighted = sigma_factor*v_sigma_weighted;
                    logsum_errors += log(v_sigma_weighted);
                    chi2 += pow(((HI_VF_boxFiltered[0].data[j0][i0] - TRmodels)/v_sigma_weighted), 2);
                    NobsPoints_available += 1;
                }
            }
            else if(TRparam[0].final_fit == 'N')
            {
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0 \
                && !isinf(HI_VF_mom0[0].data[j0][i0]) \
                && !isnan(HI_VF_mom0[0].data[j0][i0]) \
                && HI_VF_mom0[0].data[j0][i0] > 0.0) // no blank
                {
                    //v_sigma_weighted = TRparam[0].e_sigma * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[j0][i0] / HI_VF_mom0[0].data[j0][i0];
                    //v_sigma_weighted = HI_VF_boxFiltered_sigma[0].data[j0][i0] * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[j0][i0];

                    v_sigma_weighted = HI_VF_boxFiltered_sigma[0].data[j0][i0] * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[j0][i0] / HI_VF_mom0[0].data[j0][i0];
                    // put sigma_factor which is derived from the Bayesian fit
                    logsum_errors += log(v_sigma_weighted);
                    chi2 += pow(((HI_VF_boxFiltered[0].data[j0][i0] - TRmodels)/v_sigma_weighted), 2);
                    NobsPoints_available += 1;
                }
            }
        }
    }

    GLLhood0 = -(NobsPoints_available/2.0)*log(2.0*M_PI);
    slhood = GLLhood0 - logsum_errors - chi2/2.0;
    *lnew = slhood;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_trfit_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0, i0=0, j0=0;
    int n_freeParams = 0;

    double logsum_errors = 0., mean_errors=0, sum_errors=0;
    double chi2 = 0.0;
    double e_sigma, sigma_factor;
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double NobsPoints = 0.0, NobsPoints_available=0;
    double TRmodels;
    float rbm_SN, std_SN, rbm_sigma, std_sigma;
    float HI_VF_boxFiltered_sigma_temp;
    float HI_VF_boxFiltered_SN_temp;
    // uniform prior
    double _r1, _r2;

    double _r, _mu, _mu_max, _sigma, _sigma_N, _sigma_W;

    // for student-t distribution parameters
    double logsum_sigma2=0, logsum_student_chi=0, _nu, log_Likelihood_studenT;


    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. XPOS prior : xpo1 ~ xpo2
    if(TRparam[0].xpos_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0].xposF_EinastoFit_t;
        _sigma = TRparam[0].ellipse_semi_mx_boxfiltered / 5.;
        Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));

        if(Cube[n_freeParams] >= TRparam[0].nax1) Cube[n_freeParams] = TRparam[0].nax1/2.0;
        n_freeParams++;

        // uniform prior
        //Cube[n_freeParams] = TRparam[0].xpos1 + Cube[n_freeParams]*(TRparam[0].xpos2-TRparam[0].xpos1);
        //n_freeParams++;
    }
    // 2. YPOS prior : ypo1 ~ ypo2
    if(TRparam[0].ypos_fix == 'T')
    {
        _r = Cube[n_freeParams];
        _mu = TRparam[0].yposF_EinastoFit_t;
        _sigma = TRparam[0].ellipse_semi_mx_boxfiltered / 5.;
        Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));

        if(Cube[n_freeParams] >= TRparam[0].nax2) Cube[n_freeParams] = TRparam[0].nax2/2.0;
        n_freeParams++;

        // uniform prior
        //Cube[n_freeParams] = TRparam[0].ypos1 + Cube[n_freeParams]*(TRparam[0].ypos2-TRparam[0].ypos1);
        //n_freeParams++;
    }
    // 3. VSYS prior : vsys1 ~ vsys2
    if(TRparam[0].vsys_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0].vsysF_EinastoFit_t;
        _sigma = TRparam[0].LOS_hist_Gfit_sigma/1.0;
        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        //if(Cube[n_freeParams] <= TRparam[0].vlos_lower_limit) Cube[n_freeParams] = TRparam[0].vlos_lower_limit;
        //if(Cube[n_freeParams] >= TRparam[0].vlos_upper_limit) Cube[n_freeParams] = TRparam[0].vlos_upper_limit;
        n_freeParams++;

        // uniform
        //Cube[n_freeParams] = TRparam[0].vsys1 + Cube[n_freeParams]*(TRparam[0].vsys2-TRparam[0].vsys1);
        //n_freeParams++;
    }
    // 4. PA prior : pa1 ~ pa2: SEHEON
    if(TRparam[0].pa_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0].ellipse_pa_boxfiltered; // in unit
        if(_mu < 180) // take the wider wing
        {
            _sigma = (360-_mu) / 5.0;
        }
        else
        {
            _sigma = (_mu) / 5.0;
        }


        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        if(Cube[n_freeParams] > 360) Cube[n_freeParams] = Cube[n_freeParams] - (int)(Cube[n_freeParams]/360)*360;
        if(Cube[n_freeParams] < 0) Cube[n_freeParams] = Cube[n_freeParams] + (1-(int)(Cube[n_freeParams]/360))*360;
        n_freeParams++;

        // uniform
        //Cube[n_freeParams] = TRparam[0].pa1_for_TRfit*TRparam[0].PA_MAX_in_degree + Cube[n_freeParams]*(TRparam[0].pa2_for_TRfit*TRparam[0].PA_MAX_in_degree-TRparam[0].pa1_for_TRfit*TRparam[0].PA_MAX_in_degree);
        //n_freeParams++;
    }
    // 5. INCL prior : incl1 ~ incl2
    if(TRparam[0].incl_fix == 'T')
    {
        _r = Cube[n_freeParams];
        // update priors based on either the tr or the dirty fits 
        _mu = TRparam[0].ellipse_incl_boxfiltered; // in degree
        if(_mu < 45) // take the wider wing
        {
            _sigma = (90-_mu) / 5.0;
        }
        else
        {
            _sigma = (_mu) / 5.0;
        }
        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);

        if(Cube[n_freeParams] > 90 && Cube[n_freeParams] <= 180) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/180)+1)*180 - Cube[n_freeParams];

        if(Cube[n_freeParams] > 180 && Cube[n_freeParams] <= 270) Cube[n_freeParams] = 90 - (((int)(Cube[n_freeParams]/270)+1)*270 - Cube[n_freeParams]);


        if(Cube[n_freeParams] > 270 && Cube[n_freeParams] <= 360) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/360)+1)*360 - Cube[n_freeParams];

        if(Cube[n_freeParams] < 0 && Cube[n_freeParams] >= -90) Cube[n_freeParams] = fabs(Cube[n_freeParams]);
        if(Cube[n_freeParams] < -90 && Cube[n_freeParams] >= -180) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/180)+1)*180 + Cube[n_freeParams];

        if(Cube[n_freeParams] < -180 && Cube[n_freeParams] >= -270) Cube[n_freeParams] = ((int)(Cube[n_freeParams]/270)+1)*270 + Cube[n_freeParams];
        n_freeParams++;

        //Cube[n_freeParams] = TRparam[0].incl1_for_TRfit*TRparam[0].INCL_MAX_in_degree + Cube[n_freeParams]*(TRparam[0].incl2_for_TRfit*TRparam[0].INCL_MAX_in_degree-TRparam[0].incl1_for_TRfit*TRparam[0].INCL_MAX_in_degree);
        //n_freeParams++;
    }
    // 6. VROT prior : vrot1 ~ vrot2
    if(TRparam[0].vrot_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].vrot1 + Cube[n_freeParams]*(TRparam[0].vrot2-TRparam[0].vrot1);
        n_freeParams++;
    }
    // 7. VRAD prior : vrot1 ~ vrot2
    if(TRparam[0].vrad_fix == 'T')
    {
        Cube[n_freeParams] = TRparam[0].vrad1 + Cube[n_freeParams]*(TRparam[0].vrad2-TRparam[0].vrad1);
        n_freeParams++;
    }
    // 8. sigma_factor prior : sigma_factor1 ~ sigma_factor2
    //if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma >= 0)
    if(TRparam[0].e_sigma >= 0)
    {
        // update priors based on either the dirty fits 
        //_mu = 9; // initial default value (5 - 15 km/s in spirals (van der Kruit & Shostak 1982; see all Tamburro et al. 2009)
        //_mu_max = 1E3;  // maximum value of e-sigma
        //_sigma_N = _mu / 3.; // sigma for the narrow side
        //_sigma_W = (_mu_max-_mu) / 5.; // sigma for the wide side
        //Cube[n_freeParams] = gaussian_prior_esigma_skew(_r, _mu, _mu_max, _sigma_N, _sigma_W);

        _r = Cube[n_freeParams];
        _r1 = 0.1;
        _r2 = 10;
        _mu = 5; // default mean of e_sigma
        _sigma = _mu/3.;
        //Cube[n_freeParams] = fabs(gaussian_prior(_r, _mu, _sigma));
        Cube[n_freeParams] = gaussian_prior(_r, _mu, _sigma);
        //Cube[n_freeParams] = uniform_priors(_r, _r1, _r2);
        e_sigma = Cube[n_freeParams];
        n_freeParams++;
    }

    // 7. Total number of points to fit in a given ring
    NobsPoints = (double)TRparam[0].Npoints_in_tilted_ring;
    chi2 = 0.;
    GLLhood0 = 0.;
    logsum_errors = 0.;
    NobsPoints_available = 0;

    sum_errors = 0.;
    mean_errors = 0.;
    // _nu of student-T distribution : _nu = 30 for normal : _nu = 1 recommended for best removing the outliers
    _nu = TRparam[0]._nu_studenT;

    logsum_student_chi = 0;
    logsum_sigma2 = 0;


    for(i=0; i<TRparam[0].Npoints_in_tilted_ring; i++)
    {
        i0 = TRparam[0].tilted_ring[i][0];
        j0 = TRparam[0].tilted_ring[i][1];

        TiltedRingModel(Cube, i0, j0, TRparam, &TRmodels); // Update TRmodels

        //if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma >= 0 && HI_VF_weight_TRfit[0].data[j0][i0] > 0) // fitted vlos_e mode : partial fit
        if(TRparam[0].e_sigma >= 0)
        {
            if(TRparam[0].final_fit == 'Y' || TRparam[0].final_fit == 'N')
            {
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0 \
                && !isinf(HI_VF_mom0[0].data[j0][i0]) \
                && !isnan(HI_VF_mom0[0].data[j0][i0]) \
                && HI_VF_mom0[0].data[j0][i0] > 0.0) // no blank
                {
                    logsum_errors += 0.5*log(M_PI*_nu*e_sigma*e_sigma);
                    logsum_student_chi += ((1+_nu)/2.0)*log(1.0+pow((HI_VF_boxFiltered[0].data[j0][i0] - TRmodels)/e_sigma, 2)/_nu);
                    NobsPoints_available += 1;
                }
            }
        }
    }

    log_Likelihood_studenT = NobsPoints_available*(log(gsl_sf_gamma((_nu+1)/2.0)) - log(gsl_sf_gamma(_nu/2.0))) - logsum_errors - logsum_student_chi;
    *lnew = log_Likelihood_studenT;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int TiltedRingModel(double *Cube, int i, int j, TR_ringParameters *TRparam, double *Vmodel_TR) // Update TRmodels
{
    int n=0, n_freeParams=0;
    double XPOS, YPOS, VSYS, PA, INCL, VROT, VRAD, GG;
    double pixelScale=0.0;
    double rGalaxyPlane_pixel, rGalaxyPlane_arcsec, rGalaxyPlane_pc;
    double cos_theta, sin_theta;
    GG = 4.302; // kpc Msol^-1 (km/s)^2

    pixelScale = TRparam[0].pixelScale; // arcsec/pixel

    if(TRparam[0].xpos_fix == 'T')
    {
        XPOS = Cube[n_freeParams]; // xpos fit
        n_freeParams++;
    }
    else
        XPOS = TRparam[0].xposF; // xpos fixed to the derived one

    if(TRparam[0].ypos_fix == 'T')
    {
        YPOS = Cube[n_freeParams]; // ypos fit
        n_freeParams++;
    }
    else
        YPOS = TRparam[0].yposF; // ypos fixed to the derived one

    if(TRparam[0].vsys_fix == 'T')
    {
        VSYS = Cube[n_freeParams]; // vsys fit
        n_freeParams++;
    }
    else
        VSYS = TRparam[0].vsysF; // vsys fixed to the derived one

    if(TRparam[0].pa_fix == 'T')
    {
        PA = Cube[n_freeParams]; // pa fit : in degree
        n_freeParams++;
    }
    else
        PA = TRparam[0].paF; // pa (in degree) fixed to the derived one

    if(TRparam[0].incl_fix == 'T')
    {
        INCL = Cube[n_freeParams]; // incl in degree fit
        n_freeParams++;
    }
    else
        INCL = TRparam[0].inclF; // incl (in degree) fixed to the derived one

    /* Derive radii in galaxy plane in pixel, arcsec and pc units */
    rGalaxyPlane_pixel = sqrt(pow((((double)j-YPOS)*cos(PA*M_PI/180.)-((double)i-XPOS)*sin(PA*M_PI/180.)), 2) + pow((((double)i-XPOS)*cos(PA*M_PI/180.) + ((double)j-YPOS)*sin(PA*M_PI/180.))/cos(INCL*M_PI/180.), 2));

    if(TRparam[0].vrot_fix == 'T')
    {
        VROT = Cube[n_freeParams]; // vrot fit
        n_freeParams++;
    }

    if(TRparam[0].vrad_fix == 'T')
    {
        VRAD = Cube[n_freeParams]; // vrad fit
        n_freeParams++;
    }
    else
        VRAD = TRparam[0].vradF;

    cos_theta = (((double)j-YPOS)*cos(PA*M_PI/180.) - ((double)i-XPOS)*sin(PA*M_PI/180.))/rGalaxyPlane_pixel;
    sin_theta = (((double)j-YPOS)*sin(PA*M_PI/180.) + ((double)i-XPOS)*cos(PA*M_PI/180.))/(rGalaxyPlane_pixel*cos(INCL*M_PI/180.));

    /* Calculate a model line-of-sight velocity given ring parameters */
    //*Vmodel_TR = VSYS + VROT*sin(INCL*M_PI/180.)*(((double)j-YPOS)*cos(PA*M_PI/180.) - ((double)i-XPOS)*sin(PA*M_PI/180.))/rGalaxyPlane_pixel;
    *Vmodel_TR = VSYS + VROT*cos_theta*sin(INCL*M_PI/180.) + VRAD*sin_theta*sin(INCL*M_PI/180.);

    return 0;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dumper_TRfits(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam)
{
    // convert the 2D Fortran arrays to C arrays
    int i=0, j=0;
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

    if(TRparam[0].xpos_fix == 'T')
    {
        /* save the best fit for xpos */
        TRparam[0].xposF = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].xposF_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    if(TRparam[0].ypos_fix == 'T')
    {
        /* save the best fit for ypos */
        TRparam[0].yposF = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].yposF_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    if(TRparam[0].vsys_fix == 'T')
    {
        /* save the best fit for vsys */
        TRparam[0].vsysF = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].vsysF_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    if(TRparam[0].pa_fix == 'T') // in degree
    {
        /* save the best fit for pa */
        TRparam[0].paF = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].paF_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    if(TRparam[0].incl_fix == 'T') // in degree
    {
        /* save the best fit for incl */
        TRparam[0].inclF = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].inclF_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    if(TRparam[0].vrot_fix == 'T')
    {
        /* save the best fit for vrot */
        TRparam[0].vrotF = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].vrotF_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    if(TRparam[0].vrad_fix == 'T')
    {
        /* save the best fit for vrad */
        TRparam[0].vradF = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].vradF_e = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    //if(TRparam[0].sigma_factor_fix == 'T')
    //{
    //    /* save the best fit for sigma_factor */
    //    TRparam[0].sigma_factor = paramConstr[0][*nPar*2+n_freeParams];
    //    TRparam[0].sigma_factor_e = paramConstr[0][*nPar+n_freeParams];
    //    n_freeParams++;
    //}
    if(TRparam[0].e_sigma >= 0)
    {
        /* save the best fit for sigma_factor */
        TRparam[0].e_sigma_tr = paramConstr[0][*nPar*2+n_freeParams];
        TRparam[0].e_sigma_e_tr = paramConstr[0][*nPar+n_freeParams];
        n_freeParams++;
    }

    /* Save current maximum loglikelihood value, evidence and error values */
    TRparam[0].maxLogLikeF = *maxLogLike;
    //TRparam[0].logZF = *logZ;
    //TRparam[0].logZerrF = *logZerr;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void trfit_multinest_ellipsefit_rings_student(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *vrot, char vrotfix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, int find_symmetric_VF, char *finalfit, char final_fit)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0, j=0, k=0, i1, ii=0;
    int i0=0, j0=0;
    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    double efr, tol, Ztol, logZero;
    double total_Npoints_allrings=0.;
    double ri=0., ro=0., ri_tr=0, ro_tr=0, r_tr, ring;
    double xpos_iso1Dfit, ypos_iso1Dfit, pa_iso1Dfit, incl_iso1Dfit;
    double pa_temp[9999], incl_temp[9999], vrot_temp[9999], vrot_median;
    double ellipse_a, ellipse_b, xpos_prior_gab, ypos_prior_gab;

    double peri_x, peri_y, xpos_temp, ypos_temp;
    float *rmax_ellipse = malloc(sizeof(float)*(TRparam[0].xpos2-TRparam[0].xpos1)*(TRparam[0].ypos2-TRparam[0].ypos1));

    // GSL histogram statistics
    double gsl_mean_xpos, gsl_std_xpos, gsl_max_xpos, gsl_min_xpos;
    double gsl_mean_xpos_e, gsl_std_xpos_e, gsl_max_xpos_e, gsl_min_xpos_e;
    double gsl_mean_ypos, gsl_std_ypos, gsl_max_ypos, gsl_min_ypos;
    double gsl_mean_ypos_e, gsl_std_ypos_e, gsl_max_ypos_e, gsl_min_ypos_e;
    double gsl_mean_vsys, gsl_std_vsys, gsl_max_vsys, gsl_min_vsys;
    double gsl_mean_vsys_e, gsl_std_vsys_e, gsl_max_vsys_e, gsl_min_vsys_e;
    double gsl_mean_pa, gsl_std_pa, gsl_max_pa, gsl_min_pa;
    double gsl_mean_pa_e, gsl_std_pa_e, gsl_max_pa_e, gsl_min_pa_e;
    double gsl_mean_incl, gsl_std_incl, gsl_max_incl, gsl_min_incl;
    double gsl_mean_incl_e, gsl_std_incl_e, gsl_max_incl_e, gsl_min_incl_e;
    double gsl_mean_vrot, gsl_std_vrot, gsl_max_vrot, gsl_min_vrot;
    double gsl_mean_vrot_e, gsl_std_vrot_e, gsl_max_vrot_e, gsl_min_vrot_e;
    double gsl_mean_vrad, gsl_std_vrad, gsl_max_vrad, gsl_min_vrad;
    double gsl_mean_vrad_e, gsl_std_vrad_e, gsl_max_vrad_e, gsl_min_vrad_e;

    // the number of tilted-rings to fit
    TRparam[0].Nrings = (int)((TRparam[0].ring_e-TRparam[0].ring_s)/TRparam[0].ring_w)*2;

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    ndims = multinest_param[0].ndims;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive_einasto_halofit;
    if(nlive > 50) nlive = 50; // for a quick TRfit
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;

    //strcpy(root, multinest_param[0].root);
    //strncpy(root, multinest_param[0].root, strlen(multinest_param[0].root));
    strcpy(root, "trfit.ellipsefit_rings.");

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    TRparam[0].xpos_fix = xposfix;
    TRparam[0].ypos_fix = yposfix;
    TRparam[0].vsys_fix = vsysfix;
    TRparam[0].pa_fix = pafix;
    TRparam[0].incl_fix = inclfix;
    TRparam[0].vrot_fix = vrotfix;
    TRparam[0].vrad_fix = vradfix;
    TRparam[0].sigma_factor_fix = sigmafactorfix;
    TRparam[0].final_fit = final_fit;

    // dynamic total_error 1D array
    double *total_error_ring_params = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array
    double *total_error_ring_params_temp = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array
    double lower_bound_total_error_ring_params, upper_bound_total_error_ring_params;
    double hist_bin_total_error_ring_params, hist_mean_total_error_ring_params, hist_std_total_error_ring_params;

    // xpos
    double *xpos_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *xpos_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_xpos, hist_std_xpos;
    // xpos error
    double *xpos_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_xpos_err, hist_std_xpos_err;

    // ypos
    double *ypos_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *ypos_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_ypos, hist_std_ypos;
    // ypos error
    double *ypos_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_ypos_err, hist_std_ypos_err;

    // vsys
    double *vsys_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vsys_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vsys, hist_std_vsys;
    // vsys error
    double *vsys_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vsys_err, hist_std_vsys_err;

    // pa
    double *pa_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *pa_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_pa, hist_std_pa;
    // pa error
    double *pa_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_pa_err, hist_std_pa_err;

    // incl
    double *incl_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *incl_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_incl, hist_std_incl;
    // incl error
    double *incl_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_incl_err, hist_std_incl_err;

    // vrot
    double *vrot_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrot_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vrot, hist_std_vrot;
    // vrot error
    double *vrot_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vrot_err, hist_std_vrot_err;

    // vrad
    double *vrad_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrad_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vrad, hist_std_vrad;
    // vrad error
    double *vrad_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double hist_mean_vrad_err, hist_std_vrad_err;

    /* Set up ring parameters with currently derived values*/
    set_nfree_params_trfit_multinest_ellipsefit_rings_student(xposfix, yposfix, vsysfix, pafix, inclfix, vrotfix, vradfix, sigmafactorfix, TRparam);
    ndims = TRparam[0].n_freeParams;
    nPar = TRparam[0].n_freeParams;
    nClsPar = TRparam[0].n_freeParams;
    int pWrap[ndims];
    for(i = 0; i < ndims; i++) pWrap[i] = multinest_param[0].pWrap[i];


    // the number of rings for iso 1D fit (divided by 2 to increase the fit reliability)
    TRparam[0].ring_s_for_einasto1Dfit = (int)TRparam[0].ring_s;
    TRparam[0].ring_w_for_einasto1Dfit = (int)TRparam[0].ring_w;
    TRparam[0].ring_e_for_einasto1Dfit = (int)TRparam[0].ring_e;
    TRparam[0].Nrings_to_semi_mx = (int)(1.0*(TRparam[0].ellipse_semi_mx_boxfiltered-TRparam[0].ring_s_for_einasto1Dfit)/TRparam[0].ring_w_for_einasto1Dfit);

    if(TRparam[0].Nrings < 2)
    {
        printf("++ The number of rings out to the semi major axis derived from ellipse fit is less than 2. Adjust ring width etc. ++\n");

        exit(0);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /* Start tilted-ring fits */
    for(i=0; i<(int)TRparam[0].Nrings; i++) // this is to avoid the outer rotation velocities with large errors
    {
        ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
        ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
        ring = (ri+ro)/2.0;

        // To use all pixels inside the innermost radius for the first ring
        if (ri < 0.0) ri = 0.0;

        xpos_iso1Dfit = TRparam[0].xposF_EinastoFit;
        ypos_iso1Dfit = TRparam[0].yposF_EinastoFit;
        pa_iso1Dfit = TRparam[0].paF_EinastoFit;
        incl_iso1Dfit = TRparam[0].inclF_EinastoFit;

        nlive = multinest_param[0].nlive_einasto_halofit;
        if(nlive > 50) nlive = 50; // for a quick TRfit
        define_tiltedRing_ellipse_rings(TRparam[0].xpos0[i], TRparam[0].ypos0[i], TRparam[0].pa0[i], TRparam[0].incl0[i], ri, ro, TRparam, side); 

        // update xpos & ypos priors based on the current ring radius + pa0 + incl0
        ellipse_a = ri;
        xpos_prior_gab = fabs(ellipse_a*sin(TRparam[0].pa0[i]*M_PI/180.));
        ypos_prior_gab = fabs(ellipse_a*cos(TRparam[0].pa0[i]*M_PI/180.));

        if(xposfix == 'T')
        {
            TRparam[0].xpos1 = TRparam[0].xpos0[i] - xpos_prior_gab;
            TRparam[0].xpos2 = TRparam[0].xpos0[i] + xpos_prior_gab;
        }
        if(yposfix == 'T')
        {
            TRparam[0].ypos1 = TRparam[0].ypos0[i] - ypos_prior_gab;
            TRparam[0].ypos2 = TRparam[0].ypos0[i] + ypos_prior_gab;
        }

        // use the current PA model value for multinest if not fitted : in degree
        TRparam[0].paF = TRparam[0].paF_EinastoFit;

        // use the current INCL model value for multinest if not fitted : in degree
        TRparam[0].inclF = TRparam[0].inclF_EinastoFit;

        /* _ Calling multines_ */   
        MPI_Barrier(MPI_COMM_WORLD);
        run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_trfit_student, dumper_TRfits, TRparam);   
        MPI_Barrier(MPI_COMM_WORLD);

        TRparam[0].xpos[i] = TRparam[0].xposF; 
        TRparam[0].xpos0[i] = TRparam[0].xposF;
        TRparam[0].xpos_e[i] = TRparam[0].xposF_e;
        xpos_tr[i] = TRparam[0].xpos[i];
        xpos_err[i] = TRparam[0].xpos_e[i];

        TRparam[0].ypos[i] = TRparam[0].yposF;
        TRparam[0].ypos0[i] = TRparam[0].yposF;
        TRparam[0].ypos_e[i] = TRparam[0].yposF_e;
        ypos_tr[i] = TRparam[0].ypos[i];
        ypos_err[i] = TRparam[0].ypos_e[i];
        
        TRparam[0].vsys[i] = TRparam[0].vsysF;
        TRparam[0].vsys0[i] = TRparam[0].vsysF;
        TRparam[0].vsys_e[i] = TRparam[0].vsysF_e;
        vsys_tr[i] = TRparam[0].vsys[i];
        vsys_err[i] = TRparam[0].vsys_e[i];

        // in degree : for
        TRparam[0].pa[i] = TRparam[0].paF;
        TRparam[0].pa_temp[i] = TRparam[0].paF;
        TRparam[0].pa0[i] = TRparam[0].paF;
        TRparam[0].pa_e[i] = TRparam[0].paF_e;
        pa_tr[i] = TRparam[0].pa[i];
        pa_err[i] = TRparam[0].pa_e[i];

        // in degree
        TRparam[0].incl[i] = TRparam[0].inclF;
        TRparam[0].incl_temp[i] = TRparam[0].inclF;
        TRparam[0].incl0[i] = TRparam[0].inclF;
        TRparam[0].incl_e[i] = TRparam[0].inclF_e;
        incl_tr[i] = TRparam[0].incl[i];
        incl_err[i] = TRparam[0].incl_e[i];

        // VROT
        if(side == 0)
        {   
            TRparam[0].vrot[i] = TRparam[0].vrotF;
            TRparam[0].vrot0[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;
        }
        else if(side == -1)
        {   
            TRparam[0].vrot_rec[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e_rec[i] = TRparam[0].vrotF_e;
        }
        else if(side == 1)
        {   
            TRparam[0].vrot_app[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e_app[i] = TRparam[0].vrotF_e;
        }
        vrot_tr[i] = TRparam[0].vrotF;
        vrot_err[i] = TRparam[0].vrotF_e;

        // VRAD
        if(side == 0)
        {   
            TRparam[0].vrad[i] = TRparam[0].vradF;
            TRparam[0].vrad_temp[i] = TRparam[0].vradF;
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
        else if(side == -1)
        {   
            TRparam[0].vrad_rec[i] = TRparam[0].vradF;
            TRparam[0].vrad_rec_e[i] = TRparam[0].vradF_e;
        }
        else if(side == 1)
        {   
            TRparam[0].vrad_app[i] = TRparam[0].vradF;
            TRparam[0].vrad_app_e[i] = TRparam[0].vradF_e;
        }
        vrad_tr[i] = TRparam[0].vradF;
        vrad_err[i] = TRparam[0].vradF_e;

        // e_sigma
        TRparam[0].e_sigma_student_TR[i] = TRparam[0].e_sigma_tr; // student-t e_sigma  


        TRparam[0].maxLogLike[i] = TRparam[0].maxLogLikeF;
        TRparam[0].logZ[i] = TRparam[0].logZF;
        TRparam[0].logZerr[i] = TRparam[0].logZerrF;
        TRparam[0].npoints_inaring[i] = TRparam[0].Npoints_in_tilted_ring;
        TRparam[0].npoints_inaring_decim0[i] = TRparam[0].Npoints_in_tilted_ring_decim0;
        total_Npoints_allrings += TRparam[0].Npoints_in_tilted_ring;

        if(i > 0) // if the number of pixels in a ring is less then 10% of the previous one, exit!
        {   
            if(TRparam[0].npoints_inaring_decim0[i] < TRparam[0].npoints_inaring_decim0[i-1]*0.3)
            {
                ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
                ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
                ring = (ri+ro)/2.0;
                TRparam[0].ellipse_semi_mx_boxfiltered = (ri+ro)/2.0;
                TRparam[0].ring_e = (ri+ro)/2.0;
                //TRparam[0].Nrings_to_semi_mx = (int)(1.0*(TRparam[0].ellipse_semi_mx_boxfiltered-TRparam[0].ring_s)/TRparam[0].ring_w)+1;
                TRparam[0].Nrings = i;
                break;
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }

    // +++ Derive weights of individual parameters based on their normalised errors
    //FD_nbins = (int)(sqrt(TRparam[0].Nrings)+0.5);
    // XPOS
    if(xposfix == 'T')
    {
        // 0. gsl statistics
        gsl_mean_xpos_e = gsl_stats_mean(xpos_err, 1, TRparam[0].Nrings);
        gsl_std_xpos_e = sqrt(gsl_stats_variance(xpos_err, 1, TRparam[0].Nrings));
        gsl_min_xpos_e = gsl_stats_min(xpos_err, 1, TRparam[0].Nrings);
        gsl_max_xpos_e =  gsl_stats_max(xpos_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_xpos_e) || isinf(gsl_min_xpos_e) || gsl_min_xpos_e <= 0)
        if(isnan(gsl_min_xpos_e) || isinf(gsl_min_xpos_e) || gsl_min_xpos_e <= 0)
            gsl_min_xpos = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // XPOS
        robust_mean_std(xpos_tr, TRparam[0].Nrings, &hist_mean_xpos, &hist_std_xpos);
        // XPOS error
        robust_mean_std_e(xpos_err, TRparam[0].Nrings, &hist_mean_xpos_err, &hist_std_xpos_err);
    }

    // YPOS
    if(yposfix == 'T')
    {
        gsl_mean_ypos_e = gsl_stats_mean(ypos_err, 1, TRparam[0].Nrings);
        gsl_std_ypos_e = sqrt(gsl_stats_variance(ypos_err, 1, TRparam[0].Nrings));
        gsl_min_ypos_e = gsl_stats_min(ypos_err, 1, TRparam[0].Nrings);
        gsl_max_ypos_e =  gsl_stats_max(ypos_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_ypos_e) || isinf(gsl_min_ypos_e) || gsl_min_ypos_e <= 0)
            gsl_min_ypos_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // YPOS
        robust_mean_std(ypos_tr, TRparam[0].Nrings, &hist_mean_ypos, &hist_std_ypos);
        // YPOS error
        robust_mean_std_e(ypos_err, TRparam[0].Nrings, &hist_mean_ypos_err, &hist_std_ypos_err);
    }

    // VSYS
    if(vsysfix == 'T')
    {
        gsl_mean_vsys_e = gsl_stats_mean(vsys_err, 1, TRparam[0].Nrings);
        gsl_std_vsys_e = sqrt(gsl_stats_variance(vsys_err, 1, TRparam[0].Nrings));
        gsl_min_vsys_e = gsl_stats_min(vsys_err, 1, TRparam[0].Nrings);
        gsl_max_vsys_e =  gsl_stats_max(vsys_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vsys_e) || isinf(gsl_min_vsys_e) || gsl_min_vsys_e <= 0)
            gsl_min_vsys_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // VSYS
        robust_mean_std(vsys_tr, TRparam[0].Nrings, &hist_mean_vsys, &hist_std_vsys);
        // VSYS error
        robust_mean_std_e(vsys_err, TRparam[0].Nrings, &hist_mean_vsys_err, &hist_std_vsys_err);
    }
    
    // PA
    if(pafix == 'T')
    {
        gsl_mean_pa_e = gsl_stats_mean(pa_err, 1, TRparam[0].Nrings);
        gsl_std_pa_e = sqrt(gsl_stats_variance(pa_err, 1, TRparam[0].Nrings));
        gsl_min_pa_e = gsl_stats_min(pa_err, 1, TRparam[0].Nrings);
        gsl_max_pa_e = gsl_stats_max(pa_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_pa_e) || isinf(gsl_min_pa_e) || gsl_min_pa_e <= 0)
            gsl_min_pa_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // PA
        robust_mean_std(pa_tr, TRparam[0].Nrings, &hist_mean_pa, &hist_std_pa);
        // PA error
        robust_mean_std_e(pa_err, TRparam[0].Nrings, &hist_mean_pa_err, &hist_std_pa_err);
    }

    // INCL
    if(inclfix == 'T')
    {
        gsl_mean_incl_e = gsl_stats_mean(incl_err, 1, TRparam[0].Nrings);
        gsl_std_incl_e = sqrt(gsl_stats_variance(incl_err, 1, TRparam[0].Nrings));
        gsl_min_incl_e = gsl_stats_min(incl_err, 1, TRparam[0].Nrings);
        gsl_max_incl_e = gsl_stats_max(incl_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_incl_e) || isinf(gsl_min_incl_e) || gsl_min_incl_e <= 0)
            gsl_min_incl_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // INCL
        robust_mean_std(incl_tr, TRparam[0].Nrings, &hist_mean_incl, &hist_std_incl);
        // INCL error
        robust_mean_std_e(incl_err, TRparam[0].Nrings, &hist_mean_incl_err, &hist_std_incl_err);
    }

    // VROT
    if(vrotfix == 'T')
    {
        gsl_mean_vrot_e = gsl_stats_mean(vrot_err, 1, TRparam[0].Nrings);
        gsl_std_vrot_e = sqrt(gsl_stats_variance(vrot_err, 1, TRparam[0].Nrings));
        gsl_min_vrot_e = gsl_stats_min(vrot_err, 1, TRparam[0].Nrings);
        gsl_max_vrot_e = gsl_stats_max(vrot_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrot_e) || isinf(gsl_min_vrot_e) || gsl_min_vrot_e <= 0)
            gsl_min_vrot_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // VROT
        robust_mean_std(vrot_tr, TRparam[0].Nrings, &hist_mean_vrot, &hist_std_vrot);
        // VROT error
        robust_mean_std_e(vrot_err, TRparam[0].Nrings, &hist_mean_vrot_err, &hist_std_vrot_err);
    }

    // VRAD
    if(vradfix == 'T')
    {
        gsl_mean_vrad_e = gsl_stats_mean(vrad_err, 1, TRparam[0].Nrings);
        gsl_std_vrad_e = sqrt(gsl_stats_variance(vrad_err, 1, TRparam[0].Nrings));
        gsl_min_vrad_e = gsl_stats_min(vrad_err, 1, TRparam[0].Nrings);
        gsl_max_vrad_e = gsl_stats_max(vrad_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrad_e) || isinf(gsl_min_vrad_e) || gsl_min_vrad_e <= 0)
            gsl_min_vrad_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // VRAD
        robust_mean_std(vrad_tr, TRparam[0].Nrings, &hist_mean_vrad, &hist_std_vrad);
        // VRAD error
        robust_mean_std_e(vrad_err, TRparam[0].Nrings, &hist_mean_vrad_err, &hist_std_vrad_err);
    }

    // 3. Normalise the errors with respect to their minimums & derive the total errors
    // normalise errors with their mimimums to derive the total error normalised
    for(i=0; i<(int)TRparam[0].Nrings; i++) // this is to avoid the outer rotation velocities with large errors
    {
        // xpos
        if(xposfix == 'T')
        {
            if(TRparam[0].xpos_e[i] > 0)
                xpos_err[i] = TRparam[0].xpos_e[i]/gsl_min_xpos_e;
            else
                xpos_err[i] = 1E3; // large value
        }
        else
            xpos_err[i] = 0; // no error


        // ypos
        if(yposfix == 'T')
        {
            if(TRparam[0].ypos_e[i] > 0)
                ypos_err[i] = TRparam[0].ypos_e[i]/gsl_min_ypos_e;
            else
                ypos_err[i] = 1E3; // large value
        }
        else
            ypos_err[i] = 0; // no error


        // vsys
        if(vsysfix == 'T')
        {
            if(TRparam[0].vsys_e[i] > 0)
                vsys_err[i] = TRparam[0].vsys_e[i]/gsl_min_vsys_e;
            else
                vsys_err[i] = 1E3; // large value
        }
        else
            vsys_err[i] = 0; // no error

        // pa
        if(pafix == 'T')
        {
            if(TRparam[0].pa_e[i] > 0)
                pa_err[i] = TRparam[0].pa_e[i]/gsl_min_pa_e;
            else
                pa_err[i] = 1E3; // large value
        }
        else
            pa_err[i] = 0; // no error

        // incl
        if(inclfix == 'T')
        {
            if(TRparam[0].incl_e[i] > 0)
                incl_err[i] = TRparam[0].incl_e[i]/gsl_min_incl_e;
            else
                incl_err[i] = 1E3; // large value
        }
        else
            incl_err[i] = 0; // no error

        if(vrotfix == 'T')
        {
            if(TRparam[0].vrot_e[i] > 0)
                vrot_err[i] = TRparam[0].vrot_e[i]/gsl_min_vrot_e;
            else
                vrot_err[i] = 1E3;
        }
        else
            vrot_err[i] = 0; // no error

        if(vradfix == 'T')
        {
            if(TRparam[0].vrad_e[i] > 0)
                vrad_err[i] = TRparam[0].vrad_e[i]/gsl_min_vrad_e;
            else
                vrad_err[i] = 1E3;
        }
        else
            vrad_err[i] = 0; // no error


        // total errors (PA_e+INCL_e+VROT_e+VRAD_e): normalised ones to their minimums
        total_error_ring_params[i] = sqrt(pow(xpos_err[i], 2) + \
                                            pow(ypos_err[i], 2) + \
                                            pow(vsys_err[i], 2) + \
                                            pow(pa_err[i], 2) + \
                                            pow(incl_err[i], 2) + \
                                            pow(vrot_err[i], 2) + \
                                            pow(vrad_err[i], 2));

        if(isinf(total_error_ring_params[i]) || isnan(total_error_ring_params[i]) || total_error_ring_params[i] <= 0)
            total_error_ring_params[i] = 1E3; // large value

        // this is for histogram analysis
        total_error_ring_params_temp[i] = total_error_ring_params[i];
    }

    // 4. Calculate mean & std of total error histograms to match with those of ring params
    // mean & std of total error histogram 
    robust_mean_std_e(total_error_ring_params_temp, TRparam[0].Nrings, &hist_mean_total_error_ring_params, &hist_std_total_error_ring_params);

    // Re-assign histogram with scaled total errors: see above
    for(i=0; i<(int)TRparam[0].Nrings; i++) // calculate error weighted mode of the ring parameters using their histograms
    {
        // xpos
        if(xposfix == 'T')
        {
            TRparam[0].xpos0[i] = hist_mean_xpos;
            TRparam[0].xposF_EinastoFit = hist_mean_xpos; // this is for define_tiltedrings() in TRfits_multinest_using_ISOfit_rings which is run the very after this step.
        
            xpos_tr_err[i] = total_error_ring_params[i]*(hist_mean_xpos_err/hist_mean_total_error_ring_params);
        }

        // ypos
        if(yposfix == 'T')
        {
            TRparam[0].ypos0[i] = hist_mean_ypos;
            TRparam[0].yposF_EinastoFit = hist_mean_ypos; // this is for define_tiltedrings() in TRfits_multinest_using_ISOfit_rings which is run the very after this step.
            ypos_tr_err[i] = total_error_ring_params[i]*(hist_mean_ypos_err/hist_mean_total_error_ring_params);
        }

        // vsys
        if(vsysfix == 'T')
        {
            TRparam[0].vsys0[i] = hist_mean_vsys;
            TRparam[0].vsysF_EinastoFit = hist_mean_vsys; // this is for define_tiltedrings() in TRfits_multinest_using_ISOfit_rings which is run the very after this step.
            vsys_tr_err[i] = total_error_ring_params[i]*(hist_mean_vsys_err/hist_mean_total_error_ring_params);
        }

        // pa
        if(pafix == 'T')
        {
            TRparam[0].pa0[i] = hist_mean_pa;
            TRparam[0].paF_EinastoFit = hist_mean_pa; // this is for define_tiltedrings() in TRfits_multinest_using_ISOfit_rings which is run the very after this step.

            pa_tr_err[i] = total_error_ring_params[i]*(hist_mean_pa_err/hist_mean_total_error_ring_params);
        }

        // incl
        if(inclfix == 'T')
        {
            TRparam[0].incl0[i] = hist_mean_incl;
            TRparam[0].inclF_EinastoFit = hist_mean_incl; // this is for define_tiltedrings() in TRfits_multinest_using_ISOfit_rings which is run the very after this step.
            incl_tr_err[i] = total_error_ring_params[i]*(hist_mean_incl_err/hist_mean_total_error_ring_params);
        }

        // vrot
        if(vrotfix == 'T')
        {
            TRparam[0].vrot0[i] = hist_mean_vrot;
            vrot_tr_err[i] = total_error_ring_params[i]*(hist_mean_vrot_err/hist_mean_total_error_ring_params);
        }

        // vrad
        //if(vradfix == 'T')
        //{
        //    TRparam[0].vrad[i] = hist_mean_vrad;
        //    vrad_tr_err[i] = total_error_ring_params[i]*(hist_mean_vrad_err/hist_mean_total_error_ring_params);
        //}
    }

    // update priors based on the derived mode and sigma of individual ring params
    if(xposfix == 'T')
    {
        //robust_mean_std_histogram_ac(xpos_tr, xpos_tr_err, TRparam[0].Nrings, &hist_mean_xpos, &hist_std_xpos);
        robust_mean_std(xpos_tr, TRparam[0].Nrings, &hist_mean_xpos, &hist_std_xpos);
        TRparam[0].xpos1 = hist_mean_xpos - 5*hist_std_xpos;
        TRparam[0].xpos2 = hist_mean_xpos + 5*hist_std_xpos;
    }
    if(yposfix == 'T')
    {
        //robust_mean_std_histogram_ac(ypos_tr, ypos_tr_err, TRparam[0].Nrings, &hist_mean_ypos, &hist_std_ypos);
        robust_mean_std(ypos_tr, TRparam[0].Nrings, &hist_mean_ypos, &hist_std_ypos);
        TRparam[0].ypos1 = hist_mean_ypos - 5*hist_std_ypos;
        TRparam[0].ypos2 = hist_mean_ypos + 5*hist_std_ypos;
    }
    if(vsysfix == 'T')
    {
        //robust_mean_std_histogram_ac(vsys_tr, vsys_tr_err, TRparam[0].Nrings, &hist_mean_vsys, &hist_std_vsys);
        robust_mean_std(vsys_tr, TRparam[0].Nrings, &hist_mean_vsys, &hist_std_vsys);
        TRparam[0].vsys1 = hist_mean_vsys - 5*hist_std_vsys;
        TRparam[0].vsys2 = hist_mean_vsys + 5*hist_std_vsys;
    }
    if(pafix == 'T')
    {
        //robust_mean_std_histogram_ac(pa_tr, pa_tr_err, TRparam[0].Nrings, &hist_mean_pa, &hist_std_pa);
        robust_mean_std(pa_tr, TRparam[0].Nrings, &hist_mean_pa, &hist_std_pa);
        TRparam[0]._bspline_pa_hist_sigma = hist_std_pa;
    }
    if(inclfix == 'T')
    {
        //robust_mean_std_histogram_ac(incl_tr, incl_tr_err, TRparam[0].Nrings, &hist_mean_incl, &hist_std_incl);
        robust_mean_std(incl_tr, TRparam[0].Nrings, &hist_mean_incl, &hist_std_incl);
        TRparam[0]._bspline_incl_hist_sigma = hist_std_incl;
    }
    if(vradfix == 'T')
    {
        robust_mean_std(vrad_tr, TRparam[0].Nrings, &hist_mean_vrad, &hist_std_vrad);
        TRparam[0]._bspline_vrad_hist_sigma = hist_std_vrad;
    }

    /* Save the total number of points in all rings */
    TRparam[0].total_Npoints_allRings = total_Npoints_allrings;

    /* calculate median value for XPOS if fitted */
    if(TRparam[0].xpos_fix == 'T')
    {
        TRparam[0].xposF = hist_mean_xpos; // error weighted mean of xpos histogram
        TRparam[0].xposF_EinastoFit = hist_mean_xpos; // update xposF_EinastoFit initially from ellipse fit
        TRparam[0].xposF_e = hist_std_xpos; // error weighted sigma of xpos histogram
    }

    /* calculate median value for YPOS if fitted */
    if(TRparam[0].ypos_fix == 'T')
    {
        TRparam[0].yposF = hist_mean_ypos; // error weighted mean of ypos histogram
        TRparam[0].yposF_EinastoFit = hist_mean_ypos; // update yposF_EinastoFit initially from ellipse fit
        TRparam[0].yposF_e = hist_std_ypos; // error weighted sigma of ypos histogram
    }

    /* calculate median value for VSYS if fitted */
    if(TRparam[0].vsys_fix == 'T')
    {
        TRparam[0].vsysF = hist_mean_vsys; // error weighted mean of vsys histogram
        TRparam[0].vsysF_EinastoFit = hist_mean_vsys; // update vsysF_EinastoFit initially from ellipse fit
        TRparam[0].vsysF_e = hist_std_vsys; // error weighted sigma of vsys histogram
    }

    /* calculate median value for PA if fitted */
    if(TRparam[0].pa_fix == 'T')
    {
        TRparam[0].paF = hist_mean_pa; // error weighted mean of pa histogram
        TRparam[0].paF_EinastoFit = hist_mean_pa; // update paF_EinastoFit initially from ellipse fit
        TRparam[0].paF_e = hist_std_pa; // error weighted sigma of pa histogram
    }
    /* calculate median value for PA if fitted */
    if(TRparam[0].incl_fix == 'T')
    {
        TRparam[0].inclF = hist_mean_incl; // error weighted mean of incl histogram
        TRparam[0].inclF_EinastoFit = hist_mean_incl; // update inclF_EinastoFit initially from ellipse fit
        TRparam[0].inclF_e = hist_std_incl; // error weighted sigma of incl histogram
    }
    /* calculate median value for VRAD if fitted */
    if(TRparam[0].vrad_fix == 'T')
    {
        TRparam[0].vradF = hist_mean_vrad; // error weighted mean of vrad histogram
        TRparam[0].vradF_e = hist_std_vrad; // error weighted sigma of vrad histogram
    }
    /* calculate median value for VROT if fitted */
    if(TRparam[0].vrot_fix == 'T')
    {
        vrot_median = hist_mean_vrot;
    } 

    if(find_symmetric_VF == 1)
    {
        TRparam[0].xpos_fix = 'F';
        TRparam[0].ypos_fix = 'F';
        TRparam[0].vsys_fix = 'F';
        TRparam[0].pa_fix = 'F';
        TRparam[0].incl_fix = 'F';
        TRparam[0].vrot_fix = vrotfix;

        // the number of tilted-rings to fit
        //TRparam[0].Nrings = (int)((TRparam[0].ring_e-TRparam[0].ring_s)/TRparam[0].ring_w);

        // the number of rings for iso 1D fit (divided by 2 to increase the fit reliability)
        TRparam[0].ring_s_for_einasto1Dfit = (int)TRparam[0].ring_s;
        TRparam[0].ring_w_for_einasto1Dfit = (int)TRparam[0].ring_w;
        TRparam[0].ring_e_for_einasto1Dfit = (int)TRparam[0].ring_e;
        TRparam[0].Nrings_to_semi_mx = (int)(1.0*(TRparam[0].ellipse_semi_mx_boxfiltered-TRparam[0].ring_s_for_einasto1Dfit)/TRparam[0].ring_w_for_einasto1Dfit)+1;

        /* Start tilted-ring fits */
        for(i=0; i<(int)TRparam[0].Nrings_to_semi_mx; i++) // this is to avoid the outer rotation velocities with large errors
        {
            ri = TRparam[0].ring_s_for_einasto1Dfit + i*TRparam[0].ring_w_for_einasto1Dfit - 0.5*TRparam[0].ring_w_for_einasto1Dfit;
            ro = TRparam[0].ring_s_for_einasto1Dfit + i*TRparam[0].ring_w_for_einasto1Dfit + 0.5*TRparam[0].ring_w_for_einasto1Dfit;

            // to use all pixels inside the innermost radius for the first ring
            if (ri < 0.0) ri = 0.0;
            // These initial ring parameters are from ellipe fit to the largest connected area found: pa and incl are given in degree

            // fixed xpos ypos pa incl to find the optimal symmetric VF!
            xpos_iso1Dfit = TRparam[0].xposF; // from median value of TR fits
            ypos_iso1Dfit = TRparam[0].yposF; // from median value of TR fits
            pa_iso1Dfit = TRparam[0].paF; // ellipse fit PA
            incl_iso1Dfit = TRparam[0].inclF; // ellipse fit INCL

            TRparam[0].xposF_EinastoFit = TRparam[0].xposF; // from median value of TR fits
            TRparam[0].yposF_EinastoFit = TRparam[0].yposF; // from median value of TR fits
            TRparam[0].vsysF_EinastoFit = TRparam[0].vsysF; // from median value of TR fits
            TRparam[0].paF_EinastoFit = TRparam[0].paF; // from median value of TR fits
            TRparam[0].inclF_EinastoFit = TRparam[0].inclF; // from median value of TR fits

            nlive = multinest_param[0].nlive_einasto_halofit;
            if(nlive > 50) nlive = 50; // for a quick TRfit
            define_tiltedRing_ellipse_rings(xpos_iso1Dfit, ypos_iso1Dfit, pa_iso1Dfit, incl_iso1Dfit, ri, ro, TRparam, side);       
            TRparam[0].npoints_inaring[i] = TRparam[0].Npoints_in_tilted_ring;
            TRparam[0].npoints_inaring_decim0[i] = TRparam[0].Npoints_in_tilted_ring_decim0;

            total_Npoints_allrings += TRparam[0].Npoints_in_tilted_ring;
            //printf("%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", (ri+ro)/2.0, TRparam[0].xposF, TRparam[0].xposF_e, TRparam[0].yposF, TRparam[0].yposF_e, TRparam[0].vsysF, TRparam[0].vsysF_e, TRparam[0].paF, TRparam[0].paF_e, TRparam[0].inclF, TRparam[0].inclF_e, TRparam[0].vrotF, TRparam[0].vrotF_e, TRparam[0].npoints_inaring_decim0[i], TRparam[0].N_all_pixels_in_a_ring);

            if(i > 0) // if the number of pixels in a ring is less then 20% of the previous one, exit!
            {   
                if(TRparam[0].npoints_inaring_decim0[i] < TRparam[0].npoints_inaring_decim0[i-1]*0.5)
                {
                    ri = TRparam[0].ring_s_for_einasto1Dfit + (i-1)*TRparam[0].ring_w_for_einasto1Dfit - 0.5*TRparam[0].ring_w_for_einasto1Dfit;
                    ro = TRparam[0].ring_s_for_einasto1Dfit + (i-1)*TRparam[0].ring_w_for_einasto1Dfit + 0.5*TRparam[0].ring_w_for_einasto1Dfit;
                    TRparam[0].ellipse_semi_mx_boxfiltered = (ri+ro)/2.0;
                    TRparam[0].ring_e = (ri+ro)/2.0;

                    TRparam[0].ring_s_for_einasto1Dfit = (int)TRparam[0].ring_s/1.0;
                    TRparam[0].ring_w_for_einasto1Dfit = (int)TRparam[0].ring_w/1.0;
                    TRparam[0].ring_e_for_einasto1Dfit = (int)TRparam[0].ring_e;
                    TRparam[0].Nrings_to_semi_mx = (int)(1.0*(TRparam[0].ellipse_semi_mx_boxfiltered-TRparam[0].ring_s_for_einasto1Dfit)/TRparam[0].ring_w_for_einasto1Dfit)+1;
                    break;
                }
            }
        }

        peri_x = -1.0*TRparam[0].ellipse_semi_mx_boxfiltered*sin(TRparam[0].paF*M_PI/180.) + TRparam[0].xposF;
        peri_y = TRparam[0].ellipse_semi_mx_boxfiltered*cos(TRparam[0].paF*M_PI/180.) + TRparam[0].yposF; 

        k=0;
        for(i=0; i<(int)(TRparam[0].xpos2-TRparam[0].xpos1); i++)
        {
            for(j=0; j<(int)(TRparam[0].ypos2-TRparam[0].ypos1); j++)
            {
                xpos_temp = TRparam[0].xpos1 + i;
                ypos_temp = TRparam[0].ypos1 + j;
                rmax_ellipse[k] = sqrt(pow((peri_x-xpos_temp), 2) + pow((peri_y-ypos_temp), 2));
                k++;
            }
        }
        qsort(rmax_ellipse, k, sizeof(rmax_ellipse[0]), comparisonFunctionFloat);
        free(rmax_ellipse);
    }
    else
    {
        TRparam[0].xposF_EinastoFit = TRparam[0].xposF; // from median value of TR fits
        TRparam[0].yposF_EinastoFit = TRparam[0].yposF; // from median value of TR fits
        TRparam[0].vsysF_EinastoFit = TRparam[0].vsysF; // from median value of TR fits
        TRparam[0].paF_EinastoFit = TRparam[0].paF; // from median value of TR fits
        TRparam[0].inclF_EinastoFit = TRparam[0].inclF; // from median value of TR fits
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
void set_nfree_params_trfit_multinest_ellipsefit_rings_student(char xposfix, char yposfix, char vsysfix, char pafix, char inclfix, char vrotfix, char vradfix, char sigmafactorfix, TR_ringParameters *TRparam)
{
    int i=0;
    int n_freeParams = 0;
    double ring=0.;

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
        n_freeParams++;
    }   
    // 5. INCL
    TRparam[0].incl_fix = inclfix;
    if(TRparam[0].incl_fix == 'T')
    {
        n_freeParams++;
    }
    // 6. VROT
    TRparam[0].vrot_fix = vrotfix;
    if(TRparam[0].vrot_fix == 'T')
    {
        n_freeParams++;
    }
    // 7. VRAD
    TRparam[0].vrad_fix = vradfix;
    if(TRparam[0].vrad_fix == 'T')
    {
        n_freeParams++;
    }
    // 8. sigma_factor
    TRparam[0].sigma_factor_fix = sigmafactorfix;
    //if(TRparam[0].sigma_factor_fix == 'T' || TRparam[0].sigma_factor_fix == 'F')
    if(TRparam[0].e_sigma >= 0)
    {
        n_freeParams++;
    }

    // Update the total number of free tilted-ring parameters for TR fits
    TRparam[0].n_freeParams = n_freeParams;

    // Update TR ring parameters for each ring from either einasto or ellipse fit results
    for(i=0; i<TRparam[0].Nrings; i++)
    {
        ring = TRparam[0].ring_s + (double)i*TRparam[0].ring_w;
        TRparam[0].ring_radius[i] = ring;

        // 1. initial XPOS
        TRparam[0].xpos0[i] = TRparam[0].xposF_EinastoFit;

        // 2. intial YPOS
        TRparam[0].ypos0[i] = TRparam[0].yposF_EinastoFit;

        // 3. intial VSYS
        TRparam[0].vsys0[i] = TRparam[0].vsysF_EinastoFit;

        // 4. initial PA in degree (input paF_EinastoFit is given in degree)
        TRparam[0].pa0[i] = TRparam[0].paF_EinastoFit; // this is to avoid a simple use of ellipse fit without receding or approaching information.
        // 5. intial INCL in degree (input incl_ISO is given in degree)
        TRparam[0].incl0[i] = TRparam[0].inclF_EinastoFit;
    }
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void define_tiltedRing(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side)
{
    double i=0, j=0;
    int i0=0, j0=0, k0=0, k=0;
    int decimX, decimY;
    int Npoints_in_tiltedRing_without_decimals;
    int x0, y0, x1, y1, box_x, box_y;
    double theta=0., r=0., r_inner=0., r_outter=0.;
    double Npoints_in_tiltedRing=0;
    double Npoints_in_tiltedRing_total_including_blanks=0;
    double ring0_regrad, ring1_regrad;
    double sine_free_angle, costh;
    double mean_ve_temp, std_ve_temp;
    double mean_vf_e_temp_w, std_vf_e_temp_w;
    double mean_mom2_temp, std_mom2_temp;
    double mean_mom4_temp, std_mom4_temp;
    double HI_VF_weight_TRfit_max;
    double r_i0j0;
    double pa_intp, incl_intp;
    double mean_vf_e, std_vf_e;
    Ellipse_Parameter ellipse;
    ring_parameter ring;

    ellipse.xpos = xpos;
    ellipse.ypos = ypos;
    ellipse.a = ring0 + (ring1-ring0)/2.0;
    ellipse.b = ellipse.a * cos(incl*M_PI/180.);
    ellipse.e = sqrt(1-(ellipse.b*ellipse.b)/(ellipse.a*ellipse.a));

    // 0. set free angle
    sine_free_angle = fabs(sin(TRparam[0].free_angle*M_PI/180.));

    // 1. extract the region to fit with decimals starting from the centre position given
    Npoints_in_tiltedRing_total_including_blanks = 0;
    decimX = 0;
    decimY = 0;

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0;
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;

    gsl_interp_accel *acc;
    gsl_spline *spline;

    double *weight_temp = malloc(sizeof(double) * (x1-x0)*(y1-y0));
    // count the total pixels in a given ring
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            r_i0j0 = r_ij_pa_incl_given_W(i0, j0, xpos, ypos, pa, incl, TRparam);
            acc = gsl_interp_accel_alloc();
            spline = gsl_spline_alloc (gsl_interp_cspline, TRparam[0].Nrings);
            gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl, TRparam[0].Nrings);

            // interpolate incl at r_i0j0
            if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
            {
                incl_intp = TRparam[0].incl[TRparam[0].Nrings-1]; // outermost incl
            }
            else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
            {
                incl_intp = TRparam[0].incl[0]; // outermost incl
            }
            else
            {
                incl_intp = gsl_spline_eval (spline, r_i0j0, acc);
            }

            // interpolate pa at r_i0j0
            gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].pa, TRparam[0].Nrings);
            if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
            {
                pa_intp = TRparam[0].pa[TRparam[0].Nrings-1]; // outermost pa
            }
            else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
            {
                pa_intp = TRparam[0].pa[0]; // innermost pa
            }
            else
            {
                pa_intp = gsl_spline_eval (spline, r_i0j0, acc);
            }
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);


            i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa_intp*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa_intp*M_PI/180.));
            j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa_intp*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa_intp*M_PI/180.))/cos(incl_intp*M_PI/180.);
            r = sqrt(i*i+j*j); // distance from centre

            if(r < 1) // if r is smaller than one pixel
            {
                r = 1;
                theta = 0.0;
            }
            else
                theta = atan2((double)j, (double)i)*180./M_PI;

            costh = fabs(cos(theta*M_PI/180.));

            if(r >= ring0 && r < ring1 && costh > sine_free_angle)
            {
                if(side == -1 && fabs(theta) <= 90.0) // receding side
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
                else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
                else if(side == 0 || side == 999) // both sides
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }
    TRparam[0].N_all_pixels_in_a_ring = Npoints_in_tiltedRing_total_including_blanks;

    // between ring0 and ring1
    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0; 
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;

    Npoints_in_tiltedRing = 0;
    // this is for decimX=0, decimY=0
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            r_i0j0 = r_ij_pa_incl_given_W(i0, j0, xpos, ypos, pa, incl, TRparam);
            acc = gsl_interp_accel_alloc();
            spline = gsl_spline_alloc (gsl_interp_cspline, TRparam[0].Nrings);
            gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl, TRparam[0].Nrings);

            // interpolate incl at r_i0j0
            if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
            {
                incl_intp = TRparam[0].incl[TRparam[0].Nrings-1]; // outermost incl
            }
            else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
            {
                incl_intp = TRparam[0].incl[0]; // outermost incl
            }
            else
            {
                incl_intp = gsl_spline_eval (spline, r_i0j0, acc);
            }

            // interpolate pa at r_i0j0
            gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].pa, TRparam[0].Nrings);
            if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
            {
                pa_intp = TRparam[0].pa[TRparam[0].Nrings-1]; // outermost pa
            }
            else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
            {
                pa_intp = TRparam[0].pa[0]; // outermost pa
            }
            else
            {
                pa_intp = gsl_spline_eval (spline, r_i0j0, acc);
            }
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);


            i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa_intp*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa_intp*M_PI/180.));
            j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa_intp*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa_intp*M_PI/180.))/cos(incl_intp*M_PI/180.);
            r = sqrt(i*i+j*j); // distance from centre

            if(r < 1) // if r is smaller than one pixel
            {
                r = 1;
                theta = 0.0;
            }
            else
                theta = atan2((double)j, (double)i)*180./M_PI;

            costh = fabs(cos(theta*M_PI/180.));

            if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
            {
                if(side == -1 && fabs(theta) <= 90.0) // receding side
                {
                    // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                    // squared in the MLE by sigma*sigma...
                    HI_VF_weight_TRfit[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)); // weight : note that radial weight doesn't need to be applied as all points within a ring have the same radius.
                    Npoints_in_tiltedRing += 1;
                }
                else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                {
                    HI_VF_weight_TRfit[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)); // weight : note that radial weight doesn't need to be applied as all points within a ring have the same radius.
                    Npoints_in_tiltedRing += 1;
                }
                else if(side == 0 || side == 999) // both sides
                {
                    HI_VF_weight_TRfit[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)); // weight : note that radial weight doesn't need to be applied as all points within a ring have the same radius.
                    Npoints_in_tiltedRing += 1;
                }
                
                if(HI_VF_weight_TRfit[0].data[j0][i0] == 0 || isnan(HI_VF_weight_TRfit[0].data[j0][i0]) || isinf(HI_VF_weight_TRfit[0].data[j0][i0]))
                    HI_VF_weight_TRfit[0].data[j0][i0] = 1E-2; // put a small value
            }

            if(r >= ring0-TRparam[0].ring_w && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
            {
                // filling any gab
                if(HI_VF_weight_TRfit[0].data[j0][i0] == 0)
                    HI_VF_weight_TRfit[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow));
            }
        }
    }
    TRparam[0].Npoints_in_tilted_ring_decim0 = Npoints_in_tiltedRing;


    // between ring0 and ring1 : for computing TR WEIGHT : if(r > rmax) r =
    // rmax
    Npoints_in_tiltedRing = 0;
    decimX = (int)TRparam[0].decimX_trfit;
    decimY = (int)TRparam[0].decimY_trfit;

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0; 
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {
                r_i0j0 = r_ij_pa_incl_given_W(i0, j0, xpos, ypos, pa, incl, TRparam);
                acc = gsl_interp_accel_alloc();
                spline = gsl_spline_alloc (gsl_interp_cspline, TRparam[0].Nrings);
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl, TRparam[0].Nrings);

                // interpolate incl at r_i0j0
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    incl_intp = TRparam[0].incl[TRparam[0].Nrings-1]; // outermost incl
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    incl_intp = TRparam[0].incl[0]; // outermost incl
                }
                else
                {
                    incl_intp = gsl_spline_eval (spline, r_i0j0, acc);
                }

                // interpolate pa at r_i0j0
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].pa, TRparam[0].Nrings);
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    pa_intp = TRparam[0].pa[TRparam[0].Nrings-1]; // outermost pa
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    pa_intp = TRparam[0].pa[0]; // outermost pa
                }
                else
                {
                    pa_intp = gsl_spline_eval (spline, r_i0j0, acc);
                }
                gsl_spline_free(spline);
                gsl_interp_accel_free(acc);


                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa_intp*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa_intp*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa_intp*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa_intp*M_PI/180.))/cos(incl_intp*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0] && costh > sine_free_angle)
                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                        // squared in the MLE by sigma*sigma...
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_weight_TRfit[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_weight_TRfit[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_weight_TRfit[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    
                    if(HI_VF_weight_TRfit[0].data[j0][i0] == 0 || isnan(HI_VF_weight_TRfit[0].data[j0][i0]) || isinf(HI_VF_weight_TRfit[0].data[j0][i0]))
                        HI_VF_weight_TRfit[0].data[j0][i0] = 1E-2; // put a small value
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }
    TRparam[0].Npoints_in_tilted_ring = Npoints_in_tiltedRing;


    HI_VF_weight_TRfit_max = gsl_stats_max(weight_temp, 1, Npoints_in_tiltedRing);
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            // normalise (0 ~ 1)
            HI_VF_weight_TRfit[0].data[j0][i0] = HI_VF_weight_TRfit[0].data[j0][i0] / HI_VF_weight_TRfit_max;
        }
    }

    // Derive the geometry (+flux) weighted scale factor for vlos_e
    if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma >= 0) // variable vlos_e mode : partial fit
    {
        if(TRparam[0].final_fit == 'Y')
        {
            double *mom4_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *mom2_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp_w = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);

            // compute mean & std of mom4
            k=0;
            for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
            {
                i0 = TRparam[0].tilted_ring[k0][0];
                j0 = TRparam[0].tilted_ring[k0][1];

                // 0. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
                if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]) && HI_VF_mom0[0].data[j0][i0] > 0)
                {
                    mom4_temp[k] = HI_VF_mom0[0].data[j0][i0];
                    k++;
                }
            }
            robust_mean_std(mom4_temp, k, &mean_mom4_temp, &std_mom4_temp);
            TRparam[0].mean_mom4 = mean_mom4_temp;
            TRparam[0].std_mom4 = std_mom4_temp;

            for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
            {
                i0 = TRparam[0].tilted_ring[k0][0];
                j0 = TRparam[0].tilted_ring[k0][1];

                // 1. read vf_e
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0)
                    vf_e_temp[k0] = HI_VF_boxFiltered_sigma[0].data[j0][i0];
                else
                    vf_e_temp[k0] = 1; // a large value in km/s

                // 2. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
                if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]))
                {
                    if(HI_VF_mom0[0].data[j0][i0] > 1E5*mean_mom4_temp)
                    {
                        mom4_temp[k0] = mean_mom4_temp;
                    }
                    else
                    {
                        mom4_temp[k0] = HI_VF_mom0[0].data[j0][i0];
                    }
                }

                // 3. read mom2 : this is used for matching the vlos_e_w with mom2 distribution (see below)
                if(!isinf(HI_VF_mom2[0].data[j0][i0]) && !isnan(HI_VF_mom2[0].data[j0][i0]) && HI_VF_mom2[0].data[j0][i0] > 0)
                    mom2_temp[k0] = HI_VF_mom2[0].data[j0][i0];
                else
                    mom2_temp[k0] = 0; // a large value in km/s

                // 3. read HI_VF : derive dispersions of VLOS in a given ring which is used for VROT_e later
                //if(!isinf(HI_VF[0].data[j0][i0]) && !isnan(HI_VF[0].data[j0][i0]))
                //    mom2_temp[k0] = HI_VF[0].data[j0][i0];
                //else
                //    mom2_temp[k0] = 0; // a large value in km/s

                // geometry & flux-weighted vlos_e_w
                //vf_e_temp_w[k0] = vf_e_temp[k0] / (mom4_temp[k0]);
                vf_e_temp_w[k0] = vf_e_temp[k0];
            }

            // 1. compute the mean of geometry & flux-weighted variable vlos_e_w
            robust_mean_std(vf_e_temp_w, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e_temp_w, &std_vf_e_temp_w);
            // 2. compute the mean of mom2 
            robust_mean_std(mom2_temp, TRparam[0].Npoints_in_tilted_ring, &mean_mom2_temp, &std_mom2_temp);
            // compute the mean & std of vf_e
            robust_mean_std(vf_e_temp, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e, &std_vf_e);

            TRparam[0].dispersion_VLOS = mean_mom2_temp; // this will be used for the quadratic addition in the error of VROT 

            // scale factor
            //TRparam[0].scale_factor_var_vlose_w = mean_mom2_temp / mean_vf_e_temp_w;
            TRparam[0].scale_factor_var_vlose_w = mean_vf_e / mean_vf_e_temp_w;

            free(mom4_temp);
            free(mom2_temp);
            free(vf_e_temp);
            free(vf_e_temp_w);
        }
        else if(TRparam[0].final_fit == 'N')
        {
            double *mom4_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp_w = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);

            for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
            {
                i0 = TRparam[0].tilted_ring[k0][0];
                j0 = TRparam[0].tilted_ring[k0][1];

                // 1. read vf_e
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0)
                    vf_e_temp[k0] = HI_VF_boxFiltered_sigma[0].data[j0][i0];
                else
                    vf_e_temp[k0] = 1; // a large value in km/s

                // 2. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
                if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]) && HI_VF_mom0[0].data[j0][i0] > 0)
                    mom4_temp[k0] = HI_VF_mom0[0].data[j0][i0];
                else
                    mom4_temp[k0] = 1E-5; // a small value mjy/beam

                // geometry & flux-weighted vlos_e_w
                vf_e_temp_w[k0] = vf_e_temp[k0];
            }

            // 1. compute the mean of geometry & flux-weighted variable vlos_e_w
            robust_mean_std(vf_e_temp_w, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e_temp_w, &std_vf_e_temp_w);
            // compute the mean & std of vf_e
            robust_mean_std(vf_e_temp, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e, &std_vf_e);

            // scale factor
            //TRparam[0].scale_factor_var_vlose_w = TRparam[0].e_sigma / mean_vf_e_temp_w;
            TRparam[0].scale_factor_var_vlose_w = mean_vf_e / mean_vf_e_temp_w;

            free(mom4_temp);
            free(vf_e_temp);
            free(vf_e_temp_w);
        }
    }

    // Count TRparam[0].Npoints_in_tilted_ring_decim0
    // between ring0 and ring1
    // 0 - nax1
    decimX = 0;
    decimY = 0;
    Npoints_in_tiltedRing = 0;

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0; 
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    // count the total pixels in a given ring
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {
                r_i0j0 = r_ij_pa_incl_given_W(i0, j0, xpos, ypos, pa, incl, TRparam);
                acc = gsl_interp_accel_alloc();
                spline = gsl_spline_alloc (gsl_interp_cspline, TRparam[0].Nrings);
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl, TRparam[0].Nrings);

                // interpolate incl at r_i0j0
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    incl_intp = TRparam[0].incl[TRparam[0].Nrings-1]; // outermost incl
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    incl_intp = TRparam[0].incl[0]; // outermost incl
                }
                else
                {
                    incl_intp = gsl_spline_eval (spline, r_i0j0, acc);
                }

                // interpolate pa at r_i0j0
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].pa, TRparam[0].Nrings);
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    pa_intp = TRparam[0].pa[TRparam[0].Nrings-1]; // outermost pa
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    pa_intp = TRparam[0].pa[0]; // outermost pa
                }
                else
                {
                    pa_intp = gsl_spline_eval (spline, r_i0j0, acc);
                }
                gsl_spline_free(spline);
                gsl_interp_accel_free(acc);


                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa_intp*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa_intp*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa_intp*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa_intp*M_PI/180.))/cos(incl_intp*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0] && costh > sine_free_angle)
                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == -2 && theta > 0) // right side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 2 && theta < 0) // left side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }
    TRparam[0].Npoints_in_tilted_ring_decim0 = Npoints_in_tiltedRing;

    // give blank area based errors: this is to minimise the effect of
    // perimeters...
    // 0 - nax1
    //Npoints_in_tiltedRing=0;
    decimX = (int)TRparam[0].decimX_trfit;
    decimY = (int)TRparam[0].decimY_trfit;
    free(weight_temp);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void define_tiltedRing_ellipse_rings(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side)
{
    double i=0, j=0;
    int i0=0, j0=0, k0=0, k=0;
    int decimX, decimY;
    int Npoints_in_tiltedRing_without_decimals;
    int x0, y0, x1, y1, box_x, box_y;
    double theta=0., r=0., r_inner=0., r_outter=0.;
    double Npoints_in_tiltedRing=0;
    double Npoints_in_tiltedRing_total_including_blanks=0;
    double ring0_regrad, ring1_regrad;
    double sine_free_angle, costh;
    double mean_ve_temp, std_ve_temp;
    double mean_vf_e_temp_w, std_vf_e_temp_w;
    double mean_mom2_temp, std_mom2_temp;
    double mean_mom4_temp, std_mom4_temp;
    double HI_VF_weight_TRfit_max;
    double mean_vf_e, std_vf_e;
    Ellipse_Parameter ellipse;
    ring_parameter ring;

    ellipse.xpos = xpos;
    ellipse.ypos = ypos;
    ellipse.a = ring0 + (ring1-ring0)/2.0;
    ellipse.b = ellipse.a * cos(incl*M_PI/180.);
    ellipse.e = sqrt(1-(ellipse.b*ellipse.b)/(ellipse.a*ellipse.a));

    // 0. set free angle
    sine_free_angle = fabs(sin(TRparam[0].free_angle*M_PI/180.));

    // 1. extract the region to fit with decimals starting from the centre position given
    Npoints_in_tiltedRing_total_including_blanks = 0;
    decimX = 0;
    decimY = 0;

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0;
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;

    double *weight_temp = malloc(sizeof(double) * (x1-x0)*(y1-y0));

    // count the total pixels in a given ring
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            //if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            if((HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0])
            {
                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                if(r >= ring0 && r < ring1 && costh > sine_free_angle)
                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        Npoints_in_tiltedRing_total_including_blanks += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        Npoints_in_tiltedRing_total_including_blanks += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        Npoints_in_tiltedRing_total_including_blanks += 1;
                    }
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }
    TRparam[0].N_all_pixels_in_a_ring = Npoints_in_tiltedRing_total_including_blanks;

    // between ring0 and ring1
    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0; 
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    // this is for decimX=0, decimY=0
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {
                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0] && costh > sine_free_angle)
                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                        // squared in the MLE by sigma*sigma...
                        HI_VF_weight_TRfit[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)); // weight : note that radial weight doesn't need to be applied as all points within a ring have the same radius.
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        HI_VF_weight_TRfit[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)); // weight : note that radial weight doesn't need to be applied as all points within a ring have the same radius.
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        HI_VF_weight_TRfit[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)); // weight : note that radial weight doesn't need to be applied as all points within a ring have the same radius.
                    }
                    
                    if(HI_VF_weight_TRfit[0].data[j0][i0] == 0 || isnan(HI_VF_weight_TRfit[0].data[j0][i0]) || isinf(HI_VF_weight_TRfit[0].data[j0][i0]))
                        HI_VF_weight_TRfit[0].data[j0][i0] = 1E-2; // put a small value
                }
            }
        }
    }



    // between ring0 and ring1
    Npoints_in_tiltedRing = 0;
    decimX = (int)TRparam[0].decimX_trfit;
    decimY = (int)TRparam[0].decimY_trfit;

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0; 
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {
                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0] && costh > sine_free_angle)
                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                        // squared in the MLE by sigma*sigma...
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_weight_TRfit[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_weight_TRfit[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_weight_TRfit[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    
                    if(HI_VF_weight_TRfit[0].data[j0][i0] == 0 || isnan(HI_VF_weight_TRfit[0].data[j0][i0]) || isinf(HI_VF_weight_TRfit[0].data[j0][i0]))
                        HI_VF_weight_TRfit[0].data[j0][i0] = 1E-2; // put a small value
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }
    TRparam[0].Npoints_in_tilted_ring = Npoints_in_tiltedRing;


    HI_VF_weight_TRfit_max = gsl_stats_max(weight_temp, 1, Npoints_in_tiltedRing);
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            // normalise (0 ~ 1)
            HI_VF_weight_TRfit[0].data[j0][i0] = HI_VF_weight_TRfit[0].data[j0][i0] / HI_VF_weight_TRfit_max;
        }
    }

    // Derive the geometry (+flux) weighted scale factor for vlos_e
    if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma > 0) // constant vlos_e mode
    {
        double *ve_temp = malloc(sizeof(double) * Npoints_in_tiltedRing);
        // compute scale_factor_const_vlose_w
        for(i0=0; i0<(int)Npoints_in_tiltedRing; i0++)
        {
            // normalised constant one channel vlos_error with the geometry weight
            ve_temp[i0] = TRparam[0].e_sigma / (weight_temp[i0]/HI_VF_weight_TRfit_max);
        }

        // compute the mean of geometry-weighted constant vlos_e
        robust_mean_std(ve_temp, Npoints_in_tiltedRing, &mean_ve_temp, &std_ve_temp);
        // compute a scale factor for matching the geometry-weighted constant
        // vlos_e (which includes large values) with the input constant sigma_v
        // (this is to avoid large vlos_e in the bayesian fit
        TRparam[0].scale_factor_const_vlose_w = TRparam[0].e_sigma / mean_ve_temp;
        free(ve_temp);
    }
    else if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0) // variable vlos_e mode : full fit
    {
        double *mom4_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *mom2_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *vf_e_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *vf_e_temp_w = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);

        double mean_vf_e, std_vf_e;
        for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
        {
            i0 = TRparam[0].tilted_ring[k0][0];
            j0 = TRparam[0].tilted_ring[k0][1];

            // 1. read vf_e
            if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0)
                vf_e_temp[k0] = HI_VF_boxFiltered_sigma[0].data[j0][i0];
            else
                vf_e_temp[k0] = 1; // a large value in km/s

            // 2. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
            if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]) && HI_VF_mom0[0].data[j0][i0] > 0)
                mom4_temp[k0] = HI_VF_mom0[0].data[j0][i0];
            else
                mom4_temp[k0] = 1E-5; // a small value mjy/beam

            // 3. read mom2 : this is used for matching the vlos_e_w with mom2 distribution (see below)
            if(!isinf(HI_VF_mom2[0].data[j0][i0]) && !isnan(HI_VF_mom2[0].data[j0][i0]) && HI_VF_mom2[0].data[j0][i0] > 0)
                mom2_temp[k0] = HI_VF_mom2[0].data[j0][i0];
            else
                mom2_temp[k0] = 1E3; // a large value in km/s

            // geometry & flux-weighted vlos_e_w
            vf_e_temp_w[k0] = vf_e_temp[k0] / (weight_temp[k0]/HI_VF_weight_TRfit_max) / sqrt(mom4_temp[k0]);
        }

        // 1. compute the mean of geometry & flux-weighted variable vlos_e_w
        robust_mean_std(vf_e_temp_w, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e_temp_w, &std_vf_e_temp_w);
        // 2. compute the mean of mom2 
        robust_mean_std(mom2_temp, TRparam[0].Npoints_in_tilted_ring, &mean_mom2_temp, &std_mom2_temp);
        // compute mean & std of vf_e
        robust_mean_std(vf_e_temp, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e, &std_vf_e);

        // scale factor
        //TRparam[0].scale_factor_var_vlose_w = mean_mom2_temp / mean_vf_e_temp_w;
        TRparam[0].scale_factor_var_vlose_w = mean_vf_e / mean_vf_e_temp_w;

        free(mom4_temp);
        free(mom2_temp);
        free(vf_e_temp);
        free(vf_e_temp_w);
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0) // variable vlos_e mode : partial fit
    {
        if(TRparam[0].final_fit == 'Y')
        {
            double *mom4_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *mom2_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp_w = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);

            // compute mean & std of mom4
            k=0;
            for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
            {
                i0 = TRparam[0].tilted_ring[k0][0];
                j0 = TRparam[0].tilted_ring[k0][1];

                // 0. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
                if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]) && HI_VF_mom0[0].data[j0][i0] > 0)
                {
                    mom4_temp[k] = HI_VF_mom0[0].data[j0][i0];
                    k++;
                }
            }
            robust_mean_std(mom4_temp, k, &mean_mom4_temp, &std_mom4_temp);
            TRparam[0].mean_mom4 = mean_mom4_temp;
            TRparam[0].std_mom4 = std_mom4_temp;

            for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
            {
                i0 = TRparam[0].tilted_ring[k0][0];
                j0 = TRparam[0].tilted_ring[k0][1];

                // 1. read vf_e
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0)
                    vf_e_temp[k0] = HI_VF_boxFiltered_sigma[0].data[j0][i0];
                else
                    vf_e_temp[k0] = 1; // a large value in km/s

                // 2. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
                if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]))
                {
                    if(HI_VF_mom0[0].data[j0][i0] > 1E5*mean_mom4_temp)
                    {
                        mom4_temp[k0] = mean_mom4_temp;
                    }
                    else
                    {
                        mom4_temp[k0] = HI_VF_mom0[0].data[j0][i0];
                    }
                }

                // 3. read mom2 : this is used for matching the vlos_e_w with mom2 distribution (see below)
                if(!isinf(HI_VF_mom2[0].data[j0][i0]) && !isnan(HI_VF_mom2[0].data[j0][i0]) && HI_VF_mom2[0].data[j0][i0] > 0)
                    mom2_temp[k0] = HI_VF_mom2[0].data[j0][i0];
                else
                    mom2_temp[k0] = 1E3; // a large value in km/s

                // geometry & flux-weighted vlos_e_w
                //vf_e_temp_w[k0] = vf_e_temp[k0] / (weight_temp[k0]/HI_VF_weight_TRfit_max) / sqrt(mom4_temp[k0]);
                //vf_e_temp_w[k0] = vf_e_temp[k0] / (mom4_temp[k0]);
                vf_e_temp_w[k0] = vf_e_temp[k0];
            }

            // 1. compute the mean of geometry & flux-weighted variable vlos_e_w
            robust_mean_std(vf_e_temp_w, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e_temp_w, &std_vf_e_temp_w);
            // 2. compute the mean of mom2 
            robust_mean_std(mom2_temp, TRparam[0].Npoints_in_tilted_ring, &mean_mom2_temp, &std_mom2_temp);
            // compute mean & std of vf_e
            robust_mean_std(vf_e_temp, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e, &std_vf_e);

            // scale factor
            //TRparam[0].scale_factor_var_vlose_w = mean_mom2_temp / mean_vf_e_temp_w;
            TRparam[0].scale_factor_var_vlose_w = mean_vf_e / mean_vf_e_temp_w;

            free(mom4_temp);
            free(mom2_temp);
            free(vf_e_temp);
            free(vf_e_temp_w);
        }
        else if(TRparam[0].final_fit == 'N')
        {
            double *mom4_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
            double *vf_e_temp_w = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);

            for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
            {
                i0 = TRparam[0].tilted_ring[k0][0];
                j0 = TRparam[0].tilted_ring[k0][1];

                // 1. read vf_e
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0)
                    vf_e_temp[k0] = HI_VF_boxFiltered_sigma[0].data[j0][i0];
                else
                    vf_e_temp[k0] = 1; // a large value in km/s

                // 2. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
                if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]) && HI_VF_mom0[0].data[j0][i0] > 0)
                    mom4_temp[k0] = HI_VF_mom0[0].data[j0][i0];
                else
                    mom4_temp[k0] = 1E-5; // a small value mjy/beam

                // geometry & flux-weighted vlos_e_w
                //vf_e_temp_w[k0] = TRparam[0].e_sigma / (weight_temp[k0]/HI_VF_weight_TRfit_max) / mom4_temp[k0];
                //vf_e_temp_w[k0] = TRparam[0].e_sigma / (weight_temp[k0]/HI_VF_weight_TRfit_max);
                //vf_e_temp_w[k0] = vf_e_temp[k0] / (weight_temp[k0]/HI_VF_weight_TRfit_max);
                vf_e_temp_w[k0] = vf_e_temp[k0];
            }

            // 1. compute the mean of geometry & flux-weighted variable vlos_e_w
            robust_mean_std(vf_e_temp_w, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e_temp_w, &std_vf_e_temp_w);
            // compute mean & std of vf_e
            robust_mean_std(vf_e_temp, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e, &std_vf_e);

            // scale factor
            //TRparam[0].scale_factor_var_vlose_w = TRparam[0].e_sigma / mean_vf_e_temp_w;
            TRparam[0].scale_factor_var_vlose_w = mean_vf_e / mean_vf_e_temp_w;

            free(mom4_temp);
            free(vf_e_temp);
            free(vf_e_temp_w);
        }
    }

    // Count TRparam[0].Npoints_in_tilted_ring_decim0
    // between ring0 and ring1
    // 0 - nax1
    decimX = 0;
    decimY = 0;
    Npoints_in_tiltedRing = 0;

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0; 
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    // count the total pixels in a given ring
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {

                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0] && costh > sine_free_angle)
                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == -2 && theta > 0) // right side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 2 && theta < 0) // left side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }
    TRparam[0].Npoints_in_tilted_ring_decim0 = Npoints_in_tiltedRing;

    decimX = (int)TRparam[0].decimX_trfit;
    decimY = (int)TRparam[0].decimY_trfit;
    free(weight_temp);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void define_tiltedRing_ellipse_rings_extend_rings(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side)
{
    double i=0, j=0;
    int i0=0, j0=0, k0=0;
    int decimX, decimY;
    int Npoints_in_tiltedRing_without_decimals;
    int x0, y0, x1, y1, box_x, box_y;
    double theta=0., r=0., r_inner=0., r_outter=0.;
    double Npoints_in_tiltedRing=0, Npoints_in_tiltedRing_decim0=0;
    double Npoints_in_tiltedRing_total_including_blanks=0;
    double ring0_regrad, ring1_regrad;
    double sine_free_angle, costh;
    double mean_ve_temp, std_ve_temp;
    double mean_vf_e_temp_w, std_vf_e_temp_w;
    double mean_mom2_temp, std_mom2_temp;
    double HI_VF_weight_TRfit_max;
    double fract_Navail_Nall_decim0;
    Ellipse_Parameter ellipse;
    ring_parameter ring;

    ellipse.xpos = xpos;
    ellipse.ypos = ypos;
    ellipse.a = ring0 + (ring1-ring0)/2.0;
    ellipse.b = ellipse.a * cos(incl*M_PI/180.);
    ellipse.e = sqrt(1-(ellipse.b*ellipse.b)/(ellipse.a*ellipse.a));

    // 0. set free angle
    sine_free_angle = fabs(sin(TRparam[0].free_angle*M_PI/180.));

    // 1. extract the region to fit with decimals starting from the centre position given

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0;
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;

    Npoints_in_tiltedRing_total_including_blanks = 0;
    decimX = 0;
    decimY = 0;
    // count the total pixels in a given ring
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
            j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
            r = sqrt(i*i+j*j); // distance from centre

            if(r < 1) // if r is smaller than one pixel
            {
                r = 1;
                theta = 0.0;
            }
            else
                theta = atan2((double)j, (double)i)*180./M_PI;

            costh = fabs(cos(theta*M_PI/180.));

            if(r >= ring0 && r < ring1 && costh > sine_free_angle)
            {
                if(side == -1 && fabs(theta) <= 90.0) // receding side
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
                else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
                else if(side == 0 || side == 999) // both sides
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }
    //TRparam[0].N_all_pixels_in_a_ring = Npoints_in_tiltedRing_total_including_blanks;

    
    // between ring0 and ring1
    Npoints_in_tiltedRing_decim0 = 0;
    decimX = 0;
    decimY = 0;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {
                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                        // squared in the MLE by sigma*sigma...
                        Npoints_in_tiltedRing_decim0 += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        Npoints_in_tiltedRing_decim0 += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        Npoints_in_tiltedRing_decim0 += 1;
                    }
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }

    // between ring0 and ring1
    Npoints_in_tiltedRing = 0;
    //decimX = (int)TRparam[0].decimX_trfit;
    //decimY = (int)TRparam[0].decimY_trfit;
    //decimX = (int)TRparam[0].decimX;
    //decimY = (int)TRparam[0].decimY;
    decimX = 0;
    decimY = 0;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {
                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0] && costh > sine_free_angle)
                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                        // squared in the MLE by sigma*sigma...
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                }
            }
            j0 += decimY; 
        }
        i0 += decimX; 
    }


    // give blank area based errors: this is to minimise the effect of
    // perimeters...
    // 0 - nax1
    //Npoints_in_tiltedRing=0;

    // Compute fractional weight : Nall_including_blanks / Nvaialble_decim_user
    // between ring0 and ring1
    // This is used for computing combined weight for the Einasto fit
    decimX = 0;
    decimY = 0;
    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0; 
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
            j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
            r = sqrt(i*i+j*j); // distance from centre

            if(r < 1) // if r is smaller than one pixel
            {
                r = 1;
                theta = 0.0;
            }
            else
                theta = atan2((double)j, (double)i)*180./M_PI;

            costh = fabs(cos(theta*M_PI/180.));

            if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
            {
                if(side == -1 && fabs(theta) <= 90.0) // receding side
                {
                    //HI_VF_fract_navail_nall[0].data[j0][i0] = Npoints_in_tiltedRing/Npoints_in_tiltedRing_total_including_blanks;
                    fract_Navail_Nall_decim0 = Npoints_in_tiltedRing_decim0/Npoints_in_tiltedRing_total_including_blanks;
                }
                else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                {
                    //HI_VF_fract_navail_nall[0].data[j0][i0] = Npoints_in_tiltedRing/Npoints_in_tiltedRing_total_including_blanks;
                    fract_Navail_Nall_decim0 = Npoints_in_tiltedRing_decim0/Npoints_in_tiltedRing_total_including_blanks;
                }
                else if(side == 0 || side == 999) // both sides
                {
                    //HI_VF_fract_navail_nall[0].data[j0][i0] = Npoints_in_tiltedRing/Npoints_in_tiltedRing_total_including_blanks;
                    fract_Navail_Nall_decim0 = Npoints_in_tiltedRing_decim0/Npoints_in_tiltedRing_total_including_blanks;
                }
                
                //if(HI_VF_fract_navail_nall[0].data[j0][i0] == 0 || isnan(HI_VF_fract_navail_nall[0].data[j0][i0]) || isinf(HI_VF_fract_navail_nall[0].data[j0][i0]))
                if(isnan(HI_VF_fract_navail_nall[0].data[j0][i0]) || isinf(HI_VF_fract_navail_nall[0].data[j0][i0])) ;
                    // HI_VF_fract_navail_nall[0].data[j0][i0] = 1; // no weight

                //if(ring0 >= TRparam[0].ring_radius[TRparam[0].Nrings-3])
                //    HI_VF_fract_navail_nall[0].data[j0][i0] = 1; // no weight
            }

            /*
            if(r >= ring0-TRparam[0].ring_w/1.0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
            {
                if(HI_VF_fract_navail_nall[0].data[j0][i0] == 1E-3 || HI_VF_fract_navail_nall[0].data[j0][i0] == 0)
                    HI_VF_fract_navail_nall[0].data[j0][i0] = Npoints_in_tiltedRing/Npoints_in_tiltedRing_total_including_blanks;
            }
            */

            j0 += decimY; 
        }
        i0 += decimX; 
    }
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void find_Navail_Nall_pixels(double xpos, double ypos, double pa, double incl, double ring0, double ring1, TR_ringParameters *TRparam, int side, double *Nall_ring, double *Navail_ring)
{
    double i=0, j=0;
    int i0=0, j0=0, k0=0;
    int decimX, decimY;
    int Npoints_in_tiltedRing_without_decimals;
    int x0, y0, x1, y1, box_x, box_y;
    double theta=0., r=0., r_inner=0., r_outter=0.;
    double Npoints_in_tiltedRing=0, Npoints_in_tiltedRing_decim0=0;
    double Npoints_in_tiltedRing_total_including_blanks=0;
    double ring0_regrad, ring1_regrad;
    double sine_free_angle, costh;
    double mean_ve_temp, std_ve_temp;
    double mean_vf_e_temp_w, std_vf_e_temp_w;
    double mean_mom2_temp, std_mom2_temp;
    double HI_VF_weight_TRfit_max;
    double fract_Navail_Nall_decim0;
    Ellipse_Parameter ellipse;
    ring_parameter ring;

    ellipse.xpos = xpos;
    ellipse.ypos = ypos;

    // 0. set free angle
    sine_free_angle = fabs(sin(TRparam[0].free_angle*M_PI/180.));

    // 1. extract the region to fit with decimals starting from the centre position given

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0;
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;

    Npoints_in_tiltedRing_total_including_blanks = 0;
    decimX = 0;
    decimY = 0;
    //decimX = TRparam[0].decimX_einasto_halofit;
    //decimY = TRparam[0].decimY_einasto_halofit;
    // count the total pixels in a given ring
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
            j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
            r = sqrt(i*i+j*j); // distance from centre

            if(r < 1) // if r is smaller than one pixel
            {
                r = 1;
                theta = 0.0;
            }
            else
                theta = atan2((double)j, (double)i)*180./M_PI;

            costh = fabs(cos(theta*M_PI/180.));

            if(r >= ring0 && r < ring1 && costh > sine_free_angle)
            {
                if(side == -1 && fabs(theta) <= 90.0) // receding side
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
                else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
                else if(side == 0 || side == 999) // both sides
                {
                    Npoints_in_tiltedRing_total_including_blanks += 1;
                }
            }
            j0 += decimY;
        }
        i0 += decimX;
    }
    //TRparam[0].N_all_pixels_in_a_ring = Npoints_in_tiltedRing_total_including_blanks;
    *Nall_ring = Npoints_in_tiltedRing_total_including_blanks;
    // between ring0 and ring1
    Npoints_in_tiltedRing = 0;
    decimX = 0;
    decimY = 0;
    //decimX = TRparam[0].decimX_einasto_halofit;
    //decimY = TRparam[0].decimY_einasto_halofit;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered[0].data[j0][i0]))
            {
                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));

                //if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered[0].data[j0][i0]+1) > HI_VF_boxFiltered[0].data[j0][i0] && costh > sine_free_angle)
                if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                        // squared in the MLE by sigma*sigma...
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0 || side == 999) // both sides
                    {
                        Npoints_in_tiltedRing += 1;
                    }
                }
            }
            j0 += decimY;
        }
        i0 += decimX;
    }
    *Navail_ring = Npoints_in_tiltedRing;

    // give blank area based errors: this is to minimise the effect of
    // perimeters...
    // 0 - nax1
    //Npoints_in_tiltedRing=0;

    // Compute fractional weight : Nall_including_blanks / Nvaialble_decim_user
    // between ring0 and ring1
    // This is used for computing combined weight for the Einasto fit
    decimX = 0;
    decimY = 0;
    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0;
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
            j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
            r = sqrt(i*i+j*j); // distance from centre

            if(r < 1) // if r is smaller than one pixel
            {
                r = 1;
                theta = 0.0;
            }
            else
                theta = atan2((double)j, (double)i)*180./M_PI;

            costh = fabs(cos(theta*M_PI/180.));

            if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
            {
                if(side == -1 && fabs(theta) <= 90.0) // receding side
                {
                    HI_VF_fract_navail_nall[0].data[j0][i0] = Npoints_in_tiltedRing/Npoints_in_tiltedRing_total_including_blanks;
                }
                else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                {
                    HI_VF_fract_navail_nall[0].data[j0][i0] = Npoints_in_tiltedRing/Npoints_in_tiltedRing_total_including_blanks;
                }
                else if(side == 0 || side == 999) // both sides
                {
                    HI_VF_fract_navail_nall[0].data[j0][i0] = Npoints_in_tiltedRing/Npoints_in_tiltedRing_total_including_blanks;
                }

                if(isnan(HI_VF_fract_navail_nall[0].data[j0][i0]) || isinf(HI_VF_fract_navail_nall[0].data[j0][i0]))
                    HI_VF_fract_navail_nall[0].data[j0][i0] = 1E90; // no weight
            }
            j0 += decimY;
        }
        i0 += decimX;
    }
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void trfit_multinest_trfit_rings_normal(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *vrot, char vrotfix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, char *finalfit, char final_fit)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0, ii=0, j=0;
    int i0=0, j0=0;
    int x, y;
    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    int Nrings;
    double efr, tol, Ztol, logZero;
    double total_Npoints_allrings=0.;
    double ri=0., ro=0., ring;

    double A, B, del_A, del_B, costheta, costheta1, costheta2, del_ring, del_costheta1, del_costheta2, del_costheta, del_sini, sini, sigma_Vlos, Vsys, sigma_Vsys, Vhat, sigma_Vhat;

    // GSL histogram statistics
    double gsl_mean_xpos, gsl_std_xpos, gsl_max_xpos, gsl_min_xpos;
    double gsl_mean_xpos_e, gsl_std_xpos_e, gsl_max_xpos_e, gsl_min_xpos_e;
    double gsl_mean_ypos, gsl_std_ypos, gsl_max_ypos, gsl_min_ypos;
    double gsl_mean_ypos_e, gsl_std_ypos_e, gsl_max_ypos_e, gsl_min_ypos_e;
    double gsl_mean_vsys, gsl_std_vsys, gsl_max_vsys, gsl_min_vsys;
    double gsl_mean_vsys_e, gsl_std_vsys_e, gsl_max_vsys_e, gsl_min_vsys_e;
    double gsl_mean_pa, gsl_std_pa, gsl_max_pa, gsl_min_pa;
    double gsl_mean_pa_e, gsl_std_pa_e, gsl_max_pa_e, gsl_min_pa_e;
    double gsl_mean_incl, gsl_std_incl, gsl_max_incl, gsl_min_incl;
    double gsl_mean_incl_e, gsl_std_incl_e, gsl_max_incl_e, gsl_min_incl_e;
    double gsl_mean_vrot, gsl_std_vrot, gsl_max_vrot, gsl_min_vrot;
    double gsl_mean_vrot_e, gsl_std_vrot_e, gsl_max_vrot_e, gsl_min_vrot_e;
    double gsl_mean_vrad, gsl_std_vrad, gsl_max_vrad, gsl_min_vrad;
    double gsl_mean_vrad_e, gsl_std_vrad_e, gsl_max_vrad_e, gsl_min_vrad_e;

    double total_error_temp, vrot_median;

    // GSL histogram statistics

    double xpos1_0, xpos2_0, ypos1_0, ypos2_0;
    double lower_bound_Vrot_err_inaring, upper_bound_Vrot_err_inaring;
    double hist_bin_Vrot_err_inaring, hist_mean_Vrot_err_inaring, hist_std_Vrot_err_inaring;
    double hist_bin_Vrad_err_inaring, hist_mean_Vrad_err_inaring, hist_std_Vrad_err_inaring;

    double FD_h, IQR;
    double ll_bin, ul_bin;
    int FD_nbins, h_loop;

    double ellipse_a, ellipse_b, xpos_prior_gab, ypos_prior_gab;
    //double *Vrot_err_inaring = malloc(sizeof(double) * 1); // resize later using realloc in the loop below
    double Vrot_err_inaring[50000] = {0}; // resize later using realloc in the loop below
    double Vrad_err_inaring[50000] = {0}; // resize later using realloc in the loop below
    //TRparam[0].Nrings = (int)((TRparam[0].ring_e-TRparam[0].ring_s)/TRparam[0].ring_w)+1;

    // dynamic total_error 1D array
    double *total_error_ring_params = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array
    double *total_error_ring_params_temp = malloc(sizeof(double) * TRparam[0].Nrings); // 1D array
    double lower_bound_total_error_ring_params, upper_bound_total_error_ring_params;
    double hist_bin_total_error_ring_params, hist_mean_total_error_ring_params, hist_std_total_error_ring_params;

    // xpos
    double *xpos_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *xpos_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_xpos, upper_bound_xpos;
    double hist_bin_xpos, hist_mean_xpos, hist_std_xpos;
    // xpos error
    double *xpos_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_xpos_err, upper_bound_xpos_err;
    double hist_bin_xpos_err, hist_mean_xpos_err, hist_std_xpos_err;

    // ypos
    double *ypos_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *ypos_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_ypos, upper_bound_ypos;
    double hist_bin_ypos, hist_mean_ypos, hist_std_ypos;
    // ypos error
    double *ypos_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_ypos_err, upper_bound_ypos_err;
    double hist_bin_ypos_err, hist_mean_ypos_err, hist_std_ypos_err;

    // vsys
    double *vsys_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vsys_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_vsys, upper_bound_vsys;
    double hist_bin_vsys, hist_mean_vsys, hist_std_vsys;
    // vsys error
    double *vsys_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_vsys_err, upper_bound_vsys_err;
    double hist_bin_vsys_err, hist_mean_vsys_err, hist_std_vsys_err;

    // pa
    double *pa_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *pa_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_pa, upper_bound_pa;
    double hist_bin_pa, hist_mean_pa, hist_std_pa;
    // pa error
    double *pa_err = malloc(sizeof(double) * TRparam[0].Nrings);

    double lower_bound_pa_err, upper_bound_pa_err;
    double hist_bin_pa_err, hist_mean_pa_err, hist_std_pa_err;

    // incl
    double *incl_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *incl_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_incl, upper_bound_incl;
    double hist_bin_incl, hist_mean_incl, hist_std_incl;
    // incl error
    double *incl_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_incl_err, upper_bound_incl_err;
    double hist_bin_incl_err, hist_mean_incl_err, hist_std_incl_err;

    // vrot
    double *vrot_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrot_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_vrot, upper_bound_vrot;
    double hist_bin_vrot, hist_mean_vrot, hist_std_vrot;
    // vrot error
    double *vrot_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_vrot_err, upper_bound_vrot_err;
    double hist_bin_vrot_err, hist_mean_vrot_err, hist_std_vrot_err;

    // vrad
    double *vrad_tr = malloc(sizeof(double) * TRparam[0].Nrings);
    double *vrad_tr_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_vrad, upper_bound_vrad;
    double hist_bin_vrad, hist_mean_vrad, hist_std_vrad;
    // vrad error
    double *vrad_err = malloc(sizeof(double) * TRparam[0].Nrings);
    double lower_bound_vrad_err, upper_bound_vrad_err;
    double hist_bin_vrad_err, hist_mean_vrad_err, hist_std_vrad_err;

    double XPOS, YPOS, VSYS, PA, INCL, VEXP, VROT, Vlos;
    double XPOS_e, YPOS_e, VSYS_e, PA_e, INCL_e, VEXP_e, VROT_e, Vlos_e;
    double T1, T2;
    double pdVROT_pdXPOS, pdVROT_pdYPOS, pdVROT_pdVSYS, pdVROT_pdPA, pdVROT_pdINCL, pdVROT_pdVEXP, pdVROT_pdVlos;
    double dT1_dXPOS, dT2_dXPOS, dT1_dYPOS, dT2_dYPOS, dT1_dPA, dT2_dPA, dVROT_dINCL, pdVROT_pdT1, pdVROT_pdT2;

    double sigma_v;

    double _ring_w, _xpos, _ypos, _pa, _incl;

    // gsl vector for ring params : this is for calculating max/min of errors. See below.
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


    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive_einasto_halofit;
    if(nlive > 50) nlive = 50; // for a quick TRfit
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;
    maxiter = multinest_param[0].maxiter;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;
    int pWrap[ndims];

    //strcpy(root, multinest_param[0].root);
    //strncpy(root, multinest_param[0].root, strlen(multinest_param[0].root));
    strcpy(root, "trfit.einastofit_rings.");

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;

    TRparam[0].xpos_fix = xposfix;
    TRparam[0].ypos_fix = yposfix;
    TRparam[0].vsys_fix = vsysfix;
    TRparam[0].pa_fix = pafix;
    TRparam[0].incl_fix = inclfix;
    TRparam[0].vrot_fix = vrotfix;
    TRparam[0].vrad_fix = vradfix;
    TRparam[0].sigma_factor_fix = sigmafactorfix;
    TRparam[0].final_fit = final_fit;

    /* Set up ring parameters with currently derived values: note that pa0[i] and incl0[i] are given in degree */
    set_nfree_params_trfit_multinest_trfit_rings_student(TRparam);
    ndims = TRparam[0].n_freeParams;
    nPar = TRparam[0].n_freeParams;
    nClsPar = TRparam[0].n_freeParams;
    for(i = 0; i < ndims; i++) pWrap[i] = multinest_param[0].pWrap[i];

    /* Start tilted-ring fits */
    MPI_Barrier(MPI_COMM_WORLD);

    for(i=0; i<TRparam[0].Nrings; i++)
    {
        ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
        ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        nlive = multinest_param[0].nlive_einasto_halofit;
        if(nlive > 50) nlive = 50; // for a quick TRfit
        define_tiltedRing(TRparam[0].xpos0[i], TRparam[0].ypos0[i], TRparam[0].pa0[i], TRparam[0].incl0[i], ri, ro, TRparam, side);

        if(TRparam[0].xpos_fix == 'F')
        {
            // xpos for TiltedRingModel() function : see the function for more details
            TRparam[0].xposF = TRparam[0].xpos0[i]; // TRparam[0].xpos0[i] is TRparam[0].xposF_EinastoFit from Einasto halo fitting!
            TRparam[0].xposF_e = TRparam[0].xpos_e[i];
            TRparam[0].xpos[i] = TRparam[0].xpos0[i]; // Einasto halo fit result
        }

        if(TRparam[0].ypos_fix == 'F')
        {
            // ypos for TiltedRingModel() function : see the function for more details
            TRparam[0].yposF = TRparam[0].ypos0[i];
            TRparam[0].yposF_e = TRparam[0].ypos_e[i];
            TRparam[0].ypos[i] = TRparam[0].ypos0[i]; // Einasto halo fit result
        }

        if(TRparam[0].vsys_fix == 'F')
        {
            // vsys for TiltedRingModel() function : see the function for more details
            TRparam[0].vsysF = TRparam[0].vsys0[i];
            TRparam[0].vsysF_e = TRparam[0].vsys_e[i];
            TRparam[0].vsys[i] = TRparam[0].vsys0[i]; // Einasto halo fit result
        }

        if(TRparam[0].pa_fix == 'F')
        {
            // use the current PA model value for multinest if not fitted: in degree
            TRparam[0].paF = TRparam[0].pa0[i];
            TRparam[0].paF_e = TRparam[0].pa_e[i];
            TRparam[0].pa[i] = TRparam[0].pa0[i]; // in degree
            TRparam[0].pa_temp[i] = TRparam[0].pa0[i]; // in degree
        }

        if(TRparam[0].incl_fix == 'F')
        {
            // use the current INCL model value for multinest if not fitted: in degree
            TRparam[0].inclF = TRparam[0].incl0[i];
            TRparam[0].inclF_e = TRparam[0].incl_e[i];
            TRparam[0].incl[i] = TRparam[0].incl0[i]; // in degree 
            TRparam[0].incl_temp[i] = TRparam[0].incl0[i]; // in degree 
        }

        if(TRparam[0].vrad_fix == 'F')
        {
            // use the current VRAD model value for multinest if not fitted: in km/s
            TRparam[0].vradF = TRparam[0].vrad[i];
            TRparam[0].vradF_e = TRparam[0].vrad_e[i];
            TRparam[0].vrad[i] = TRparam[0].vrad[i];
            TRparam[0].vrad_temp[i] = TRparam[0].vrad[i];
        }

        // update xpos & ypos priors based on the current ring radius: This is not necessary for this fitting as xpos & ypos are fixed here
        // update xpos & ypos priors based on the current ring radius + pa0 + incl0
        ellipse_a = ri;
        xpos_prior_gab = fabs(ellipse_a*sin(TRparam[0].pa0[i]*M_PI/180.));
        ypos_prior_gab = fabs(ellipse_a*cos(TRparam[0].pa0[i]*M_PI/180.));

        if(TRparam[0].xpos_fix == 'T')
        {
            TRparam[0].xpos1 = TRparam[0].xpos0[i] - xpos_prior_gab;
            TRparam[0].xpos2 = TRparam[0].xpos0[i] + xpos_prior_gab;
        }
        if(TRparam[0].ypos_fix == 'T')
        {
            TRparam[0].ypos1 = TRparam[0].ypos0[i] - ypos_prior_gab;
            TRparam[0].ypos2 = TRparam[0].ypos0[i] + ypos_prior_gab;
        }
        if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0)// partial fit
        {
            TRparam[0].sigma_factor1 = 0;
            TRparam[0].sigma_factor2 = 200;
        }

        /* Calling multinest */ 
        MPI_Barrier(MPI_COMM_WORLD);
        if(TRparam[0].Npoints_in_tilted_ring > 5)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_trfit, dumper_TRfits, TRparam);   
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // XPOS
        TRparam[0].xpos[i] = TRparam[0].xposF; // Einasto halo fit result
        TRparam[0].xpos0[i] = TRparam[0].xposF; // Einasto halo fit result
        TRparam[0].xpos_e[i] = TRparam[0].xposF_e;
        xpos_tr[i] = TRparam[0].xpos[i];
        xpos_err[i] = TRparam[0].xpos_e[i];


        // YPOS
        TRparam[0].ypos[i] = TRparam[0].yposF; // Einasto halo fit result
        TRparam[0].ypos0[i] = TRparam[0].yposF; // Einasto halo fit result
        TRparam[0].ypos_e[i] = TRparam[0].yposF_e;
        ypos_tr[i] = TRparam[0].ypos[i];
        ypos_err[i] = TRparam[0].ypos_e[i];

        
        // VSYS
        TRparam[0].vsys[i] = TRparam[0].vsysF; // Einasto halo fit result
        TRparam[0].vsys0[i] = TRparam[0].vsysF; // Einasto halo fit result
        TRparam[0].vsys_e[i] = TRparam[0].vsysF_e;
        vsys_tr[i] = TRparam[0].vsys[i];
        vsys_err[i] = TRparam[0].vsys_e[i];

        // PA
        TRparam[0].pa[i] = TRparam[0].paF; // Einasto halo fit result
        TRparam[0].pa_temp[i] = TRparam[0].paF; // Einasto halo fit result
        TRparam[0].pa0[i] = TRparam[0].paF; // Einasto halo fit result
        TRparam[0].pa_e[i] = TRparam[0].paF_e;
        pa_tr[i] = TRparam[0].pa[i];
        pa_err[i] = TRparam[0].pa_e[i];

        // INCL
        TRparam[0].incl[i] = TRparam[0].inclF; // Einasto halo fit result
        TRparam[0].incl_temp[i] = TRparam[0].inclF; // Einasto halo fit result
        TRparam[0].incl0[i] = TRparam[0].inclF; // Einasto halo fit result
        TRparam[0].incl_e[i] = TRparam[0].inclF_e;
        incl_tr[i] = TRparam[0].incl[i];
        incl_err[i] = TRparam[0].incl_e[i];

        // VROT
        if(side == 999) // error propagation not used. Use vrotF_e from multinest fitting
        {   
            TRparam[0].vrot[i] = TRparam[0].vrotF;
            TRparam[0].vrot0[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrot_temp[i] = TRparam[0].vrotF;
            TRparam[0].vrot_temp_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad[i] = TRparam[0].vradF;
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
        else if(side == 0 || side == -1 || side == 1)
        {
            for(j=0; j<TRparam[0].Npoints_in_tilted_ring; j++)
            {
                x = TRparam[0].tilted_ring[j][0];
                y = TRparam[0].tilted_ring[j][1];

                if(x-TRparam[0].xposF == 0) x += 1; // to avoid zero!
                if(y-TRparam[0].yposF == 0) y += 1; // to avoid zero!

                XPOS = TRparam[0].xposF;
                XPOS_e = TRparam[0].xposF_e;
                if(isinf(TRparam[0].xposF_e) || isnan(TRparam[0].xposF_e) || TRparam[0].xposF_e == 0) TRparam[0].xposF_e = 1;

                YPOS = TRparam[0].yposF;
                YPOS_e = TRparam[0].yposF_e;
                if(isinf(TRparam[0].yposF_e) || isnan(TRparam[0].yposF_e) || TRparam[0].yposF_e == 0) TRparam[0].yposF_e = 1;

                VSYS = TRparam[0].vsysF;
                VSYS_e = TRparam[0].vsysF_e;
                if(isinf(TRparam[0].vsysF_e) || isnan(TRparam[0].vsysF_e) || TRparam[0].vsysF_e == 0) TRparam[0].vsysF_e = 1;

                if(isinf(TRparam[0].paF_e) || isnan(TRparam[0].paF_e) || TRparam[0].paF_e == 0) TRparam[0].paF_e = 1;
                PA = TRparam[0].paF*M_PI/180.;
                PA_e = TRparam[0].paF_e*M_PI/180.; // in radian

                if(isinf(TRparam[0].inclF_e) || isnan(TRparam[0].inclF_e) || TRparam[0].inclF_e == 0) TRparam[0].inclF_e = 1;
                INCL = TRparam[0].inclF*M_PI/180.;
                INCL_e = TRparam[0].inclF_e*M_PI/180.; // in radian

                VEXP = TRparam[0].vradF;
                VEXP_e = TRparam[0].vradF_e;

                Vlos = HI_VF_boxFiltered[0].data[y][x];
                Vlos_e = HI_VF_boxFiltered_sigma_e_norm[0].data[y][x];

                T1 = -(x-XPOS)*sin(PA) + (y-YPOS)*cos(PA);
                T2 =  (x-XPOS)*cos(PA) + (y-YPOS)*sin(PA);

                dT1_dXPOS = sin(PA);
                dT2_dXPOS = -cos(PA);

                dT1_dYPOS = -cos(PA);
                dT2_dYPOS = -sin(PA);

                dT1_dPA = -(x-XPOS)*cos(PA) - (y-YPOS)*sin(PA);
                dT2_dPA = -(x-XPOS)*sin(PA) + (y-YPOS)*cos(PA);

                dVROT_dINCL = (Vlos/T1)*(-(cos(INCL)/pow(sin(INCL),2))*sqrt(T1*T1+T2*T2/pow(cos(INCL),2))
                            + (1./sin(INCL))*(T2*T2*sin(INCL)/pow(cos(INCL),3))/sqrt(T1*T1+T2*T2/pow(cos(INCL),2)))
                            - (VSYS/T1)*(-(cos(INCL)/pow(sin(INCL),2))*sqrt(T1*T1+T2*T2/pow(cos(INCL),2))
                            + (1./sin(INCL))*(T2*T2*sin(INCL)/pow(cos(INCL),3))/sqrt(T1*T1+T2*T2/pow(cos(INCL),2)))
                            + (T2/T1)*VEXP*sin(INCL)/pow(cos(INCL),2);

                pdVROT_pdT1 = (Vlos/sin(INCL))*(1./(T1*T1))*T1/sqrt(T1*T1+T2*T2/pow(cos(INCL),2))
                            + (VSYS/sin(INCL))*(1./(T1*T1))*T1/sqrt(T1*T1+T2*T2/pow(cos(INCL),2))
                            + (VEXP/cos(INCL))*(-T2/(T1*T1));

                pdVROT_pdT2 = (Vlos/sin(INCL))*(1./T1)*T2/(sqrt(T1*T1+T2*T2/pow(cos(INCL),2))*pow(cos(INCL),2))
                            - (VSYS/sin(INCL))*(1./T1)*T2/(sqrt(T1*T1+T2*T2/pow(cos(INCL),2))*pow(cos(INCL),2))
                            + (VEXP/cos(INCL))*(1./T1);

                pdVROT_pdXPOS = pdVROT_pdT1*dT1_dXPOS + pdVROT_pdT2*dT2_dXPOS;
                pdVROT_pdYPOS = pdVROT_pdT1*dT1_dYPOS + pdVROT_pdT2*dT2_dYPOS;
                pdVROT_pdVSYS = -(1./T1)*(1./sin(INCL))*sqrt(T1*T1+T2*T2/pow(cos(INCL),2)) + (T2/T1)*(VEXP/cos(INCL));
                pdVROT_pdPA = pdVROT_pdT1*dT1_dPA + pdVROT_pdT2*dT2_dPA;
                pdVROT_pdINCL = dVROT_dINCL;
                pdVROT_pdVEXP = (T2/T1)*(1/cos(INCL));
                pdVROT_pdVlos = (1./T1)*(1./sin(INCL))*sqrt(T1*T1+T2*T2/pow(cos(INCL),2));

                Vrot_err_inaring[j] = sqrt(pow(pdVROT_pdXPOS*XPOS_e, 2)
                                          +pow(pdVROT_pdYPOS*YPOS_e, 2)
                                          +pow(pdVROT_pdVSYS*VSYS_e, 2)
                                          +pow(pdVROT_pdPA*PA_e, 2)
                                          +pow(pdVROT_pdINCL*INCL_e, 2)
                                          +pow(pdVROT_pdVEXP*VEXP_e, 2)
                                          +pow(pdVROT_pdVlos*Vlos_e, 2));

                if(isnan(Vrot_err_inaring[j]) || isinf(Vrot_err_inaring[j]))
                    Vrot_err_inaring[j] = 999;
            }

            // Calculate mean & std of VROT errors
            if(TRparam[0].Npoints_in_tilted_ring == 0)
            {
                hist_mean_Vrot_err_inaring = 999;
                hist_std_Vrot_err_inaring = 999;
            }
            else
            {
                // Calculate robust mean & std of VROT errors based on histogram
                robust_mean_std_e(Vrot_err_inaring, TRparam[0].Npoints_in_tilted_ring, &hist_mean_Vrot_err_inaring, &hist_std_Vrot_err_inaring);
                //robust_mean_std_e(Vrot_err_inaring, TRparam[0].Npoints_in_tilted_ring, &hist_mean_Vrot_err_inaring, &hist_std_Vrot_err_inaring);
                //robust_mean_std_e(Vrad_err_inaring, TRparam[0].Npoints_in_tilted_ring, &hist_mean_Vrad_err_inaring, &hist_std_Vrad_err_inaring);
            }
        }

        if(side == 0)
        {   
            TRparam[0].vrot[i] = TRparam[0].vrotF;
            TRparam[0].vrot0[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e[i] = hist_mean_Vrot_err_inaring;
            //TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad[i] = TRparam[0].vradF;
            TRparam[0].vrad_temp[i] = TRparam[0].vradF;
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
        else if(side == -1)
        {   
            TRparam[0].vrot_rec[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e_rec[i] = hist_mean_Vrot_err_inaring;
            //TRparam[0].vrot_e_rec[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad_rec[i] = TRparam[0].vradF;
            TRparam[0].vrad_rec_e[i] = TRparam[0].vradF_e;

            // this is for printing : TRfits_multinest_using_ISOfit_rings with both sides (0) should be done lastely
            TRparam[0].vrot[i] = TRparam[0].vrot_rec[i];
            TRparam[0].vrot_e[i] = hist_mean_Vrot_err_inaring;
            //TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad[i] = TRparam[0].vrad_rec[i];
            TRparam[0].vrad_e[i] = TRparam[0].vrad_rec_e[i];
        }
        else if(side == 1)
        {   
            TRparam[0].vrot_app[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e_app[i] = hist_mean_Vrot_err_inaring;
            //TRparam[0].vrot_e_app[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad_app[i] = TRparam[0].vradF;
            TRparam[0].vrad_app_e[i] = TRparam[0].vradF_e;

            // this is for printing : TRfits_multinest_using_ISOfit_rings with both sides (0) should be done lastely
            TRparam[0].vrot[i] = TRparam[0].vrot_app[i];
            TRparam[0].vrot_e[i] = hist_mean_Vrot_err_inaring;
            //TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;
            TRparam[0].vrad[i] = TRparam[0].vrad_app[i];
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
        vrot_tr[i] = TRparam[0].vrotF;
        vrot_err[i] = TRparam[0].vrotF_e;

        vrad_tr[i] = TRparam[0].vradF;
        vrad_err[i] = TRparam[0].vradF_e;

        // sigma Factor
        TRparam[0].e_sigma_student_TR[i] = TRparam[0].e_sigma_tr; // Einasto halo fit result

        // copy ring params to gsl vectors for finding min/max of errors 
        gsl_vector_set(gsl_xpos, i, TRparam[0].xposF);
        gsl_vector_set(gsl_xpos_e, i, TRparam[0].xposF_e);
        gsl_vector_set(gsl_ypos, i, TRparam[0].yposF);
        gsl_vector_set(gsl_ypos_e, i, TRparam[0].yposF_e);
        gsl_vector_set(gsl_vsys, i, TRparam[0].vsysF);
        gsl_vector_set(gsl_vsys_e, i, TRparam[0].vsysF_e);
        gsl_vector_set(gsl_pa, i, TRparam[0].paF);
        gsl_vector_set(gsl_pa_e, i, TRparam[0].paF_e);
        gsl_vector_set(gsl_incl, i, TRparam[0].inclF);
        gsl_vector_set(gsl_incl_e, i, TRparam[0].inclF_e);
        gsl_vector_set(gsl_vrot, i, TRparam[0].vrotF);
        gsl_vector_set(gsl_vrot_e, i, TRparam[0].vrotF_e);
        gsl_vector_set(gsl_vrad, i, TRparam[0].vradF);
        gsl_vector_set(gsl_vrad_e, i, TRparam[0].vradF_e);

        TRparam[0].maxLogLike[i] = TRparam[0].maxLogLikeF;
        TRparam[0].logZ[i] = TRparam[0].logZF;
        TRparam[0].logZerr[i] = TRparam[0].logZerrF;

        TRparam[0].npoints_inaring[i] = TRparam[0].Npoints_in_tilted_ring;
        TRparam[0].npoints_inaring_decim0[i] = TRparam[0].Npoints_in_tilted_ring_decim0;

        total_Npoints_allrings += TRparam[0].Npoints_in_tilted_ring;


        // Update geometry+flux-weighted error map : ring by ring
        // for the pixes in the current ring
        for(j=0; j<TRparam[0].Npoints_in_tilted_ring; j++)
        {
            x = TRparam[0].tilted_ring[j][0];
            y = TRparam[0].tilted_ring[j][1];

            if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma > 0) // if sigma_factor is not fitted and sigma_v is provided
            {
                HI_VF_boxFiltered_sigma_e_norm[0].data[y][x] = TRparam[0].e_sigma * TRparam[0].scale_factor_const_vlose_w / HI_VF_weight_TRfit[0].data[y][x];
            }
            else if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0)// full fit
            {
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[y][x]) \
                && !isnan(HI_VF_boxFiltered_sigma[0].data[y][x]) \
                && HI_VF_boxFiltered_sigma[0].data[y][x] > 0 \
                && !isinf(HI_VF_mom0[0].data[y][x]) \
                && !isnan(HI_VF_mom0[0].data[y][x]) \
                && HI_VF_mom0[0].data[y][x] > 0.0) // no blank
                {
                    HI_VF_boxFiltered_sigma_e_norm[0].data[y][x] = TRparam[0].sigma_factor * HI_VF_boxFiltered_sigma[0].data[y][x] * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[y][x] / HI_VF_mom0[0].data[y][x];
                }
            }
            else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0)// partial fit
            {
                if(TRparam[0].final_fit == 'Y')
                {
                    if(!isinf(HI_VF_mom0[0].data[y][x]) \
                    && !isnan(HI_VF_mom0[0].data[y][x]) \
                    && HI_VF_mom0[0].data[y][x] > 0.0 \
                    && !isnan(HI_VF_weight_TRfit[0].data[y][x]) \
                    && !isinf(HI_VF_weight_TRfit[0].data[y][x]))
                    {
                        //HI_VF_boxFiltered_sigma_e_norm[0].data[y][x] = TRparam[0].sigma_factor * HI_VF_boxFiltered_sigma[0].data[y][x] * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[y][x] / HI_VF_mom0[0].data[y][x];
                        HI_VF_boxFiltered_sigma_e_norm[0].data[y][x] = TRparam[0].sigma_factor * HI_VF_boxFiltered_sigma[0].data[y][x] * TRparam[0].scale_factor_var_vlose_w / HI_VF_mom0[0].data[y][x];
                    }
                }
                else if(TRparam[0].final_fit == 'N')
                {
                    if(!isinf(HI_VF_mom0[0].data[y][x]) \
                    && !isnan(HI_VF_mom0[0].data[y][x]) \
                    && !isnan(HI_VF_weight_TRfit[0].data[y][x]) \
                    && !isinf(HI_VF_weight_TRfit[0].data[y][x]) \
                    && HI_VF_mom0[0].data[y][x] > 0.0) // no blank
                    {
                        sigma_v = 0.1;
                        //HI_VF_boxFiltered_sigma_e_norm[0].data[y][x] = TRparam[0].sigma_factor * sigma_v * TRparam[0].scale_factor_var_vlose_w / HI_VF_weight_TRfit[0].data[y][x];
                        HI_VF_boxFiltered_sigma_e_norm[0].data[y][x] = TRparam[0].sigma_factor * sigma_v * TRparam[0].scale_factor_var_vlose_w;
                    }
                }
            }
        }
    }
    //free(Vrot_err_inaring);

    /* Save the total number of points in all rings */
    // update the total number of pixels available
    TRparam[0].total_Npoints_allRings = total_Npoints_allrings;

    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0; i<2*TRparam[0].Nrings; i++)
    {
        ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
        ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        if(i >= TRparam[0].Nrings)
        {
            TRparam[0].xpos0[i] = TRparam[0].xpos0[TRparam[0].Nrings-1];
            TRparam[0].ypos0[i] = TRparam[0].ypos0[TRparam[0].Nrings-1];
            TRparam[0].pa0[i] = TRparam[0].pa0[TRparam[0].Nrings-1];
            TRparam[0].incl0[i] = TRparam[0].incl0[TRparam[0].Nrings-1];
        }

        //define_tiltedRing_ellipse_rings_extend_rings(TRparam[0].xpos0[i], TRparam[0].ypos0[i], TRparam[0].pa0[i], TRparam[0].incl0[i], ri, ro, TRparam, side);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*
    for(i=0; i<200*TRparam[0].Nrings; i++)
    {
        _ring_w = 0.1;
        ri = 0 + i*_ring_w - 0.5*_ring_w;
        ro = 0 + i*_ring_w + 0.5*_ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        if(ring >= TRparam[0].ring_radius[TRparam[0].Nrings-1])
        {
            _xpos = TRparam[0].xpos0[TRparam[0].Nrings-1];
            _ypos = TRparam[0].ypos0[TRparam[0].Nrings-1];
            _pa = TRparam[0].pa0[TRparam[0].Nrings-1];
            _incl = TRparam[0].incl0[TRparam[0].Nrings-1];
        }

        find_Navail_Nall_pixels(_xpos, _ypos, _pa, _incl, ri, ro, TRparam, side);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    */

    // +++ Derive weights of individual parameters based on their normalised errors
    //FD_nbins = (int)(sqrt(TRparam[0].Nrings)+0.5);
    // XPOS
    if(xposfix == 'T')
    {
        // 0. gsl statistics
        gsl_mean_xpos_e = gsl_stats_mean(xpos_err, 1, TRparam[0].Nrings);
        gsl_std_xpos_e = sqrt(gsl_stats_variance(xpos_err, 1, TRparam[0].Nrings));
        gsl_min_xpos_e = gsl_stats_min(xpos_err, 1, TRparam[0].Nrings);
        gsl_max_xpos_e =  gsl_stats_max(xpos_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_xpos_e) || isinf(gsl_min_xpos_e) || gsl_min_xpos_e <= 0)
        if(isnan(gsl_min_xpos_e) || isinf(gsl_min_xpos_e) || gsl_min_xpos_e <= 0)
            gsl_min_xpos = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // XPOS
        robust_mean_std(xpos_tr, TRparam[0].Nrings, &hist_mean_xpos, &hist_std_xpos);
        // XPOS error
        robust_mean_std_e(xpos_err, TRparam[0].Nrings, &hist_mean_xpos_err, &hist_std_xpos_err);
    }

    // YPOS
    if(yposfix == 'T')
    {
        gsl_mean_ypos_e = gsl_stats_mean(ypos_err, 1, TRparam[0].Nrings);
        gsl_std_ypos_e = sqrt(gsl_stats_variance(ypos_err, 1, TRparam[0].Nrings));
        gsl_min_ypos_e = gsl_stats_min(ypos_err, 1, TRparam[0].Nrings);
        gsl_max_ypos_e =  gsl_stats_max(ypos_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_ypos_e) || isinf(gsl_min_ypos_e) || gsl_min_ypos_e <= 0)
            gsl_min_ypos_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // YPOS
        robust_mean_std(ypos_tr, TRparam[0].Nrings, &hist_mean_ypos, &hist_std_ypos);
        // YPOS error
        robust_mean_std_e(ypos_err, TRparam[0].Nrings, &hist_mean_ypos_err, &hist_std_ypos_err);
    }

    // VSYS
    if(vsysfix == 'T')
    {
        gsl_mean_vsys_e = gsl_stats_mean(vsys_err, 1, TRparam[0].Nrings);
        gsl_std_vsys_e = sqrt(gsl_stats_variance(vsys_err, 1, TRparam[0].Nrings));
        gsl_min_vsys_e = gsl_stats_min(vsys_err, 1, TRparam[0].Nrings);
        gsl_max_vsys_e =  gsl_stats_max(vsys_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vsys_e) || isinf(gsl_min_vsys_e) || gsl_min_vsys_e <= 0)
            gsl_min_vsys_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // VSYS
        robust_mean_std(vsys_tr, TRparam[0].Nrings, &hist_mean_vsys, &hist_std_vsys);
        // VSYS error
        robust_mean_std_e(vsys_err, TRparam[0].Nrings, &hist_mean_vsys_err, &hist_std_vsys_err);
    }

    // PA
    if(pafix == 'T')
    {
        gsl_mean_pa_e = gsl_stats_mean(pa_err, 1, TRparam[0].Nrings);
        gsl_std_pa_e = sqrt(gsl_stats_variance(pa_err, 1, TRparam[0].Nrings));
        gsl_min_pa_e = gsl_stats_min(pa_err, 1, TRparam[0].Nrings);
        gsl_max_pa_e = gsl_stats_max(pa_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_pa_e) || isinf(gsl_min_pa_e) || gsl_min_pa_e <= 0)
            gsl_min_pa_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // PA
        robust_mean_std(pa_tr, TRparam[0].Nrings, &hist_mean_pa, &hist_std_pa);
        // PA error
        robust_mean_std_e(pa_err, TRparam[0].Nrings, &hist_mean_pa_err, &hist_std_pa_err);
    }

    // INCL
    if(inclfix == 'T')
    {
        gsl_mean_incl_e = gsl_stats_mean(incl_err, 1, TRparam[0].Nrings);
        gsl_std_incl_e = sqrt(gsl_stats_variance(incl_err, 1, TRparam[0].Nrings));
        gsl_min_incl_e = gsl_stats_min(incl_err, 1, TRparam[0].Nrings);
        gsl_max_incl_e = gsl_stats_max(incl_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_incl_e) || isinf(gsl_min_incl_e) || gsl_min_incl_e <= 0)
            gsl_min_incl_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // INCL
        robust_mean_std(incl_tr, TRparam[0].Nrings, &hist_mean_incl, &hist_std_incl);
        // INCL error
        robust_mean_std_e(incl_err, TRparam[0].Nrings, &hist_mean_incl_err, &hist_std_incl_err);
    }

    // VROT
    if(vrotfix == 'T')
    {
        gsl_mean_vrot_e = gsl_stats_mean(vrot_err, 1, TRparam[0].Nrings);
        gsl_std_vrot_e = sqrt(gsl_stats_variance(vrot_err, 1, TRparam[0].Nrings));
        gsl_min_vrot_e = gsl_stats_min(vrot_err, 1, TRparam[0].Nrings);
        gsl_max_vrot_e = gsl_stats_max(vrot_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrot_e) || isinf(gsl_min_vrot_e) || gsl_min_vrot_e <= 0)
            gsl_min_vrot_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // VROT
        robust_mean_std(vrot_tr, TRparam[0].Nrings, &hist_mean_vrot, &hist_std_vrot);
        // VROT error
        robust_mean_std_e(vrot_err, TRparam[0].Nrings, &hist_mean_vrot_err, &hist_std_vrot_err);
    }

    // VRAD
    if(vradfix == 'T')
    {
        gsl_mean_vrad_e = gsl_stats_mean(vrad_err, 1, TRparam[0].Nrings);
        gsl_std_vrad_e = sqrt(gsl_stats_variance(vrad_err, 1, TRparam[0].Nrings));
        gsl_min_vrad_e = gsl_stats_min(vrad_err, 1, TRparam[0].Nrings);
        gsl_max_vrad_e = gsl_stats_max(vrad_err, 1, TRparam[0].Nrings);
        if(isnan(gsl_min_vrad_e) || isinf(gsl_min_vrad_e) || gsl_min_vrad_e <= 0)
            gsl_min_vrad_e = 1E3; // small value

        // 1. Derive individual errors : robust mode values based on histograms
        // VRAD
        robust_mean_std(vrad_tr, TRparam[0].Nrings, &hist_mean_vrad, &hist_std_vrad);
        // VRAD error
        robust_mean_std_e(vrad_err, TRparam[0].Nrings, &hist_mean_vrad_err, &hist_std_vrad_err);
    }

    // 3. Normalise the errors with respect to their minimums & derive the total errors
    // normalise errors with their mimimums to derive the total error normalised
    for(i=0; i<(int)TRparam[0].Nrings; i++) // this is to avoid the outer rotation velocities with large errors
    {
        // xpos
        if(xposfix == 'T')
        {
            if(TRparam[0].xpos_e[i] > 0)
                xpos_err[i] = TRparam[0].xpos_e[i]/gsl_min_xpos_e;
            else
                xpos_err[i] = 1E3; // large value
        }
        else
            xpos_err[i] = 0; // no error


        // ypos
        if(yposfix == 'T')
        {
            if(TRparam[0].ypos_e[i] > 0)
                ypos_err[i] = TRparam[0].ypos_e[i]/gsl_min_ypos_e;
            else
                ypos_err[i] = 1E3; // large value
        }
        else
            ypos_err[i] = 0; // no error


        // vsys
        if(vsysfix == 'T')
        {
            if(TRparam[0].vsys_e[i] > 0)
                vsys_err[i] = TRparam[0].vsys_e[i]/gsl_min_vsys_e;
            else
                vsys_err[i] = 1E3; // large value
        }
        else
            vsys_err[i] = 0; // no error

        // pa
        if(pafix == 'T')
        {
            if(TRparam[0].pa_e[i] > 0)
                pa_err[i] = TRparam[0].pa_e[i]/gsl_min_pa_e;
            else
                pa_err[i] = 1E3; // large value
        }
        else
            pa_err[i] = 0; // no error

        // incl
        if(inclfix == 'T')
        {
            if(TRparam[0].incl_e[i] > 0)
                incl_err[i] = TRparam[0].incl_e[i]/gsl_min_incl_e;
            else
                incl_err[i] = 1E3; // large value
        }
        else
            incl_err[i] = 0; // no error

        if(vrotfix == 'T')
        {
            if(TRparam[0].vrot_e[i] > 0)
                vrot_err[i] = TRparam[0].vrot_e[i]/gsl_min_vrot_e;
            else
                vrot_err[i] = 1E3;
        }
        else
            vrot_err[i] = 0; // no error

        if(vradfix == 'T')
        {
            if(TRparam[0].vrad_e[i] > 0)
                vrad_err[i] = TRparam[0].vrad_e[i]/gsl_min_vrad_e;
            else
                vrad_err[i] = 1E3;
        }
        else
            vrad_err[i] = 0; // no error


        // total errors (PA_e+INCL_e+VROT_e): normalised ones to their minimums
        total_error_ring_params[i] = sqrt(pow(xpos_err[i], 2) + \
                                            pow(ypos_err[i], 2) + \
                                            pow(vsys_err[i], 2) + \
                                            pow(pa_err[i], 2) + \
                                            pow(incl_err[i], 2) + \
                                            pow(vrot_err[i], 2) + \
                                            pow(vrad_err[i], 2));

        if(isinf(total_error_ring_params[i]) || isnan(total_error_ring_params[i]) || total_error_ring_params[i] <= 0)
            total_error_ring_params[i] = 1E3; // large value

        // this is for histogram analysis
        total_error_ring_params_temp[i] = total_error_ring_params[i];
    }

    // 4. Calculate mean & std of total error histograms to match with those of ring params
    // mean & std of total error histogram 
    robust_mean_std_e(total_error_ring_params_temp, TRparam[0].Nrings, &hist_mean_total_error_ring_params, &hist_std_total_error_ring_params);

    // Re-assign histogram with scaled total errors: see above
    for(i=0; i<(int)TRparam[0].Nrings; i++) // calculate error weighted mode of the ring parameters using their histograms
    {
        // xpos
        if(xposfix == 'T')
        {
            xpos_tr_err[i] = total_error_ring_params[i]*(hist_mean_xpos_err/hist_mean_total_error_ring_params);
        }

        // ypos
        if(yposfix == 'T')
        {
            ypos_tr_err[i] = total_error_ring_params[i]*(hist_mean_ypos_err/hist_mean_total_error_ring_params);
        }

        // vsys
        if(vsysfix == 'T')
        {
            vsys_tr_err[i] = total_error_ring_params[i]*(hist_mean_vsys_err/hist_mean_total_error_ring_params);
        }

        // pa
        if(pafix == 'T')
        {
            pa_tr_err[i] = total_error_ring_params[i]*(hist_mean_pa_err/hist_mean_total_error_ring_params);
        }

        // incl
        if(inclfix == 'T')
        {
            incl_tr_err[i] = total_error_ring_params[i]*(hist_mean_incl_err/hist_mean_total_error_ring_params);
        }

        // vrot
        if(vrotfix == 'T')
        {
            vrot_tr_err[i] = total_error_ring_params[i]*(hist_mean_vrot_err/hist_mean_total_error_ring_params);
        }

        // vrad
        if(vradfix == 'T')
        {
            vrad_tr_err[i] = total_error_ring_params[i]*(hist_mean_vrad_err/hist_mean_total_error_ring_params);
        }
    }


    // update priors based on the derived mode and sigma of individual ring params
    if(xposfix == 'T')
    {
        //robust_mean_std_histogram_ac(xpos_tr, xpos_tr_err, TRparam[0].Nrings, &hist_mean_xpos, &hist_std_xpos);
        robust_mean_std(xpos_tr, TRparam[0].Nrings, &hist_mean_xpos, &hist_std_xpos);
        TRparam[0].xpos1 = hist_mean_xpos - 10*hist_std_xpos;
        TRparam[0].xpos2 = hist_mean_xpos + 10*hist_std_xpos;
    }
    if(yposfix == 'T')
    {
        //robust_mean_std_histogram_ac(ypos_tr, ypos_tr_err, TRparam[0].Nrings, &hist_mean_ypos, &hist_std_ypos);
        robust_mean_std(ypos_tr, TRparam[0].Nrings, &hist_mean_ypos, &hist_std_ypos);
        TRparam[0].ypos1 = hist_mean_ypos - 10*hist_std_ypos;
        TRparam[0].ypos2 = hist_mean_ypos + 10*hist_std_ypos;
    }
    if(vsysfix == 'T')
    {
        //robust_mean_std_histogram_ac(vsys_tr, vsys_tr_err, TRparam[0].Nrings, &hist_mean_vsys, &hist_std_vsys);
        robust_mean_std(vsys_tr, TRparam[0].Nrings, &hist_mean_vsys, &hist_std_vsys);
        TRparam[0].vsys1 = hist_mean_vsys - 10*hist_std_vsys;
        TRparam[0].vsys2 = hist_mean_vsys + 10*hist_std_vsys;
    }
    if(pafix == 'T')
    {
        //robust_mean_std_histogram_ac(pa_tr, pa_tr_err, TRparam[0].Nrings, &hist_mean_pa, &hist_std_pa);
        robust_mean_std(pa_tr, TRparam[0].Nrings, &hist_mean_pa, &hist_std_pa);
        TRparam[0]._bspline_pa_hist_sigma = hist_std_pa;
    }
    if(inclfix == 'T')
    {
        //robust_mean_std_histogram_ac(incl_tr, incl_tr_err, TRparam[0].Nrings, &hist_mean_incl, &hist_std_incl);
        robust_mean_std(incl_tr, TRparam[0].Nrings, &hist_mean_incl, &hist_std_incl);
        TRparam[0]._bspline_incl_hist_sigma = hist_std_incl;
    }
    if(vradfix == 'T')
    {
        robust_mean_std(vrad_tr, TRparam[0].Nrings, &hist_mean_vrad, &hist_std_vrad);
        TRparam[0]._bspline_vrad_hist_sigma = hist_std_vrad;
    }


    /* calculate median value for XPOS if fitted */
    if(TRparam[0].xpos_fix == 'T')
    {
        TRparam[0].xposF = hist_mean_xpos; // error weighted mean of xpos histogram
        TRparam[0].xposF_EinastoFit = hist_mean_xpos; // update xposF_EinastoFit initially from ellipse fit
        TRparam[0].xposF_e = hist_std_xpos; // error weighted sigma of xpos histogram
    }

    /* calculate median value for YPOS if fitted */
    if(TRparam[0].ypos_fix == 'T')
    {
        TRparam[0].yposF = hist_mean_ypos; // error weighted mean of ypos histogram
        TRparam[0].yposF_EinastoFit = hist_mean_ypos; // update yposF_EinastoFit initially from ellipse fit
        TRparam[0].yposF_e = hist_std_ypos;
    }

    /* calculate median value for VSYS if fitted */
    if(TRparam[0].vsys_fix == 'T')
    {
        TRparam[0].vsysF = hist_mean_vsys; // error weighted mean of vsys histogram
        TRparam[0].vsysF_EinastoFit = hist_mean_vsys; // update vsysF_EinastoFit initially from ellipse fit
        TRparam[0].vsysF_e = hist_std_vsys;
    }

    /* calculate median value for PA if fitted */
    if(TRparam[0].pa_fix == 'T')
    {
        TRparam[0].paF = hist_mean_pa; // error weighted mean of pa histogram
        TRparam[0].paF_EinastoFit = hist_mean_pa; // update paF_EinastoFit initially from ellipse fit
        TRparam[0].paF_e = hist_std_pa;
    }
    /* calculate median value for INCL if fitted */
    if(TRparam[0].incl_fix == 'T')
    {
        TRparam[0].inclF = hist_mean_incl; // error weighted mean of incl histogram
        TRparam[0].inclF_EinastoFit = hist_mean_incl; // update inclF_EinastoFit initially from ellipse fit
        TRparam[0].inclF_e = hist_std_incl;
    }
    /* calculate median value for VRAD if fitted */
    if(TRparam[0].vrad_fix == 'T')
    {
        TRparam[0].vradF = hist_mean_vrad; // error weighted mean of vrad histogram
        TRparam[0].vradF_e = hist_std_vrad;
    }
    /* calculate median value for VROT if fitted */
    if(TRparam[0].vrot_fix == 'T')
    {
        vrot_median = hist_mean_vrot;
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

    gsl_vector_free(gsl_xpos);
    gsl_vector_free(gsl_xpos_e);
    gsl_vector_free(gsl_ypos);
    gsl_vector_free(gsl_ypos_e);
    gsl_vector_free(gsl_vsys);
    gsl_vector_free(gsl_vsys_e);
    gsl_vector_free(gsl_pa);
    gsl_vector_free(gsl_pa_e);
    gsl_vector_free(gsl_incl);
    gsl_vector_free(gsl_incl_e);
    gsl_vector_free(gsl_vrot);
    gsl_vector_free(gsl_vrot_e);
    gsl_vector_free(gsl_vrad);
    gsl_vector_free(gsl_vrad_e);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void trfit_multinest_trfit_rings_student(char *xpos, char xposfix, char *ypos, char yposfix, char *vsys, char vsysfix, char *pa, char pafix, char *incl, char inclfix, char *vrot, char vrotfix, char *vrad, char vradfix, char *sigmafactor, char sigmafactorfix, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, char *finalfit, char final_fit)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0, ii=0, j=0;
    int i0=0, j0=0;
    int x, y;
    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    int Nrings;
    double efr, tol, Ztol, logZero;
    double total_Npoints_allrings=0.;
    double ri=0., ro=0., ring;

    double ellipse_a, ellipse_b, xpos_prior_gab, ypos_prior_gab;
    //TRparam[0].Nrings = (int)((TRparam[0].ring_e-TRparam[0].ring_s)/TRparam[0].ring_w)+1;

    // dynamic total_error 1D array
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

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive_einasto_halofit;
    //if(nlive > 50) nlive = 50; // for a quick TRfit
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;
    maxiter = multinest_param[0].maxiter;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;
    int pWrap[ndims];

    //strcpy(root, multinest_param[0].root);
    //strncpy(root, multinest_param[0].root, strlen(multinest_param[0].root));
    strcpy(root, "trfit.trfit_rings.");

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    outfile = multinest_param[0].outfile;
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;

    TRparam[0].xpos_fix = xposfix;
    TRparam[0].ypos_fix = yposfix;
    TRparam[0].vsys_fix = vsysfix;
    TRparam[0].pa_fix = pafix;
    TRparam[0].incl_fix = inclfix;
    TRparam[0].vrot_fix = vrotfix;
    TRparam[0].vrad_fix = vradfix;
    TRparam[0].sigma_factor_fix = sigmafactorfix;
    TRparam[0].final_fit = final_fit;

    /* Set up ring parameters with currently derived values: note that pa0[i] and incl0[i] are given in degree */
    set_nfree_params_trfit_multinest_trfit_rings_student(TRparam);
    ndims = TRparam[0].n_freeParams;
    nPar = TRparam[0].n_freeParams;
    nClsPar = TRparam[0].n_freeParams;
    for(i = 0; i < ndims; i++) pWrap[i] = multinest_param[0].pWrap[i];

    /* Start tilted-ring fits */
    MPI_Barrier(MPI_COMM_WORLD);

    for(i=0; i<TRparam[0].Nrings; i++)
    {
        ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
        ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;
        //printf(" ");
        nlive = multinest_param[0].nlive_einasto_halofit;
        //if(nlive > 50) nlive = 50; // for a quick TRfit
        define_tiltedRing(TRparam[0].xpos0[i], TRparam[0].ypos0[i], TRparam[0].pa0[i], TRparam[0].incl0[i], ri, ro, TRparam, side);

        if(TRparam[0].xpos_fix == 'F')
        {
            // xpos for TiltedRingModel() function : see the function for more details
            TRparam[0].xposF = TRparam[0].xpos0[i]; // TRparam[0].xpos0[i] is TRparam[0].xposF_EinastoFit from Einasto halo fitting!
            TRparam[0].xposF_e = TRparam[0].xpos_e[i];
            TRparam[0].xpos[i] = TRparam[0].xpos0[i]; // Einasto halo fit result
        }

        if(TRparam[0].ypos_fix == 'F')
        {
            // ypos for TiltedRingModel() function : see the function for more details
            TRparam[0].yposF = TRparam[0].ypos0[i];
            TRparam[0].yposF_e = TRparam[0].ypos_e[i];
            TRparam[0].ypos[i] = TRparam[0].ypos0[i]; // Einasto halo fit result
        }

        if(TRparam[0].vsys_fix == 'F')
        {
            // vsys for TiltedRingModel() function : see the function for more details
            TRparam[0].vsysF = TRparam[0].vsys0[i];
            TRparam[0].vsysF_e = TRparam[0].vsys_e[i];
            TRparam[0].vsys[i] = TRparam[0].vsys0[i]; // Einasto halo fit result
        }

        if(TRparam[0].pa_fix == 'F')
        {
            // use the current PA model value for multinest if not fitted: in degree
            TRparam[0].paF = TRparam[0].pa0[i];
            TRparam[0].paF_e = TRparam[0].pa_e[i];
            TRparam[0].pa[i] = TRparam[0].pa0[i]; // in degree
            TRparam[0].pa_temp[i] = TRparam[0].pa0[i]; // in degree
        }

        if(TRparam[0].incl_fix == 'F')
        {
            // use the current INCL model value for multinest if not fitted: in degree
            TRparam[0].inclF = TRparam[0].incl0[i];
            TRparam[0].inclF_e = TRparam[0].incl_e[i];
            TRparam[0].incl[i] = TRparam[0].incl0[i]; // in degree 
            TRparam[0].incl_temp[i] = TRparam[0].incl0[i]; // in degree 
        }

        if(TRparam[0].vrad_fix == 'F')
        {
            // use the current VRAD model value for multinest if not fitted: in km/s
            TRparam[0].vradF = TRparam[0].vrad[i];
            TRparam[0].vradF_e = TRparam[0].vrad_e[i];
            TRparam[0].vrad[i] = TRparam[0].vrad[i];
            TRparam[0].vrad_temp[i] = TRparam[0].vrad[i];
        }

        // update xpos & ypos priors based on the current ring radius + pa0 + incl0
        ellipse_a = ri;
        xpos_prior_gab = fabs(ellipse_a*sin(TRparam[0].pa0[i]*M_PI/180.));
        ypos_prior_gab = fabs(ellipse_a*cos(TRparam[0].pa0[i]*M_PI/180.));

        if(TRparam[0].xpos_fix == 'T')
        {
            TRparam[0].xpos1 = TRparam[0].xpos0[i] - xpos_prior_gab;
            TRparam[0].xpos2 = TRparam[0].xpos0[i] + xpos_prior_gab;
        }
        if(TRparam[0].ypos_fix == 'T')
        {
            TRparam[0].ypos1 = TRparam[0].ypos0[i] - ypos_prior_gab;
            TRparam[0].ypos2 = TRparam[0].ypos0[i] + ypos_prior_gab;
        }
        if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0)// not neccessary for this version
        {
            TRparam[0].sigma_factor1 = 0;
            TRparam[0].sigma_factor2 = 200;
        }

        /* Calling multinest */ 
        MPI_Barrier(MPI_COMM_WORLD);
        run(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_trfit_student, dumper_TRfits, TRparam);   
        MPI_Barrier(MPI_COMM_WORLD);

        // XPOS
        TRparam[0].xpos[i] = TRparam[0].xposF; // Einasto halo fit result
        TRparam[0].xpos0[i] = TRparam[0].xposF; // Einasto halo fit result
        TRparam[0].xpos_e[i] = TRparam[0].xposF_e;
        xpos_tr[i] = TRparam[0].xpos[i];
        xpos_err[i] = TRparam[0].xpos_e[i];


        // YPOS
        TRparam[0].ypos[i] = TRparam[0].yposF; // Einasto halo fit result
        TRparam[0].ypos0[i] = TRparam[0].yposF; // Einasto halo fit result
        TRparam[0].ypos_e[i] = TRparam[0].yposF_e;
        ypos_tr[i] = TRparam[0].ypos[i];
        ypos_err[i] = TRparam[0].ypos_e[i];

        
        // VSYS
        TRparam[0].vsys[i] = TRparam[0].vsysF; // Einasto halo fit result
        TRparam[0].vsys0[i] = TRparam[0].vsysF; // Einasto halo fit result
        TRparam[0].vsys_e[i] = TRparam[0].vsysF_e;
        vsys_tr[i] = TRparam[0].vsys[i];
        vsys_err[i] = TRparam[0].vsys_e[i];

        // PA
        TRparam[0].pa[i] = TRparam[0].paF; // Einasto halo fit result
        TRparam[0].pa_temp[i] = TRparam[0].paF; // Einasto halo fit result
        TRparam[0].pa0[i] = TRparam[0].paF; // Einasto halo fit result
        TRparam[0].pa_e[i] = TRparam[0].paF_e;
        pa_tr[i] = TRparam[0].pa[i];
        pa_err[i] = TRparam[0].pa_e[i];

        // INCL
        TRparam[0].incl[i] = TRparam[0].inclF; // Einasto halo fit result
        TRparam[0].incl_temp[i] = TRparam[0].inclF; // Einasto halo fit result
        TRparam[0].incl0[i] = TRparam[0].inclF; // Einasto halo fit result
        TRparam[0].incl_e[i] = TRparam[0].inclF_e;
        incl_tr[i] = TRparam[0].incl[i];
        incl_err[i] = TRparam[0].incl_e[i];

        // VROT
        if(side == 999) // error propagation not used. Use vrotF_e from multinest fitting
        {   
            TRparam[0].vrot[i] = TRparam[0].vrotF;
            TRparam[0].vrot0[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrot_temp[i] = TRparam[0].vrotF;
            TRparam[0].vrot_temp_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad[i] = TRparam[0].vradF;
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
        else if(side == 0)
        {   
            TRparam[0].vrot[i] = TRparam[0].vrotF;
            TRparam[0].vrot0[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad[i] = TRparam[0].vradF;
            TRparam[0].vrad_temp[i] = TRparam[0].vradF;
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
        else if(side == -1)
        {   
            TRparam[0].vrot_rec[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e_rec[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad_rec[i] = TRparam[0].vradF;
            TRparam[0].vrad_rec_e[i] = TRparam[0].vradF_e;

            // this is for printing : TRfits_multinest_using_ISOfit_rings with both sides (0) should be done lastely
            TRparam[0].vrot[i] = TRparam[0].vrot_rec[i];
            TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad[i] = TRparam[0].vrad_rec[i];
            TRparam[0].vrad_e[i] = TRparam[0].vrad_rec_e[i];
        }
        else if(side == 1)
        {   
            TRparam[0].vrot_app[i] = TRparam[0].vrotF;
            TRparam[0].vrot_e_app[i] = TRparam[0].vrotF_e;

            TRparam[0].vrad_app[i] = TRparam[0].vradF;
            TRparam[0].vrad_app_e[i] = TRparam[0].vradF_e;

            // this is for printing : TRfits_multinest_using_ISOfit_rings with both sides (0) should be done lastely
            TRparam[0].vrot[i] = TRparam[0].vrot_app[i];
            TRparam[0].vrot_e[i] = TRparam[0].vrotF_e;
            TRparam[0].vrad[i] = TRparam[0].vrad_app[i];
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
        vrot_tr[i] = TRparam[0].vrot[i];
        vrot_err[i] = TRparam[0].vrot_e[i];

        vrad_tr[i] = TRparam[0].vrad[i];
        vrad_err[i] = TRparam[0].vrad_e[i];

        // e_sigma 
        TRparam[0].e_sigma_student_TR[i] = TRparam[0].e_sigma_tr; // e_sigma from TR fit
    }

    /* Save the total number of points in all rings */
    // update the total number of pixels available
    TRparam[0].total_Npoints_allRings = total_Npoints_allrings;

    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0; i<2*TRparam[0].Nrings; i++)
    {
        ri = TRparam[0].ring_s + i*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
        ro = TRparam[0].ring_s + i*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        if(i >= TRparam[0].Nrings)
        {
            TRparam[0].xpos0[i] = TRparam[0].xpos0[TRparam[0].Nrings-1];
            TRparam[0].ypos0[i] = TRparam[0].ypos0[TRparam[0].Nrings-1];
            TRparam[0].pa0[i] = TRparam[0].pa0[TRparam[0].Nrings-1];
            TRparam[0].incl0[i] = TRparam[0].incl0[TRparam[0].Nrings-1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


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
void set_nfree_params_trfit_multinest_trfit_rings_student(TR_ringParameters *TRparam)
{
    int i=0, j=0, k=0, p=0;
    int pa_nS=0;
    int incl_nS=0;
    int vrad_nS=0;
    int n_freeParams = 0;
    int bspline_order;
    double ring=0., ring_max=0.;
    double _ring_temp, _value_temp, _e_temp;

    int ncoeffs, n;
    int nbreak;
    double xi, yi, yerr, Bj;
    gsl_bspline_workspace *_bw_pa, *_bw_incl, *_bw_vrad;
    gsl_vector *_B_pa, *_B_incl, *_B_vrad;
    gsl_vector *_c_pa, *_c_incl, *_c_vrad, *_w_pa, *_w_incl, *_w_vrad;
    gsl_vector *_x_pa, *_x_incl, *_x_vrad, *_y_pa, *_y_incl, *_y_vrad;
    gsl_matrix *_X_pa, *_X_incl, *_X_vrad, *_cov_pa, *_cov_incl, *_cov_vrad;
    gsl_multifit_linear_workspace *_mw_pa, *_mw_incl, *_mw_vrad;

    n = TRparam[0].Nrings;
    // dynamic 1D array
    double *x_dat = malloc(sizeof(double) * n);
    double *y_dat = malloc(sizeof(double) * n);
    double *e_dat = malloc(sizeof(double) * n);

    double pa_e_upper_bound;
    double pa_e_lower_bound;
    double incl_e_upper_bound;
    double incl_e_lower_bound;
    double vrad_e_upper_bound;
    double vrad_e_lower_bound;

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

    // SAVE FIT RESULTS FROM EINASTO FIT WHICH ARE USED FOR UPDATING GAUSSIAN PRIORS 

    // XPOS
    TRparam[0].xposF_EinastoFit_t = TRparam[0].xposF_EinastoFit;
    // YPOS
    TRparam[0].yposF_EinastoFit_t = TRparam[0].yposF_EinastoFit;
    // VSYS
    TRparam[0].vsysF_EinastoFit_t = TRparam[0].vsysF_EinastoFit;
    // PA
    for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++) // bspline
    {
        TRparam[0]._p_bs_t[pa_nS] = TRparam[0]._p_bs[pa_nS];
        TRparam[0]._p_bs_e_t[pa_nS] = TRparam[0]._p_bs_e[pa_nS];
    }
    // INCL
    for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++) // bspline
    {
        TRparam[0]._i_bs_t[incl_nS] = TRparam[0]._i_bs[incl_nS];
        TRparam[0]._i_bs_e_t[incl_nS] = TRparam[0]._i_bs_e[incl_nS];
    }
    // _n
    TRparam[0]._n_t = TRparam[0]._n;
    TRparam[0]._ne_t = TRparam[0]._ne;

    // r_2
    TRparam[0]._r_2_t = TRparam[0].r_2;
    TRparam[0]._r_2e_t = TRparam[0].r_2e;

    // rho_2
    TRparam[0]._rho_2_t = TRparam[0].rho_2;
    TRparam[0]._rho_2e_t = TRparam[0].rho_2e;

    // VRAD
    for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++) // bspline
    {
        TRparam[0]._vr_bs_t[vrad_nS] = TRparam[0]._vr_bs[vrad_nS];
        TRparam[0]._vr_bs_e_t[vrad_nS] = TRparam[0]._vr_bs_e[vrad_nS];
    }


    // 1. XPOS
    //TRparam[0].xpos_fix = xposfix; 
    if(TRparam[0].xpos_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }
    // 2. YPOS
    //TRparam[0].ypos_fix = yposfix;
    if(TRparam[0].ypos_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }
    // 3. VSYS
    //TRparam[0].vsys_fix = vsysfix;
    if(TRparam[0].vsys_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }
    // 4. PA
    //TRparam[0].pa_fix = pafix;
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++) // bspline
                n_freeParams++;
        }
        else if(strcmp(TRparam[0].pa_function, "const") == 0)
        {
            n_freeParams++;
        }
    }   

    // 5. INCL
    //TRparam[0].incl_fix = inclfix;
    if(TRparam[0].incl_fix == 'T')
    {
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++) // bspline
                n_freeParams++;
        }
        else if(strcmp(TRparam[0].incl_function, "const") == 0)
        {
            n_freeParams++;
        }
    }   

    // 6. VROT
    //TRparam[0].vrot_fix = vrotfix;
    if(TRparam[0].vrot_fix == 'T')
    {
        n_freeParams++;
    }

    // 7. VRAD
    //TRparam[0].vrad_fix = vradfix;
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

    // 8. sigma_factor
    //TRparam[0].sigma_factor_fix = sigmafactorfix;
    /*
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].final_fit == 'Y')
    {
        n_freeParams++;
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && TRparam[0].final_fit == 'Y')
    {
        n_freeParams++;
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && TRparam[0].final_fit == 'N')
    {
        n_freeParams++;
    }
    */

    if(TRparam[0].e_sigma >= 0)
    {
        n_freeParams++;
    }

    // Update the total number of free tilted-ring parameters for TR fits
    TRparam[0].n_freeParams = n_freeParams;

    // set vector(x, y) for PA b-spline
    if(strcmp(TRparam[0].pa_function, "bspline") == 0)
    {
        for(i=0; i<n; i++)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].pa_temp[i]/TRparam[0].PA_MAX_in_degree;
            _e_temp = TRparam[0].pa_temp_e[i]/TRparam[0].PA_MAX_in_degree;

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = _e_temp;

            gsl_vector_set (_x_pa, i, x_dat[i]);
            gsl_vector_set (_y_pa, i, y_dat[i]);
            gsl_vector_set (_w_pa, i, 1/e_dat[i]);
        }
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_pa);
    }

    // set vector(x, y) for INCL b-spline
    if(strcmp(TRparam[0].incl_function, "bspline") == 0)
    {
        for(i=0; i<n; i++)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].incl_temp[i]/TRparam[0].INCL_MAX_in_degree;
            _e_temp = TRparam[0].incl_temp_e[i]/TRparam[0].INCL_MAX_in_degree;

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = _e_temp;
            gsl_vector_set (_x_incl, i, x_dat[i]);
            gsl_vector_set (_y_incl, i, y_dat[i]);
            gsl_vector_set (_w_incl, i, 1/e_dat[i]);
        }
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_incl);
    }

    // set vector(x, y) for VRAD b-spline
    if(strcmp(TRparam[0].vrad_function, "bspline") == 0 && TRparam[0].vrad_nbreak_bspline > 1 && TRparam[0].vrad_order_bspline >= 0)
    {
        for(i=0; i<n; i++)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].vrad_temp[i]/TRparam[0].vrad_max;
            _e_temp = TRparam[0].vrad_temp_e[i]/TRparam[0].vrad_max;

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = _e_temp;
            gsl_vector_set (_x_vrad, i, x_dat[i]);
            gsl_vector_set (_y_vrad, i, y_dat[i]);
            gsl_vector_set (_w_vrad, i, 1/e_dat[i]);
        }
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_vrad);
    }

    for(i=0; i<TRparam[0].Nrings; i++)
    {
        ring = TRparam[0].ring_s + (double)i*TRparam[0].ring_w;

        // 1. XPOS
        TRparam[0].xpos0[i] = TRparam[0].xposF_EinastoFit;
        TRparam[0].xpos_e[i] = TRparam[0].xposF_EinastoFit_e;
        TRparam[0].xposF_EinastoFit_t = TRparam[0].xposF_EinastoFit; // this is used for updating priors later

        // 2. YPOS
        TRparam[0].ypos0[i] = TRparam[0].yposF_EinastoFit;
        TRparam[0].ypos_e[i] = TRparam[0].yposF_EinastoFit_e;
        TRparam[0].yposF_EinastoFit_t = TRparam[0].yposF_EinastoFit; // this is used for updating priors later

        // 3. VSYS
        TRparam[0].vsys0[i] = TRparam[0].vsysF_EinastoFit;
        TRparam[0].vsys_e[i] = TRparam[0].vsysF_EinastoFit_e;
        TRparam[0].vsysF_EinastoFit_t = TRparam[0].vsysF_EinastoFit; // this is used for updating priors later

        // 4. PA
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_pa);
            /* construct the fit matrix _X */
            for (k=0; k<n; k++)
            {
                xi = gsl_vector_get(_x_pa, k);
                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_pa, _bw_pa);

                /* fill in row i of _X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_pa; j++)
                {
                    //Bj = gsl_vector_get(_B, j);
                    Bj = 0.01;
                    gsl_matrix_set(_X_pa, k, j, Bj);
                }
            }
            /* construct the fit matrix _cov */
            for (k=0; k<TRparam[0].n_coeffs_bspline_pa; k++)
            {
                xi = gsl_vector_get(_x_pa, k);
                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_pa, _bw_pa);

                /* fill in row i of X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_pa; j++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B, _j);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_pa, k, j, Bj);
                }
            }

            // PA
            xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _c_pa->data[pa_nS] = TRparam[0]._p_bs[pa_nS];
            }
            gsl_bspline_eval(xi, _B_pa, _bw_pa);
            gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            TRparam[0].pa0[i] = yi*TRparam[0].PA_MAX_in_degree;

            // 1. PA_e : upper bound
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _c_pa->data[pa_nS] = TRparam[0]._p_bs[pa_nS] + TRparam[0]._p_bs_e[pa_nS];
            }
            gsl_bspline_eval(xi, _B_pa, _bw_pa);
            gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            pa_e_upper_bound = yi*TRparam[0].PA_MAX_in_degree;

            // 2. PA_e : lower bound
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _c_pa->data[pa_nS] = TRparam[0]._p_bs[pa_nS] - TRparam[0]._p_bs_e[pa_nS];
            }
            gsl_bspline_eval(xi, _B_pa, _bw_pa);
            gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            pa_e_lower_bound = yi*TRparam[0].PA_MAX_in_degree;
            TRparam[0].pa_e[i] = (pa_e_upper_bound - pa_e_lower_bound)/2.0;

            // PA_e
            //xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            //if(xi < x_dat[0]) xi = x_dat[0];
            //if(xi > x_dat[n-1]) xi = x_dat[n-1];
            //for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            //{
            //    _c_pa->data[pa_nS] = TRparam[0]._p_bs_e[pa_nS];
            //}
            //gsl_bspline_eval(xi, _B_pa, _bw_pa);
            //gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            //TRparam[0].pa_e[i] = yi*TRparam[0].PA_MAX_in_degree;
        }
        else if(strcmp(TRparam[0].pa_function, "const") == 0)
        {
            TRparam[0].pa0[i] = TRparam[0].paF_EinastoFit; // this is to avoid a simple use of ellipse fit without receding or approaching information.
            TRparam[0].pa_e[i] = TRparam[0].paF_EinastoFit_e; // this is to avoid a simple use of ellipse fit without receding or approaching information.
        }


        // INCL
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_incl);
            /* construct the fit matrix _X */
            for (k=0; k<n; k++)
            {
                xi = gsl_vector_get(_x_incl, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_incl, _bw_incl);

                /* fill in row i of _X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_incl; j++)
                {
                    //Bj = gsl_vector_get(_B, j);
                    Bj = 0.01;
                    gsl_matrix_set(_X_incl, k, j, Bj);
                }
            }

            /* construct the fit matrix _cov */
            for (k=0; k<TRparam[0].n_coeffs_bspline_incl; k++)
            {
                xi = gsl_vector_get(_x_incl, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_incl, _bw_incl);

                /* fill in row i of X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_incl; j++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B, _j);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_incl, k, j, Bj);
                }
            }

            xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];
            // INCL
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _c_incl->data[incl_nS] = TRparam[0]._i_bs[incl_nS];
            }
            gsl_bspline_eval(xi, _B_incl, _bw_incl);
            gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
            TRparam[0].incl0[i] = yi*TRparam[0].INCL_MAX_in_degree;

            // 1. INCL_e : upper bound
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _c_incl->data[incl_nS] = TRparam[0]._i_bs[incl_nS] + TRparam[0]._i_bs_e[incl_nS];
            }
            gsl_bspline_eval(xi, _B_incl, _bw_incl);
            gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            incl_e_upper_bound = yi*TRparam[0].INCL_MAX_in_degree;

            // 2. INCL_e : lower bound
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _c_incl->data[incl_nS] = TRparam[0]._i_bs[incl_nS] - TRparam[0]._i_bs_e[incl_nS];
            }
            gsl_bspline_eval(xi, _B_incl, _bw_incl);
            gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
            // compute the lower bound of uncertainty from the fit
            incl_e_lower_bound = yi*TRparam[0].INCL_MAX_in_degree;
            TRparam[0].incl_e[i] = (incl_e_upper_bound - incl_e_lower_bound)/2.0;
        }
        else if(strcmp(TRparam[0].incl_function, "const") == 0)
        {
            TRparam[0].incl0[i] = TRparam[0].inclF_EinastoFit; // this is to avoid a simple use of ellipse fit without receding or approaching information.
            TRparam[0].incl_e[i] = TRparam[0].inclF_EinastoFit_e; // this is to avoid a simple use of ellipse fit without receding or approaching information.
        }


        // VRAD
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0 && TRparam[0].vrad_fix == 'T' && TRparam[0].vrad_nbreak_bspline > 1 && TRparam[0].vrad_order_bspline >= 0)
        {
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_vrad);
            /* construct the fit matrix _X */
            for (k=0; k<n; k++)
            {
                xi = gsl_vector_get(_x_vrad, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_vrad, _bw_vrad);

                /* fill in row i of _X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_vrad; j++)
                {
                    //Bj = gsl_vector_get(_B, j);
                    Bj = 0.01;
                    gsl_matrix_set(_X_vrad, k, j, Bj);
                }
            }

            /* construct the fit matrix _cov */
            for (k=0; k<TRparam[0].n_coeffs_bspline_vrad; k++)
            {
                xi = gsl_vector_get(_x_vrad, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_vrad, _bw_vrad);

                /* fill in row i of X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_vrad; j++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B, _j);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_vrad, k, j, Bj);
                }
            }

            xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];
            // VRAD
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _c_vrad->data[vrad_nS] = TRparam[0]._vr_bs[vrad_nS];
            }
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);
            gsl_multifit_linear_est(_B_vrad, _c_vrad, _cov_vrad, &yi, &yerr);
            TRparam[0].vrad_einastofit_bs[i] = yi*TRparam[0].vrad_max;

            // 1. VRAD_e : upper bound
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _c_vrad->data[vrad_nS] = TRparam[0]._vr_bs[vrad_nS] + TRparam[0]._vr_bs_e[vrad_nS];
            }
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);
            gsl_multifit_linear_est(_B_vrad, _c_vrad, _cov_vrad, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            vrad_e_upper_bound = yi*TRparam[0].vrad_max;

            // 2. VRAD_e : lower bound
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _c_vrad->data[vrad_nS] = TRparam[0]._vr_bs[vrad_nS] - TRparam[0]._vr_bs_e[vrad_nS];
            }
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);
            gsl_multifit_linear_est(_B_vrad, _c_vrad, _cov_vrad, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            vrad_e_lower_bound = yi*TRparam[0].vrad_max;
            TRparam[0].vrad_einastofit_bs_e[i] = (vrad_e_upper_bound - vrad_e_lower_bound)/2.0;
        }
        else if(strcmp(TRparam[0].vrad_function, "const") == 0)
        {
            TRparam[0].vrad[i] = TRparam[0].vradF;
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
    }
    free(x_dat);
    free(y_dat);
    free(e_dat);

    if(strcmp(TRparam[0].pa_function, "bspline") == 0)
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

    if(strcmp(TRparam[0].incl_function, "bspline") == 0)
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

    if(strcmp(TRparam[0].vrad_function, "bspline") == 0 && TRparam[0].vrad_nbreak_bspline > 1 && TRparam[0].vrad_order_bspline >= 0)
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
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void set_nfree_params_trfit_multinest_trfit_rings_student_test(TR_ringParameters *TRparam)
{
    int i=0, j=0, k=0, p=0;
    int pa_nS=0;
    int incl_nS=0;
    int vrad_nS=0;
    int n_freeParams = 0;
    int bspline_order;
    double ring=0., ring_max=0.;
    double _ring_temp, _value_temp, _e_temp;

    int ncoeffs, n;
    int nbreak;
    double xi, yi, yerr, Bj;
    gsl_bspline_workspace *_bw_pa, *_bw_incl, *_bw_vrad;
    gsl_vector *_B_pa, *_B_incl, *_B_vrad;
    gsl_vector *_c_pa, *_c_incl, *_c_vrad, *_w_pa, *_w_incl, *_w_vrad;
    gsl_vector *_x_pa, *_x_incl, *_x_vrad, *_y_pa, *_y_incl, *_y_vrad;
    gsl_matrix *_X_pa, *_X_incl, *_X_vrad, *_cov_pa, *_cov_incl, *_cov_vrad;
    gsl_multifit_linear_workspace *_mw_pa, *_mw_incl, *_mw_vrad;

    n = TRparam[0].Nrings;
    // dynamic 1D array
    double *x_dat = malloc(sizeof(double) * n);
    double *y_dat = malloc(sizeof(double) * n);
    double *e_dat = malloc(sizeof(double) * n);

    double pa_e_upper_bound;
    double pa_e_lower_bound;
    double incl_e_upper_bound;
    double incl_e_lower_bound;
    double vrad_e_upper_bound;
    double vrad_e_lower_bound;

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
    if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
    {
        TRparam[0].vrad_nbreak_bspline = 2;
        TRparam[0].vrad_order_bspline = 0;
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

    // SAVE FIT RESULTS FROM EINASTO FIT WHICH ARE USED FOR UPDATING GAUSSIAN PRIORS 

    // XPOS
    TRparam[0].xposF_EinastoFit_t = TRparam[0].xposF_EinastoFit;
    // YPOS
    TRparam[0].yposF_EinastoFit_t = TRparam[0].yposF_EinastoFit;
    // VSYS
    TRparam[0].vsysF_EinastoFit_t = TRparam[0].vsysF_EinastoFit;
    // PA
    for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++) // bspline
    {
        TRparam[0]._p_bs_t[pa_nS] = TRparam[0]._p_bs[pa_nS];
        TRparam[0]._p_bs_e_t[pa_nS] = TRparam[0]._p_bs_e[pa_nS];
    }
    // INCL
    for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++) // bspline
    {
        TRparam[0]._i_bs_t[incl_nS] = TRparam[0]._i_bs[incl_nS];
        TRparam[0]._i_bs_e_t[incl_nS] = TRparam[0]._i_bs_e[incl_nS];
    }
    // _n
    TRparam[0]._n_t = TRparam[0]._n;
    TRparam[0]._ne_t = TRparam[0]._ne;

    // r_2
    TRparam[0]._r_2_t = TRparam[0].r_2;
    TRparam[0]._r_2e_t = TRparam[0].r_2e;

    // rho_2
    TRparam[0]._rho_2_t = TRparam[0].rho_2;
    TRparam[0]._rho_2e_t = TRparam[0].rho_2e;

    // VRAD
    for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++) // bspline
    {
        TRparam[0]._vr_bs_t[vrad_nS] = TRparam[0]._vr_bs[vrad_nS];
        TRparam[0]._vr_bs_e_t[vrad_nS] = TRparam[0]._vr_bs_e[vrad_nS];
    }


    // 1. XPOS
    //TRparam[0].xpos_fix = xposfix; 
    if(TRparam[0].xpos_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }
    // 2. YPOS
    //TRparam[0].ypos_fix = yposfix;
    if(TRparam[0].ypos_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }
    // 3. VSYS
    //TRparam[0].vsys_fix = vsysfix;
    if(TRparam[0].vsys_fix == 'T') // update the pre value with the current one only if newly fitted
    {
        n_freeParams++;
    }
    // 4. PA
    //TRparam[0].pa_fix = pafix;
    if(TRparam[0].pa_fix == 'T')
    {
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++) // bspline
                n_freeParams++;
        }
        else if(strcmp(TRparam[0].pa_function, "const") == 0)
        {
            n_freeParams++;
        }
    }   

    // 5. INCL
    //TRparam[0].incl_fix = inclfix;
    if(TRparam[0].incl_fix == 'T')
    {
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++) // bspline
                n_freeParams++;
        }
        else if(strcmp(TRparam[0].incl_function, "const") == 0)
        {
            n_freeParams++;
        }
    }   

    // 6. VROT
    //TRparam[0].vrot_fix = vrotfix;
    if(TRparam[0].vrot_fix == 'T')
    {
        n_freeParams++;
    }

    // 7. VRAD
    //TRparam[0].vrad_fix = vradfix;
    if(TRparam[0].vrad_fix == 'F')
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

    // 8. sigma_factor
    //TRparam[0].sigma_factor_fix = sigmafactorfix;
    /*
    if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].final_fit == 'Y')
    {
        n_freeParams++;
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && TRparam[0].final_fit == 'Y')
    {
        n_freeParams++;
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && TRparam[0].final_fit == 'N')
    {
        n_freeParams++;
    }
    */

    if(TRparam[0].e_sigma >= 0)
    {
        n_freeParams++;
    }

    // Update the total number of free tilted-ring parameters for TR fits
    TRparam[0].n_freeParams = n_freeParams;

    // set vector(x, y) for PA b-spline
    if(strcmp(TRparam[0].pa_function, "bspline") == 0)
    {
        for(i=0; i<n; i++)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].pa_temp[i]/TRparam[0].PA_MAX_in_degree;
            _e_temp = TRparam[0].pa_temp_e[i]/TRparam[0].PA_MAX_in_degree;

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = _e_temp;

            gsl_vector_set (_x_pa, i, x_dat[i]);
            gsl_vector_set (_y_pa, i, y_dat[i]);
            gsl_vector_set (_w_pa, i, 1/e_dat[i]);
        }
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_pa);
    }

    // set vector(x, y) for INCL b-spline
    if(strcmp(TRparam[0].incl_function, "bspline") == 0)
    {
        for(i=0; i<n; i++)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].incl_temp[i]/TRparam[0].INCL_MAX_in_degree;
            _e_temp = TRparam[0].incl_temp_e[i]/TRparam[0].INCL_MAX_in_degree;

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = _e_temp;
            gsl_vector_set (_x_incl, i, x_dat[i]);
            gsl_vector_set (_y_incl, i, y_dat[i]);
            gsl_vector_set (_w_incl, i, 1/e_dat[i]);
        }
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_incl);
    }

    // set vector(x, y) for VRAD b-spline
    if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
    {
        for(i=0; i<n; i++)
        {
            _ring_temp = TRparam[0].ring_radius[i]/TRparam[0].rGalaxyPlane_pixel_max; // from TR rings
            _value_temp = TRparam[0].vrad_temp[i]/TRparam[0].vrad_max;
            _e_temp = TRparam[0].vrad_temp_e[i]/TRparam[0].vrad_max;

            x_dat[i] = _ring_temp;
            y_dat[i] = _value_temp;
            e_dat[i] = _e_temp;
            gsl_vector_set (_x_vrad, i, x_dat[i]);
            gsl_vector_set (_y_vrad, i, y_dat[i]);
            gsl_vector_set (_w_vrad, i, 1/e_dat[i]);
        }
        gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_vrad);
    }

    for(i=0; i<TRparam[0].Nrings; i++)
    {
        ring = TRparam[0].ring_s + (double)i*TRparam[0].ring_w;

        // 1. XPOS
        TRparam[0].xpos0[i] = TRparam[0].xposF_EinastoFit;
        TRparam[0].xpos_e[i] = TRparam[0].xposF_EinastoFit_e;
        TRparam[0].xposF_EinastoFit_t = TRparam[0].xposF_EinastoFit; // this is used for updating priors later

        // 2. YPOS
        TRparam[0].ypos0[i] = TRparam[0].yposF_EinastoFit;
        TRparam[0].ypos_e[i] = TRparam[0].yposF_EinastoFit_e;
        TRparam[0].yposF_EinastoFit_t = TRparam[0].yposF_EinastoFit; // this is used for updating priors later

        // 3. VSYS
        TRparam[0].vsys0[i] = TRparam[0].vsysF_EinastoFit;
        TRparam[0].vsys_e[i] = TRparam[0].vsysF_EinastoFit_e;
        TRparam[0].vsysF_EinastoFit_t = TRparam[0].vsysF_EinastoFit; // this is used for updating priors later

        // 4. PA
        if(strcmp(TRparam[0].pa_function, "bspline") == 0)
        {
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_pa);
            /* construct the fit matrix _X */
            for (k=0; k<n; k++)
            {
                xi = gsl_vector_get(_x_pa, k);
                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_pa, _bw_pa);

                /* fill in row i of _X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_pa; j++)
                {
                    //Bj = gsl_vector_get(_B, j);
                    Bj = 0.01;
                    gsl_matrix_set(_X_pa, k, j, Bj);
                }
            }
            /* construct the fit matrix _cov */
            for (k=0; k<TRparam[0].n_coeffs_bspline_pa; k++)
            {
                xi = gsl_vector_get(_x_pa, k);
                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_pa, _bw_pa);

                /* fill in row i of X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_pa; j++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B, _j);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_pa, k, j, Bj);
                }
            }

            // PA
            xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _c_pa->data[pa_nS] = TRparam[0]._p_bs[pa_nS];
            }
            gsl_bspline_eval(xi, _B_pa, _bw_pa);
            gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            TRparam[0].pa0[i] = yi*TRparam[0].PA_MAX_in_degree;

            // 1. PA_e : upper bound
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _c_pa->data[pa_nS] = TRparam[0]._p_bs[pa_nS] + TRparam[0]._p_bs_e[pa_nS];
            }
            gsl_bspline_eval(xi, _B_pa, _bw_pa);
            gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            pa_e_upper_bound = yi*TRparam[0].PA_MAX_in_degree;

            // 2. PA_e : lower bound
            for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            {
                _c_pa->data[pa_nS] = TRparam[0]._p_bs[pa_nS] - TRparam[0]._p_bs_e[pa_nS];
            }
            gsl_bspline_eval(xi, _B_pa, _bw_pa);
            gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            pa_e_lower_bound = yi*TRparam[0].PA_MAX_in_degree;
            TRparam[0].pa_e[i] = (pa_e_upper_bound - pa_e_lower_bound)/2.0;

            // PA_e
            //xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            //if(xi < x_dat[0]) xi = x_dat[0];
            //if(xi > x_dat[n-1]) xi = x_dat[n-1];
            //for(pa_nS=0; pa_nS<TRparam[0].n_coeffs_bspline_pa; pa_nS++)
            //{
            //    _c_pa->data[pa_nS] = TRparam[0]._p_bs_e[pa_nS];
            //}
            //gsl_bspline_eval(xi, _B_pa, _bw_pa);
            //gsl_multifit_linear_est(_B_pa, _c_pa, _cov_pa, &yi, &yerr);
            //TRparam[0].pa_e[i] = yi*TRparam[0].PA_MAX_in_degree;
        }
        else if(strcmp(TRparam[0].pa_function, "const") == 0)
        {
            TRparam[0].pa0[i] = TRparam[0].paF_EinastoFit; // this is to avoid a simple use of ellipse fit without receding or approaching information.
            TRparam[0].pa_e[i] = TRparam[0].paF_EinastoFit_e; // this is to avoid a simple use of ellipse fit without receding or approaching information.
        }


        // INCL
        if(strcmp(TRparam[0].incl_function, "bspline") == 0)
        {
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_incl);
            /* construct the fit matrix _X */
            for (k=0; k<n; k++)
            {
                xi = gsl_vector_get(_x_incl, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_incl, _bw_incl);

                /* fill in row i of _X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_incl; j++)
                {
                    //Bj = gsl_vector_get(_B, j);
                    Bj = 0.01;
                    gsl_matrix_set(_X_incl, k, j, Bj);
                }
            }

            /* construct the fit matrix _cov */
            for (k=0; k<TRparam[0].n_coeffs_bspline_incl; k++)
            {
                xi = gsl_vector_get(_x_incl, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_incl, _bw_incl);

                /* fill in row i of X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_incl; j++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B, _j);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_incl, k, j, Bj);
                }
            }

            xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];
            // INCL
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _c_incl->data[incl_nS] = TRparam[0]._i_bs[incl_nS];
            }
            gsl_bspline_eval(xi, _B_incl, _bw_incl);
            gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
            TRparam[0].incl0[i] = yi*TRparam[0].INCL_MAX_in_degree;

            // 1. INCL_e : upper bound
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _c_incl->data[incl_nS] = TRparam[0]._i_bs[incl_nS] + TRparam[0]._i_bs_e[incl_nS];
            }
            gsl_bspline_eval(xi, _B_incl, _bw_incl);
            gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            incl_e_upper_bound = yi*TRparam[0].INCL_MAX_in_degree;

            // 2. INCL_e : lower bound
            for(incl_nS=0; incl_nS<TRparam[0].n_coeffs_bspline_incl; incl_nS++)
            {
                _c_incl->data[incl_nS] = TRparam[0]._i_bs[incl_nS] - TRparam[0]._i_bs_e[incl_nS];
            }
            gsl_bspline_eval(xi, _B_incl, _bw_incl);
            gsl_multifit_linear_est(_B_incl, _c_incl, _cov_incl, &yi, &yerr);
            // compute the lower bound of uncertainty from the fit
            incl_e_lower_bound = yi*TRparam[0].INCL_MAX_in_degree;
            TRparam[0].incl_e[i] = (incl_e_upper_bound - incl_e_lower_bound)/2.0;
        }
        else if(strcmp(TRparam[0].incl_function, "const") == 0)
        {
            TRparam[0].incl0[i] = TRparam[0].inclF_EinastoFit; // this is to avoid a simple use of ellipse fit without receding or approaching information.
            TRparam[0].incl_e[i] = TRparam[0].inclF_EinastoFit_e; // this is to avoid a simple use of ellipse fit without receding or approaching information.
        }


        // VRAD
        if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
        {
            gsl_bspline_knots_uniform(x_dat[0], x_dat[n-1], _bw_vrad);
            /* construct the fit matrix _X */
            for (k=0; k<n; k++)
            {
                xi = gsl_vector_get(_x_vrad, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_vrad, _bw_vrad);

                /* fill in row i of _X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_vrad; j++)
                {
                    //Bj = gsl_vector_get(_B, j);
                    Bj = 0.01;
                    gsl_matrix_set(_X_vrad, k, j, Bj);
                }
            }

            /* construct the fit matrix _cov */
            for (k=0; k<TRparam[0].n_coeffs_bspline_vrad; k++)
            {
                xi = gsl_vector_get(_x_vrad, k);

                /* compute B_j(xi) for all j */
                gsl_bspline_eval(xi, _B_vrad, _bw_vrad);

                /* fill in row i of X */
                for (j=0; j<TRparam[0].n_coeffs_bspline_vrad; j++)
                {
                    // Bj actually doesn't matter much for now as _cov is saying out errors in coefficients _c
                    //Bj = gsl_vector_get(_B, _j);
                    Bj = 0.01;
                    gsl_matrix_set(_cov_vrad, k, j, Bj);
                }
            }

            xi = ring/TRparam[0].rGalaxyPlane_pixel_max;
            if(xi < x_dat[0]) xi = x_dat[0];
            if(xi > x_dat[n-1]) xi = x_dat[n-1];
            // VRAD
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _c_vrad->data[vrad_nS] = TRparam[0]._vr_bs[vrad_nS];
            }
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);
            gsl_multifit_linear_est(_B_vrad, _c_vrad, _cov_vrad, &yi, &yerr);
            TRparam[0].vrad_einastofit_bs[i] = yi*TRparam[0].vrad_max;

            // 1. VRAD_e : upper bound
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _c_vrad->data[vrad_nS] = TRparam[0]._vr_bs[vrad_nS] + TRparam[0]._vr_bs_e[vrad_nS];
            }
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);
            gsl_multifit_linear_est(_B_vrad, _c_vrad, _cov_vrad, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            vrad_e_upper_bound = yi*TRparam[0].vrad_max;

            // 2. VRAD_e : lower bound
            for(vrad_nS=0; vrad_nS<TRparam[0].n_coeffs_bspline_vrad; vrad_nS++)
            {
                _c_vrad->data[vrad_nS] = TRparam[0]._vr_bs[vrad_nS] - TRparam[0]._vr_bs_e[vrad_nS];
            }
            gsl_bspline_eval(xi, _B_vrad, _bw_vrad);
            gsl_multifit_linear_est(_B_vrad, _c_vrad, _cov_vrad, &yi, &yerr);
            // compute the upper bound of uncertainty from the fit
            vrad_e_lower_bound = yi*TRparam[0].vrad_max;
            TRparam[0].vrad_einastofit_bs_e[i] = (vrad_e_upper_bound - vrad_e_lower_bound)/2.0;
        }
        else if(strcmp(TRparam[0].vrad_function, "const") == 0)
        {
            TRparam[0].vrad[i] = TRparam[0].vradF;
            TRparam[0].vrad_e[i] = TRparam[0].vradF_e;
        }
    }
    free(x_dat);
    free(y_dat);
    free(e_dat);

    if(strcmp(TRparam[0].pa_function, "bspline") == 0)
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

    if(strcmp(TRparam[0].incl_function, "bspline") == 0)
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

    if(strcmp(TRparam[0].vrad_function, "bspline") == 0)
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
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* Define a tilted ring with the given ring parameters before TR fit*/
void define_area_tofit(double ring0, double ring1, int decimX, int decimY, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side)
{
    int i0=0, j0=0;
    int x0, y0, x1, y1;

    double i=0, j=0, id=0;
    double pa=0, incl=0, theta=0., r=0.;
    double Npoints_in_tiltedRing=0;
    double sine_free_angle, costh;
    double machep;
    Ellipse_Parameter ellipse;


    // 1. extract the region to fit with decimals starting from the centre position given
    Npoints_in_tiltedRing = 0;
    decimX = TRparam[0].decimX_einasto_halofit;
    decimY = TRparam[0].decimY_einasto_halofit;

    machep = r8_epsilon();
    pa = TRparam[0].paF_EinastoFit; // in degree
    incl = TRparam[0].inclF_EinastoFit; // in degree

    ellipse.xpos = TRparam[0].xposF_EinastoFit;
    ellipse.ypos = TRparam[0].yposF_EinastoFit;
    ellipse.a = ring0 + (ring1-ring0)/2.0;
    ellipse.b = ellipse.a * cos(incl*M_PI/180.);
    ellipse.e = sqrt(1-(ellipse.b*ellipse.b)/(ellipse.a*ellipse.a));

    if(ellipse.xpos == 0 && ellipse.ypos == 0 && pa == 0 && incl == 0)
    {
        ellipse.xpos = TRparam[0].nax1/2.0;
        ellipse.ypos = TRparam[0].nax2/2.0;
        pa = 45;
        incl = 45;
    }

    // set free angle
    sine_free_angle = fabs(sin(TRparam[0].free_angle*M_PI/180.));
    // between ring0 and ring1
    if(ring0 < 0) ring0 = 0;
    if(ring1 > TRparam[0].nax1) ring1 = TRparam[0].nax1;
    if(ring1 > TRparam[0].nax2) ring1 = TRparam[0].nax2;
    // 0 - nax1
    Npoints_in_tiltedRing = 0;

    x0 = 0;
    y0 = 0;
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;


    // Reset the connected area to fit
    for(id=0; id<4024*4024; id++)
    {
        TRparam[0].tilted_ring[(int)id][0] = 0;
        TRparam[0].tilted_ring[(int)id][1] = 0;
    }

    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
            j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
            r = sqrt(i*i+j*j); // distance from centre

            if(r < 1) // if r is smaller than one pixel
            {
                r = 1;
                theta = 0.0;
            }
            else
                theta = atan2((double)j, (double)i)*180./M_PI;

            costh = fabs(cos(theta*M_PI/180.));

            if(r >= ring0 && r < ring1 && (HI_VF_boxFiltered_decim0[0].data[j0][i0]+1) > HI_VF_boxFiltered_decim0[0].data[j0][i0] && costh > sine_free_angle)
            {
                if(side == -1 && fabs(theta) <= 90.0) // receding side
                {
                    TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                    TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                    Npoints_in_tiltedRing += 1;
                }
                else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                {
                    TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                    TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                    Npoints_in_tiltedRing += 1;
                }
                else if(side == 0) // both sides
                {
                    TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                    TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                    Npoints_in_tiltedRing += 1;
//printf("%d %d\n", i0, j0);
                }
            }
            j0 += decimY;
        }
        i0 += decimX;
    }
    TRparam[0].Npoints_in_tilted_ring = Npoints_in_tiltedRing;
    return;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void define_area_tofit_set_geo_weight(double ring0, double ring1, int decimX, int decimY, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int side, int histogram)
{
    double i=0, j=0;
    int i0=0, j0=0, k0=0, mn=0, mi=0, mj=0, k=0;
    int pa_n=0, incl_n=0;
    int Npoints_in_tiltedRing_without_decimals;
    int num_bins, hist_pts, Npoints_in_tilted_ring_filtered;

    float *LOS_velocities, interval, amp_bin_min, amp_bin_max, LOS_vel_hist_rbm, LOS_vel_hist_std, x_shift;
    float a0, v0, sig0, bg0, a_left, a_right, x_left, x_right;
    double n_freeParams_sGfit;
    double pa=0, incl=0, theta=0., r=0., r_inner=0., r_outter=0.;
    double Npoints_in_tiltedRing=0;
    double ring0_regrad, ring1_regrad;
    double sine_free_angle, costh;
    double rGalaxyPlane_pixel, machep, machine_epsilon, r_i0j0;
    Ellipse_Parameter ellipse;
    double a, b, perimeter, perimeter_max, perimeter_r1beam;
    int n_inner_avoid;

    int box_x, box_y, x0, y0, x1, y1;
    double FD_h, IQR;
    double ll_bin, ul_bin;
    int FD_nbins, h_loop;
    int bi=0, bn=0;

    double *_filterbox, *weight_temp, *ve_temp, mean_ve_temp, std_ve_temp;
    double mean_vf_e_temp_w, std_vf_e_temp_w;
    double mean_mom2_temp, std_mom2_temp;
    double mean_mom4_temp, std_mom4_temp;
    double HI_VF_geo_radial_angle_w_max;
    double lower_bound_LOS_velocities, upper_bound_LOS_velocities;
    double hist_mean_LOS, hist_std_LOS;

    double lower_bound_filterbox, upper_bound_filterbox;
    double hist_mean_filterbox, hist_std_filterbox;
    double mean_vf_e, std_vf_e;
    double v_sigma_weighted;

    sine_free_angle = fabs(sin(TRparam[0].free_angle*M_PI/180.));

    gsl_interp_accel *acc;
    gsl_spline *spline;

    // 1. extract the region to fit with decimals starting from the centre position given
    Npoints_in_tiltedRing = 0;
    //decimX = TRparam[0].decimX_einasto_halofit;
    //decimY = TRparam[0].decimY_einasto_halofit;

    machep = r8_epsilon();
    machine_epsilon = machep;

    for(i0=0; i0<TRparam[0].nax1; i0++)
    {
        for(j0=0; j0<TRparam[0].nax2; j0++)
        {
            // initialise
            HI_VF_geo_radial_angle_w[0].data[j0][i0] = 1; // default weighting
        }
    }

    //pa = TRparam[0].paF_EinastoFit * TRparam[0].PA_MAX_in_degree;
    pa = TRparam[0].pa0[TRparam[0].Nrings-1]; // in degree

    incl = TRparam[0].incl0[TRparam[0].Nrings-1]; // in degree

    ellipse.xpos = TRparam[0].xpos[0];
    ellipse.ypos = TRparam[0].ypos[0];

    // calculate rms of the residual between the input VF and the ISO model derived in the very previous step: for filtering out outliers
    for(i0=0; i0<TRparam[0].nax1; i0++)
    {
        for(j0=0; j0<TRparam[0].nax2; j0++)
        {
            if(!isinf(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && !isinf(HI_VF_einasto_halomodel[0].data[j0][i0]) && !isnan(HI_VF_einasto_halomodel[0].data[j0][i0]))
            {
                HI_VF_res[0].data[j0][i0] = HI_VF_boxFiltered_decim0[0].data[j0][i0]-HI_VF_einasto_halomodel[0].data[j0][i0];
                HI_VF_geo_radial_angle_w[0].data[j0][i0] = 1.0;
            }
            else
            {
                HI_VF_res[0].data[j0][i0] = 1E90; // blank
                HI_VF_geo_radial_angle_w[0].data[j0][i0] = 1.0;
            }
        }
    }

    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    x0 = 0;
    y0 = 0;
    x1 = TRparam[0].nax1;
    y1 = TRparam[0].nax2;
    _filterbox = malloc(sizeof(double) * box_x*box_y);
    weight_temp = malloc(sizeof(double) * (x1-x0)*(y1-y0));
    // decimX = 0
    // decimY = 0
    // put HI_VF_geo_radial_angle_w values to all the available pixels in the
    // velocity field (this is for calculating VROT_e in the last step where
    // vlos_e_w is needed.
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_res[0].data[j0][i0]) && !isnan(HI_VF_res[0].data[j0][i0])) // no blank
            {
                for (bi=0; bi<box_x*box_y; bi++)
                {
                    _filterbox[bi] = 1E9;
                }

                // extract _boxfilter
                bn = 0;
                for(mi=-(box_x-1)/2; mi<(box_x+1)/2; mi++)
                {
                    for(mj=-(box_y-1)/2; mj<(box_y+1)/2; mj++)
                    {
                        if(i0+mi < 0 || i0+mi >= TRparam[0].nax1 || j0+mj < 0 || j0+mj >= TRparam[0].nax2) continue;
                        if(!isinf(HI_VF_res[0].data[j0+mj][i0+mi]) && !isnan(HI_VF_res[0].data[j0+mj][i0+mi]) && HI_VF_res[0].data[j0+mj][i0+mi] != 0.0) // if no blank
                        {
                            _filterbox[bn] = HI_VF_res[0].data[j0+mj][i0+mi];
                            bn++;
                        }
                    }
                }
                // derive a robust mean of data in the _boxfilter based on their histogram
                robust_mean_std(_filterbox, bn, &hist_mean_filterbox, &hist_std_filterbox);

                // for each pixel
                // derive r(i0,j0) with the B-splines of PA+INCL
                r_i0j0 = r_ij_pa_incl_bspline_W(i0, j0, TRparam);

                acc = gsl_interp_accel_alloc();
                spline = gsl_spline_alloc (gsl_interp_cspline, TRparam[0].Nrings);
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl0, TRparam[0].Nrings);

                // interpolate incl at r_i0j0
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    incl = TRparam[0].incl0[TRparam[0].Nrings-1]; // outermost incl
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    incl = TRparam[0].incl0[0]; // outermost incl
                }
                else
                {
                    incl = gsl_spline_eval (spline, r_i0j0, acc);
                }

                // interpolate pa at r_i0j0
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].pa0, TRparam[0].Nrings);
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    pa = TRparam[0].pa0[TRparam[0].Nrings-1]; // outermost pa
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    pa = TRparam[0].pa0[0]; // outermost pa
                }
                else
                {
                    pa = gsl_spline_eval (spline, r_i0j0, acc);
                }
                gsl_spline_free(spline);
                gsl_interp_accel_free(acc);

                ellipse.a = ring0 + (ring1-ring0)/2.0;
                ellipse.b = ellipse.a * cos(incl*M_PI/180.);
                ellipse.e = sqrt(1-(ellipse.b*ellipse.b)/(ellipse.a*ellipse.a));

                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                //if(r < TRparam[0].ring_radius[0]) // if r is smaller than one pixel
                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));


                // perimeter
                a = 1*TRparam[0].rGalaxyPlane_pixel_max;
                b = a * cos(incl*M_PI/180.);
                perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
       
                n_inner_avoid = (int)0.0*TRparam[0].Nrings; // inner 0%
                a = TRparam[0].ring_radius[n_inner_avoid];
                b = a * cos(incl*M_PI/180.);
                perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
                perimeter_r1beam = perimeter_max/perimeter;
                //perimeter_r1beam = -1; // negative : outlied in the Einasto fitting

                a = r; 
                b = a * cos(incl*M_PI/180.);
                perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
                //TRparam[0].perimeter = perimeter/perimeter_max;
                TRparam[0].perimeter = perimeter_max/perimeter;

                // filter out any spike-like outliers
                //if(!isinf(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && costh > sine_free_angle && (HI_VF_res[0].data[j0][i0] >= hist_mean_filterbox - 3*hist_std_filterbox) && (HI_VF_res[0].data[j0][i0] < hist_mean_filterbox + 3*hist_std_filterbox))
                if(!isinf(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && costh > sine_free_angle)
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        // sqrt(pow(costh, TRparam[0].wpow)) is going to be
                        // squared in the MLE by sigma*sigma..
                        //HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) / (TRparam[0].perimeter*HI_VF_fract_navail_nall[0].data[j0][i0]); // weight

                        //if(r > TRparam[0].ring_radius[1])
                        if(r > 0.1*TRparam[0].ring_radius[n_inner_avoid])
                        {
                            HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) * (TRparam[0].perimeter); // weight
                        }
                        else
                        {
                            //HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) * (perimeter_r1beam); // weight
                            HI_VF_geo_radial_angle_w[0].data[j0][i0] = -1; // pixels in the region within one ring are not used...
                        }


                        //TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        //TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        //Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        //HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) / (TRparam[0].perimeter*HI_VF_fract_navail_nall[0].data[j0][i0]); // weight
                        //if(r > TRparam[0].ring_radius[1])
                        if(r > TRparam[0].ring_radius[n_inner_avoid])
                        {
                            HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) * (TRparam[0].perimeter); // weight
                        }
                        else
                        {
                            HI_VF_geo_radial_angle_w[0].data[j0][i0] = -1; // pixels in the region within one ring are not used...
                        }

                        //TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        //TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        //Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0) // both sides
                    {
                        //HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) / (TRparam[0].perimeter*HI_VF_fract_navail_nall[0].data[j0][i0]); // weight
                        if(r > TRparam[0].ring_radius[n_inner_avoid])
                        {
                            HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) * (TRparam[0].perimeter); // weight
                        }
                        else
                        {
                            //HI_VF_geo_radial_angle_w[0].data[j0][i0] = sqrt(pow(costh, TRparam[0].wpow)) * pow(TRparam[0].rGalaxyPlane_pixel_max/r, TRparam[0].rwpow) * (perimeter_r1beam); // weight
                            HI_VF_geo_radial_angle_w[0].data[j0][i0] = -1; // pixels in the region within one ring are not used...
                        }
                        //TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        //TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        //Npoints_in_tiltedRing += 1;
                    }

                    if(HI_VF_geo_radial_angle_w[0].data[j0][i0] == 0 || isnan(HI_VF_geo_radial_angle_w[0].data[j0][i0]) || isinf(HI_VF_geo_radial_angle_w[0].data[j0][i0]))
                        HI_VF_geo_radial_angle_w[0].data[j0][i0] = 1;
                }
                else
                    HI_VF_geo_radial_angle_w[0].data[j0][i0] = 1;
                //else // for connected area issue
                //    HI_VF_boxFiltered[0].data[j0][i0] = 1E90;
            }
        }
    }

    // given Einasto decim x, y
    Npoints_in_tiltedRing = 0;
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            if(!isinf(HI_VF_res[0].data[j0][i0]) && !isnan(HI_VF_res[0].data[j0][i0])) // no blank
            {

                for (bi=0; bi<box_x*box_y; bi++)
                {
                    _filterbox[bi] = 1E9;
                }

                // extract _boxfilter
                bn = 0;
                for(mi=-(box_x-1)/2; mi<(box_x+1)/2; mi++)
                {
                    for(mj=-(box_y-1)/2; mj<(box_y+1)/2; mj++)
                    {
                        if(i0+mi < 0 || i0+mi >= TRparam[0].nax1 || j0+mj < 0 || j0+mj >= TRparam[0].nax2) continue;
                        if(!isinf(HI_VF_res[0].data[j0+mj][i0+mi]) && !isnan(HI_VF_res[0].data[j0+mj][i0+mi]) && HI_VF_res[0].data[j0+mj][i0+mi] != 0.0) // if no blank
                        {
                            _filterbox[bn] = HI_VF_res[0].data[j0+mj][i0+mi];
                            bn++;
                        }
                    }
                }
                // derive a robust mean of data in the _boxfilter based on their histogram
                robust_mean_std(_filterbox, bn, &hist_mean_filterbox, &hist_std_filterbox);

                // for each pixel
                // derive r(i0,j0) with the B-splines of PA+INCL
                r_i0j0 = r_ij_pa_incl_bspline_W(i0, j0, TRparam);

                acc = gsl_interp_accel_alloc();
                spline = gsl_spline_alloc (gsl_interp_cspline, TRparam[0].Nrings);
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].incl0, TRparam[0].Nrings);

                // interpolate incl at r_i0j0
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    incl = TRparam[0].incl0[TRparam[0].Nrings-1]; // outermost incl
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    incl = TRparam[0].incl0[0]; // outermost incl
                }
                else
                {
                    incl = gsl_spline_eval (spline, r_i0j0, acc);
                }

                // interpolate pa at r_i0j0
                gsl_spline_init(spline, TRparam[0].ring_radius, TRparam[0].pa0, TRparam[0].Nrings);
                if(r_i0j0 >= TRparam[0].ring_radius[TRparam[0].Nrings-1]) // if a is outside the ring range where interpolation can be done
                {
                    pa = TRparam[0].pa0[TRparam[0].Nrings-1]; // outermost pa
                }
                else if(r_i0j0 <= TRparam[0].ring_radius[0]) // if a is outside the ring range where interpolation can be done
                {
                    pa = TRparam[0].pa0[0]; // outermost pa
                }
                else
                {
                    pa = gsl_spline_eval (spline, r_i0j0, acc);
                }
                gsl_spline_free(spline);
                gsl_interp_accel_free(acc);

                ellipse.a = ring0 + (ring1-ring0)/2.0;
                ellipse.b = ellipse.a * cos(incl*M_PI/180.);
                ellipse.e = sqrt(1-(ellipse.b*ellipse.b)/(ellipse.a*ellipse.a));


                i = (double)( - ( (double)i0 - ellipse.xpos ) * sin(pa*M_PI/180.) + ( (double)j0 - ellipse.ypos ) * cos(pa*M_PI/180.));
                j = (double)( - ( (double)i0 - ellipse.xpos ) * cos(pa*M_PI/180.) - ( (double)j0 - ellipse.ypos ) * sin(pa*M_PI/180.))/cos(incl*M_PI/180.);
                r = sqrt(i*i+j*j); // distance from centre

                if(r < 1) // if r is smaller than one pixel
                {
                    r = 1;
                    theta = 0.0;
                }
                else
                    theta = atan2((double)j, (double)i)*180./M_PI;

                costh = fabs(cos(theta*M_PI/180.));


                // perimeter
                a = TRparam[0].rGalaxyPlane_pixel_max;
                b = a * cos(incl*M_PI/180.);
                perimeter_max = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));

                a = r; 
                b = a * cos(incl*M_PI/180.);
                perimeter = M_PI*(3.0*(a+b) - sqrt((3.0*a+b)*(a+3.0*b)));
                TRparam[0].perimeter = perimeter_max/perimeter;

                // filter out any spike-like outliers
                //if(!isinf(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && costh > sine_free_angle && (HI_VF_res[0].data[j0][i0] >= hist_mean_filterbox - 3*hist_std_filterbox) && (HI_VF_res[0].data[j0][i0] < hist_mean_filterbox + 3*hist_std_filterbox))
                if(!isinf(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_decim0[0].data[j0][i0]) && costh > sine_free_angle && r > TRparam[0].ring_radius[n_inner_avoid])
                {
                    if(side == -1 && fabs(theta) <= 90.0) // receding side
                    {
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_geo_radial_angle_w[0].data[j0][i0];
                        HI_VF_boxFiltered[0].data[j0][i0] = HI_VF_boxFiltered_decim0[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 1 && fabs(theta) >= 90.0) // approaching side
                    {
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_geo_radial_angle_w[0].data[j0][i0];
                        HI_VF_boxFiltered[0].data[j0][i0] = HI_VF_boxFiltered_decim0[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }
                    else if(side == 0) // both sides
                    {
                        weight_temp[(int)Npoints_in_tiltedRing] = HI_VF_geo_radial_angle_w[0].data[j0][i0];
                        HI_VF_boxFiltered[0].data[j0][i0] = HI_VF_boxFiltered_decim0[0].data[j0][i0];
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][0] = i0;
                        TRparam[0].tilted_ring[(int)Npoints_in_tiltedRing][1] = j0;
                        Npoints_in_tiltedRing += 1;
                    }

                    if(HI_VF_geo_radial_angle_w[0].data[j0][i0] == 0 || isnan(HI_VF_geo_radial_angle_w[0].data[j0][i0]) || isinf(HI_VF_geo_radial_angle_w[0].data[j0][i0]))
                        HI_VF_geo_radial_angle_w[0].data[j0][i0] = 1;
                }
            }
            j0 += decimY;
        }
        i0 += decimX;
    }
    TRparam[0].Npoints_in_tilted_ring = Npoints_in_tiltedRing;

    HI_VF_geo_radial_angle_w_max = gsl_stats_max(weight_temp, 1, Npoints_in_tiltedRing);
    for(i0=x0; i0<x1; i0++)
    {
        for(j0=y0; j0<y1; j0++)
        {
            // normalise (0 ~ 1)
            //HI_VF_geo_radial_angle_w[0].data[j0][i0] = HI_VF_geo_radial_angle_w[0].data[j0][i0] / HI_VF_geo_radial_angle_w_max;
            HI_VF_geo_radial_angle_w[0].data[j0][i0] = HI_VF_geo_radial_angle_w[0].data[j0][i0];
        }
    }

    /*
    double _ring_w, _xpos, _ypos, _pa, _incl, ri, ro, ring;
    for(i=0; i<TRparam[0].nax1; i++)
    {
        _ring_w = 1;
        ri = 0 + i*_ring_w - 0.5*_ring_w;
        ro = 0 + i*_ring_w + 0.5*_ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        if(ring >= TRparam[0].ring_radius[TRparam[0].Nrings-1])
        {
            _xpos = TRparam[0].xpos0[TRparam[0].Nrings-1];
            _ypos = TRparam[0].ypos0[TRparam[0].Nrings-1];
            _pa = TRparam[0].pa0[TRparam[0].Nrings-1];
            _incl = TRparam[0].incl0[TRparam[0].Nrings-1];
        }

        find_Navail_Nall_pixels(_xpos, _ypos, _pa, _incl, ri, ro, TRparam, side);
    }
    */


    // Derive the geometry (+flux) weighted scale factor for vlos_e
    if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma > 0) // constant vlos_e mode
    {
        ve_temp = malloc(sizeof(double) * Npoints_in_tiltedRing);
        // compute scale_factor_const_vlose_w
        for(k0=0; k0<(int)Npoints_in_tiltedRing; k0++)
        {
        
            i0 = TRparam[0].tilted_ring[k0][0];
            j0 = TRparam[0].tilted_ring[k0][1];
            // normalised constant one channel vlos_error with the geometry weight
            ve_temp[k0] = TRparam[0].e_sigma / (weight_temp[k0]/HI_VF_geo_radial_angle_w_max);
        }

        // compute the mean of geometry-weighted constant vlos_e
        robust_mean_std(ve_temp, Npoints_in_tiltedRing, &mean_ve_temp, &std_ve_temp);
        // compute a scale factor for matching the geometry-weighted constant
        // vlos_e (which includes large values) with the input constant sigma_v
        // (this is to avoid large vlos_e in the bayesian fit
        TRparam[0].scale_factor_const_vlose_w = TRparam[0].e_sigma / mean_ve_temp;
        free(ve_temp);
    }
    else if(TRparam[0].sigma_factor_fix == 'T' && TRparam[0].e_sigma == 0) // variable vlos_e mode : full fit
    {
        double *mom4_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *mom2_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *vf_e_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *vf_e_temp_w = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);

        double mean_vf_e, std_vf_e;
        for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++) 
        {
            i0 = TRparam[0].tilted_ring[k0][0];
            j0 = TRparam[0].tilted_ring[k0][1];

            // 1. read vf_e
            if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0)
                vf_e_temp[k0] = HI_VF_boxFiltered_sigma[0].data[j0][i0];
            else
                vf_e_temp[k0] = 1; // a large value in km/s

            // 2. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
            if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]) && HI_VF_mom0[0].data[j0][i0] > 0)
                mom4_temp[k0] = HI_VF_mom0[0].data[j0][i0];
            else
                mom4_temp[k0] = 1E-5; // a small value mjy/beam

            // 3. read mom2 : this is used for matching the vlos_e_w with mom2 distribution (see below)
            if(!isinf(HI_VF_mom2[0].data[j0][i0]) && !isnan(HI_VF_mom2[0].data[j0][i0]) && HI_VF_mom2[0].data[j0][i0] > 0)
                mom2_temp[k0] = HI_VF_mom2[0].data[j0][i0];
            else
                mom2_temp[k0] = 1E3; // a large value in km/s

            // geometry & flux-weighted vlos_e_w
            vf_e_temp_w[k0] = vf_e_temp[k0] / (weight_temp[k0]/HI_VF_geo_radial_angle_w_max) / sqrt(mom4_temp[k0]);
        }

        // 1. compute the mean of geometry & flux-weighted variable vlos_e_w
        robust_mean_std(vf_e_temp_w, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e_temp_w, &std_vf_e_temp_w);

        // 2. compute the mean of mom2 
        robust_mean_std(mom2_temp, TRparam[0].Npoints_in_tilted_ring, &mean_mom2_temp, &std_mom2_temp);
        robust_mean_std(vf_e_temp, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e, &std_vf_e);

        // scale factor
        TRparam[0].scale_factor_var_vlose_w = mean_vf_e / mean_vf_e_temp_w;


        free(mom4_temp);
        free(mom2_temp);
        free(vf_e_temp);
        free(vf_e_temp_w);
    }
    else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0) // variable vlos_e mode : partial fit
    {
        double *mom4_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *vf_e_temp = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);
        double *vf_e_temp_w = malloc(sizeof(double) * TRparam[0].Npoints_in_tilted_ring);

        // compute mean & std of mom4
        k=0;
        for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++)
        {
            i0 = TRparam[0].tilted_ring[k0][0];
            j0 = TRparam[0].tilted_ring[k0][1];

            // 0. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
            if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]) && HI_VF_mom0[0].data[j0][i0] > 0)
            {
                mom4_temp[k] = HI_VF_mom0[0].data[j0][i0];
                k++;
            }
        }
        robust_mean_std(mom4_temp, k, &mean_mom4_temp, &std_mom4_temp);
        TRparam[0].mean_mom4 = mean_mom4_temp;
        TRparam[0].std_mom4 = std_mom4_temp;

        for(k0=0; k0<TRparam[0].Npoints_in_tilted_ring; k0++) 
        {
            i0 = TRparam[0].tilted_ring[k0][0];
            j0 = TRparam[0].tilted_ring[k0][1];

            // 1. read vf_e
            if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0)
                vf_e_temp[k0] = HI_VF_boxFiltered_sigma[0].data[j0][i0];
            else
                vf_e_temp[k0] = 1; // a large value in km/s

            // 2. read mom4 (peak) : HI_VF_mom0[0].data is actually mom4 : the parameter name will be changed later (tbd)
            if(!isinf(HI_VF_mom0[0].data[j0][i0]) && !isnan(HI_VF_mom0[0].data[j0][i0]))
            {
                if(HI_VF_mom0[0].data[j0][i0] > 1E5*mean_mom4_temp)
                {
                    mom4_temp[k0] = mean_mom4_temp;
                }
                else
                {
                    mom4_temp[k0] = HI_VF_mom0[0].data[j0][i0];
                }
            }


            // geometry & flux-weighted vlos_e_w
            //vf_e_temp_w[k0] = vf_e_temp[k0] / (weight_temp[k0]/HI_VF_geo_radial_angle_w_max) / sqrt(mom4_temp[k0]);
            //vf_e_temp_w[k0] = vf_e_temp[k0] / (mom4_temp[k0]);
            vf_e_temp_w[k0] = vf_e_temp[k0];
        }

        // 1. compute the mean of geometry & flux-weighted variable vlos_e_w
        robust_mean_std(vf_e_temp_w, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e_temp_w, &std_vf_e_temp_w);
        // compute mean and std of vf_e
        robust_mean_std(vf_e_temp, TRparam[0].Npoints_in_tilted_ring, &mean_vf_e, &std_vf_e);

        // scale factor
        TRparam[0].scale_factor_var_vlose_w = mean_vf_e / mean_vf_e_temp_w;


        free(mom4_temp);
        free(vf_e_temp);
        free(vf_e_temp_w);
    }

    // save the geometrical + flux weighted error map
    for(i0=0; i0<TRparam[0].nax1; i0++)
    {
        for(j0=0; j0<TRparam[0].nax2; j0++)
        {
            HI_VF_sigma_geo_flux_weighted[0].data[j0][i0] = 1E90;
        }
    }

    for(i0=0; i0<TRparam[0].nax1; i0++)
    {
        for(j0=0; j0<TRparam[0].nax2; j0++)
        {

            if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma > 0 && HI_VF_geo_radial_angle_w[0].data[j0][i0] > 0) // constant vlos_e mode
            { 
                // same vlos_e weight 
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0) 
                {
                    HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] = 1;
                }
            }
            else if(TRparam[0].sigma_factor_fix == 'F' && TRparam[0].e_sigma == 0 && HI_VF_geo_radial_angle_w[0].data[j0][i0] > 0) // constant vlos_e mode : partial fit
            { 
                // geometry (perimeter + cos(theta)) weighted constant vlos_e
                // HI_VF_geo_radial_angle_w[0].data[j0][i0] is weighted by perimeter + cos(theta) and the normalised (0 ~ 1) 
                if(!isinf(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && !isnan(HI_VF_boxFiltered_sigma[0].data[j0][i0]) \
                && HI_VF_boxFiltered_sigma[0].data[j0][i0] > 0) 
                {
                    // normalised & inverse weighted
                    HI_VF_boxFiltered_vlos_ew[0].data[j0][i0] = 1.0/(HI_VF_boxFiltered_sigma[0].data[j0][i0]/mean_vf_e);
                }
            }
        }
    }

    free(_filterbox);
    free(weight_temp);
    return;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
                            char *fullfit, char fullfit_fix)
{
    int i=0;
    int decimX, decimY;
    int x0, y0, x1, y1, box_x, box_y;
    double ring, _ring_w, _xpos, _ypos, _pa, _incl, ri, ro;

    TRparam[0].xpos_fix = xposfix;
    TRparam[0].ypos_fix = yposfix;
    TRparam[0].vsys_fix = vsysfix;
    TRparam[0].pa_fix = pafix;
    TRparam[0].incl_fix = inclfix;
    TRparam[0]._n_fix = _nfix;
    TRparam[0].r_2_fix = _r_2fix;
    TRparam[0].rho_2_fix = _rho_2fix;
    TRparam[0].fullFit = fullfit_fix;

    TRparam[0].sigma_factor1 = 0;
    TRparam[0].sigma_factor2 = 10;

    reset_HI_VF_fract_navail_nall(TRparam);
    for(i=0; i<2*TRparam[0].Nrings; i++)
    {
        ri = 0 + i*_ring_w - _ring_w;
        ro = 0 + i*_ring_w + 0.5*_ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        if(ring > TRparam[0].ring_radius[TRparam[0].Nrings-1])
        {
            _ring_w = TRparam[0].ring_w;
        }
        else
        {
            _ring_w = TRparam[0].ring_w;
        }

        ri = 0 + i*_ring_w - _ring_w;
        ro = 0 + i*_ring_w + 0.5*_ring_w;
        if ( ri < 0.0) ri = 0.0;
        ring = (ri+ro)/2.0;

        if(ring <= TRparam[0].ring_radius[0])
        {
            _xpos = TRparam[0].xpos0[0];
            _ypos = TRparam[0].ypos0[0];
            _pa = TRparam[0].pa0[0];
            _incl = TRparam[0].incl0[0];
        }
        if(ring > TRparam[0].ring_radius[0] && ring < TRparam[0].ring_radius[TRparam[0].Nrings-1])
        {
            _xpos = spline_intp_ringparam(TRparam, "XPOS", ring);
            _ypos = spline_intp_ringparam(TRparam, "YPOS", ring);
            _pa = spline_intp_ringparam(TRparam, "PA", ring);
            _incl = spline_intp_ringparam(TRparam, "INCL", ring);
        }
        else
        {
            _xpos = TRparam[0].xpos0[TRparam[0].Nrings-1];
            _ypos = TRparam[0].ypos0[TRparam[0].Nrings-1];
            _pa = TRparam[0].pa0[TRparam[0].Nrings-1];
            _incl = TRparam[0].incl0[TRparam[0].Nrings-1];
        }
        find_Navail_Nall_pixels(_xpos, _ypos, _pa, _incl, ri, ro, TRparam, 0, Nall_ring, Navail_ring);
    }

    return;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// --- End of line
