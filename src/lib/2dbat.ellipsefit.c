#include "2dbat.ellipsefit.h"

// 2DBAT user defined functions
// ellipse fit related


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* Ellipse equation in rectangular coordinate */
double ellipse_function_rect(double a, double b, double x)
{
    double FF;
    FF = b*sqrt(1. - (x*x)/(a*a));
    return FF;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* Ellipse equation in polar coordinate */
double ellipse_function_polar(double a, double e, double theta)
{
    double FF;
    FF = a*sqrt((1.-e*e)/(1.-e*e*cos(theta*M_PI/180.)*cos(theta*M_PI/180.)));
    return FF;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void do_ellipseFit_update_uniformPriors_ringParam(TR_ringParameters *TRparam)
{
    // Define an area to fit from ellipse fit
    // Tilted-ring fit using the rings defined with the ellipse fit : this is
    // for deriving a first-approximate rotation velocity which will be used for
    // fitting ISO rotation velocity. From the ISO fit, initial estimates for rc
    // and rho are derived, which are used then for guessing their uniform
    // priors: Note the later Bayesian ISO model fitting will be sensitive to
    // these uniform priors.

    int i, j, k, i0, j0;
    // declare variables
    int     m[5];
    float   p[5];
    float   p_pre[5];
    float   e[5];
    float delt_X;
    float delt_Y;
    float GG = 4.302;
    double peri_x=0, peri_y=0, xpos_temp, ypos_temp;
    double vel_geo_rec_side, vel_geo_app_side, R_semi_mx, init_Vrot_max, init_ISO_rF;
    float *x_temp = malloc(sizeof(float)*TRparam[0].Npoints_in_tilted_ring);
    float *y_temp = malloc(sizeof(float)*TRparam[0].Npoints_in_tilted_ring);

    // 2-0. Initialise parametres for ellipse fit
    for (k=0; k <5; k++)
    {
        p[k] = 0.0;
        p_pre[k] = 0.0;
        e[k] = 0.0;
        m[k] = 1;
    }

    // 2-1. Set the defined area for ellipse fit
    for(i=0; i<TRparam[0].Npoints_in_tilted_ring; i++)
    {
        x_temp[i] = TRparam[0].tilted_ring[i][0];
        y_temp[i] = TRparam[0].tilted_ring[i][1];
    }

    // 2-2. Perform ellipse fit
    printf("Perform ellipse fit to the defined area: ");
    k = ellipse1_c(&TRparam[0].Npoints_in_tilted_ring, x_temp, y_temp, p);
    k = ellipse2_c(&TRparam[0].Npoints_in_tilted_ring, x_temp, y_temp, p, e, m);
    printf("[Done] \n\n");

    // 2-3. Save the ellipse fit results for the optimal blank pixel filtering : note pa and incl in degree
    TRparam[0].ellipse_xpos_boxfiltered = p[2];
    TRparam[0].ellipse_ypos_boxfiltered = p[3];
    TRparam[0].ellipse_pa_boxfiltered = p[4];
    TRparam[0].pa_EllipseFit_e = e[4];

    TRparam[0].ellipse_incl_boxfiltered = p[1];
    TRparam[0].incl_EllipseFit_e = e[1];
    TRparam[0].ellipse_semi_mx_boxfiltered = 1.0*p[0]; // p[0] semi-major axis

//printf("incl:%f xpos:%f ypos:%f mx:%f nax1:%d nax2:%d, pa:%f\n", TRparam[0].ellipse_incl_boxfiltered, TRparam[0].ellipse_xpos_boxfiltered, TRparam[0].ellipse_ypos_boxfiltered, TRparam[0].ellipse_semi_mx_boxfiltered, TRparam[0].nax1, TRparam[0].nax2, TRparam[0].ellipse_pa_boxfiltered);

    if(TRparam[0].ellipse_xpos_boxfiltered < 0 || TRparam[0].ellipse_xpos_boxfiltered > TRparam[0].nax1 || \
       TRparam[0].ellipse_ypos_boxfiltered < 0 || TRparam[0].ellipse_ypos_boxfiltered > TRparam[0].nax2 || \
       TRparam[0].ellipse_semi_mx_boxfiltered < 0 || TRparam[0].ellipse_semi_mx_boxfiltered > TRparam[0].nax2 || \
       TRparam[0].ellipse_pa_boxfiltered < 0 || TRparam[0].ellipse_pa_boxfiltered > 360)
    {
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("Ellipse fitting seems weired!!!\n");
        printf("Check the binning used\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        //exit(1);
    }

    //TRparam[0].ellipse_xpos_boxfiltered = 46;
    //TRparam[0].ellipse_ypos_boxfiltered = 30;
    //TRparam[0].ellipse_incl_boxfiltered = 50;
    //TRparam[0].ellipse_pa_boxfiltered = 150;
    //TRparam[0].ellipse_semi_mx_boxfiltered = 20;

    // 2-4. Check velocity of the geometric receding side
    printf("!+++ DETERMINE THE RECEDING AND APPROACHING SIDES: ");
    if((TRparam[0].ellipse_pa_boxfiltered >= 0 && TRparam[0].ellipse_pa_boxfiltered < 90) || (TRparam[0].ellipse_pa_boxfiltered >= 180 && TRparam[0].ellipse_pa_boxfiltered < 270) )
    {
        i0 = (int)(TRparam[0].ellipse_xpos_boxfiltered - 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*sin(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));
        j0 = (int)(TRparam[0].ellipse_ypos_boxfiltered + 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*cos(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));

        while(1)
        {
            if(isinf(HI_VF[0].data[j0][i0]) != 1 && isnan(HI_VF[0].data[j0][i0]) != 1)
            {
                vel_geo_rec_side = HI_VF[0].data[j0][i0];
                break;
            }
            else
            {
                i0 += 1;
                j0 -= 1;
            }
        }
        // Check velocity of the geometric approaching side
        i0 = (int)(TRparam[0].ellipse_xpos_boxfiltered + 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*sin(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));
        j0 = (int)(TRparam[0].ellipse_ypos_boxfiltered - 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*cos(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));
        while(1)
        {
            if(isinf(HI_VF[0].data[j0][i0]) != 1 && isnan(HI_VF[0].data[j0][i0]) != 1)
            {
                vel_geo_app_side = HI_VF[0].data[j0][i0];
                break;
            }
            else
            {
                i0 -= 1;
                j0 += 1;
            }
        }
    }

    if((TRparam[0].ellipse_pa_boxfiltered >= 90 && TRparam[0].ellipse_pa_boxfiltered < 180) || (TRparam[0].ellipse_pa_boxfiltered >= 270 && TRparam[0].ellipse_pa_boxfiltered < 360))
    {
        i0 = (int)(TRparam[0].ellipse_xpos_boxfiltered - 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*sin(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));
        j0 = (int)(TRparam[0].ellipse_ypos_boxfiltered + 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*cos(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));

        while(1)
        {
            if(isinf(HI_VF[0].data[j0][i0]) != 1 && isnan(HI_VF[0].data[j0][i0]) != 1)
            {
                vel_geo_rec_side = HI_VF[0].data[j0][i0];
                break;
            }
            else
            {
                i0 += 1;
                j0 += 1;
            }
        }
        // Check velocity of the geometric approaching side
        i0 = (int)(TRparam[0].ellipse_xpos_boxfiltered + 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*sin(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));
        j0 = (int)(TRparam[0].ellipse_ypos_boxfiltered - 0.5*TRparam[0].ellipse_semi_mx_boxfiltered*cos(TRparam[0].ellipse_pa_boxfiltered*M_PI/180.));
        while(1)
        {
            if(isinf(HI_VF[0].data[j0][i0]) != 1 && isnan(HI_VF[0].data[j0][i0]) != 1)
            {
                vel_geo_app_side = HI_VF[0].data[j0][i0];
                break;
            }
            else
            {
                i0 -= 1;
                j0 -= 1;
            }
        }
    }

    // 2-5. Update the velocities of receding and approaching sides
    TRparam[0].vel_geo_app_side = vel_geo_app_side;
    TRparam[0].vel_geo_rec_side = vel_geo_rec_side;

    // 2-6. Initialise TR ring parametres with the ellipse fit results 
    // save the final ellipse fit results: note pa and incl in units
    // initialise ring parametres + their uniform priors based on the first ellipse fit
    printf("!+++ SET GAUSSIAN (XPOS, YPOS, VSYS, PA, INCL) + UNIFORM (VROT) PRIORS BASED ON THE ELLIPSE FIT: ");
    printf("[Done] \n\n");
    // A. xpos
    TRparam[0].xposF_EinastoFit = p[2];
    TRparam[0].xposF_EinastoFit_t = p[2];
    TRparam[0].xposF = p[2];
    TRparam[0].xpos1 = TRparam[0].xposF_EinastoFit - (TRparam[0].ellipse_semi_mx_boxfiltered/2.0); // use a wider range as based on the initial ellipse fit
    TRparam[0].xpos2 = TRparam[0].xposF_EinastoFit + (TRparam[0].ellipse_semi_mx_boxfiltered/2.0);
    TRparam[0].xpos1_from_ellipsefit = TRparam[0].xposF_EinastoFit - (TRparam[0].ellipse_semi_mx_boxfiltered/2.0); // use as a hard lower limit 
    TRparam[0].xpos2_from_ellipsefit = TRparam[0].xposF_EinastoFit + (TRparam[0].ellipse_semi_mx_boxfiltered/2.0); // use as a hard upper limit 


    // B. ypos
    TRparam[0].yposF_EinastoFit = p[3];
    TRparam[0].yposF_EinastoFit_t = p[3];
    TRparam[0].yposF = p[3];
    TRparam[0].ypos1 = TRparam[0].yposF_EinastoFit - (TRparam[0].ellipse_semi_mx_boxfiltered/2.0);
    TRparam[0].ypos2 = TRparam[0].yposF_EinastoFit + (TRparam[0].ellipse_semi_mx_boxfiltered/2.0);
    TRparam[0].ypos1_from_ellipsefit = TRparam[0].yposF_EinastoFit - (TRparam[0].ellipse_semi_mx_boxfiltered/2.0); // use as a hard lower limit
    TRparam[0].ypos2_from_ellipsefit = TRparam[0].yposF_EinastoFit + (TRparam[0].ellipse_semi_mx_boxfiltered/2.0); // use as a hard upper limit

    // C. pa
    if(vel_geo_rec_side >= vel_geo_app_side)
    {
        //TRparam[0].paF_EinastoFit = p[4]/TRparam[0].PA_MAX_in_degree;
        //TRparam[0].paF = p[4]/TRparam[0].PA_MAX_in_degree;
        TRparam[0].paF_EinastoFit = p[4];
        TRparam[0].paF = p[4];
        TRparam[0].ellipse_pa_boxfiltered = TRparam[0].paF;
    }
    else
    {
        //TRparam[0].paF_EinastoFit = (p[4]+180)/TRparam[0].PA_MAX_in_degree;
        //TRparam[0].paF = (p[4]+180)/TRparam[0].PA_MAX_in_degree;
        TRparam[0].paF_EinastoFit = p[4]+180;
        TRparam[0].paF = p[4]+180;
        TRparam[0].ellipse_pa_boxfiltered = TRparam[0].paF;
    }

    TRparam[0].pa1 = (TRparam[0].paF_EinastoFit - 150.)/TRparam[0].PA_MAX_in_degree; // pa_ISO - 150 degree
    TRparam[0].pa1_for_TRfit = (TRparam[0].paF_EinastoFit - 150.)/TRparam[0].PA_MAX_in_degree; // pa_ISO - 150 degree
    TRparam[0].pa2 = (TRparam[0].paF_EinastoFit + 150.)/TRparam[0].PA_MAX_in_degree; // pa_ISO - 150 degree
    TRparam[0].pa2_for_TRfit = (TRparam[0].paF_EinastoFit + 150.)/TRparam[0].PA_MAX_in_degree; // pa_ISO - 150 degree

/*
    if(TRparam[0].pa1 < 0) // in this case, we need to cover all pa range
    {
        TRparam[0].pa1 = 0.0;
        TRparam[0].pa1_for_TRfit = 0.0;
    }
    if(TRparam[0].pa2 > 1)
    {
        TRparam[0].pa2 = 1.0;
        TRparam[0].pa2_for_TRfit = 1.0;
    }
    TRparam[0].pa1 = 0.0;
    TRparam[0].pa1_for_TRfit = 0.0;
    TRparam[0].pa2 = 1.0;
    TRparam[0].pa2_for_TRfit = 1.0;
*/

    // Find the coordinates of the edges of the major axis to derive more correct receding and approaching side velocities at the edges of the major axis
    delt_X = -1.0*TRparam[0].ellipse_semi_mx_boxfiltered * sin(TRparam[0].paF_EinastoFit*M_PI/180.);
    delt_Y = TRparam[0].ellipse_semi_mx_boxfiltered * cos(TRparam[0].paF_EinastoFit*M_PI/180.);

    // D. VSYS: combine the information from Gfit of LOS histogram and RBM estimates 
    TRparam[0].vsysF_EinastoFit = TRparam[0].LOS_hist_Gfit_V0;
    TRparam[0].vsysF_EinastoFit_t = TRparam[0].LOS_hist_Gfit_V0;
    TRparam[0].vsysF = TRparam[0].LOS_hist_Gfit_V0;

    TRparam[0].vsys1 = TRparam[0].LOS_vel_hist_rbm - 3*TRparam[0].LOS_vel_hist_std;
    TRparam[0].vsys2 = TRparam[0].LOS_vel_hist_rbm + 3*TRparam[0].LOS_vel_hist_std;
    TRparam[0].vsys1_from_vlosfit = TRparam[0].LOS_vel_hist_rbm - 3*TRparam[0].LOS_vel_hist_std; // use as a hard lower limit
    TRparam[0].vsys2_from_vlosfit = TRparam[0].LOS_vel_hist_rbm + 3*TRparam[0].LOS_vel_hist_std; // use as a hard upper limit

    // E. incl
    TRparam[0].inclF_EinastoFit = p[1];
    TRparam[0].inclF = p[1];
    TRparam[0].ellipse_incl_boxfiltered = TRparam[0].inclF;

    TRparam[0].incl1 = (TRparam[0].inclF_EinastoFit - 40)/TRparam[0].INCL_MAX_in_degree; // incl_ISO - 50 degree: use a bit wider uniform as it is based on the initial ellipse fit
    TRparam[0].incl1_for_TRfit = (TRparam[0].inclF_EinastoFit - 40)/TRparam[0].INCL_MAX_in_degree; // incl_ISO - 50 degree: use a bit wider uniform as it is based on the initial ellipse fit

/*
    if(TRparam[0].incl1_for_TRfit < 0)
        TRparam[0].incl1_for_TRfit = 0;

    if(TRparam[0].incl1 < 0)
    {
        TRparam[0].incl1 = 0.0;
        TRparam[0].incl1_for_TRfit = 0.0;
    }

    TRparam[0].incl2 = (TRparam[0].inclF_EinastoFit + 40)/TRparam[0].INCL_MAX_in_degree; // incl_ISO + 50 degree
    TRparam[0].incl2_for_TRfit = (TRparam[0].inclF_EinastoFit + 40)/TRparam[0].INCL_MAX_in_degree; // incl_ISO + 50 degree
    if(TRparam[0].incl2_for_TRfit > 0)
        TRparam[0].incl2_for_TRfit = 1;
    if(TRparam[0].incl2 > 1)
    {
        TRparam[0].incl2 = 1.0;
        TRparam[0].incl2_for_TRfit = 1.0;
    }

    TRparam[0].incl1 = 0.0;
    TRparam[0].incl1_for_TRfit = 0.0;
    TRparam[0].incl2 = 1.0;
    TRparam[0].incl2_for_TRfit = 1.0;
*/
    // INCL: Initialise sersic-polynomial function coefficients with initial INCL as given above: i.e., TRparam[0].inclF
    //for(i=0; i<TRparam[0].SersicPoly_order+1; i++)
    //{
    //  if(i==0)
    //      TRparam[0].a_sersicpoly[i] = TRparam[0].inclF_EinastoFit/TRparam[0].INCL_MAX_in_degree;
    //      //TRparam[0].a_sersicpoly[i] = TRparam[0].inclF_EinastoFit*TRparam[0].INCL_MAX_in_degree;
    //  else
    //      TRparam[0].a_sersicpoly[i] = 0.;
    //}

    // F. vrot
    vel_geo_app_side = TRparam[0].vlos_lower_limit; // update with the Gfit of LOS histogram : more robust
    vel_geo_rec_side = TRparam[0].vlos_upper_limit;
    TRparam[0].vrot1 = 0;
    TRparam[0].vrot2 = 3.0*(fabs(vel_geo_rec_side-vel_geo_app_side)/2.0)/sin(TRparam[0].inclF_EinastoFit*M_PI/180.);

    // G. vrad
    TRparam[0].vrad1 = -0.5*fabs(TRparam[0].vrot2);
    TRparam[0].vrad2 = +0.5*fabs(TRparam[0].vrot2);

    R_semi_mx = TRparam[0].ellipse_semi_mx_boxfiltered;
    init_ISO_rF = TRparam[0].ellipse_semi_mx_boxfiltered/10.;
    init_Vrot_max = TRparam[0].vrot2;


    printf("\n\t+++ ELLIPSE FIT +++\n");
    printf( "\t- X position     = %10.4f (+- %10.4f)\n", p[2], e[2]);
    printf( "\t- Y position     = %10.4f (+- %10.4f)\n", p[3], e[3]);
    printf( "\t- Position Angle = %10.4f (+- %10.4f)\n", TRparam[0].paF_EinastoFit, e[4]);
    printf( "\t- Inclination    = %10.4f (+- %10.4f)\n", p[1], e[1]);
    printf( "\t- Major axis     = %10.4f (+- %10.4f)\n\n", p[0]*2, e[0]*2);

    printf("\t+++ UPDATED INITIAL UNIFORM PRIORS OF TILTED-RING PARAMETRES +++\n");
    printf("\t- XPOS: (%.2f - %.2f)\n", TRparam[0].xpos1, TRparam[0].xpos2);
    printf("\t- YPOS: (%.2f - %.2f)\n", TRparam[0].ypos1, TRparam[0].ypos2);
    printf("\t- VSYS: (%.2f - %.2f)\n", TRparam[0].vsys1, TRparam[0].vsys2);
    printf("\t- VROT: (%.2f - %.2f)\n", TRparam[0].vrot1, TRparam[0].vrot2);
    printf("\t- VRAD: (%.2f - %.2f)\n", TRparam[0].vrad1, TRparam[0].vrad2);
    printf("\t- _n  : (%.2e - %.2e)\n", TRparam[0]._n1, TRparam[0]._n2);
    printf("\t- r_2 : (%.2e - %.2e)\n", TRparam[0].r_21, TRparam[0].r_22);
    printf("\t- rho_2 : (%.2e - %.2e)\n", TRparam[0].rho_21, TRparam[0].rho_22);
    printf("\t- sigma_factor : (%.2f - %.2f)\n\n", TRparam[0].sigma_factor1, TRparam[0].sigma_factor2);


    float *rmax_ellipse = malloc(sizeof(float)*(TRparam[0].xpos2-TRparam[0].xpos1)*(TRparam[0].ypos2-TRparam[0].ypos1));
    peri_x = -1.0*p[0]*sin(TRparam[0].paF_EinastoFit*M_PI/180.) + p[2];
    peri_y = p[0]*cos(TRparam[0].paF_EinastoFit*M_PI/180.) + p[3];

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

    free(x_temp);
    free(y_temp);
    free(rmax_ellipse);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GIPSY ellipse code
//invmat calculates the inverse of matrix. The algorithm used is the
//Gauss-Jordan algorithm described in Stoer, Numerische matematik, 1 Teil.
int invmat( double matrix[PARS][PARS], int nfree )
{
    double even;
    double hv[PARS];
    double mjk;
    double rowmax;
    int   evin;
    int   i;
    int   j;
    int   k;
    int   per[PARS];
    int   row;

    for (i = 0; i < nfree; i++) per[i] = i; /* set permutation array */
    for (j = 0; j < nfree; j++)
    {       /* in j-th column, ... */
        rowmax = fabs( matrix[j][j] );      /* determine row with ... */
        row = j;                    /* largest element. */
        for (i = j + 1; i < nfree; i++)
        {
            if (fabs( matrix[i][j] ) > rowmax)
            {
                rowmax = fabs( matrix[i][j] );
                row = i;
            }
        }
        if (matrix[row][j] == 0.0) return( -6 );    /* determinant is zero! */
        if (row > j)
        {               /* if largest element not ... */
            for (k = 0; k < nfree; k++)
            {       /* on diagonal, then ... */
                even = matrix[j][k];        /* permutate rows. */
                matrix[j][k] = matrix[row][k];
                matrix[row][k] = even;
            }
            evin = per[j];              /* keep track of permutation */
            per[j] = per[row];
            per[row] = evin;
        }
        even = 1.0 / matrix[j][j];      /* modify column */
        for (i = 0; i < nfree; i++) matrix[i][j] *= even;
        matrix[j][j] = even;
        for (k = 0; k < j; k++)
        {
            mjk = matrix[j][k];
            for (i = 0; i < j; i++) matrix[i][k] -= matrix[i][j] * mjk;
            for (i = j + 1; i < nfree; i++) matrix[i][k] -= matrix[i][j] * mjk;
            matrix[j][k] = -even * mjk;
        }
        for (k = j + 1; k < nfree; k++)
        {
            mjk = matrix[j][k];
            for (i = 0; i < j; i++) matrix[i][k] -= matrix[i][j] * mjk;
            for (i = j + 1; i < nfree; i++) matrix[i][k] -= matrix[i][j] * mjk;
            matrix[j][k] = -even * mjk;
        }
    }
    for (i = 0; i < nfree; i++)
    {       /* finally, repermute the ... */
        for (k = 0; k < nfree; k++)
        {       /* columns. */
            hv[per[k]] = matrix[i][k];
        }
        for (k = 0; k < nfree; k++)
        {
            matrix[i][k] = hv[k];
        }
    }
    return( 0 );                    /* all is well */
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//inimat sets up the matrix to be inverted in inivec. 
//inimat returns the reduced chi-squared.
double  inimat( double  s[PARS][PARS] ,
                        double  rl[PARS] ,
                        float   *x ,
                        float   *y ,
                        int n ,
                        float   *p ,
                        float   *e ,
                        int ip[PARS] ,
                        int nfree )
{
    double  chi = 0.0;          /* return value */
    double  cosi, cosp, sini, sinp;     /* sines and cosines */
    int     i, j, k;            /* counters */

    /*
     * initialize the matrix and the vector
     */
    for (j = 0; j < nfree; j++)
    {   
        rl[j] = 0.0;
        for (k = 0; k <= j; k++)
        {
            s[j][k] = 0.0;
        }
    }
    cosp = cos( F * p[4] );         /* cosine of p.a. */
    sinp = sin( F * p[4] );         /* sine of p.a. */
    cosi = cos( F * p[1] );         /* cosine of inclination */
    sini = sin( F * p[1] );         /* sine of inclination */
    for (i = 0; i < n; i++)
    {
        double  cost, sint, u;

         //Calculate the rotated x and y coordinates
        cost=( -( x[i] - p[2] ) * sinp + ( y[i] - p[3] ) * cosp ) / p[0];
        sint=(- ( x[i] - p[2] ) * cosp - ( y[i] - p[3] ) * sinp ) / p[0] / cosi;
        u = 1.0 - cost * cost - sint * sint;    /* difference with model */

        //Now calculate the partial derivatives
        e[0] = -2.0 * ( cost * cost + sint * sint ) / p[0];
        e[1] =  2.0 * F * sint * sint * sini / cosi;
        e[2] =  2.0 * ( cost * sinp + sint * cosp / cosi ) / p[0];
        e[3] = -2.0 * ( cost * cosp - sint * sinp / cosi ) / p[0];
        e[4] = -2.0 * F * sint * cost * sini * sini / cosi;
        chi = chi + u * u;          /* add to reduced chi-squared */

        //Now we fill the matrix and the vector
        for (j = 0; j < nfree; j++)
        {
            rl[j] += u * e[ip[j]];
            for (k = 0; k <= j; k++)
            {
                s[j][k] += e[ip[j]] * e[ip[k]];
            }
        }
    }
    return( chi );              /* return chi-squared */
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//inivec calculates the correction vector.
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
                        int nfree )
{
    double  cosi, cosp, sinp;       /* sines and cosines */
    int     i, j, k;            /* counters */

    //First we modify and scale the matrix.
    for (j = 0; j < nfree; j++)
    {
        double  sjj = s[j][j];          /* diagonal */

        if (sjj == 0.0) { return( -3 ); }       /* error */
        for (k = 0; k < j; k++)
        {
            double  sjk = s[j][k] / sqrt( sjj * s[k][k] );
            s1[j][k] = s1[k][j] = sjk;
        }
        s1[j][j] = 1.0 + labda;         /* new value on diagonal */
    }
    if (invmat( s1, nfree )) return( -4 );  /* error inverting matrix */
    for (i = 0; i < PARS; e[i++] = 0.0);        /* zero difference vector */
    for (j = 0; j < nfree; j++)
    {
        double  sjj = s[j][j];

        for (k = 0; k < nfree; k++)
        {
            e[ip[j]] = e[ip[j]] + rl[k] * s1[j][k] / sqrt( sjj * s[k][k] );
        }
    }
    for (i = 0; i < PARS; i++)
    {
        e[i] += p[i];
    }
    (*q) = 0.0;
    cosp = cos( F * e[4] );         /* cosine of p.a. */
    sinp = sin( F * e[4] );         /* sine of p.a. */
    cosi = cos( F * e[1] );         /* cosine of inclination */
    for (i = 0; i < n; i++)
    {
        double  cost, sint, u;

        //Calculate the rotated x and y coordinates
        cost=( -( x[i] - e[2] ) * sinp + ( y[i] - e[3] ) * cosp ) / e[0];
        sint=( -( x[i] - e[2] ) * cosp - ( y[i] - e[3] ) * sinp ) / e[0] / cosi;
        u = 1.0 - cost * cost - sint * sint;    /* difference with model */
        (*q) += u * u;              /* add to chi-squared */
    }
    return( 0 );                    /* return to caller */
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
              ellipse1.dc2
Purpose:      Routine which gives an initial estimate for an ellipse to N
              points in the array X,Y.
Use:          INTEGER ELLIPSE1( N ,         Input      INTEGER
                                X ,         Input      REAL ARRAY
                                Y ,         Input      REAL ARRAY
                                P )         Output     REAL ARRAY

              ELLIPSE1   Returns 0 when successfull, and 1 on error 
              N          Number of points in X and Y.
              X          Array with X-coordinates.
              Y          Array with Y-coordinates.
              P          Estimated ellipse parameters:
                         P(1) = radius
                         P(2) = inclination (0..90 degrees)
                         P(3) = X0 (centre)
                         P(4) = Y0 (centre)
                         P(5) = position angle of major axis
                                (0..180 degrees) w.r.t. Y axis.
*/
int ellipse1_c(double *n ,      /* number of points */
                    float *x ,      /* X coordinates */
                    float *y ,      /* Y coordinates */
                    float *p )      /* ellipse parameters */
{
    double  ellips[PARS];           /* ellipse parameters */
    double  matrix[PARS][PARS];     /* the matrix */
    double  vector[PARS];           /* the vector */
    double  xc, yc;             /* centre of ellipse */
    int     i, j, k;            /* loop counters */

    //Due to stability problems with ellipses
    //where the origin is far away from the centre
    //of the ellipse we first have to find an approximate

    for (xc = yc = 0.0, k = 0; k < (*n); k++)
    {   /* loop over all positions */
        xc += x[k];             /* sum x coordinates */
        yc += y[k];             /* sum y coordinates */
    }
    xc /= (double) (*n);                /* mean x position */
    yc /= (double) (*n);                /* mean y position */

    //Next: initialize the matrix to be inverted by invmat and
    //the vector
    for (i = 0; i < PARS; i++)
    {
        for (j = 0; j < PARS; j++)
        {
            matrix[i][j] = 0.0;
        }
        vector[i] = 0.0;
    }

    //Then: fill matrix
    for (k = 0; k < (*n); k++)
    {
        double  xs = x[k] - xc;         /* x coord. w.r.t. centre */
        double  ys = y[k] - yc;         /* y coord. w.r.t. centre */
        double  xx;             /* x * x */
        double  xy;             /* cross product */
        double  yy;             /* y * y */

        xx = xs * xs;               /* x * x */
        xy = xs * ys;               /* cross product */
        yy = ys * ys;               /* y * y */
        matrix[0][0] += xx * xx;
        matrix[0][1] += 2.0 * xx * xy;
        matrix[0][2] += xy * xy;
        matrix[0][3] += xx * xs;
        matrix[0][4] += xy * xs;
        matrix[1][1] += 4.0 * xy * xy;
        matrix[1][2] += 2.0 * xy * yy;
        matrix[1][3] += 2.0 * xy * xs;
        matrix[1][4] += 2.0 * xy * ys;
        matrix[2][2] += yy * yy;
        matrix[2][3] += xy * ys;
        matrix[2][4] += yy * ys;
        matrix[3][3] += xx;
        matrix[3][4] += xy;
        matrix[4][4] += yy;
        vector[0] += xx;
        vector[1] += 2.0 * xy;
        vector[2] += yy;
        vector[3] += xs;
        vector[4] += ys;
    }

    //make the matrix symmetric since this saves a lot of typing

    for (i = 0; i < PARS; i++)
    {
        for (j = 0; j < i; j++)
        {
            matrix[i][j] = matrix[j][i];
        }
    }

    //invert the matrix
    if (invmat( matrix, PARS ))
    {
        return( 1 );                /* cannot invert matrix */
    }
    
    //solve MATRIX * ELLIPS = VECTOR
    for (i = 0; i < PARS; i++)
    {
        ellips[i] = 0.0;                /* reset this ellipse parm. */
        for (j = 0; j < PARS; j++)
        {
            ellips[i] += matrix[i][j] * vector[j];
        }
    }
     //NOTE: The ellipse equation taken is
     //AA.x^2 + 2.BB.x.y + CC.y^2 + DD.x + EE.y = 1
     //Where AA..EE were now solved for
     //from which now the ellipse-parameters are derived
    {
        double  aa = ellips[0];
        double  bb = ellips[1];
        double  cc = ellips[2];
        double  dd = ellips[3];
        double  ee = ellips[4];
        double  pa, pp;
        double  cospa, sinpa, sinpp;
        double  s1, s2, s3, y1, y2, y3;
        double  x0, y0;
        double  ab, al, r;

        pp = atan( 2.0 * bb / ( aa - cc ) );    /* estimate of position angle */
        pa = 0.5 * pp;              /* p.a. of an (UNDETERMINED) axis */
        cospa = cos( pa );          /* cosine */
        sinpp = sin( pp );          /* sine of double angle */
        sinpa = sin( pa );          /* sine */
        al = 2.0 * bb / sinpp / ( aa + cc );    /* auxiliary */
        r = sqrt( ( 1.0 + al ) / ( 1.0 - al ) );    /* axial ratio (OR ITS RECIPROCAL) */
        s1 = bb * ee - cc * dd;         /* three other auxiliaries */
        s2 = bb * dd - aa * ee;
        s3 = aa * cc - bb * bb;
        x0 = s1 / s3;               /* X-centre of ellipse */
        y0 = s2 / s3;               /* Y-centre of ellipse */
        y1 = sinpa * sinpa + r * r * cospa * cospa;
        y2 = x0 * y0 * ( r * r - 1.0 ) * sinpp;
        y3 = y0 * y0 * ( cospa * cospa + r * r * sinpa * sinpa );
                        /* length of (yet undetermined) axis (A or B) */
        ab = sqrt( y1 * ( x0 * x0 + 1.0 / aa ) + y2 + y3 );
        pa = pa * 45.0 / atan( 1.0 );       /* convert to degrees */
    
        //determination which axis is the long axis
        if (r < 1.0)
        {               /* the other is */
            ab /= r;                /* this one is the major axis */
            pa -= 90.0;             /* so change position angle */
        } else
        {                   /* the right one is */
            r = 1.0 / r;                /* so change axial ratio */
        }
        if (pa < 0.0) pa += 180.0;      /* position angle in range 0.....180 */
        p[0] = ab;              /* radius */
        p[1] = acos( r ) * 45.0 / atan( 1.0 );  /* inclination assuming projected circle */
        p[2] = x0 + xc;             /* new x-position of centre */
        p[3] = y0 + yc;             /* new y-position of centre */
        p[4] = pa;              /* position angle major axis */
    }
    return( 0 );                    /* return to caller */
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
#>            ellipse2.dc2
Function:     ELLIPSE2
Purpose:      Fits an ellipse to a set of X and Y positions.
Use:          INTEGER ELLIPSE2( N ,            Input       INTEGER
                                X ,            Input       REAL ARRAY
                                Y ,            Input       REAL ARRAY
                                P ,        Input/Output    REAL ARRAY
                                E ,            Output      REAL ARRAY
                                M )            Input       INTEGER ARRAY

              ELLIPSE2     Returns number of iterations on success, else
                           -1: No free parameters
                           -2: Exceede iteration limit (50)
                           -3: Diagonal of matrix has zeroes
                           -4: Matrix could not be inverted
              N            Number of positions.
              X            Array with X coordinates (in units).
              Y            Array with Y coordinates (in units).
              P            On input contains the initial estimates
                           of the ellipse parameters. On output the
                           fitted parameters. P contains:
                           P(1) = radius (in units)
                           P(2) = inclination (degrees)
                           P(3) = X centre of ellipse (units)
                           P(4) = Y centre of ellipse (units)
                           P(5) = Position angle (degrees)
              E            Contains the errors in the fitted parameters.
              M            Mask for free (1) or fixed (0) parameters.
#<
*/
int ellipse2_c(double *n ,      /* number of coordinates */
                    float *x ,      /* X coordinates */
                    float *y ,      /* Y coordinates */
                    float *p ,      /* ellipse parameters */
                    float *e ,      /* errors in ellipse parms. */
                    int *m)     /* fixed/free mask */
{
    double       chi;                            /* red. chi-squared */
    double  labda;              /* mixing parameter */
    double       q;                              /* red. chi-squared */
    double  rl[PARS];           /* vector */
    double  s[PARS][PARS], s1[PARS][PARS];  /* matrices */
    int     h = 0;              /* iteration counter */
    int         i;                              /* counter */
    int     ip[PARS];           /* permutation array */
    int     nfree = 0;          /* number of free parameters */
    int     r = 0;              /* return value */

    labda = LAB * FAC;              /* start value */
    for (i = 0; i < PARS; i++)
    {
        if (m[i]) ip[nfree++] = i;      /* fit this parameter */
    }
    if (nfree == 0 || nfree >= (*n)) return( -1 );
    do{                     /* iteration loop */
        if (++h > T) { r = -2; break; }     /* too many iterations */
        chi = inimat( s, rl, x, y, (*n), p, e, ip, nfree );
        if (labda > LABMIN) labda /= FAC;       /* new labda */
        r = inivec( s, s1, rl, labda, &q, p, e, x, y, (*n), ip, nfree );
        if (r) break;                   /* error from inivec */
        while (q >= chi)
        {               /* interpolation loop */
           if (labda > LABMAX) break;       /* leave loop */
           labda *= FAC;                /* new labda */
           r = inivec( s, s1, rl, labda, &q, p, e, x, y, (*n), ip, nfree );
           if (r) break;                /* error from inivec */
        }
        if (labda <= LABMAX)
        {
            for (i = 0; i < PARS; i++) p[i] = e[i];
        }
        if (fabs( chi - q ) / q < TOL || labda > LABMAX)
        {
            labda = 0.0;                /* for Taylor solution */
            chi = inimat( s, rl, x, y, (*n), p, e, ip, nfree );
            r = inivec( s, s1, rl, labda, &q, p, e, x, y, (*n), ip, nfree );
            if (!r)
            {
                r = h;              /* number of iterations */
                for (i = 0; i < PARS; i++)
                {
                    p[i] = e[i];            /* the parameters */
                    e[i] = 0.0;         /* reset errors */
                }
                q = sqrt( q / (double) ( (*n) - nfree ) );
                for (i = 0; i < nfree; i++)
                {
                    e[ip[i]] = q * sqrt( s1[i][i] / s[i][i] );
                }
            }
        }
    } while(!r);                /* until error or finished */
    return(r);                  /* return to caller */
}


// --- End of line



