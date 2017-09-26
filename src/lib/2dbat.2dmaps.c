#include "2dbat.2dmaps.h"

// 2DBAT user defined functions
// 2D MAPS related

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// find the largest connected blob
void find_the_largest_connected_area(TR_ringParameters *TRparam, int x0, int y0, int x1, int y1, int histogram, int use_allPixels)
{
    int i, j, k, blob_index_largest_area, ib, jb;
    int mi, mj, mn, mm, n;
    int decimX, decimY;
    int box_x, box_y;
    int bi=0, bn=0;
    int *blob_index_marked, *blob_index_marked_for_qsort;

    double label_nn, label_ww, blob_index=0, N_nonblanks=0, N_pixels_in_connected_blob=0, id, vlos_allpixels;
    double *LOS_velocities, *LOS_velocities_err, LOS_vel_hist_rbm, LOS_vel_hist_std;
    double *_filterbox;
    double lower_bound_LOS_velocities, upper_bound_LOS_velocities;
    double hist_mean_LOS, hist_std_LOS;
    double hist_mean_filterbox, hist_std_filterbox;
    double gsl_mean_vlos, gsl_std_vlos, gsl_max_vlos, gsl_min_vlos;

// 1. Initialise with blanks
    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            HI_VF_boxFiltered[0].data[j][i] = 1E90;
            HI_VF_boxFiltered_sigma[0].data[j][i] = 1E90;
            HI_VF_boxFiltered_SN[0].data[j][i] = 1E90;
            HI_VF_boxFiltered_decim0[0].data[j][i] = 1E90;
            HI_VF_boxFiltered_decim_user[0].data[j][i] = 1E90;
            HI_VF_temp[0].data[j][i] = 1E90;
        }
    }

    // median filtering
    decimX = 0;
    decimY = 0;

    n = 0;
    box_x = TRparam[0].box_x;
    box_y = TRparam[0].box_y;
    _filterbox = malloc(sizeof(double) * box_x*box_y);
    x0 = 0;
    y0 = 0;
    if(TRparam[0].nax1 < TRparam[0].nax2) // IF nax1 != nax2
    {
        x1 = TRparam[0].nax1;
        y1 = TRparam[0].nax1;
    }
    else if(TRparam[0].nax1 > TRparam[0].nax2) // IF nax1 != nax2
    {
        x1 = TRparam[0].nax2;
        y1 = TRparam[0].nax2;
    }
    else // IF nax1 == nax2
    {
        x1 = TRparam[0].nax1;
        y1 = TRparam[0].nax2;
    }

    // Reset the connected area to fit
    for(id=0; id<4024*4024; id++)
    {
        TRparam[0].tilted_ring[(long int)id][0] = 0;
        TRparam[0].tilted_ring[(long int)id][1] = 0;
    }
    // Calculate vlos_mean & vlos_std
    id = 0;
    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            if(!isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i]) && fabs(HI_VF[0].data[j][i]) > 1E-7) // no blank
            {
                HI_VF_temp[0].data[j][i] = HI_VF[0].data[j][i];
                TRparam[0].tilted_ring[(long int)id][0] = i;
                TRparam[0].tilted_ring[(long int)id][1] = j;
                id += 1;
            }
            else
            {
                HI_VF_temp[0].data[j][i] = 1E90;
            }
        }
    }
    vlos_allpixels = id; // save the number of all the available pixels for vlos statistics
               
    if(histogram == 1) // this is for deriving VLOS statistics
    {
        LOS_velocities = malloc(sizeof(double)*(long int)id);
        LOS_velocities_err = malloc(sizeof(double)*(long int)id);
        for(id=0; id<vlos_allpixels; id++)
        {
            LOS_velocities[(long int)id] = HI_VF_temp[0].data[TRparam[0].tilted_ring[(long int)id][1]][TRparam[0].tilted_ring[(long int)id][0]];
            LOS_velocities_err[(long int)id] = HI_VF_sigma[0].data[TRparam[0].tilted_ring[(long int)id][1]][TRparam[0].tilted_ring[(long int)id][0]];
        }
        robust_mean_std(LOS_velocities, vlos_allpixels, &hist_mean_LOS, &hist_std_LOS);
        robust_mean_std_histogram_ac(LOS_velocities, LOS_velocities_err, (long int)vlos_allpixels, &gsl_mean_vlos, &gsl_std_vlos);

        TRparam[0].LOS_vel_hist_rbm = gsl_mean_vlos;
        TRparam[0].LOS_vel_hist_std = gsl_std_vlos;
        TRparam[0].vlos_lower = gsl_mean_vlos - 5*gsl_std_vlos;
        TRparam[0].vlos_upper = gsl_mean_vlos + 5*gsl_std_vlos;

        TRparam[0].vlos_lower_limit = gsl_mean_vlos - 5*gsl_std_vlos;
        TRparam[0].vlos_upper_limit = gsl_mean_vlos + 5*gsl_std_vlos;

        TRparam[0].LOS_hist_Gfit_V0 = gsl_mean_vlos;
        TRparam[0].LOS_hist_Gfit_sigma = gsl_std_vlos;

        free(LOS_velocities);
        free(LOS_velocities_err);

        printf("\n\t++Found lower and upper values of LOS velocities++\n");
        printf("\t- Vlos_lower  = %10.4f\tVlos_upper = %f\n\n", TRparam[0].vlos_lower, TRparam[0].vlos_upper);
        printf("\t++Robust mean & std of LOS velocities++\n");
        printf("\t- R-mean      = %10.4f\tSTD        = %f\n\n", TRparam[0].LOS_vel_hist_rbm, TRparam[0].LOS_vel_hist_std);
        //printf("\t++Gaussian fit of the histogram of LOS velocities++\n");
        //printf("\t- Vlos_centre = %10.4f\tVlos_sigma = %f\n\n", TRparam[0].LOS_hist_Gfit_V0, TRparam[0].LOS_hist_Gfit_sigma);
        printf("\t- Vlos_centre_gsl = %10.4f\tVlos_sigma_gsl = %f\n\n", gsl_mean_vlos, gsl_std_vlos);
    }


    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            // 1. filtering input VF
            if(!isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i]) && HI_VF[0].data[j][i] > TRparam[0].vlos_lower && HI_VF[0].data[j][i] < TRparam[0].vlos_upper) // no blank
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
                        if(i+mi < 0 || i+mi >= TRparam[0].nax1 || j+mj < 0 || j+mj >= TRparam[0].nax2) continue;
                        if(!isinf(HI_VF[0].data[j+mj][i+mi]) && !isnan(HI_VF[0].data[j+mj][i+mi]) && HI_VF[0].data[j+mj][i+mi] != 0.0) // if no blank
                        {
                            _filterbox[bn] = HI_VF[0].data[j+mj][i+mi];
                            bn++;
                        }
                    }
                }
                robust_mean_std(_filterbox, bn, &hist_mean_filterbox, &hist_std_filterbox);

                if(HI_VF[0].data[j][i] > hist_mean_filterbox-3*hist_std_filterbox && HI_VF[0].data[j][i] < hist_mean_filterbox+3*hist_std_filterbox)
                {
                    HI_VF_boxFiltered[0].data[j][i] = HI_VF[0].data[j][i];
                }
                else
                {
                    //HI_VF_boxFiltered[0].data[j][i] = hist_mean_filterbox; // replace with the mean value derived
                    HI_VF_boxFiltered[0].data[j][i] = 1E90; // replace with the mean value derived
                }
            }
            else
            {
                HI_VF_boxFiltered_decim0[0].data[j][i] = -1.0;
                //HI_VF_boxFiltered[0].data[j][i] = 1E90;
            }

            // 2. filtering sigma VF
            if(!isinf(HI_VF_sigma[0].data[j][i]) && !isnan(HI_VF_sigma[0].data[j][i])) // no blank
            {
                // extract _boxfilter
                bn = 0;
                for(mi=-(box_x-1)/2; mi<(box_x+1)/2; mi++)
                {
                    for(mj=-(box_y-1)/2; mj<(box_y+1)/2; mj++)
                    {
                        if(i+mi < 0 || i+mi >= TRparam[0].nax1 || j+mj < 0 || j+mj >= TRparam[0].nax2) continue;
                        if(!isinf(HI_VF_sigma[0].data[j+mj][i+mi]) && !isnan(HI_VF_sigma[0].data[j+mj][i+mi]) && HI_VF_sigma[0].data[j+mj][i+mi] != 0.0) // if no blank
                        {
                            _filterbox[bn] = HI_VF_sigma[0].data[j+mj][i+mi];
                            bn++;
                        }
                    }
                }
                robust_mean_std(_filterbox, bn, &hist_mean_filterbox, &hist_std_filterbox);


                if(HI_VF_sigma[0].data[j][i] > hist_mean_filterbox-3*hist_std_filterbox && HI_VF_sigma[0].data[j][i] < hist_mean_filterbox+3*hist_std_filterbox)
                {
                    HI_VF_boxFiltered_sigma[0].data[j][i] = HI_VF_sigma[0].data[j][i];
                }
                else
                {
                    //HI_VF_boxFiltered_sigma[0].data[j][i] = hist_mean_filterbox; // replace with the mean value derived
                    HI_VF_boxFiltered_sigma[0].data[j][i] = 1E90; // replace with the mean value derived
                }
            }
            else
            {
                HI_VF_boxFiltered_sigma[0].data[j][i] = 1E90;
            }

            // 3. filtering flux (or s/n) VF
            if(!isinf(HI_VF_sn[0].data[j][i]) && !isnan(HI_VF_sn[0].data[j][i])) // no blank
            {
                // extract _boxfilter
                bn = 0;
                for(mi=-(box_x-1)/2; mi<(box_x+1)/2; mi++)
                {
                    for(mj=-(box_y-1)/2; mj<(box_y+1)/2; mj++)
                    {
                        if(i+mi < 0 || i+mi >= TRparam[0].nax1 || j+mj < 0 || j+mj >= TRparam[0].nax2) continue;
                        if(!isinf(HI_VF_sn[0].data[j+mj][i+mi]) && !isnan(HI_VF_sn[0].data[j+mj][i+mi]) && HI_VF_sn[0].data[j+mj][i+mi] != 0.0) // if no blank
                        {
                            _filterbox[bn] = HI_VF_sn[0].data[j+mj][i+mi];
                            bn++;
                        }
                    }
                }
                robust_mean_std(_filterbox, bn, &hist_mean_filterbox, &hist_std_filterbox);


                if(HI_VF_sn[0].data[j][i] > hist_mean_filterbox-3*hist_std_filterbox && HI_VF_sn[0].data[j][i] < hist_mean_filterbox+3*hist_std_filterbox)
                {
                    HI_VF_boxFiltered_SN[0].data[j][i] = HI_VF_sn[0].data[j][i];
                }
                else
                {
                    //HI_VF_boxFiltered_SN[0].data[j][i] = hist_mean_filterbox; // replace with the mean value derived
                    HI_VF_boxFiltered_SN[0].data[j][i] = 1E90; // replace with the mean value derived
                }
            }
            else
            {
                HI_VF_boxFiltered_SN[0].data[j][i] = 1E90; 
            }
            j += decimY;
        }
        i += decimX; 
    }

// 2. make a boxfiltered map (-1.0 is blank!)
    decimX = 0;
    decimY = 0;
    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            if(!isinf(HI_VF_boxFiltered[0].data[j][i]) && !isnan(HI_VF_boxFiltered[0].data[j][i]) && HI_VF[0].data[j][i] > TRparam[0].vlos_lower && HI_VF[0].data[j][i] < TRparam[0].vlos_upper) // no blank
            {
                HI_VF_boxFiltered_decim0[0].data[j][i] = 999999999.0; // no blank
            }
            else
            {
                HI_VF_boxFiltered_decim0[0].data[j][i] = -1.0; // blank
            }
            j += decimY;
        }
        i += decimX; 
    }


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // I. decim user
    // find the largest connected area
    blob_index=0;
    decimX = TRparam[0].decimX_einasto_halofit;
    decimY = TRparam[0].decimY_einasto_halofit;
    for(j=y1-1-TRparam[0].decim_y0; j>=y0; j--)
    {
        for(i=x0+TRparam[0].decim_x0; i<x1; i++)
        {
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i-1-decimX >= x0)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j+1+decimY < y1)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = label_ww;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = label_nn;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                }
            }
            i += decimX;
        }
        j -= decimY; 
    }

    for(j=y0+TRparam[0].decim_y0; j<y1; j++)
    {
        for(i=x1-1-TRparam[0].decim_x0; i>=x0; i--)
        {
            if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i+1+decimX < x1)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j-1-decimY >= y0)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = label_ww;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = label_nn;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                }
            }
            i -= decimX; 
        }
        j += decimY;
    }

    for(j=y0+TRparam[0].decim_y0; j<y1; j++)
    {
        for(i=x0+TRparam[0].decim_x0; i<x1; i++)
        {
            if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i-1-decimX >= x0)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j-1-decimY >= y0)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = label_ww;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = label_nn;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                }
            }
            i += decimX;
        }
        j += decimY;
    }

    for(j=y1-1-TRparam[0].decim_y0; j>=y0; j--)
    {
        for(i=x1-1-TRparam[0].decim_x0; i>=x0; i--)
        {
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i+1+decimX < x1)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j+1+decimY < y1)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX])
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    }
                    else
                    {
                        HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                    }
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = label_ww;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = label_nn;

                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];

                        }
                    }
                }
            }
            i -= decimX; 
        }
        j -= decimY; 
    }

    blob_index_marked = malloc(sizeof(int)*blob_index);
    blob_index_marked_for_qsort = malloc(sizeof(int)*blob_index);
    for(k=0; k<blob_index; k++)
    {
        blob_index_marked[k] = 0;
        blob_index_marked_for_qsort[k] = 0;
    }

    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            for(k=0; k<blob_index; k++)
            {
                if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
                {
                    if(fabs(HI_VF_boxFiltered_decim0[0].data[j][i] - k) < 0.01)
                    {
                        blob_index_marked[k] += 1;
                        blob_index_marked_for_qsort[k] += 1;

                    }
                }
            }
        }
    }

    qsort(blob_index_marked_for_qsort, blob_index, sizeof(blob_index_marked_for_qsort[0]), comparisonFunctionFloat);
    for(k=0; k<blob_index; k++)
    {
        if(fabs(blob_index_marked[k] - blob_index_marked_for_qsort[(int)blob_index-1]) < 0.01)
        {
            blob_index_largest_area = k;
        }
    }

    // The largest connected area found
    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            if(fabs(HI_VF_boxFiltered_decim0[0].data[j][i] - blob_index_largest_area) < 0.01) // the largest connected area
            {
                HI_VF_boxFiltered_decim_user[0].data[j][i] = HI_VF_boxFiltered[0].data[j][i];
                // save the coordinates of the largest connected blob temporarily for the histogram analysis below
            }
            else
            {
                HI_VF_boxFiltered_decim_user[0].data[j][i] = 1E90;
            }
        }
    }


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // II. decim user
    // find the largest connected area
    blob_index=0;
    N_pixels_in_connected_blob = 0;
    decimX = 0;
    decimY = 0;
    for(j=y1-1; j>=y0; j--)
    {
        for(i=x0; i<x1; i++)
        {
            if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i-1-decimX >= x0)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j+1+decimY < y1)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = label_ww;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = label_nn;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                }
            }
            i += decimX;
        }
        j -= decimY; 
    }

    for(j=y0; j<y1; j++)
    {
        for(i=x1-1; i>=x0; i--)
        {
            if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i+1+decimX < x1)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j-1-decimY >= y0)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = label_ww;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = label_nn;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                }
            }
            i -= decimX; 
        }
        j += decimY;
    }

    for(j=y0; j<y1; j++)
    {
        for(i=x0; i<x1; i++)
        {
            if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i-1-decimX >= x0)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j-1-decimY >= y0)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = label_ww;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = label_nn;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i-1-decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j-1-decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                }
            }
            i += decimX;
        }
        j += decimY;
    }

    for(j=y1-1; j>=y0; j--)
    {
        for(i=x1-1; i>=x0; i--)
        {
            if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
            if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
            {
                // 1. the left to the current pixel
                if(i+1+decimX < x1)
                    label_ww = HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX]; 
                else
                    label_ww = -1.0;

                // 2. the upper to the current pixel
                if(j+1+decimY < y1)
                    label_nn = HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i]; 
                else
                    label_nn = -1.0;

                if(label_ww == -1.0 && label_nn == -1.0)
                {
                    blob_index += 1.0;
                    HI_VF_boxFiltered_decim0[0].data[j][i] = blob_index;
                }
                else if(label_ww == -1.0 && label_nn != -1.0)
                {
                    if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn == -1.0)
                {
                    if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX])
                        HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                    else
                        HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                }
                else if(label_ww != -1.0 && label_nn != -1.0)
                {
                    if(label_ww <= label_nn)
                    {
                        if(label_ww <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_ww;
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = label_ww;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                    else
                    {
                        if(label_nn <= HI_VF_boxFiltered_decim0[0].data[j][i])
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i] = label_nn;
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = label_nn;
                        }
                        else
                        {
                            HI_VF_boxFiltered_decim0[0].data[j][i+1+decimX] = HI_VF_boxFiltered_decim0[0].data[j][i];
                            HI_VF_boxFiltered_decim0[0].data[j+1+decimY][i] = HI_VF_boxFiltered_decim0[0].data[j][i];
                        }
                    }
                }
            }
            i -= decimX; 
        }
        j -= decimY; 
    }

    blob_index_marked = malloc(sizeof(int)*blob_index);
    blob_index_marked_for_qsort = malloc(sizeof(int)*blob_index);
    for(k=0; k<blob_index; k++)
    {
        blob_index_marked[k] = 0;
        blob_index_marked_for_qsort[k] = 0;
    }

    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            for(k=0; k<blob_index; k++)
            {
                if(HI_VF_boxFiltered_decim0[0].data[j][i] != -1.0) // no blank
                {
                    if(fabs(HI_VF_boxFiltered_decim0[0].data[j][i] - k) < 0.01)
                    {
                        blob_index_marked[k] += 1;
                        blob_index_marked_for_qsort[k] += 1;

                    }
                }
            }
        }
    }

    qsort(blob_index_marked_for_qsort, blob_index, sizeof(blob_index_marked_for_qsort[0]), comparisonFunctionFloat);
    for(k=0; k<blob_index; k++)
    {
        if(fabs(blob_index_marked[k] - blob_index_marked_for_qsort[(int)blob_index-1]) < 0.01)
        {
            blob_index_largest_area = k;
        }
    }

    // Reset the connected area to fit
    N_pixels_in_connected_blob=0;
    for(id=0; id<4024*4024; id++)
    {
        TRparam[0].tilted_ring[(long int)id][0] = 0;
        TRparam[0].tilted_ring[(long int)id][1] = 0;
    }

    // The largest connected area found
    N_pixels_in_connected_blob=0;
    for(i=x0; i<x1; i++)
    {
        for(j=y0; j<y1; j++)
        {
            if(fabs(HI_VF_boxFiltered_decim0[0].data[j][i] - blob_index_largest_area) < 0.01) // the largest connected area
            {
                HI_VF_boxFiltered_decim0[0].data[j][i] = HI_VF_boxFiltered[0].data[j][i];
                // save the coordinates of the largest connected blob temporarily for the histogram analysis below
                TRparam[0].tilted_ring[(long int)N_pixels_in_connected_blob][0] = i;
                TRparam[0].tilted_ring[(long int)N_pixels_in_connected_blob][1] = j;
                N_pixels_in_connected_blob += 1;
            }
            else
            {
                HI_VF_boxFiltered_decim0[0].data[j][i] = 1E90;
                HI_VF_boxFiltered[0].data[j][i] = 1E90;
                HI_VF_boxFiltered_sigma[0].data[j][i] = 1E90;
                HI_VF_boxFiltered_SN[0].data[j][i] = 1E90;
            }
        }
    }



    if(TRparam[0].use_allPixels == 1) // pixels with VLOS outside the range between vlos_lower and vlos_upper derived above are being filtered out
    {
        // Reset the connected area to fit
        for(id=0; id<4024*4024; id++)
        {
            //TRparam[0].tilted_ring[(int)N_pixels_in_connected_blob][0] = 0;
            //TRparam[0].tilted_ring[(int)N_pixels_in_connected_blob][1] = 0;
            //N_pixels_in_connected_blob += 1;
            TRparam[0].tilted_ring[(int)id][0] = 0;
            TRparam[0].tilted_ring[(int)id][1] = 0;
        }

        N_pixels_in_connected_blob=0;
        // decimX=0 && decimY=0
        for(i=x0; i<x1; i++)
        {
            for(j=y0; j<y1; j++)
            {
                if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
                // HI_VF_temp[0].data[j][i] : the largest connected area found as above
                if(!isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i]) && HI_VF[0].data[j][i] > TRparam[0].vlos_lower && HI_VF[0].data[j][i] < TRparam[0].vlos_upper) // no blank
                {
                    HI_VF_boxFiltered_decim0[0].data[j][i] = HI_VF[0].data[j][i];
                    TRparam[0].tilted_ring[(int)N_pixels_in_connected_blob][0] = i;
                    TRparam[0].tilted_ring[(int)N_pixels_in_connected_blob][1] = j;
                    N_pixels_in_connected_blob += 1;
                }
                else
                {
                    HI_VF_boxFiltered_decim0[0].data[j][i] = 1E90;
                    HI_VF_boxFiltered[0].data[j][i] = 1E90;
                    HI_VF_boxFiltered_sigma[0].data[j][i] = 1E90;
                    HI_VF_boxFiltered_SN[0].data[j][i] = 1E90;
                }
            }
        }

        // decimX=user && decimY=user
        decimX = TRparam[0].decimX_einasto_halofit;
        decimY = TRparam[0].decimY_einasto_halofit;

        // 1. Initialise with blanks
        for(i=x0; i<x1; i++)
        {
            for(j=y0; j<y1; j++)
            {
                HI_VF_boxFiltered_decim_user[0].data[j][i] = 1E90;
            }
        }
        for(i=x0; i<x1; i++)
        {
            for(j=y0; j<y1; j++)
            {
                if(j<0 || j>=TRparam[0].nax1 || i<0 || i >=TRparam[0].nax1) continue;
                // HI_VF_temp[0].data[j][i] : the largest connected area found as above
                if(!isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i]) && HI_VF[0].data[j][i] > TRparam[0].vlos_lower && HI_VF[0].data[j][i] < TRparam[0].vlos_upper) // no blank
                {
                    HI_VF_boxFiltered_decim_user[0].data[j][i] = HI_VF[0].data[j][i];
                }
                else
                {
                    HI_VF_boxFiltered_decim_user[0].data[j][i] = 1E90;
                }
                j += decimY; 
            }
            i += decimX; 
        }
    }

    free(_filterbox);
    free(blob_index_marked);
    free(blob_index_marked_for_qsort);
    return;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void reset_HI_VF_weight_trfit(TR_ringParameters *TRparam)
{
    int i0=0, j0=0;

    // reset
    for(i0=0; i0<TRparam[0].nax1; i0++)
    {
        for(j0=0; j0<TRparam[0].nax2; j0++)
        {
            HI_VF_weight_TRfit[0].data[j0][i0] = 0.0;
        }
    }
    return;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void reset_HI_VF_fract_navail_nall(TR_ringParameters *TRparam)
{
    int i0=0, j0=0;

    // reset
    for(i0=0; i0<TRparam[0].nax1; i0++)
    {
        for(j0=0; j0<TRparam[0].nax2; j0++)
        {
            HI_VF_fract_navail_nall[0].data[j0][i0] = 1E90; // no weight
        }
    }
    return;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int save_2dmapsfits(fitsfile *fptr1, char *inputfile, char *outputfile, int nax1, int nax2, velocity_field *twodmap_array)
{
    /* perfect single gaussian */
    int group=1;
    int status=0;
    int nelements=0;
    long fpixel=1;
    FILE *file;
    fitsfile *fptr_out; 

    group=1; status=0;
    nelements = nax1*nax2;

    // remove the output fits if it exists.
    if((file = fopen(outputfile, "r")) != NULL)
    {
        remove(outputfile);
    }

    // open input vf fits
    fits_open_file(&fptr1, inputfile, READONLY, &status);
    // create output fits file
    fits_create_file(&fptr_out, outputfile, &status);
    // copy header
    fits_copy_header(fptr1, fptr_out, &status);
    fits_close_file(fptr_out, &status);

    // write results
    fits_open_file(&fptr_out, outputfile, READWRITE, &status);
    fits_write_img(fptr_out, TFLOAT, fpixel, nelements, twodmap_array[0].data[0], &status);

    fits_close_file(fptr1, &status);
    fits_close_file(fptr_out, &status);
    return 0;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int make_model_vfs(TR_ringParameters *TRparam, fitsfile *fptr1, char *inputfile, int nax1, int nax2, filename_2dbat *fname)
{
    int i=0, j=0, pa_n=0, incl_n=0, n=0, k=0;
    int ring_i=0;
    int status=0;
    double r1, r2, GG, e;
    double _n, r_2, rho_2, _x, iGamma;
    int n_res_input_TR=0;
    double mean_res_input_TR, std_res_input_TR;

    GG = 4.302*pow(10., -3); // pc Msol^-1 (km/s)^2
    e = 2.71828;

    double x_pixel_from_xpos, y_pixel_from_ypos, r_pixel_from_centre;
    double d1, d2, w1, w2, ri, ro, ring00, ring01, ring_outermost, ring_semi_mx;
    double incl01, posa01, rads01, vrad01, vrot01, vsys01, cost, sint;
    double vrot_Einasto_fit00, vrot_Einasto_fit01, vrot01_einasto;
    double ring00_pc, ring01_pc;

    // set outermost radius : 3 x semi-mx
    ri = TRparam[0].ring_s + (3*TRparam[0].Nrings-1)*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
    ro = TRparam[0].ring_s + (3*TRparam[0].Nrings-1)*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
    ring_outermost = (ri+ro)/2.0;

    for(i=0; i<TRparam[0].nax1; i++)
        for(j=0; j<TRparam[0].nax2; j++)
            HI_VF_res_input_minus_trfit[0].data[j][i] = 1E90; // blank


    for(i=0; i<TRparam[0].nax1; i++)
    {
        for(j=0; j<TRparam[0].nax2; j++)
        {
            x_pixel_from_xpos = (i - TRparam[0].xposF_EinastoFit);
            y_pixel_from_ypos = (j - TRparam[0].yposF_EinastoFit);
            r_pixel_from_centre = sqrt(x_pixel_from_xpos*x_pixel_from_xpos + y_pixel_from_ypos*y_pixel_from_ypos);
            if(r_pixel_from_centre > ring_outermost) // if radius is outside the outermost ring
            {
                HI_VF_tr_model[0].data[j][i] = 1E90; // blank
                HI_VF_einasto_halomodel[0].data[j][i] = 1E90; // blank
                HI_VF_res_input_minus_trfit[0].data[j][i] = 1E90; // blank
                HI_VF_res_input_minus_einastofit[0].data[j][i] = 1E90; // blank
                HI_VF_res_trfit_minus_einastofit[0].data[j][i] = 1E90; // blank
            }
            else
            {
                d1 = radius_galaxy_plane(0, x_pixel_from_xpos, y_pixel_from_ypos, TRparam[0].pa[0], TRparam[0].incl[0]);
                d2 = radius_galaxy_plane(ring_outermost, x_pixel_from_xpos, y_pixel_from_ypos, TRparam[0].pa[TRparam[0].Nrings-1], TRparam[0].incl[TRparam[0].Nrings-1]);
                if(d1*d2 > 0.0)
                {
                    HI_VF_tr_model[0].data[j][i] = 1E90; // blank
                    HI_VF_einasto_halomodel[0].data[j][i] = 1E90; // blank
                    //HI_VF_res_input_minus_trfit[0].data[j][i] = 1E90; // blank
                    HI_VF_res_input_minus_einastofit[0].data[j][i] = 1E90; // blank
                    HI_VF_res_trfit_minus_einastofit[0].data[j][i] = 1E90; // blank
                }
                else
                {
                    for(n=1; n<3*TRparam[0].Nrings; n++)
                    {
                        if(n >= TRparam[0].Nrings)
                        {
                            TRparam[0].pa[n] = TRparam[0].pa[TRparam[0].Nrings-1];
                            TRparam[0].incl[n] = TRparam[0].incl[TRparam[0].Nrings-1];
                            TRparam[0].vrad[n] = TRparam[0].vrad[TRparam[0].Nrings-1];
                            TRparam[0].vsys[n] = TRparam[0].vsys[TRparam[0].Nrings-1];
                            TRparam[0].vrot[n] = TRparam[0].vrot[TRparam[0].Nrings-1];
                        }
                        ri = TRparam[0].ring_s + (n-0)*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
                        ro = TRparam[0].ring_s + (n-0)*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
                        if(ri < 0) ri = 0;
                        d2 = radius_galaxy_plane((ri+ro)/2.0, x_pixel_from_xpos, y_pixel_from_ypos, TRparam[0].pa[n], TRparam[0].incl[n]);
                        if(d1*d2 > 0.0)
                        {
                            d1 = d2;
                        }
                        else
                        {
                            break;
                        }
                    }

                    ri = TRparam[0].ring_s + (n-1)*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
                    ro = TRparam[0].ring_s + (n-1)*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
                    if(ri < 0) ri = 0;
                    ring00 = (ri+ro)/2.0;

                    ri = TRparam[0].ring_s + (n-0)*TRparam[0].ring_w - 0.5*TRparam[0].ring_w;
                    ro = TRparam[0].ring_s + (n-0)*TRparam[0].ring_w + 0.5*TRparam[0].ring_w;
                    if(ri < 0) ri = 0;
                    ring01 = (ri+ro)/2.0;

                    _n = TRparam[0]._n;
                    r_2 = TRparam[0].r_2;
                    rho_2 = TRparam[0].rho_2;
                    ring00_pc = ring00 * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265;
                    _x = ring00_pc/r_2;
                    iGamma =  gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);
                    vrot_Einasto_fit00 = sqrt((4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/ring00_pc)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma));

                    ring01_pc = ring01 * TRparam[0].pixelScale * 1 * pow(10, 3) / 206.265;
                    _x = ring01_pc/r_2;
                    iGamma =  gsl_sf_gamma(3.0*_n) - gsl_sf_gamma_inc(3.0*_n, _x);
                    vrot_Einasto_fit01 = sqrt((4.0*M_PI*GG*_n*rho_2*pow(r_2, 3)/ring01_pc)*(pow(e, 2.0*_n)*pow(2*_n, -3.0*_n) *iGamma));

                    rads01 = ( ring01 * d1 - ring00 * d2 ) / ( d1 - d2 );
                    w1 = ( ring00 - rads01 ) / ( ring00 - ring01 );
                    w2 = 1.0 - w1;

                    if(n < 100*TRparam[0].Nrings)
                    {
                        posa01 = w1 * TRparam[0].pa[n] + w2 * TRparam[0].pa[n-1];
                        incl01 = w1 * TRparam[0].incl[n] + w2 * TRparam[0].incl[n-1];
                        vrot01 = w1 * TRparam[0].vrot[n] + w2 * TRparam[0].vrot[n-1];
                        vrad01 = w1 * TRparam[0].vrad[n] + w2 * TRparam[0].vrad[n-1];
                        vsys01 = w1 * TRparam[0].vsys[n] + w2 * TRparam[0].vsys[n-1];
                    }
                    /*
                    else
                    {
                        posa01 = w1 * TRparam[0].pa[TRparam[0].Nrings-1] + w2 * TRparam[0].pa[TRparam[0].Nrings-1];
                        incl01 = w1 * TRparam[0].incl[TRparam[0].Nrings-1] + w2 * TRparam[0].incl[TRparam[0].Nrings-1];
                        vrot01 = w1 * TRparam[0].vrot[TRparam[0].Nrings-1] + w2 * TRparam[0].vrot[TRparam[0].Nrings-1];
                        vrad01 = w1 * TRparam[0].vrad[TRparam[0].Nrings-1] + w2 * TRparam[0].vrad[TRparam[0].Nrings-1];
                        vsys01 = w1 * TRparam[0].vsys[TRparam[0].Nrings-1] + w2 * TRparam[0].vsys[TRparam[0].Nrings-1];
                    }
                    */
                    vrot01_einasto = w1 * vrot_Einasto_fit01 + w2 * vrot_Einasto_fit00;

                    if ( rads01 > FLT_EPSILON )
                    {
                       cost = ( - x_pixel_from_xpos * sind( posa01 ) + y_pixel_from_ypos * cosd( posa01 ) ) / rads01;
                       sint = ( - x_pixel_from_xpos * cosd( posa01 ) - y_pixel_from_ypos * sind( posa01 ) ) / rads01 / cosd( incl01 );
                    }
                    else
                    {
                       cost = 1.0;
                       sint = 0.0;
                    }

                    //if(!isinf(HI_VF_boxFiltered_decim0[0].data[j][i]) && !isnan(HI_VF_boxFiltered_decim0[0].data[j][i])) // no blank
                    if(!isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i])) // no blank
                    {
                        HI_VF_tr_model[0].data[j][i] = vsys01 + ( vrot01 * cost + vrad01 * sint ) * sind( incl01 );
                        HI_VF_einasto_halomodel[0].data[j][i] = vsys01 + ( vrot01_einasto * cost + vrad01 * sint ) * sind( incl01 );

                        HI_VF_res_input_minus_trfit[0].data[j][i] = HI_VF[0].data[j][i] - HI_VF_tr_model[0].data[j][i];
                        HI_VF_res_input_minus_einastofit[0].data[j][i] = HI_VF[0].data[j][i] - HI_VF_einasto_halomodel[0].data[j][i];
                        HI_VF_res_trfit_minus_einastofit[0].data[j][i] = HI_VF_tr_model[0].data[j][i] - HI_VF_einasto_halomodel[0].data[j][i];
                        n_res_input_TR++;
                    }
                    else
                    {
                        HI_VF_tr_model[0].data[j][i] = 1E90;
                        HI_VF_einasto_halomodel[0].data[j][i] = 1E90;

                        //HI_VF_res_input_minus_trfit[0].data[j][i] = 1E90;
                        HI_VF_res_input_minus_einastofit[0].data[j][i] = 1E90;
                        HI_VF_res_trfit_minus_einastofit[0].data[j][i] = 1E90;
                    }
                }
            }
        }
    }


    double *res_input_TR = (double*)malloc(sizeof(double) * n_res_input_TR);
    double *res_e_input_TR = (double*)malloc(sizeof(double) * n_res_input_TR);
    k=0;
    for(i=0; i<TRparam[0].nax1; i++)
    {
        for(j=0; j<TRparam[0].nax2; j++)
        {
            //if(!isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i])) // no blank
            if(!isinf(HI_VF_res_input_minus_trfit[0].data[j][i]) && !isnan(HI_VF_res_input_minus_trfit[0].data[j][i])) // no blank
            {
                res_input_TR[k] = HI_VF_res_input_minus_trfit[0].data[j][i];
                res_e_input_TR[k] = 1; // default
                k++;
            }
        }
    }
    robust_mean_std(res_input_TR, k, &mean_res_input_TR, &std_res_input_TR);

    //robust_mean_std_histogram_ac(res_input_TR, res_e_input_TR, k, &mean_res_input_TR, &std_res_input_TR);

    // filtering outlying pixels based on the TR model
    if(TRparam[0]._nfilter > 1)
    {
        for(i=0; i<TRparam[0].nax1; i++)
        {
            for(j=0; j<TRparam[0].nax2; j++)
            {
                if(!isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i])) // no blank
                {
                    if(fabs(HI_VF_res_input_minus_trfit[0].data[j][i]) > mean_res_input_TR + 3*std_res_input_TR)
                    {
                        HI_VF[0].data[j][i] = 1E90;
                    }
                }
            }
        }
    }

    // save the model VF in fits
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_trfit_model, nax1, nax2, HI_VF_tr_model);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_einasto_modelvf, nax1, nax2, HI_VF_einasto_halomodel);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_res_input_minus_trfit, nax1, nax2, HI_VF_res_input_minus_trfit);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_res_input_minus_einastofit, nax1, nax2, HI_VF_res_input_minus_einastofit);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_res_trfit_minus_einastofit, nax1, nax2, HI_VF_res_trfit_minus_einastofit);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_trfit_w, nax1, nax2, HI_VF_weight_TRfit);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_einasto_halofit_geo_radial_w, nax1, nax2, HI_VF_geo_radial_angle_w);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_fract_Navail_Nall, nax1, nax2, HI_VF_fract_navail_nall);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_2Dinput_VF_boxfiltered_sigma_ew, nax1, nax2, HI_VF_boxFiltered_vlos_ew);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_2Dinput_VF_boxfiltered_decim0, nax1, nax2, HI_VF_boxFiltered_decim0);
    save_2dmapsfits(fptr1, inputfile, fname[0].fitsfile_2Dinput_VF_boxfiltered_decim_user, nax1, nax2, HI_VF_boxFiltered_decim_user);

    free(res_input_TR);
    free(res_e_input_TR);
    return 0;
}


int read_2dmaps(TR_ringParameters *TRparam, filename_2dbat *fname, char *card)
{
    int nkeys;
    int status=0;
    int anynul=0;
    int i=0, j=0;

    int vf_e_const_flag=0;
    int mom2_const_flag=0;

    FILE *file_exist;
    fitsfile *fptr1, *fptr2, *fptr3, *fptr4;

    // --------------------------------------------------------------------------------------------- //
    // +++ C. READ VELOCITY FIELDS +++
    // --------------------------------------------------------------------------------------------- //
    // Check if the iput FITS exists...
    if((file_exist = fopen(fname[0].fitsfile_2Dinput_VF, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_2Dinput_VF, TRparam[0].wdir);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        return 0;
    }
    if((file_exist = fopen(fname[0].fitsfile_2Dinput_VF_e, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_2Dinput_VF_e, TRparam[0].wdir);
        printf("\nInstead, the user supplied constnat VF_error (km/s) will be used: %.2f\n", TRparam[0].vf_e_user);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        // copy the input VF name for initializing the VF_e's 2D array below
        vf_e_const_flag = 1; // constant vf_e used
        sprintf(fname[0].fitsfile_2Dinput_VF_e, "%s", fname[0].fitsfile_2Dinput_VF);
    }
    if((file_exist = fopen(fname[0].fitsfile_mom0, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_mom0, TRparam[0].wdir);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        return 0;
    }
    if((file_exist = fopen(fname[0].fitsfile_mom2, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_mom2, TRparam[0].wdir);
        printf("\nInstead, the user supplied constnat Vdisp (km/s) will be used: %.2f\n", TRparam[0].vdisp_user);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        // copy the input VF name for initializing the mom2's 2D array below
        mom2_const_flag = 1; // constant mom2 used
        sprintf(fname[0].fitsfile_mom2, "%s", fname[0].fitsfile_2Dinput_VF);
    }


    int fits_read_key_status=0;
    ffopen(&fptr1, fname[0].fitsfile_2Dinput_VF, READONLY, &status); // open fits
    fits_get_hdrspace(fptr1, &nkeys, NULL, &status);
    // READ NAXIS1
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS1", &TRparam[0].nax1, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS2", &TRparam[0].nax2, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CDELT1", &TRparam[0].pixelScale, card, &status);
    if(fits_read_key_status == 202) // if CDELT1 keyword is not in the header, check if CD1_1 is
    {
        status = 0;
        fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CD1_1", &TRparam[0].pixelScale, NULL, &status);
        if(fits_read_key_status == 202) // if CD1_1 keyword is not in the header, put a null value
        {
            TRparam[0].pixelScale = 0.999;
        }
    }
    TRparam[0].pixelScale = fabs(3600*TRparam[0].pixelScale); // in arcsec

    if(status)
    {
        printf("CHECK FITS HEADER: CDELT1 or CD1_1 KEYWORDS DOESN'T EXIST...\n");
        fits_report_error(stderr, status);
        return(status);
    }

    // C-2. Read input velocity field
    ffg2de(fptr1, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr1, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_boxFiltered[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr1, &status); // close fits

    // C-3. Read input sigma field (or mom2)
    ffopen(&fptr2, fname[0].fitsfile_2Dinput_VF_e, READONLY, &status); // open fits
    ffg2de(fptr2, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_sigma[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr2, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_boxFiltered_sigma[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr2, &status); // close fits

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // if constant vf_e is used (e.g., 1 channel resolution in km/s)
    if(vf_e_const_flag == 1)
    {
        for(i=0; i<TRparam[0].nax1; i++)
        {
            for(j=0; j<TRparam[0].nax2; j++)
            {
                HI_VF_sigma[0].data[j][i] = TRparam[0].vf_e_user; // constant vf_e in km/s
            }   
        }
    }

    // C-4. Read input mom0 field
    ffopen(&fptr3, fname[0].fitsfile_mom0, READONLY, &status); // open fits
    ffg2de(fptr3, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_mom0[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr3, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_sn[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr3, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_boxFiltered_SN[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr3, &status); // close fits

    // C-5. Read input mom2 field
    ffopen(&fptr4, fname[0].fitsfile_mom2, READONLY, &status); // open fits
    ffg2de(fptr4, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_mom2[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr4, &status); // close fits

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // if constant vdisp is used (e.g., xxx in km/s)
    if(mom2_const_flag == 1)
    {
        for(i=0; i<TRparam[0].nax1; i++)
        {
            for(j=0; j<TRparam[0].nax2; j++)
            {
                HI_VF_mom2[0].data[j][i] = TRparam[0].vdisp_user; // constant vdisp in km/s
            }   
        }
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // put a null value to the blank pixels
    for(i=0; i<TRparam[0].nax1; i++)
    {
        for(j=0; j<TRparam[0].nax2; j++)
        {
            // Zero value is often used as a null value in some fits images like SAMI
            if(isnan(HI_VF[0].data[j][i]) || isinf(HI_VF[0].data[j][i]) || fabs(HI_VF[0].data[j][i])<1E-10)
            {
                HI_VF[0].data[j][i] = 1E90;
                HI_VF_boxFiltered[0].data[j][i] = 1E90;
                HI_VF_boxFiltered_SN[0].data[j][i] = 1E90;
                HI_VF_sn[0].data[j][i] = 1E90;
                HI_VF_mom0[0].data[j][i] = 1E90;
                HI_VF_mom2[0].data[j][i] = 1E90;
                HI_VF_sigma[0].data[j][i] = 1E90;
                HI_VF_boxFiltered_sigma[0].data[j][i] = 1E90;
            }
        }
    }

    return 0;
}

// --- End of line

