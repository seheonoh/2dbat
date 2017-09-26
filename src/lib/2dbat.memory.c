#include "2dbat.memory.h"

// 2DBAT user defined functions
// memory related

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void free_2dmap_arrays(TR_ringParameters *TRparam, multinest_paramters *multinest_param, float *fits_pointer)
{
    free(HI_VF[0].data);
    free(HI_VF_mom0[0].data);
    free(HI_VF_mom2[0].data);
    free(HI_VF_median[0].data);
    free(HI_VF_res[0].data);
    free(HI_VF_sn[0].data);
    free(HI_VF_sn_median[0].data);
    free(HI_VF_sn_res[0].data);
    free(HI_VF_sigma[0].data);
    free(HI_VF_boxFiltered_vlos_ew[0].data);
    free(HI_VF_sigma_res[0].data);
    free(HI_VF_boxFiltered[0].data);
    free(HI_VF_boxFiltered_sigma[0].data);
    free(HI_VF_boxFiltered_sigma_e_norm[0].data);
    free(HI_VF_sigma_geo_flux_weighted[0].data);
    free(HI_VF_boxFiltered_SN[0].data);
    free(HI_VF_boxFiltered_decim0[0].data);
    free(HI_VF_boxFiltered_decim_user[0].data);
    free(HI_VF_geo_radial_angle_w[0].data);
    free(HI_VF_weight_TRfit[0].data);
    free(HI_VF_fract_navail_nall[0].data);
    free(HI_VF_tr_model[0].data);
    free(HI_VF_einasto_halomodel[0].data);
    free(HI_VF_res_input_minus_trfit[0].data);
    free(HI_VF_res_input_minus_einastofit[0].data);
    free(HI_VF_res_trfit_minus_einastofit[0].data);
    free(HI_VF_temp[0].data);
    free(TRparam);
    free(multinest_param);
    free(fits_pointer);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void malloc_2dmap_arrays(int nax1, int nax2, float *fits_pointer)
{
    int i=0;
    // --------------------------------------------------------------------------------------------- //
    // +++ B. INITIALISATION +++
    // --------------------------------------------------------------------------------------------- //
    // --- B-1. Initialise parameters for velocity field ---    
    // 1. Velocity field area

    // --- A-4. Arrays for 2D velocity fields ---   
    // Dynamic allocation
    HI_VF[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF[0].data[i] = fits_pointer + i * nax2;

    HI_VF_mom0[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_mom0[0].data[i] = fits_pointer + i * nax2;

    HI_VF_mom2[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_mom2[0].data[i] = fits_pointer + i * nax2;

    HI_VF_median[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_median[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sigma[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sigma[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_vlos_ew[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_vlos_ew[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sigma_res[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sigma_res[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sn[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sn[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sn_median[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sn_median[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sn_res[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sn_res[0].data[i] = fits_pointer + i * nax2;

    HI_VF_geo_radial_angle_w[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_geo_radial_angle_w[0].data[i] = fits_pointer + i * nax2;

    HI_VF_weight_TRfit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_weight_TRfit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_fract_navail_nall[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_fract_navail_nall[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_sigma[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_sigma[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_sigma_e_norm[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_sigma_e_norm[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_SN[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_SN[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_decim0[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_decim0[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_decim_user[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_decim_user[0].data[i] = fits_pointer + i * nax2;

    HI_VF_tr_model[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_tr_model[0].data[i] = fits_pointer + i * nax2;

    HI_VF_einasto_halomodel[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_einasto_halomodel[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res_input_minus_trfit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res_input_minus_trfit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res_input_minus_einastofit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res_input_minus_einastofit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res_trfit_minus_einastofit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res_trfit_minus_einastofit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sigma_geo_flux_weighted[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sigma_geo_flux_weighted[0].data[i] = fits_pointer + i * nax2;

    HI_VF_temp[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_temp[0].data[i] = fits_pointer + i * nax2;
}



// --- End of line



