
#include "2dbat.main.h"

// --- Main program --- //
int main(int argc, char *argv[])
{
    // ++ : Copy pa/incl user functions
    if(argc != 48)
    {
        usage_2dbat();
        exit(0);
    }

    // --------------------------------------------------------------------------------------------- //
    // +++ DECLARE PARAMETRES +++
    // --------------------------------------------------------------------------------------------- //

    // File names for FITS  
    int nkeys;
    char card[FLEN_CARD];

    // Dimensions   
    int i=0, j=0, i0=0, j0=0, k=0;
    int status=0;
    int anynul=0;
    int rank=0;
    int loop_iteration=0;

    FILE *file_exist, *ftemp;
    fitsfile *fptr1, *fptr2;
    float *fits_pointer;
    clock_t start_time = clock();

    // multinest parameters temp
    // Einasto rotation curve
    double _n_t=9999, _r_2_t=9999, _rho_2_t=9999;
    // fraction
    double Nall_ring, Navail_ring;

    // ETC params
    int ein_dirty=0;
    FILE *fp;
    int conv_i=0;
    int en = 0;
    int nf=0;

    // MPI variables
    int flag=0;
    int n_node;

    // TR param struct
    TRparam = (TR_ringParameters *)malloc(sizeof (TR_ringParameters) * 1);              // not required by MultiNest, any additional information user wants to pass
    // multinest param struct
    multinest_param = (multinest_paramters *)malloc(sizeof (multinest_paramters) * 1);
    // filename param struct
    fname = (filename_2dbat *)malloc(sizeof (filename_2dbat) * 1);

    // read user input
    //read_user_input_init_params(TRparam, multinest_param, argc, argv, fitsfile_2Dinput_VF, fitsfile_2Dinput_VF_e, fitsfile_mom0, fitsfile_mom2);
    read_user_input_init_params(TRparam, multinest_param, fname, argc, argv);

    // --------------------------------------------------------------------------------------------- //
    // +++ A. READ INPUT VF KEYWORDS ++W.+
    // --------------------------------------------------------------------------------------------- //
    int fits_read_key_status=0;
    ffopen(&fptr1, fname[0].fitsfile_2Dinput_VF, READONLY, &status); // open fits
    fits_get_hdrspace(fptr1, &nkeys, NULL, &status);
    // READ NAXIS1
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS1", &TRparam[0].nax1, NULL, &status);
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS2", &TRparam[0].nax2, NULL, &status);
    fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CDELT1", &TRparam[0].pixelScale, NULL, &status);

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

    // declare 2D arrays for 2d maps
    malloc_2dmap_arrays(TRparam[0].nax1, TRparam[0].nax2, fits_pointer);

    // read 2D maps
    read_2dmaps(TRparam, fname, card);


    // --- B-2. Initialise parameters
    TRparam[0].ellipse_semi_mx_boxfiltered = TRparam[0].nax1;
    TRparam[0].sigma_factor = 1.0;
    // pa
    TRparam[0].pa1 = 0; // 
    TRparam[0].pa2 = 1; // 
    TRparam[0].PA_MAX_in_degree = 360.; // maximum value for PA
    // incl
    TRparam[0].incl1 = 0; // 
    TRparam[0].incl2 = 1; // 
    TRparam[0].INCL_MAX_in_degree = 90.; // maximum value for INCL


    MPI_Init(&argc, &argv); /* initiate MPI ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    MPI_Initialized(&flag);
    MPI_Comm_size(MPI_COMM_WORLD, &n_node);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //to get MPI ranks

    MPI_Request send_req[4096], recv_req[4096]; // maximum number of cores : 4096
    MPI_Status mpistatus[4096];

    MPI_Datatype TRparam_mpi;
    TRparam_mpi = allocate_mpi_dataset(TRparam, type, blocklen);    
    MPI_Type_commit(&TRparam_mpi);

    // ------------------------------------------------------------ //
    // ------------------------------------------------------------ //
    // SINGLE CORE MODE
    // ------------------------------------------------------------ //
    // ------------------------------------------------------------ //
    if(n_node == 1) // single node
    {
        // ++ : Copy pa/incl user functions
        if(argc != 48)
        {
            usage_2dbat();
            exit(0);
        }
        strcpy(fname[0].pa_function_user, TRparam[0].pa_function); // "bspline"
        strcpy(fname[0].incl_function_user, TRparam[0].incl_function); // "bspline"
        strcpy(fname[0].vrad_function_user, TRparam[0].vrad_function); // "bspline"

        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // ++ : nfilter: # of iteration to filter out any outlying pixels
        // ++ : set to none (nfilter=1) (optional) 
        nf = 0;
        while(nf < TRparam[0]._nfilter)
        {
            // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
            // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
            // -- PART A
            // A-1. Find the largest connected area
            // A-2. Perform ellipse fit
            nf++;
            // ******************************************************
            // --- a-1. Find the largest connected area to fit based on connected-area labelling algorithm ---
            printf("\n!+++ A-%d-1. FIND THE LARGEST CONNECTED AREA: ", nf);
            TRparam[0].use_allPixels = 0;
            find_the_largest_connected_area(TRparam, 0, 0, TRparam[0].nax1, TRparam[0].nax2, 1, TRparam[0].use_allPixels);
            printf("[Done] \n\n");

            // ******************************************************
            // --- a-2. Define the found connected area for ellipse fit: both sides ---
            printf("!+++ A-%d-2. DEFINE THE AREA FOR ELLIPSE FIT: ", nf);
            define_area_tofit(0, TRparam[0].nax1, 0, 0, multinest_param, TRparam, 0);
            printf("[Done] \n\n");


            // ******************************************************
            // --- a-3. Perform ellipse fit and update uniform priors of ring parametres ---
            printf("!+++ A-%d-3. PERFORM ELLIPSE FIT AND UPDATE PRIORS OF RING PARAMETERS: ", nf);
            do_ellipseFit_update_uniformPriors_ringParam(TRparam);
            define_area_tofit(0, TRparam[0].ellipse_semi_mx_boxfiltered*5.0, 0, 0, multinest_param, TRparam, 0);
            TRparam[0].rGalaxyPlane_pixel_max = TRparam[0].ellipse_semi_mx_boxfiltered*1.2; 
            printf("[Done] \n\n");

            // update Nrings based on filling factor
            set_Nrings(TRparam, &Nall_ring, &Navail_ring);

            // VRAD is fitted?
            if(TRparam[0].vrad_nbreak_bspline == 1 && TRparam[0].vrad_order_bspline == 0) // not fitted
                TRparam[0].vrad_fix = 'F';
            else
                TRparam[0].vrad_fix = 'T';

            printf("!+++ A-%d-4. TILTED-RING FIT WITH ELLIPSE FIT RINGS: ", nf);
            printf("[Done]\n");


            // ******************************************************
            // A-3. Perform TR fits with (XPOS, YPOS, VSYS, PA, INCL, VROT, sigmafactor) free + ellipse rings defined by the previous ellipse fit
            trfit_multinest_ellipsefit_rings_student("xpos", 'T', "ypos", 'T', "vsys", 'T', "pa", 'T', "incl", 'T', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, 0, "finalfit", 'N'); // both sides

            // NULL VALUE for n, r_2, and rho_2 : TR fit results with all ring params free
            // write the TR fit results with all ring parameters free: this is for a reference check
            TRparam[0]._n = 0.1;
            TRparam[0].r_2 = TRparam[0].rGalaxyPlane_pixel_max;
            TRparam[0].rho_2 = 0.5;
            write_trfit_results(TRparam, fname[0].bayesFit_outputfile_allfree);

            // set PA/INCL/VRAD constant
            strcpy(TRparam[0].pa_function, "const");
            strcpy(TRparam[0].incl_function, "const");
            strcpy(TRparam[0].vrad_function, "const");
            reset_HI_VF_weight_trfit(TRparam);
            print_trfit_results(TRparam);
            printf("!+++ A-%d-5. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
            printf("[Done]\n");

            // ******************************************************
            // A-4. Perform TR fitting with (PA, INCL, VROT, sigmafactor) free
            trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, "finalfit", 'N');

            if(nf != TRparam[0]._nfilter) // no filter out
            {
                // make a model VF to filter out any outliers : currently
                make_model_vfs(TRparam, fptr1, fname[0].fitsfile_2Dinput_VF, TRparam[0].nax1, TRparam[0].nax2, fname); 
            }

            // ******************************************************
            // A-5. Restore the PA/INCL functions supplied by the user
            strcpy(TRparam[0].pa_function, fname[0].pa_function_user);
            strcpy(TRparam[0].incl_function, fname[0].incl_function_user);
            strcpy(TRparam[0].vrad_function, "const");
            reset_HI_VF_weight_trfit(TRparam);
            print_trfit_results(TRparam);

            // ******************************************************
            // A-6. Filter outliers of PA & INCL out
            // A-7. Perform B-spline fit and update priors of PA+INCL B-spline priors
            if(strcmp(TRparam[0].pa_function, "bspline") == 0)
            {
                bsplinefit_set_unipriors("PA", TRparam);
            }

            if(strcmp(TRparam[0].incl_function, "bspline") == 0)
            {
                bsplinefit_set_unipriors("INCL", TRparam);
            }
            printf("!+++ A-%d-6. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
            printf("[Done]\n");

            // ******************************************************
            // A-8. Perform TR fit with (XPOS, YPOS, VSYS, VROT, sigmafactor) free
            trfit_multinest_trfit_rings_student("xpos", 'T', "ypos", 'T', "vsys", 'T', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, "finalfit", 'N');
            reset_HI_VF_weight_trfit(TRparam);
            print_trfit_results(TRparam);

            // ******************************************************
            // A-9. Filter outliers of (XPOS, YPOS, VSYS) out & update the priors of (XPOS, YPOS, VSYS)
            // Filter outliers out of XPOS, YPOS & VSYS, and update priors
            bsplinefit_set_unipriors("XPOS", TRparam);
            bsplinefit_set_unipriors("YPOS", TRparam);
            bsplinefit_set_unipriors("VSYS", TRparam);
            printf("!+++ A-%d-7. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
            printf("[Done]\n");

            // ******************************************************
            // A-10. Perform TR fit with (VROT, sigmafactor) free
            trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 999, "finalfit", 'N');

            print_trfit_results(TRparam);
            reset_HI_VF_weight_trfit(TRparam);
            // reset HI_VF_fract_navail_nall before the final TR fit :
            // HI_VF_fract_navail_nall from the last TR fit is used for Einasto fit
            reset_HI_VF_fract_navail_nall(TRparam);
            // Filter out any outliers, derive initial estimates of Rc & rho0, and setup their priors
            if(TRparam[0].vrad_fix == 'T') // if VRAD is fitted
            {
                strcpy(TRparam[0].vrad_function, fname[0].vrad_function_user);
                bsplinefit_set_unipriors("VRAD", TRparam); // change the order : liner, quadrature or cubic 
            }
            // set vrad_fix = 'F' for VROT fitting only. This is for deriving robust initial values of rho and rc
            TRparam[0].vrad_fix = 'F';
            printf("!+++ A-%d-8. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
            printf("[Done]\n");
        
            // ******************************************************
            // A-10-1. Perform TR fit with (VROT) free + VRAD fixed: if VRAD is intially fitted
            if(TRparam[0].vrad_nbreak_bspline != 1 || TRparam[0].vrad_order_bspline != 0) // if fitted: re-run TRfit with the fitted value of VRAD from bsplinefit_set_unipriors() to derive more stable(?) VROT. This is for deriving more robust n, r2, and rho2 below.
            {
                trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 999, "finalfit", 'N');
            }

            // ******************************************************
            // A-11. Filter outliers in VROT out, and update priors
            print_trfit_results(TRparam);
            reset_HI_VF_weight_trfit(TRparam);

            // restore TRparam[0].vrad_fix
            if(TRparam[0].vrad_nbreak_bspline == 1 && TRparam[0].vrad_order_bspline == 0) // not fitted
                TRparam[0].vrad_fix = 'F';
            else
                TRparam[0].vrad_fix = 'T';

            // Filter outliers out of VROT, and update priors
            // Filter out any outliers, derive initial estimates of Rc & rho0, and setup their priors
            bsplinefit_set_unipriors("VROT-einasto_halofit_rings", TRparam); // change the order : liner, quadrature or cubic 

            // initialise einasto halo params and update vrot_intp
            init_einastohalo_vrot_intp(TRparam);

            // make model vfs
            make_model_vfs(TRparam, fptr1, fname[0].fitsfile_2Dinput_VF, TRparam[0].nax1, TRparam[0].nax2, fname); 
        }

        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // -- PART B
        // B-1. Perform vEinasto 1D fit to update the priors of (n, r2, rho2)
        TRparam[0].loop_check = 0;
        int _tag_t=999;
        TRparam[0].n_hist_post = 0;
        while(TRparam[0].loop_check == 0)
        {
            // ******************************************************
            // B-2. First Einasto halo rotation curve 1D fit : dirty fit with the given ranges of priors
            v_einasto_1d_multinest4p(multinest_param, TRparam);

            // Update Einasto halo parameter priors
            conv_i = 0;
            TRparam[0].n_hist_post = 999;
            v_einasto_1d_multinest4p_remove_outliers(TRparam);
            //update_Einasto_params_priors(TRparam, 2);
            TRparam[0].N_reliable_rings = TRparam[0].Nrings;
            update_vrot_prior(TRparam);

            if(fabs((_n_t-TRparam[0]._n)/TRparam[0]._n) < 0.2) conv_i++;
            if(fabs((_r_2_t-TRparam[0].r_2)/TRparam[0].r_2) < 0.2) conv_i++;
            if(fabs((_rho_2_t-TRparam[0].rho_2)/TRparam[0].rho_2) < 0.2) conv_i++;

            if(conv_i == 3)
            {
                TRparam[0].loop_check = 1; // exit loop
                TRparam[0]._n_t = TRparam[0]._n;
                TRparam[0]._r_2_t = TRparam[0].r_2;
                TRparam[0]._rho_2_t = TRparam[0].rho_2;
                TRparam[0].e_sigma_fitted_t = TRparam[0].e_sigma_fitted;
            }
            else
            {
                TRparam[0].loop_check = 0; // continue loop
                // update temp einasto params
                _n_t = TRparam[0]._n;
                _r_2_t = TRparam[0].r_2;
                _rho_2_t = TRparam[0].rho_2;

                TRparam[0]._n_t = TRparam[0]._n;
                TRparam[0]._r_2_t = TRparam[0].r_2;
                TRparam[0]._rho_2_t = TRparam[0].rho_2;
                TRparam[0].e_sigma_fitted_t = TRparam[0].e_sigma_fitted;
            }

            if(loop_iteration > 50)
            {
                printf("Einasto 1D fit was not successful! Please check...\n");
                TRparam[0].loop_check = 1; // exit loop
            }
        }
        
        // ******************************************************
        // B-3. Reset with the multinest params for quick dirty Einasto fits
        update_multinest_params_dirtyfit(multinest_param);

        // ******************************************************
        // B-4. Reset the parameters for the dirty fit
        define_area_tofit_set_geo_weight(0, TRparam[0].ellipse_semi_mx_boxfiltered*5.0, TRparam[0].decimX_einasto_halofit_d, TRparam[0].decimY_einasto_halofit_d, multinest_param, TRparam, 0, 0);

        // update fract_Navail_Nall + set free params to fit
        //set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "n", 'T', "r_2", 'T', "rho_2", 'T', "fullfit", 'N');
        set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'T', "ypos", 'T', "vsys", 'T', "pa", 'T', "incl", 'T', "n", 'T', "r_2", 'T', "rho_2", 'T', "fullfit", 'N');

        // ******************************************************
        // B-5. DIRTY FIT : initial 2D einasto fit with weakly informative priors from the above 1D einasto fit
        einasto_halofit_multinest_student("xpos", TRparam[0].xpos_fix, "ypos", TRparam[0].ypos_fix, "vsys", TRparam[0].vsys_fix, "pa", TRparam[0].pa_fix, "incl", TRparam[0].incl_fix, "_n_Einasto", TRparam[0]._n_fix, "r_2_Einasto", TRparam[0].r_2_fix, "rho_2_Einasto", TRparam[0].rho_2_fix, "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, rank, "einastoFit.dirty.");

        // reset the parameters for the final fit
        printf("!+++ B-%d-1. DIRTY 2D EINASTO FIT : GAUSSIAN+UNIFORM PRIORS BASED ON the 1D EINASTO FIT +++\n", nf);
        print_priors_info("DIRTY", TRparam, n_node, ein_dirty+1);

        // Write down the MAP vaules of the parameters into the first line : this is for plotting later
        sprintf(fname[0].txtfile_multinest_output, "%s", "einastoFit.dirty..txt");
        sprintf(TRparam[0].txtfile_multinest_output, "%s", "einastoFit.dirty..txt");
        add_MAP_line_multinest_output(fname[0].txtfile_multinest_output, TRparam);

        // define are to fit & set geometrical weight
        define_area_tofit_set_geo_weight(0, TRparam[0].ellipse_semi_mx_boxfiltered*5.0, TRparam[0].decimX_einasto_halofit, TRparam[0].decimY_einasto_halofit, multinest_param, TRparam, 0, 0);
        // update fract_Navail_Nall + set free params to fit
        //set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "n", 'T', "r_2", 'T', "rho_2", 'T', "fullfit", 'Y');
        set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "n", 'T', "r_2", 'T', "rho_2", 'T', "fullfit", 'F');
        set_nfree_params_trfit_multinest_trfit_rings_student(TRparam);  

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-6. Full fit : final einasto fit with a larger Nlive & multinset paramers supplied by the user + updated priors based on the dirty fit
        update_multinest_params_fullfit(multinest_param);

        einasto_halofit_multinest_student("xpos", TRparam[0].xpos_fix, "ypos", TRparam[0].ypos_fix, "vsys", TRparam[0].vsys_fix, "pa", TRparam[0].pa_fix, "incl", TRparam[0].incl_fix, "_n_Einasto", TRparam[0]._n_fix, "r_2_Einasto", TRparam[0].r_2_fix, "rho_2_Einasto", TRparam[0].rho_2_fix, "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, rank, "einastoFit.");

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-7. Update the status of TR ring params
        printf("!+++ B-%d-2. FULL 2D EINASTO FIT : GAUSSIAN+UNIFORM PRIORS BASED ON the DIRTY 2D EINASTO FIT +++\n", nf);
        print_priors_info("FULL", TRparam, n_node, ein_dirty+1);

        // update fract_Navail_Nall

        // Write down the MAP vaules of the parameters into the first line : this is for plotting later
        sprintf(fname[0].txtfile_multinest_output, "%s", "einastoFit..txt");
        sprintf(TRparam[0].txtfile_multinest_output, "%s", "einastoFit..txt");
        add_MAP_line_multinest_output(fname[0].txtfile_multinest_output, TRparam);

        // set vrad_function 'const' to make fits ring-by-ring : not a model based...
        strcpy(TRparam[0].vrad_function, "const");

        // set the final TR fit : n_free params : VROT = T
        //set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "n", 'F', "r_2", 'F', "rho_2", 'F', "fullfit", 'Y');
        set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "n", 'F', "r_2", 'F', "rho_2", 'F', "fullfit", 'Y');
        set_nfree_params_trfit_multinest_trfit_rings_student(TRparam);  

        printf("!+++ B-%d-3. UPDATE PRIORS OF TILTED-RING PARAMETRES: ", nf);
        printf("[Done]\n\n");

        // Update VROT for TR fitting based on the Einasto fit
        update_vrot_prior(TRparam);

        // reset HI_VF_TRfit to zero
        reset_HI_VF_weight_trfit(TRparam);
        //TRparam[0].Nrings += 1;

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-8. Perform TR fitting with (VROT) free : RECEDING SIDE
        trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, -1, "finalfit", 'Y'); // receding side

        printf("!+++ B-%d-4. PERFORM THE FINAL TILTED-RING FIT (RECEDING SIDE): ", nf);
        printf("[Done]\n");
        reset_HI_VF_weight_trfit(TRparam);
        print_trfit_results(TRparam);

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-9. Perform TR fitting with (VROT) free : APPROACHING SIDE
        trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 1, "finalfit", 'Y'); // approaching side

        printf("!+++ B-%d-5. PERFORM THE FINAL TILTED-RING FIT (APPROACHING SIDE): ", nf);
        printf("[Done]\n");
        reset_HI_VF_weight_trfit(TRparam);
        print_trfit_results(TRparam);

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-10. Perform TR fitting with (VROT) free : BOTH SIDE
        trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, "finalfit", 'Y'); // both sides
        printf("!+++ B-%d-6. PERFORM THE FINAL TILTED-RING FIT (BOTH SIDES) ", nf);
        printf("[Done]\n");
        estimate_vrot_error(multinest_param, TRparam, 0, "finalfit", 'Y'); // both sides
        print_trfit_results(TRparam);

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-11. Save the results and free memories
        printf("!+++ B-%d-7. SAVE THE FIT RESULTS AND MAKE MODEL VELOCITY FIELDS: ", nf);
        // write the TR fit results into a file
        write_trfit_results(TRparam, fname[0].bayesFit_outputfile);
        // make model 2dmaps
        make_model_vfs(TRparam, fptr1, fname[0].fitsfile_2Dinput_VF, TRparam[0].nax1, TRparam[0].nax2, fname); 
        printf("[Done]\n\n");

        printf("\t++ CHECK THE FIT RESULTS IN THE FOLLOWING DIRECTORY: ++\n\n");
        printf("\t%s/2dbat_output\n\n", TRparam[0].wdir);

        clock_t stop_time = clock();
        double elapsed = (double)(stop_time - start_time) / CLOCKS_PER_SEC;
        //printf("sigma_factor: %f\n", TRparam[0].sigma_factor);
        printf("THE EXECUTION TIME ELAPSED: %f SECS\n\n", elapsed);

        // free dynamically allocated memories
        free_2dmap_arrays(TRparam, multinest_param, fits_pointer);
        return 0;
    }
    // ------------------------------------------------------------ //
    // ------------------------------------------------------------ //
    // MULTI CORE MODE
    // ------------------------------------------------------------ //
    // ------------------------------------------------------------ //
    else // multi-node in MPI
    {
        // ++ : Copy pa/incl user functions
        if(rank == 0)
        {
            if(argc != 48)
            {
                usage_2dbat();
                exit(0);
            }
            strcpy(fname[0].pa_function_user, TRparam[0].pa_function); // "bspline"
            strcpy(fname[0].incl_function_user, TRparam[0].incl_function); // "bspline"
            strcpy(fname[0].vrad_function_user, TRparam[0].vrad_function); // "bspline"
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // ++ : nfilter: # of iteration to filter out any outlying pixels
        // ++ : set to none (nfilter=1) (optional) 
        nf = 0;
        while(nf < TRparam[0]._nfilter)
        {
            // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
            // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
            // -- PART A
            // A-1. Find the largest connected area
            // A-2. Perform ellipse fit
            nf++;
            if(rank == 0)
            {
                // ******************************************************
                // --- a-1. Find the largest connected area to fit based on connected-area labelling algorithm ---
                printf("\n!+++ A-%d-1. FIND THE LARGEST CONNECTED AREA: ", nf);
                TRparam[0].use_allPixels = 0;
                find_the_largest_connected_area(TRparam, 0, 0, TRparam[0].nax1, TRparam[0].nax2, 1, TRparam[0].use_allPixels);
                printf("[Done] \n\n");

                // ******************************************************
                // --- a-2. Define the found connected area for ellipse fit: both sides ---
                printf("!+++ A-%d-2. DEFINE THE AREA FOR ELLIPSE FIT: ", nf);
                define_area_tofit(0, TRparam[0].nax1, 0, 0, multinest_param, TRparam, 0);
                printf("[Done] \n\n");

                // ******************************************************
                // --- a-3. Perform ellipse fit and update uniform priors of ring parametres ---
                printf("!+++ A-%d-3. PERFORM ELLIPSE FIT AND UPDATE PRIORS OF RING PARAMETERS: ", nf);
                do_ellipseFit_update_uniformPriors_ringParam(TRparam);
                define_area_tofit(0, TRparam[0].ellipse_semi_mx_boxfiltered*5.0, 0, 0, multinest_param, TRparam, 0);
                TRparam[0].rGalaxyPlane_pixel_max = TRparam[0].ellipse_semi_mx_boxfiltered*1.2; 
                printf("[Done] \n\n");

                // update Nrings based on filling factor
                set_Nrings(TRparam, &Nall_ring, &Navail_ring);

                // VRAD is fitted?
                if(TRparam[0].vrad_nbreak_bspline == 1 && TRparam[0].vrad_order_bspline == 0) // not fitted
                    TRparam[0].vrad_fix = 'F';
                else
                    TRparam[0].vrad_fix = 'T';

                printf("!+++ A-%d-4. TILTED-RING FIT WITH ELLIPSE FIT RINGS: ", nf);
                printf("[Done]\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);


            // ******************************************************
            // A-3. Perform TR fits with (XPOS, YPOS, VSYS, PA, INCL, VROT, sigmafactor) free + ellipse rings defined by the previous ellipse fit
            trfit_multinest_ellipsefit_rings_student("xpos", 'T', "ypos", 'T', "vsys", 'T', "pa", 'T', "incl", 'T', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, 0, "finalfit", 'N'); // both sides
            MPI_Barrier(MPI_COMM_WORLD);

            if(rank == 0)
            {
                // NULL VALUE for n, r_2, and rho_2 : TR fit results with all ring params free
                // write the TR fit results with all ring parameters free: this is for a reference check
                TRparam[0]._n = 0.1;
                TRparam[0].r_2 = TRparam[0].rGalaxyPlane_pixel_max;
                TRparam[0].rho_2 = 0.5;
                write_trfit_results(TRparam, fname[0].bayesFit_outputfile_allfree);

                // set PA/INCL/VRAD constant
                strcpy(TRparam[0].pa_function, "const");
                strcpy(TRparam[0].incl_function, "const");
                strcpy(TRparam[0].vrad_function, "const");
                reset_HI_VF_weight_trfit(TRparam);
                print_trfit_results(TRparam);
                printf("!+++ A-%d-5. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
                printf("[Done]\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

            // ******************************************************
            // A-4. Perform TR fitting with (PA, INCL, VROT, sigmafactor) free
            trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, "finalfit", 'N');
            MPI_Barrier(MPI_COMM_WORLD);


            if(rank == 0)
            {
                if(nf != TRparam[0]._nfilter) // no filter out
                {
                    // make a model VF to filter out any outliers : currently
                    //make_model_vfs(TRparam, fptr1, fname[0].fitsfile_2Dinput_VF, TRparam[0].nax1, TRparam[0].nax2, fname); 
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

            // ******************************************************
            // A-5. Restore the PA/INCL functions supplied by the user
            if(rank == 0)
            {
                strcpy(TRparam[0].pa_function, fname[0].pa_function_user);
                strcpy(TRparam[0].incl_function, fname[0].incl_function_user);
                strcpy(TRparam[0].vrad_function, "const");
                reset_HI_VF_weight_trfit(TRparam);
                print_trfit_results(TRparam);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // ******************************************************
            // A-6. Filter outliers of PA & INCL out
            // A-7. Perform B-spline fit and update priors of PA+INCL B-spline priors
            if(rank == 0)
            {
                if(strcmp(TRparam[0].pa_function, "bspline") == 0)
                {
                    bsplinefit_set_unipriors("PA", TRparam);
                }

                if(strcmp(TRparam[0].incl_function, "bspline") == 0)
                {
                    bsplinefit_set_unipriors("INCL", TRparam);
                }
                printf("!+++ A-%d-6. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
                printf("[Done]\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

            // ******************************************************
            // A-8. Perform TR fit with (XPOS, YPOS, VSYS, VROT, sigmafactor) free
            trfit_multinest_trfit_rings_student("xpos", 'T', "ypos", 'T', "vsys", 'T', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, "finalfit", 'N');
            MPI_Barrier(MPI_COMM_WORLD);

            if(rank == 0)
            {
                reset_HI_VF_weight_trfit(TRparam);
                print_trfit_results(TRparam);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // ******************************************************
            // A-9. Filter outliers of (XPOS, YPOS, VSYS) out & update the priors of (XPOS, YPOS, VSYS)
            if(rank == 0)
            {
                // Filter outliers out of XPOS, YPOS & VSYS, and update priors
                bsplinefit_set_unipriors("XPOS", TRparam);
                bsplinefit_set_unipriors("YPOS", TRparam);
                bsplinefit_set_unipriors("VSYS", TRparam);
                printf("!+++ A-%d-7. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
                printf("[Done]\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

            // ******************************************************
            // A-10. Perform TR fit with (VROT, sigmafactor) free
            MPI_Barrier(MPI_COMM_WORLD);
            trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 999, "finalfit", 'N');
            MPI_Barrier(MPI_COMM_WORLD);

            if(rank == 0)
            {
                print_trfit_results(TRparam);
                reset_HI_VF_weight_trfit(TRparam);
                // reset HI_VF_fract_navail_nall before the final TR fit :
                // HI_VF_fract_navail_nall from the last TR fit is used for Einasto fit
                reset_HI_VF_fract_navail_nall(TRparam);
                // Filter out any outliers, derive initial estimates of Rc & rho0, and setup their priors
                if(TRparam[0].vrad_fix == 'T') // if VRAD is fitted
                {
                    strcpy(TRparam[0].vrad_function, fname[0].vrad_function_user);
                    bsplinefit_set_unipriors("VRAD", TRparam); // change the order : liner, quadrature or cubic 
                }
                // set vrad_fix = 'F' for VROT fitting only. This is for deriving robust initial values of rho and rc
                TRparam[0].vrad_fix = 'F';
                printf("!+++ A-%d-8. TILTED-RING FIT WITH TRFIT RINGS: ", nf);
                printf("[Done]\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

            // ******************************************************
            // A-10-1. Perform TR fit with (VROT) free + VRAD fixed: if VRAD is intially fitted
            if(TRparam[0].vrad_nbreak_bspline != 1 || TRparam[0].vrad_order_bspline != 0) // if fitted: re-run TRfit with the fitted value of VRAD from bsplinefit_set_unipriors() to derive more stable(?) VROT. This is for deriving more robust n, r2, and rho2 below.
            {
                trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 999, "finalfit", 'N');
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

            // ******************************************************
            // A-11. Filter outliers in VROT out, and update priors
            if(rank == 0)
            {
                print_trfit_results(TRparam);
                reset_HI_VF_weight_trfit(TRparam);

                // restore TRparam[0].vrad_fix
                if(TRparam[0].vrad_nbreak_bspline == 1 && TRparam[0].vrad_order_bspline == 0) // not fitted
                    TRparam[0].vrad_fix = 'F';
                else
                    TRparam[0].vrad_fix = 'T';

                // Filter outliers out of VROT, and update priors
                // Filter out any outliers, derive initial estimates of Rc & rho0, and setup their priors
                bsplinefit_set_unipriors("VROT-einasto_halofit_rings", TRparam); // change the order : liner, quadrature or cubic 

                // initialise einasto halo params and update vrot_intp
                init_einastohalo_vrot_intp(TRparam);

                // make model vfs
                make_model_vfs(TRparam, fptr1, fname[0].fitsfile_2Dinput_VF, TRparam[0].nax1, TRparam[0].nax2, fname); 
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);
         }

        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // --------------------------------------------------------------------------------------------------------------------------------------------------------- //
        // -- PART B
        // B-1. Perform vEinasto 1D fit to update the priors of (n, r2, rho2)
        TRparam[0].loop_check = 0;
        int _tag_t=999;
        TRparam[0].n_hist_post = 0;
        while(TRparam[0].loop_check == 0)
        {
            // ******************************************************
            // B-2. First Einasto halo rotation curve 1D fit : dirty fit with the given ranges of priors
            v_einasto_1d_multinest4p(multinest_param, TRparam);
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

            // Update Einasto halo parameter priors
            if(rank == 0)
            {
                conv_i = 0;
                loop_iteration++;
                TRparam[0].n_hist_post = 999;
                v_einasto_1d_multinest4p_remove_outliers(TRparam);
                //update_Einasto_params_priors(TRparam, loop_iteration);
                TRparam[0].N_reliable_rings = TRparam[0].Nrings;
                update_vrot_prior(TRparam);

                if(fabs((_n_t-TRparam[0]._n)/TRparam[0]._n) < 0.2) conv_i++;
                if(fabs((_r_2_t-TRparam[0].r_2)/TRparam[0].r_2) < 0.2) conv_i++;
                if(fabs((_rho_2_t-TRparam[0].rho_2)/TRparam[0].rho_2) < 0.2) conv_i++;

                if(conv_i == 3) 
                {
                    TRparam[0].loop_check = 1; // exit loop
                    TRparam[0]._n_t = TRparam[0]._n;
                    TRparam[0]._r_2_t = TRparam[0].r_2;
                    TRparam[0]._rho_2_t = TRparam[0].rho_2;
                    TRparam[0].e_sigma_fitted_t = TRparam[0].e_sigma_fitted;

                    //printf("%f %f %f\n", TRparam[0]._n, TRparam[0].r_2, TRparam[0].rho_2);
                    for (k=1; k<n_node; k++)
                    {
                        MPI_Isend(&TRparam[0].loop_check, 1, MPI_INT, k, _tag_t, MPI_COMM_WORLD, &send_req[k]);
                        //MPI_Wait(&send_req[k], &mpistatus[k]);
                    }
                }
                else
                {
                    TRparam[0].loop_check = 0; // continue loop
                    // update temp einasto params
                    _n_t = TRparam[0]._n;
                    _r_2_t = TRparam[0].r_2;
                    _rho_2_t = TRparam[0].rho_2;

                    TRparam[0]._n_t = TRparam[0]._n;
                    TRparam[0]._r_2_t = TRparam[0].r_2;
                    TRparam[0]._rho_2_t = TRparam[0].rho_2;
                    TRparam[0].e_sigma_fitted_t = TRparam[0].e_sigma_fitted;
                    TRparam[0]._rho_2_t = TRparam[0].rho_2;

                    for (k=1; k<n_node; k++)
                    {
                        MPI_Isend(&TRparam[0].loop_check, 1, MPI_INT, k, _tag_t, MPI_COMM_WORLD, &send_req[k]);
                        //MPI_Wait(&send_req[k], &mpistatus[k]);
                    }
                }

                if(loop_iteration > 50)
                {
                    printf("Einasto 1D fit was not successful! Please check...\n");
                    TRparam[0].loop_check = 1; // exit loop
                }
            }
            else
            {
                MPI_Irecv(&TRparam[0].loop_check, 1, MPI_INT, 0, _tag_t, MPI_COMM_WORLD, &recv_req[rank]);
                MPI_Wait(&recv_req[rank], &mpistatus[rank]);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);
        }
        
        // ******************************************************
        // B-3. Reset with the multinest params for quick dirty Einasto fits
        update_multinest_params_dirtyfit(multinest_param);
        MPI_Barrier(MPI_COMM_WORLD);

        // ******************************************************
        // B-4. Reset the parameters for the dirty fit
        if(rank == 0)
        {
            define_area_tofit_set_geo_weight(0, TRparam[0].ellipse_semi_mx_boxfiltered*5.0, TRparam[0].decimX_einasto_halofit_d, TRparam[0].decimY_einasto_halofit_d, multinest_param, TRparam, 0, 0);

            // update fract_Navail_Nall + set free params to fit
            //set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "n", 'T', "r_2", 'T', "rho_2", 'T', "fullfit", 'N');
            set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'T', "ypos", 'T', "vsys", 'T', "pa", 'T', "incl", 'T', "n", 'T', "r_2", 'T', "rho_2", 'T', "fullfit", 'N');
        }
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

        // ******************************************************
        // B-5. DIRTY FIT : initial 2D einasto fit with weakly informative priors from the above 1D einasto fit
        einasto_halofit_multinest_student("xpos", TRparam[0].xpos_fix, "ypos", TRparam[0].ypos_fix, "vsys", TRparam[0].vsys_fix, "pa", TRparam[0].pa_fix, "incl", TRparam[0].incl_fix, "_n_Einasto", TRparam[0]._n_fix, "r_2_Einasto", TRparam[0].r_2_fix, "rho_2_Einasto", TRparam[0].rho_2_fix, "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, rank, "einastoFit.dirty.");

        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

        // reset the parameters for the final fit
        if(rank == 0)
        {
            printf("!+++ B-%d-1. DIRTY 2D EINASTO FIT : GAUSSIAN+UNIFORM PRIORS BASED ON the 1D EINASTO FIT +++\n", nf);
            // print the results and compute the BIC value from the einasto 2D fit
            print_priors_info("DIRTY", TRparam, n_node, ein_dirty+1);

            // Write down the MAP vaules of the parameters into the first line : this is for plotting later
            sprintf(fname[0].txtfile_multinest_output, "%s", "einastoFit.dirty..txt");
            sprintf(TRparam[0].txtfile_multinest_output, "%s", "einastoFit.dirty..txt");
            add_MAP_line_multinest_output(fname[0].txtfile_multinest_output, TRparam);

            // define are to fit & set geometrical weight
            define_area_tofit_set_geo_weight(0, TRparam[0].ellipse_semi_mx_boxfiltered*5.0, TRparam[0].decimX_einasto_halofit, TRparam[0].decimY_einasto_halofit, multinest_param, TRparam, 0, 0);
            set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'T', "incl", 'T', "n", 'T', "r_2", 'T', "rho_2", 'T', "fullfit", 'F'); // set full einast fit flag for using the fitted e-sigma value from the dirty fit (see the esigma prior in einasto_halofit_multinest_student() for more details)
            set_nfree_params_trfit_multinest_trfit_rings_student(TRparam);  
        }
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-6. Full fit : final einasto fit with a larger Nlive & multinset paramers supplied by the user + updated priors based on the dirty fit
        update_multinest_params_fullfit(multinest_param);
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

        einasto_halofit_multinest_student("xpos", TRparam[0].xpos_fix, "ypos", TRparam[0].ypos_fix, "vsys", TRparam[0].vsys_fix, "pa", TRparam[0].pa_fix, "incl", TRparam[0].incl_fix, "_n_Einasto", TRparam[0]._n_fix, "r_2_Einasto", TRparam[0].r_2_fix, "rho_2_Einasto", TRparam[0].rho_2_fix, "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, rank, "einastoFit.");

        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

        //printf("e_sigma_fitted : %f\n", TRparam[0].e_sigma_fitted);
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-7. Update the status of TR ring params
        if(rank == 0)
        {
            printf("!+++ B-%d-2. FULL 2D EINASTO FIT : GAUSSIAN+UNIFORM PRIORS BASED ON the DIRTY 2D EINASTO FIT +++\n", nf);
            print_priors_info("FULL", TRparam, n_node, ein_dirty+1);

            // update fract_Navail_Nall

            // Write down the MAP vaules of the parameters into the first line : this is for plotting later
            sprintf(fname[0].txtfile_multinest_output, "%s", "einastoFit..txt");
            sprintf(TRparam[0].txtfile_multinest_output, "%s", "einastoFit..txt");
            add_MAP_line_multinest_output(fname[0].txtfile_multinest_output, TRparam);

            // set vrad_function 'const' to make fits ring-by-ring : not a model based...
            strcpy(TRparam[0].vrad_function, "const");

            // set the final TR fit : n_free params : VROT = T
            //set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "n", 'F', "r_2", 'F', "rho_2", 'F', "fullfit", 'Y');
            set_navail_nall_pixels(TRparam, &Nall_ring, &Navail_ring, "xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "n", 'F', "r_2", 'F', "rho_2", 'F', "fullfit", 'Y');
            //set_nfree_params_trfit_multinest_trfit_rings_student(TRparam);  
            set_nfree_params_trfit_multinest_trfit_rings_student_test(TRparam);
            printf("!+++ B-%d-3. UPDATE PRIORS OF TILTED-RING PARAMETRES: ", nf);
            printf("[Done]\n\n");

            // Update VROT for TR fitting based on the Einasto fit
            update_vrot_prior(TRparam);

            // reset HI_VF_TRfit to zero
            reset_HI_VF_weight_trfit(TRparam);
        }

        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-8. Perform TR fitting with (VROT) free : RECEDING SIDE
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);
        //trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, -1, "finalfit", 'Y'); // receding side
        trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", 'T', "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, -1, "finalfit", 'Y'); // receding side
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);
        if(rank == 0)
        {
            printf("!+++ B-%d-4. PERFORM THE FINAL TILTED-RING FIT (RECEDING SIDE): ", nf);
            printf("[Done]\n");
            reset_HI_VF_weight_trfit(TRparam);
            print_trfit_results(TRparam);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);


        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-9. Perform TR fitting with (VROT) free : APPROACHING SIDE
        //trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 1, "finalfit", 'Y'); // approaching side
        trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", 'T', "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 1, "finalfit", 'Y'); // approaching side
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

        if(rank == 0)
        {
            printf("!+++ B-%d-5. PERFORM THE FINAL TILTED-RING FIT (APPROACHING SIDE): ", nf);
            printf("[Done]\n");
            reset_HI_VF_weight_trfit(TRparam);
            print_trfit_results(TRparam);
        }
        MPI_Barrier(MPI_COMM_WORLD);


        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // B-10. Perform TR fitting with (VROT) free : BOTH SIDE
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);
        //trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", TRparam[0].vrad_fix, "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, "finalfit", 'N'); // both sides
        trfit_multinest_trfit_rings_student("xpos", 'F', "ypos", 'F', "vsys", 'F', "pa", 'F', "incl", 'F', "vrot", 'T', "vrad", 'T', "sigmafactor", TRparam[0].sigma_factor_fix, multinest_param, TRparam, 0, "finalfit", 'N'); // both sides
        MPI_Barrier(MPI_COMM_WORLD);
        send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);
        if(rank == 0)
        {
            printf("!+++ B-%d-6. PERFORM THE FINAL TILTED-RING FIT (BOTH SIDES) ", nf);
            printf("[Done]\n");
            estimate_vrot_error(multinest_param, TRparam, 0, "finalfit", 'Y'); // both sides
            print_trfit_results(TRparam);

            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // B-11. Save the results and free memories
            printf("!+++ B-%d-7. SAVE THE FIT RESULTS AND MAKE MODEL VELOCITY FIELDS: ", nf);
            // write the TR fit results into a file
            write_trfit_results(TRparam, fname[0].bayesFit_outputfile);
            // make model 2dmaps
            make_model_vfs(TRparam, fptr1, fname[0].fitsfile_2Dinput_VF, TRparam[0].nax1, TRparam[0].nax2, fname); 
            printf("[Done]\n\n");

            printf("\t++ CHECK THE FIT RESULTS IN THE FOLLOWING DIRECTORY: ++\n\n");
            printf("\t%s/2dbat_output\n\n", TRparam[0].wdir);

            clock_t stop_time = clock();
            double elapsed = (double)(stop_time - start_time) / CLOCKS_PER_SEC;
            //printf("sigma_factor: %f\n", TRparam[0].sigma_factor);
            printf("THE EXECUTION TIME ELAPSED: %f SECS\n\n", elapsed);

            // free dynamically allocated memories
            free_2dmap_arrays(TRparam, multinest_param, fits_pointer);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalized(&flag);
    MPI_Finalize();
    return 0;
}

// --- End of line


