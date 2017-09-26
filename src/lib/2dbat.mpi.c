#include "2dbat.mpi.h"

// 2DBAT user defined functions
// MPI related

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// sending mpi data to other processes
void send_mpi_data(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi,  MPI_Status *mpistatus)
{
    int i=0;

    int a = 0, b=0;
    int tag1 = 1;
    int tag2 = 2;
    int tag3 = 3;
    int tag4 = 4;
    int tag5 = 5;
    int tag6 = 6;
    int tag7 = 7;
    int tag8 = 8;
    int tag9 = 9;
    int tag10 = 10;
    int tag11 = 11;
    int tag12 = 12;
    int tag13 = 13;

    if(rank == 0)
    {
        // --- D-8. Non-blocking sending TRparam_mpi to other processes
        // Wait until all messages have been sent
        // Note that we use requests and statuses starting from index 1
        for (i=1; i<n_node; i++)
        {
            MPI_Isend(TRparam, 1, TRparam_mpi, i, tag1, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_boxFiltered[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag2, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_boxFiltered_sigma[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag3, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_boxFiltered_vlos_ew[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag13, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_boxFiltered_sigma_e_norm[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag4, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_mom2[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag5, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_boxFiltered_decim0[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag6, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag7, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_geo_radial_angle_w[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag8, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_weight_TRfit[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag9, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_fract_navail_nall[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag10, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_einasto_halomodel[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag11, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }

        for (i=1; i<n_node; i++)
        {
            MPI_Isend(&(HI_VF_sigma_geo_flux_weighted[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag12, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus[i]);
        }
    }
    else
    {
        MPI_Irecv(TRparam, 1, TRparam_mpi, 0, tag1, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_boxFiltered[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag2, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_boxFiltered_sigma[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag3, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_boxFiltered_vlos_ew[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag13, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_boxFiltered_sigma_e_norm[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag4, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_mom2[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag5, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_boxFiltered_decim0[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag6, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag7, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_geo_radial_angle_w[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag8, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_weight_TRfit[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag9, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_fract_navail_nall[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag10, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_einasto_halomodel[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag11, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);

        MPI_Irecv(&(HI_VF_sigma_geo_flux_weighted[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag12, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus[rank]);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        //printf("MPI successfully sent/received\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// receiving mpi data from other processes
void receive_mpiData(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi, velocity_field *HI_VF_boxFiltered, velocity_field *HI_VF_boxFiltered_sigma, velocity_field *HI_VF_boxFiltered_sigma_e_norm, velocity_field *HI_VF_mom2, velocity_field *HI_VF_boxFiltered_decim0, velocity_field *HI_VF, MPI_Status *mpistatus)
{
    int i=0;

    int a = 0, b=0;
    int tag1 = 1;
    int tag2 = 2;
    int tag3 = 3;
    int tag4 = 4;
    int tag5 = 5;
    int tag6 = 6;
    int tag7 = 7;
    int tag8 = 8;
    int tag9 = 9;
    int tag10 = 10;

    if(rank == 0)
    {
        // --- D-8. Non-blocking sending TRparam_mpi to other processes
        // Wait until all messages have been sent
        // Note that we use requests and statuses starting from index 1
        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(TRparam, 1, TRparam_mpi, i, tag1, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);

        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(&(HI_VF_boxFiltered[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag2, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);

        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(&(HI_VF_boxFiltered_sigma[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag3, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);

        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(&(HI_VF_boxFiltered_sigma_e_norm[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag4, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);

        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(&(HI_VF_mom2[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag5, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);

        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(&(HI_VF_boxFiltered_decim0[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag6, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);

        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(&(HI_VF[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag7, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);
        printf("\n");
    }
    else
    {
        MPI_Isend(TRparam, 1, TRparam_mpi, 0, tag1, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);

        MPI_Isend(&(HI_VF_boxFiltered[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag2, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);

        MPI_Isend(&(HI_VF_boxFiltered_sigma[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag3, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);

        MPI_Isend(&(HI_VF_boxFiltered_sigma_e_norm[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag4, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);

        MPI_Isend(&(HI_VF_mom2[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag5, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);

        MPI_Isend(&(HI_VF_boxFiltered_decim0[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag6, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);

        MPI_Isend(&(HI_VF[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, tag7, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);
    }
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// receiving mpi data from other processes
void receive_mpiData_TRparam(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi, velocity_field *HI_VF_boxFiltered, velocity_field *HI_VF_boxFiltered_sigma, velocity_field *HI_VF_boxFiltered_sigma_e_norm, velocity_field *HI_VF_mom2, velocity_field *HI_VF_boxFiltered_decim0, velocity_field *HI_VF, MPI_Status *mpistatus)
{
    int i=0;

    int a = 0, b=0;
    int tag1 = 1;
    int tag2 = 2;
    int tag3 = 3;
    int tag4 = 4;
    int tag5 = 5;
    int tag6 = 6;
    int tag7 = 7;
    int tag8 = 8;
    int tag9 = 9;
    int tag10 = 10;

    if(rank == 0)
    {
        // --- D-8. Non-blocking sending TRparam_mpi to other processes
        // Wait until all messages have been sent
        // Note that we use requests and statuses starting from index 1
        for (i=1; i<n_node; i++)
        {
            MPI_Irecv(TRparam, 1, TRparam_mpi, i, tag1, MPI_COMM_WORLD, &recv_req[i]);
        }
        MPI_Waitall(n_node-1, &recv_req[1], &mpistatus[1]);
        printf("\n");
    }
    else
    {
        MPI_Isend(TRparam, 1, TRparam_mpi, 0, tag1, MPI_COMM_WORLD, &send_req[0]);
        MPI_Wait(&send_req[0], &mpistatus[0]);
    }
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MPI_Datatype allocate_mpi_dataset(TR_ringParameters *TRparam, MPI_Datatype *type, int *blocklen)
{
    // --------------------------------------------------------------------------------------------- //
    // +++ D. START PROCESS +++
    // --------------------------------------------------------------------------------------------- //
    // --- D-0. Declare MPI struct ---

    MPI_Datatype TRparam_mpi;

    MPI_Aint disp[302];
    MPI_Aint address1, address2;
    MPI_Get_address(TRparam, &address1);

    // 1. wdir
    MPI_Get_address(&TRparam[0].wdir, &address2);
    disp[0] = address2 - address1;
    // 2. txtfile_multinest_output
    MPI_Get_address(&TRparam[0].txtfile_multinest_output, &address2);
    disp[1] = address2 - address1;
    // 3. loop_check
    MPI_Get_address(&TRparam[0].loop_check, &address2);
    disp[2] = address2 - address1;
    // 4. n_freeParams
    MPI_Get_address(&TRparam[0].n_freeParams, &address2);
    disp[3] = address2 - address1;
    // 5. e_sigma
    MPI_Get_address(&TRparam[0].e_sigma, &address2);
    disp[4] = address2 - address1;
    // 6. e_sigma_fitted
    MPI_Get_address(&TRparam[0].e_sigma_fitted, &address2);
    disp[5] = address2 - address1;
    // 7. e_sigma_fitted_t
    MPI_Get_address(&TRparam[0].e_sigma_fitted_t, &address2);
    disp[6] = address2 - address1;
    // 8. maxLogLike
    MPI_Get_address(&TRparam[0].maxLogLike, &address2);
    disp[7] = address2 - address1;
    // 9. logZ
    MPI_Get_address(&TRparam[0].logZ, &address2);
    disp[8] = address2 - address1;
    // 10. logZerr
    MPI_Get_address(&TRparam[0].logZerr, &address2);
    disp[9] = address2 - address1;
    // 11. maxLogLikeF
    MPI_Get_address(&TRparam[0].maxLogLikeF, &address2);
    disp[10] = address2 - address1;
    // 12. logZF
    MPI_Get_address(&TRparam[0].logZF, &address2);
    disp[11] = address2 - address1;
    // 13. logZerrF
    MPI_Get_address(&TRparam[0].logZerrF, &address2);
    disp[12] = address2 - address1;
    // 14. nax1
    MPI_Get_address(&TRparam[0].nax1, &address2);
    disp[13] = address2 - address1;
    // 15. nax2
    MPI_Get_address(&TRparam[0].nax2, &address2);
    disp[14] = address2 - address1;
    // 16. decim_x0
    MPI_Get_address(&TRparam[0].decim_x0, &address2);
    disp[15] = address2 - address1;
    // 17. decim_y0
    MPI_Get_address(&TRparam[0].decim_y0, &address2);
    disp[16] = address2 - address1;
    // 18. decimX
    MPI_Get_address(&TRparam[0].decimX, &address2);
    disp[17] = address2 - address1;
    // 19. decimY
    MPI_Get_address(&TRparam[0].decimY, &address2);
    disp[18] = address2 - address1;
    // 20. decimX_einasto_halofit
    MPI_Get_address(&TRparam[0].decimX_einasto_halofit, &address2);
    disp[19] = address2 - address1;
    // 21. decimY_einasto_halofit
    MPI_Get_address(&TRparam[0].decimY_einasto_halofit, &address2);
    disp[20] = address2 - address1;
    // 22. decimX_trfit
    MPI_Get_address(&TRparam[0].decimX_trfit, &address2);
    disp[21] = address2 - address1;
    // 23. decimY_trfit
    MPI_Get_address(&TRparam[0].decimY_trfit, &address2);
    disp[22] = address2 - address1;
    // 24. pixelScale
    MPI_Get_address(&TRparam[0].pixelScale, &address2);
    disp[23] = address2 - address1;
    // 25. N_all_pixels_in_a_ring
    MPI_Get_address(&TRparam[0].N_all_pixels_in_a_ring, &address2);
    disp[24] = address2 - address1;
    // 26. ellipse_xpos_boxfiltered
    MPI_Get_address(&TRparam[0].ellipse_xpos_boxfiltered, &address2);
    disp[25] = address2 - address1;
    // 27. ellipse_ypos_boxfiltered
    MPI_Get_address(&TRparam[0].ellipse_ypos_boxfiltered, &address2);
    disp[26] = address2 - address1;
    // 28. ellipse_pa_boxfiltered
    MPI_Get_address(&TRparam[0].ellipse_pa_boxfiltered, &address2);
    disp[27] = address2 - address1;
    // 29. pa_EllipseFit_e
    MPI_Get_address(&TRparam[0].pa_EllipseFit_e, &address2);
    disp[28] = address2 - address1;
    // 30. ellipse_incl_boxfiltered
    MPI_Get_address(&TRparam[0].ellipse_incl_boxfiltered, &address2);
    disp[29] = address2 - address1;
    // 31. incl_EllipseFit_e
    MPI_Get_address(&TRparam[0].incl_EllipseFit_e, &address2);
    disp[30] = address2 - address1;
    // 32. ellipse_semi_mx_boxfiltered
    MPI_Get_address(&TRparam[0].ellipse_semi_mx_boxfiltered, &address2);
    disp[31] = address2 - address1;
    // 33. xpos1_from_ellipsefit
    MPI_Get_address(&TRparam[0].xpos1_from_ellipsefit, &address2);
    disp[32] = address2 - address1;
    // 34. xpos2_from_ellipsefit
    MPI_Get_address(&TRparam[0].xpos2_from_ellipsefit, &address2);
    disp[33] = address2 - address1;
    // 35. ypos1_from_ellipsefit
    MPI_Get_address(&TRparam[0].ypos1_from_ellipsefit, &address2);
    disp[34] = address2 - address1;
    // 36. ypos2_from_ellipsefit
    MPI_Get_address(&TRparam[0].ypos2_from_ellipsefit, &address2);
    disp[35] = address2 - address1;
    // 37. vsys1_from_vlosfit 
    MPI_Get_address(&TRparam[0].vsys1_from_vlosfit, &address2);
    disp[36] = address2 - address1;
    // 38. vsys2_from_vlosfit
    MPI_Get_address(&TRparam[0].vsys2_from_vlosfit, &address2);
    disp[37] = address2 - address1;
    // 39. perimeter
    MPI_Get_address(&TRparam[0].perimeter, &address2);
    disp[38] = address2 - address1;
    // 40. final_fit
    MPI_Get_address(&TRparam[0].final_fit, &address2);
    disp[39] = address2 - address1;
    // 41. fullFit
    MPI_Get_address(&TRparam[0].fullFit, &address2);
    disp[40] = address2 - address1;
    // 42. Nrings
    MPI_Get_address(&TRparam[0].Nrings, &address2);
    disp[41] = address2 - address1;
    // 43. Nrings_intp
    MPI_Get_address(&TRparam[0].Nrings_intp, &address2);
    disp[42] = address2 - address1;
    // 44. Nrings_to_semi_mx
    MPI_Get_address(&TRparam[0].Nrings_to_semi_mx, &address2);
    disp[43] = address2 - address1;
    // 45. N_reliable_rings
    MPI_Get_address(&TRparam[0].N_reliable_rings, &address2);
    disp[44] = address2 - address1;
    // 46. tilted_ring
    MPI_Get_address(&TRparam[0].tilted_ring, &address2);
    disp[45] = address2 - address1;
    // 47. Npoints_in_tilted_ring
    MPI_Get_address(&TRparam[0].Npoints_in_tilted_ring, &address2);
    disp[46] = address2 - address1;
    // 48. total_Npoints_allRings
    MPI_Get_address(&TRparam[0].total_Npoints_allRings, &address2);
    disp[47] = address2 - address1;
    // 49. ring_radius
    MPI_Get_address(&TRparam[0].ring_radius, &address2);
    disp[48] = address2 - address1;
    // 50. ring_intp
    MPI_Get_address(&TRparam[0].ring_intp, &address2);
    disp[49] = address2 - address1;
    // 51. ring_s
    MPI_Get_address(&TRparam[0].ring_s, &address2);
    disp[50] = address2 - address1;
    // 52. ring_e
    MPI_Get_address(&TRparam[0].ring_e, &address2);
    disp[51] = address2 - address1;
    // 53. ring_w
    MPI_Get_address(&TRparam[0].ring_w, &address2);
    disp[52] = address2 - address1;
    // 54. ring_s_for_einasto1Dfit
    MPI_Get_address(&TRparam[0].ring_s_for_einasto1Dfit, &address2);
    disp[53] = address2 - address1;
    // 55. ring_e_for_einasto1Dfit
    MPI_Get_address(&TRparam[0].ring_e_for_einasto1Dfit, &address2);
    disp[54] = address2 - address1;
    // 56. ring_w_for_einasto1Dfit
    MPI_Get_address(&TRparam[0].ring_w_for_einasto1Dfit, &address2);
    disp[55] = address2 - address1;
    // 57. Npoints_in_tilted_ring_decim0
    MPI_Get_address(&TRparam[0].Npoints_in_tilted_ring_decim0, &address2);
    disp[56] = address2 - address1;
    // 58. npoints_inaring
    MPI_Get_address(&TRparam[0].npoints_inaring, &address2);
    disp[57] = address2 - address1;
    // 59. npoints_inaring_decim0
    MPI_Get_address(&TRparam[0].npoints_inaring_decim0, &address2);
    disp[58] = address2 - address1;
    // 60. rGalaxyPlane_pixel
    MPI_Get_address(&TRparam[0].rGalaxyPlane_pixel, &address2);
    disp[59] = address2 - address1;
    // 61. rGalaxyPlane_pixel_max
    MPI_Get_address(&TRparam[0].rGalaxyPlane_pixel_max, &address2);
    disp[60] = address2 - address1;
    // 62. free_angle
    MPI_Get_address(&TRparam[0].free_angle, &address2);
    disp[61] = address2 - address1;
    // 63. wpow
    MPI_Get_address(&TRparam[0].wpow, &address2);
    disp[62] = address2 - address1;
    // 64. rwpow
    MPI_Get_address(&TRparam[0].rwpow, &address2);
    disp[63] = address2 - address1;
    // 65. xpos_fix
    MPI_Get_address(&TRparam[0].xpos_fix, &address2);
    disp[64] = address2 - address1;
    // 66. xposF
    MPI_Get_address(&TRparam[0].xposF, &address2);
    disp[65] = address2 - address1;
    // 67. xposF_e
    MPI_Get_address(&TRparam[0].xposF_e, &address2);
    disp[66] = address2 - address1;
    // 68. xpos0
    MPI_Get_address(&TRparam[0].xpos0, &address2);
    disp[67] = address2 - address1;
    // 69. xpos
    MPI_Get_address(&TRparam[0].xpos, &address2);
    disp[68] = address2 - address1;
    // 70. xpos_e
    MPI_Get_address(&TRparam[0].xpos_e, &address2);
    disp[69] = address2 - address1;
    // 71. xpos_temp
    MPI_Get_address(&TRparam[0].xpos_temp, &address2);
    disp[70] = address2 - address1;
    // 72. xpos_mode
    MPI_Get_address(&TRparam[0].xpos_mode, &address2);
    disp[71] = address2 - address1;
    // 73. xpos_std
    MPI_Get_address(&TRparam[0].xpos_std, &address2);
    disp[72] = address2 - address1;
    // 74. xpos1
    MPI_Get_address(&TRparam[0].xpos1, &address2);
    disp[73] = address2 - address1;
    // 75. xpos2
    MPI_Get_address(&TRparam[0].xpos2, &address2);
    disp[74] = address2 - address1;
    // 76. xposF_EinastoFit
    MPI_Get_address(&TRparam[0].xposF_EinastoFit, &address2);
    disp[75] = address2 - address1;
    // 77. xposF_EinastoFit_e
    MPI_Get_address(&TRparam[0].xposF_EinastoFit_e, &address2);
    disp[76] = address2 - address1;
    // 78. xposF_EinastoFit_t
    MPI_Get_address(&TRparam[0].xposF_EinastoFit_t, &address2);
    disp[77] = address2 - address1;
    // 79. xposF_EinastoFit_e_t
    MPI_Get_address(&TRparam[0].xposF_EinastoFit_e_t, &address2);
    disp[78] = address2 - address1;
    // 80. ypos_fix
    MPI_Get_address(&TRparam[0].ypos_fix, &address2);
    disp[79] = address2 - address1;
    // 81. yposF
    MPI_Get_address(&TRparam[0].yposF, &address2);
    disp[80] = address2 - address1;
    // 82. yposF_e
    MPI_Get_address(&TRparam[0].yposF_e, &address2);
    disp[81] = address2 - address1;
    // 83. ypos0
    MPI_Get_address(&TRparam[0].ypos0, &address2);
    disp[82] = address2 - address1;
    // 84. ypos
    MPI_Get_address(&TRparam[0].ypos, &address2);
    disp[83] = address2 - address1;
    // 85. ypos_e
    MPI_Get_address(&TRparam[0].ypos_e, &address2);
    disp[84] = address2 - address1;
    // 86. ypos_temp
    MPI_Get_address(&TRparam[0].ypos_temp, &address2);
    disp[85] = address2 - address1;
    // 87. ypos_mode
    MPI_Get_address(&TRparam[0].ypos_mode, &address2);
    disp[86] = address2 - address1;
    // 88. ypos_std
    MPI_Get_address(&TRparam[0].ypos_std, &address2);
    disp[87] = address2 - address1;
    // 89. ypos1
    MPI_Get_address(&TRparam[0].ypos1, &address2);
    disp[88] = address2 - address1;
    // 90. ypos2
    MPI_Get_address(&TRparam[0].ypos2, &address2);
    disp[89] = address2 - address1;
    // 91. yposF_EinastoFit
    MPI_Get_address(&TRparam[0].yposF_EinastoFit, &address2);
    disp[90] = address2 - address1;
    // 92. yposF_EinastoFit_e
    MPI_Get_address(&TRparam[0].yposF_EinastoFit_e, &address2);
    disp[91] = address2 - address1;
    // 93. yposF_EinastoFit_t
    MPI_Get_address(&TRparam[0].yposF_EinastoFit_t, &address2);
    disp[92] = address2 - address1;
    // 94. yposF_EinastoFit_e_t
    MPI_Get_address(&TRparam[0].yposF_EinastoFit_e_t, &address2);
    disp[93] = address2 - address1;
    // 95. vsys_fix
    MPI_Get_address(&TRparam[0].vsys_fix, &address2);
    disp[94] = address2 - address1;
    // 96. vsysF
    MPI_Get_address(&TRparam[0].vsysF, &address2);
    disp[95] = address2 - address1;
    // 97. vsysF_e
    MPI_Get_address(&TRparam[0].vsysF_e, &address2);
    disp[96] = address2 - address1;
    // 98. vsys0
    MPI_Get_address(&TRparam[0].vsys0, &address2);
    disp[97] = address2 - address1;
    // 99. vsys
    MPI_Get_address(&TRparam[0].vsys, &address2);
    disp[98] = address2 - address1;
    // 100. vsys_e
    MPI_Get_address(&TRparam[0].vsys_e, &address2);
    disp[99] = address2 - address1;
    // 101. vsys_temp
    MPI_Get_address(&TRparam[0].vsys_temp, &address2);
    disp[100] = address2 - address1;
    // 102. vsys_mode
    MPI_Get_address(&TRparam[0].vsys_mode, &address2);
    disp[101] = address2 - address1;
    // 103. vsys_std
    MPI_Get_address(&TRparam[0].vsys_std, &address2);
    disp[102] = address2 - address1;
    // 104. vsys1
    MPI_Get_address(&TRparam[0].vsys1, &address2);
    disp[103] = address2 - address1;
    // 105. vsys2
    MPI_Get_address(&TRparam[0].vsys2, &address2);
    disp[104] = address2 - address1;
    // 106. vsysF_EinastoFit
    MPI_Get_address(&TRparam[0].vsysF_EinastoFit, &address2);
    disp[105] = address2 - address1;
    // 107. vsysF_EinastoFit_e
    MPI_Get_address(&TRparam[0].vsysF_EinastoFit_e, &address2);
    disp[106] = address2 - address1;
    // 108. vsysF_EinastoFit_t
    MPI_Get_address(&TRparam[0].vsysF_EinastoFit_t, &address2);
    disp[107] = address2 - address1;
    // 109. vsysF_EinastoFit_e_t
    MPI_Get_address(&TRparam[0].vsysF_EinastoFit_e_t, &address2);
    disp[108] = address2 - address1;
    // 110. pa_fix
    MPI_Get_address(&TRparam[0].pa_fix, &address2);
    disp[109] = address2 - address1;
    // 111. paF
    MPI_Get_address(&TRparam[0].paF, &address2);
    disp[110] = address2 - address1;
    // 112. paF_e
    MPI_Get_address(&TRparam[0].paF_e, &address2);
    disp[111] = address2 - address1;
    // 113. pa0
    MPI_Get_address(&TRparam[0].pa0, &address2);
    disp[112] = address2 - address1;
    // 114. pa
    MPI_Get_address(&TRparam[0].pa, &address2);
    disp[113] = address2 - address1;
    // 115. pa_e
    MPI_Get_address(&TRparam[0].pa_e, &address2);
    disp[114] = address2 - address1;
    // 116. pa_temp
    MPI_Get_address(&TRparam[0].pa_temp, &address2);
    disp[115] = address2 - address1;
    // 117. pa_temp_e
    MPI_Get_address(&TRparam[0].pa_temp_e, &address2);
    disp[116] = address2 - address1;
    // 118. pa1
    MPI_Get_address(&TRparam[0].pa1, &address2);
    disp[117] = address2 - address1;
    // 119. pa2
    MPI_Get_address(&TRparam[0].pa2, &address2);
    disp[118] = address2 - address1;
    // 120. pa1_for_TRfit
    MPI_Get_address(&TRparam[0].pa1_for_TRfit, &address2);
    disp[119] = address2 - address1;
    // 121..pa2_for_TRfit
    MPI_Get_address(&TRparam[0].pa2_for_TRfit, &address2);
    disp[120] = address2 - address1;
    // 122. _p1_tr
    MPI_Get_address(&TRparam[0]._p1_tr, &address2);
    disp[121] = address2 - address1;
    // 123. _p2_tr
    MPI_Get_address(&TRparam[0]._p2_tr, &address2);
    disp[122] = address2 - address1;
    // 124. PA_MAX_in_degree
    MPI_Get_address(&TRparam[0].PA_MAX_in_degree, &address2);
    disp[123] = address2 - address1;
    // 125. paF_EinastoFit
    MPI_Get_address(&TRparam[0].paF_EinastoFit, &address2);
    disp[124] = address2 - address1;
    // 126. paF_EinastoFit_e
    MPI_Get_address(&TRparam[0].paF_EinastoFit_e, &address2);
    disp[125] = address2 - address1;
    // 127. pa_function
    MPI_Get_address(&TRparam[0].pa_function, &address2);
    disp[126] = address2 - address1;
    // 128. n_coeffs_bspline_pa
    MPI_Get_address(&TRparam[0].n_coeffs_bspline_pa, &address2);
    disp[127] = address2 - address1;
    // 129. pa_nbreak_bspline
    MPI_Get_address(&TRparam[0].pa_nbreak_bspline, &address2);
    disp[128] = address2 - address1;
    // 130. pa_order_bspline
    MPI_Get_address(&TRparam[0].pa_order_bspline, &address2);
    disp[129] = address2 - address1;
    // 131. _p_bs
    MPI_Get_address(&TRparam[0]._p_bs, &address2);
    disp[130] = address2 - address1;
    // 132. _p_bs_e
    MPI_Get_address(&TRparam[0]._p_bs_e, &address2);
    disp[131] = address2 - address1;
    // 133. _p_bs_t
    MPI_Get_address(&TRparam[0]._p_bs_t, &address2);
    disp[132] = address2 - address1;
    // 134. _p_bs_e_t
    MPI_Get_address(&TRparam[0]._p_bs_e_t, &address2);
    disp[133] = address2 - address1;
    // 135. bspline1pa
    MPI_Get_address(&TRparam[0].bspline1pa, &address2);
    disp[134] = address2 - address1;
    // 136. bspline2pa
    MPI_Get_address(&TRparam[0].bspline2pa, &address2);
    disp[135] = address2 - address1;
    // 137. _p_bs_tr
    MPI_Get_address(&TRparam[0]._p_bs_tr, &address2);
    disp[136] = address2 - address1;
    // 138. _bspline_pa_hist_sigma
    MPI_Get_address(&TRparam[0]._bspline_pa_hist_sigma, &address2);
    disp[137] = address2 - address1;
    // 139. incl_fix
    MPI_Get_address(&TRparam[0].incl_fix, &address2);
    disp[138] = address2 - address1;
    // 140. inclF
    MPI_Get_address(&TRparam[0].inclF, &address2);
    disp[139] = address2 - address1;
    // 141. inclF_e
    MPI_Get_address(&TRparam[0].inclF_e, &address2);
    disp[140] = address2 - address1;
    // 142. incl0
    MPI_Get_address(&TRparam[0].incl0, &address2);
    disp[141] = address2 - address1;
    // 143. incl
    MPI_Get_address(&TRparam[0].incl, &address2);
    disp[142] = address2 - address1;
    // 144. incl_e
    MPI_Get_address(&TRparam[0].incl_e, &address2);
    disp[143] = address2 - address1;
    // 145. incl_temp
    MPI_Get_address(&TRparam[0].incl_temp, &address2);
    disp[144] = address2 - address1;
    // 146. incl_temp_e
    MPI_Get_address(&TRparam[0].incl_temp_e, &address2);
    disp[145] = address2 - address1;
    // 147. incl1
    MPI_Get_address(&TRparam[0].incl1, &address2);
    disp[146] = address2 - address1;
    // 148. incl2
    MPI_Get_address(&TRparam[0].incl2, &address2);
    disp[147] = address2 - address1;
    // 149. incl1_for_TRfit
    MPI_Get_address(&TRparam[0].incl1_for_TRfit, &address2);
    disp[148] = address2 - address1;
    // 150. incl2_for_TRfit
    MPI_Get_address(&TRparam[0].incl2_for_TRfit, &address2);
    disp[149] = address2 - address1;
    // 151. _i1_tr
    MPI_Get_address(&TRparam[0]._i1_tr, &address2);
    disp[150] = address2 - address1;
    // 152. _i2_tr
    MPI_Get_address(&TRparam[0]._i2_tr, &address2);
    disp[151] = address2 - address1;
    // 153. INCL_MAX_in_degree
    MPI_Get_address(&TRparam[0].INCL_MAX_in_degree, &address2);
    disp[152] = address2 - address1;
    // 154. inclF_EinastoFit
    MPI_Get_address(&TRparam[0].inclF_EinastoFit, &address2);
    disp[153] = address2 - address1;
    // 155. inclF_EinastoFit_e
    MPI_Get_address(&TRparam[0].inclF_EinastoFit_e, &address2);
    disp[154] = address2 - address1;
    // 156. incl_function
    MPI_Get_address(&TRparam[0].incl_function, &address2);
    disp[155] = address2 - address1;
    // 157. n_coeffs_bspline_incl
    MPI_Get_address(&TRparam[0].n_coeffs_bspline_incl, &address2);
    disp[156] = address2 - address1;
    // 158. incl_nbreak_bspline
    MPI_Get_address(&TRparam[0].incl_nbreak_bspline, &address2);
    disp[157] = address2 - address1;
    // 159. incl_order_bspline
    MPI_Get_address(&TRparam[0].incl_order_bspline, &address2);
    disp[158] = address2 - address1;
    // 160. _i_bs
    MPI_Get_address(&TRparam[0]._i_bs, &address2);
    disp[159] = address2 - address1;
    // 161. _i_bs_e
    MPI_Get_address(&TRparam[0]._i_bs_e, &address2);
    disp[160] = address2 - address1;
    // 162. _i_bs_t
    MPI_Get_address(&TRparam[0]._i_bs_t, &address2);
    disp[161] = address2 - address1;
    // 163. _i_bs_e_t
    MPI_Get_address(&TRparam[0]._i_bs_e_t, &address2);
    disp[162] = address2 - address1;
    // 164. bspline1incl
    MPI_Get_address(&TRparam[0].bspline1incl, &address2);
    disp[163] = address2 - address1;
    // 165. bspline2incl
    MPI_Get_address(&TRparam[0].bspline2incl, &address2);
    disp[164] = address2 - address1;
    // 166. _i_bs_tr
    MPI_Get_address(&TRparam[0]._i_bs_tr, &address2);
    disp[165] = address2 - address1;
    // 167. _bspline_incl_hist_sigma
    MPI_Get_address(&TRparam[0]._bspline_incl_hist_sigma, &address2);
    disp[166] = address2 - address1;
    // 168. vrot_fix
    MPI_Get_address(&TRparam[0].vrot_fix, &address2);
    disp[167] = address2 - address1;
    // 169. vrotF
    MPI_Get_address(&TRparam[0].vrotF, &address2);
    disp[168] = address2 - address1;
    // 170. vrotF_e
    MPI_Get_address(&TRparam[0].vrotF_e, &address2);
    disp[169] = address2 - address1;
    // 171. vrot1
    MPI_Get_address(&TRparam[0].vrot1, &address2);
    disp[170] = address2 - address1;
    // 172. vrot2
    MPI_Get_address(&TRparam[0].vrot2, &address2);
    disp[171] = address2 - address1;
    // 173. vrot0
    MPI_Get_address(&TRparam[0].vrot0, &address2);
    disp[172] = address2 - address1;
    // 174. vrot 
    MPI_Get_address(&TRparam[0].vrot, &address2);
    disp[173] = address2 - address1;
    // 175. vrot_e
    MPI_Get_address(&TRparam[0].vrot_e, &address2);
    disp[174] = address2 - address1;
    // 176. vrot_rec
    MPI_Get_address(&TRparam[0].vrot_rec, &address2);
    disp[175] = address2 - address1;
    // 177. vrot_e_rec
    MPI_Get_address(&TRparam[0].vrot_e_rec, &address2);
    disp[176] = address2 - address1;
    // 178. vrot_app
    MPI_Get_address(&TRparam[0].vrot_app, &address2);
    disp[177] = address2 - address1;
    // 179. vrot_e_app
    MPI_Get_address(&TRparam[0].vrot_e_app, &address2);
    disp[178] = address2 - address1;
    // 180. vrot_temp
    MPI_Get_address(&TRparam[0].vrot_temp, &address2);
    disp[179] = address2 - address1;
    // 181. vrot_temp_e
    MPI_Get_address(&TRparam[0].vrot_temp_e, &address2);
    disp[180] = address2 - address1;
    // 182. vrot_intp
    MPI_Get_address(&TRparam[0].vrot_intp, &address2);
    disp[181] = address2 - address1;
    // 183. vrot_e_intp 
    MPI_Get_address(&TRparam[0].vrot_e_intp, &address2);
    disp[182] = address2 - address1;
    // 184. vrad_fix
    MPI_Get_address(&TRparam[0].vrad_fix, &address2);
    disp[183] = address2 - address1;
    // 185. vrad1
    MPI_Get_address(&TRparam[0].vrad1, &address2);
    disp[184] = address2 - address1;
    // 186. vrad2
    MPI_Get_address(&TRparam[0].vrad2, &address2);
    disp[185] = address2 - address1;
    // 187. vradF
    MPI_Get_address(&TRparam[0].vradF, &address2);
    disp[186] = address2 - address1;
    // 188. vradF_e
    MPI_Get_address(&TRparam[0].vradF_e, &address2);
    disp[187] = address2 - address1;
    // 189. vrad_max
    MPI_Get_address(&TRparam[0].vrad_max, &address2);
    disp[188] = address2 - address1;
    // 190. vrad
    MPI_Get_address(&TRparam[0].vrad, &address2);
    disp[189] = address2 - address1;
    // 191. vrad_e
    MPI_Get_address(&TRparam[0].vrad_e, &address2);
    disp[190] = address2 - address1;
    // 192. vrad_rec
    MPI_Get_address(&TRparam[0].vrad_rec, &address2);
    disp[191] = address2 - address1;
    // 193. vrad_rec_e
    MPI_Get_address(&TRparam[0].vrad_rec_e, &address2);
    disp[192] = address2 - address1;
    // 194. vrad_app
    MPI_Get_address(&TRparam[0].vrad_app, &address2);
    disp[193] = address2 - address1;
    // 195. vrad_app_e
    MPI_Get_address(&TRparam[0].vrad_app_e, &address2);
    disp[194] = address2 - address1;
    // 196. vrad_function
    MPI_Get_address(&TRparam[0].vrad_function, &address2);
    disp[195] = address2 - address1;
    // 197. n_coeffs_bspline_vrad
    MPI_Get_address(&TRparam[0].n_coeffs_bspline_vrad, &address2);
    disp[196] = address2 - address1;
    // 198. vrad_nbreak_bspline
    MPI_Get_address(&TRparam[0].vrad_nbreak_bspline, &address2);
    disp[197] = address2 - address1;
    // 199. vrad_order_bspline
    MPI_Get_address(&TRparam[0].vrad_order_bspline, &address2);
    disp[198] = address2 - address1;
    // 200. bspline1vrad
    MPI_Get_address(&TRparam[0].bspline1vrad, &address2);
    disp[199] = address2 - address1;
    // 201. bspline2vrad
    MPI_Get_address(&TRparam[0].bspline2vrad, &address2);
    disp[200] = address2 - address1;
    // 202. _vr_bs
    MPI_Get_address(&TRparam[0]._vr_bs, &address2);
    disp[201] = address2 - address1;
    // 203. _vr_bs_e
    MPI_Get_address(&TRparam[0]._vr_bs_e, &address2);
    disp[202] = address2 - address1;
    // 204. _vr_bs_t
    MPI_Get_address(&TRparam[0]._vr_bs_t, &address2);
    disp[203] = address2 - address1;
    // 205. _vr_bs_e_t
    MPI_Get_address(&TRparam[0]._vr_bs_e_t, &address2);
    disp[204] = address2 - address1;
    // 206. _vr_bs_tr
    MPI_Get_address(&TRparam[0]._vr_bs_tr, &address2);
    disp[205] = address2 - address1;
    // 207. vrad_temp
    MPI_Get_address(&TRparam[0].vrad_temp, &address2);
    disp[206] = address2 - address1;
    // 208. vrad_temp_e
    MPI_Get_address(&TRparam[0].vrad_temp_e, &address2);
    disp[207] = address2 - address1;
    // 209. vrad_einastofit_bs
    MPI_Get_address(&TRparam[0].vrad_einastofit_bs, &address2);
    disp[208] = address2 - address1;
    // 210. vrad_einastofit_bs_e  
    MPI_Get_address(&TRparam[0].vrad_einastofit_bs_e, &address2);
    disp[209] = address2 - address1;
    // 211. _bspline_vrad_hist_sigma
    MPI_Get_address(&TRparam[0]._bspline_vrad_hist_sigma, &address2);
    disp[210] = address2 - address1;
    // 212. n_gauss
    MPI_Get_address(&TRparam[0].n_gauss, &address2);
    disp[211] = address2 - address1;
    // 213. g_param
    MPI_Get_address(&TRparam[0].g_param, &address2);
    disp[212] = address2 - address1;
    // 214. g01
    MPI_Get_address(&TRparam[0].g01, &address2);
    disp[213] = address2 - address1;
    // 215. g02
    MPI_Get_address(&TRparam[0].g02, &address2);
    disp[214] = address2 - address1;
    // 216. gA1
    MPI_Get_address(&TRparam[0].gA1, &address2);
    disp[215] = address2 - address1;
    // 217. gA2
    MPI_Get_address(&TRparam[0].gA2, &address2);
    disp[216] = address2 - address1;
    // 218. gS1
    MPI_Get_address(&TRparam[0].gS1, &address2);
    disp[217] = address2 - address1;
    // 219. gS2
    MPI_Get_address(&TRparam[0].gS2, &address2);
    disp[218] = address2 - address1;
    // 220. gX1
    MPI_Get_address(&TRparam[0].gX1, &address2);
    disp[219] = address2 - address1;
    // 221. gX2
    MPI_Get_address(&TRparam[0].gX2, &address2);
    disp[220] = address2 - address1;
    // 222. LOS_hist_Gfit_V0
    MPI_Get_address(&TRparam[0].LOS_hist_Gfit_V0, &address2);
    disp[221] = address2 - address1;
    // 223. LOS_hist_Gfit_sigma
    MPI_Get_address(&TRparam[0].LOS_hist_Gfit_sigma, &address2);
    disp[222] = address2 - address1;
    // 224. LOS_vel_hist_rbm
    MPI_Get_address(&TRparam[0].LOS_vel_hist_rbm, &address2);
    disp[223] = address2 - address1;
    // 225. LOS_vel_hist_std
    MPI_Get_address(&TRparam[0].LOS_vel_hist_std, &address2);
    disp[224] = address2 - address1;
    // 226. vlos_lower
    MPI_Get_address(&TRparam[0].vlos_lower, &address2);
    disp[225] = address2 - address1;
    // 227. vlos_upper
    MPI_Get_address(&TRparam[0].vlos_upper, &address2);
    disp[226] = address2 - address1;
    // 228. vlos_lower_limit
    MPI_Get_address(&TRparam[0].vlos_lower_limit, &address2);
    disp[227] = address2 - address1;
    // 229. vlos_upper_limit
    MPI_Get_address(&TRparam[0].vlos_upper_limit, &address2);
    disp[228] = address2 - address1;
    // 230. vel_geo_app_side
    MPI_Get_address(&TRparam[0].vel_geo_app_side, &address2);
    disp[229] = address2 - address1;
    // 231. vel_geo_rec_side
    MPI_Get_address(&TRparam[0].vel_geo_rec_side, &address2);
    disp[230] = address2 - address1;
    // 232. dispersion_VLOS
    MPI_Get_address(&TRparam[0].dispersion_VLOS, &address2);
    disp[231] = address2 - address1;
    // 233. sigma_factor_fix
    MPI_Get_address(&TRparam[0].sigma_factor_fix, &address2);
    disp[232] = address2 - address1;
    // 234. sigma_factor
    MPI_Get_address(&TRparam[0].sigma_factor, &address2);
    disp[233] = address2 - address1;
    // 235. sigma_factor_e
    MPI_Get_address(&TRparam[0].sigma_factor_e, &address2);
    disp[234] = address2 - address1;
    // 236. sigma_factor1
    MPI_Get_address(&TRparam[0].sigma_factor1, &address2);
    disp[235] = address2 - address1;
    // 237. sigma_factor2
    MPI_Get_address(&TRparam[0].sigma_factor2, &address2);
    disp[236] = address2 - address1;
    // 238. mean_mom4
    MPI_Get_address(&TRparam[0].mean_mom4, &address2);
    disp[237] = address2 - address1;
    // 239. std_mom4
    MPI_Get_address(&TRparam[0].std_mom4, &address2);
    disp[238] = address2 - address1;
    // 240. sigma_factor_mode
    MPI_Get_address(&TRparam[0].sigma_factor_mode, &address2);
    disp[239] = address2 - address1;
    // 241. sigma_factor_std
    MPI_Get_address(&TRparam[0].sigma_factor_std, &address2);
    disp[240] = address2 - address1;
    // 242. e_sigma_student_TR 
    MPI_Get_address(&TRparam[0].e_sigma_student_TR, &address2);
    disp[241] = address2 - address1;
    // 243. scale_factor_const_vlose_w
    MPI_Get_address(&TRparam[0].scale_factor_const_vlose_w, &address2);
    disp[242] = address2 - address1;
    // 244. scale_factor_var_vlose_w
    MPI_Get_address(&TRparam[0].scale_factor_var_vlose_w, &address2);
    disp[243] = address2 - address1;
    // 245. scale_factor_mom0weightedvf_e_to_mom2
    MPI_Get_address(&TRparam[0].scale_factor_mom0weightedvf_e_to_mom2, &address2);
    disp[244] = address2 - address1;
    // 246. hist_mean_vf_e_mom0_weighted_scaled
    MPI_Get_address(&TRparam[0].hist_mean_vf_e_mom0_weighted_scaled, &address2);
    disp[245] = address2 - address1;
    // 247. hist_std_vf_e_mom0_weighted_scaled
    MPI_Get_address(&TRparam[0].hist_std_vf_e_mom0_weighted_scaled, &address2);
    disp[246] = address2 - address1;
    // 248. einastofit_BIC
    MPI_Get_address(&TRparam[0].einastofit_BIC, &address2);
    disp[247] = address2 - address1;
    // 249. box_x
    MPI_Get_address(&TRparam[0].box_x, &address2);
    disp[248] = address2 - address1;
    // 250. box_y
    MPI_Get_address(&TRparam[0].box_y, &address2);
    disp[249] = address2 - address1;
    // 251. _n_fix
    MPI_Get_address(&TRparam[0]._n_fix, &address2);
    disp[250] = address2 - address1;
    // 252. _n
    MPI_Get_address(&TRparam[0]._n, &address2);
    disp[251] = address2 - address1;
    // 253. _ne
    MPI_Get_address(&TRparam[0]._ne, &address2);
    disp[252] = address2 - address1;
    // 254. _n1
    MPI_Get_address(&TRparam[0]._n1, &address2);
    disp[253] = address2 - address1;
    // 255. _n2
    MPI_Get_address(&TRparam[0]._n2, &address2);
    disp[254] = address2 - address1;
    // 256. r_2_fix
    MPI_Get_address(&TRparam[0].r_2_fix, &address2);
    disp[255] = address2 - address1;
    // 257. r_2
    MPI_Get_address(&TRparam[0].r_2, &address2);
    disp[256] = address2 - address1;
    // 258. r_2e
    MPI_Get_address(&TRparam[0].r_2e, &address2);
    disp[257] = address2 - address1;
    // 259. r_21
    MPI_Get_address(&TRparam[0].r_21, &address2);
    disp[258] = address2 - address1;
    // 260. r_22
    MPI_Get_address(&TRparam[0].r_22, &address2);
    disp[259] = address2 - address1;
    // 261. rho_2_fix
    MPI_Get_address(&TRparam[0].rho_2_fix, &address2);
    disp[260] = address2 - address1;
    // 262. rho_2
    MPI_Get_address(&TRparam[0].rho_2, &address2);
    disp[261] = address2 - address1;
    // 263. rho_2e
    MPI_Get_address(&TRparam[0].rho_2e, &address2);
    disp[262] = address2 - address1;
    // 264. rho_21
    MPI_Get_address(&TRparam[0].rho_21, &address2);
    disp[263] = address2 - address1;
    // 265. rho_22
    MPI_Get_address(&TRparam[0].rho_22, &address2);
    disp[264] = address2 - address1;
    // 266. _n_t
    MPI_Get_address(&TRparam[0]._n_t, &address2);
    disp[265] = address2 - address1;
    // 267. _r_2_t
    MPI_Get_address(&TRparam[0]._r_2_t, &address2);
    disp[266] = address2 - address1;
    // 268. _rho_2_t
    MPI_Get_address(&TRparam[0]._rho_2_t, &address2);
    disp[267] = address2 - address1;
    // 269. _ne_t
    MPI_Get_address(&TRparam[0]._ne_t, &address2);
    disp[268] = address2 - address1;
    // 270. _r_2e_t
    MPI_Get_address(&TRparam[0]._r_2e_t, &address2);
    disp[269] = address2 - address1;
    // 271. _rho_2e_t
    MPI_Get_address(&TRparam[0]._rho_2e_t, &address2);
    disp[270] = address2 - address1;
    // 272. einasto1D_ne
    MPI_Get_address(&TRparam[0].einasto1D_ne, &address2);
    disp[271] = address2 - address1;
    // 273. einasto1D_r2e
    MPI_Get_address(&TRparam[0].einasto1D_r2e, &address2);
    disp[272] = address2 - address1;
    // 274. einasto1D_rho2e
    MPI_Get_address(&TRparam[0].einasto1D_rho2e, &address2);
    disp[273] = address2 - address1;
    // 275. _n1_t
    MPI_Get_address(&TRparam[0]._n1_t, &address2);
    disp[274] = address2 - address1;
    // 276. _n2_t
    MPI_Get_address(&TRparam[0]._n2_t, &address2);
    disp[275] = address2 - address1;
    // 277. r_21_t
    MPI_Get_address(&TRparam[0].r_21_t, &address2);
    disp[276] = address2 - address1;
    // 278. r_22_t
    MPI_Get_address(&TRparam[0].r_22_t, &address2);
    disp[277] = address2 - address1;
    // 279. rho_21_t
    MPI_Get_address(&TRparam[0].rho_21_t, &address2);
    disp[278] = address2 - address1;
    // 280. rho_22_t
    MPI_Get_address(&TRparam[0].rho_22_t, &address2);
    disp[279] = address2 - address1;
    // 281. Ein_n_min
    MPI_Get_address(&TRparam[0].Ein_n_min, &address2);
    disp[280] = address2 - address1;
    // 282. Ein_n_max
    MPI_Get_address(&TRparam[0].Ein_n_max, &address2);
    disp[281] = address2 - address1;
    // 283. Ein_r_2_min
    MPI_Get_address(&TRparam[0].Ein_r_2_min, &address2);
    disp[282] = address2 - address1;
    // 284. Ein_r_2_max
    MPI_Get_address(&TRparam[0].Ein_r_2_max, &address2);
    disp[283] = address2 - address1;
    // 285. Ein_rho_2_min
    MPI_Get_address(&TRparam[0].Ein_rho_2_min, &address2);
    disp[284] = address2 - address1;
    // 286. Ein_rho_2_max
    MPI_Get_address(&TRparam[0].Ein_rho_2_max, &address2);
    disp[285] = address2 - address1;
    // 287. n_hist_post
    MPI_Get_address(&TRparam[0].n_hist_post, &address2);
    disp[286] = address2 - address1;
    // 288. hist_x
    MPI_Get_address(&TRparam[0].hist_x, &address2);
    disp[287] = address2 - address1;
    // 90. hist_y
    MPI_Get_address(&TRparam[0].hist_y, &address2);
    disp[288] = address2 - address1;
    // 290. hist_ye
    MPI_Get_address(&TRparam[0].hist_ye, &address2);
    disp[289] = address2 - address1;
    // 291. decimX_einasto_halofit_d 
    MPI_Get_address(&TRparam[0].decimX_einasto_halofit_d, &address2);
    disp[290] = address2 - address1;
    // 292. decimY_einasto_halofit_d 
    MPI_Get_address(&TRparam[0].decimY_einasto_halofit_d, &address2);
    disp[291] = address2 - address1;
    // 293. use_allPixels 
    MPI_Get_address(&TRparam[0].use_allPixels, &address2);
    disp[292] = address2 - address1;
    // 294. vrot_einasto_error 
    MPI_Get_address(&TRparam[0].vrot_einasto_error, &address2);
    disp[293] = address2 - address1;
    // 295. vrot_asymmetry_error 
    MPI_Get_address(&TRparam[0].vrot_asymmetry_error, &address2);
    disp[294] = address2 - address1;
    // 296. vrot_dispersion_error 
    MPI_Get_address(&TRparam[0].vrot_dispersion_error, &address2);
    disp[295] = address2 - address1;
    // 297. _nu_studenT 
    MPI_Get_address(&TRparam[0]._nu_studenT, &address2);
    disp[296] = address2 - address1;
    // 298. e_sigma_tr 
    MPI_Get_address(&TRparam[0].e_sigma_tr, &address2);
    disp[297] = address2 - address1;
    // 299. e_sigma_e_tr 
    MPI_Get_address(&TRparam[0].e_sigma_e_tr, &address2);
    disp[298] = address2 - address1;
    // 300. vf_e_user 
    MPI_Get_address(&TRparam[0].vf_e_user, &address2);
    disp[299] = address2 - address1;
    // 301. vdisp_user 
    MPI_Get_address(&TRparam[0].vdisp_user, &address2);
    disp[300] = address2 - address1;
    // 302. _nfilter 
    MPI_Get_address(&TRparam[0]._nfilter, &address2);
    disp[301] = address2 - address1;

    //MPI_Datatype TRparam_mpi;
    MPI_Type_create_struct(302, blocklen, disp, type, &TRparam_mpi);
    return TRparam_mpi;
    //MPI_Type_commit(&TRparam_mpi);
}

// --- End of line


