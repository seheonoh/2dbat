#ifndef __2DBAT_MPI_DATATYPE_H__
#define __2DBAT_MPI_DATATYPE_H__

// declare 2DBAT MPI dataype 
// MPI datatype related

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MPI_Datatype type[302] = {/* Process + einasto output file + Parameters for multinest */
                        MPI_CHAR, MPI_CHAR, MPI_INT, MPI_INT, // 1:char wdir, 2:char txtfile_multinest_output, 3:int loop_check, 4:int n_freeParams
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 5:double e_sigma, 6:double e_sigma_fitted, 7:double e_sigma_fitted_t
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 8:double maxLogLike, 9:double logZ, 10:double logZerr
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 11:double maxLogLikeF, 12:double logZF, 13:double logZerrF
                        // VF fields
                        MPI_INT, MPI_INT, MPI_INT, MPI_INT, // 14:int nax1, 15:int nax2, 16:int decim_x0, 17:int decim_y0
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 18:double decimX, 19:double decimY, 20:double decimX_einasto_halofit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 21:double decimY_einasto_halofit, 22:double decimX_trfit, 23:double decimY_trfit
                        MPI_FLOAT, // 24:float pixelScale
                        // simple ellipse fit parametres
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 25:double N_all_pixels_in_a_ring, 25:double ellipse_xpos_boxfiltered, 26:double ellipse_ypos_boxfiltered
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 28:double ellipse_pa_boxfiltered, 29:double pa_EllipseFit_e, 30:double ellipse_incl_boxfiltered
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 31:double incl_EllipseFit_e, 32:double ellipse_semi_mx_boxfiltered, 33:double xpos1_from_ellipsefit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 34:double xpos2_from_ellipsefit, 35:double ypos1_from_ellipsefit, 36:double ypos2_from_ellipsefit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 37:double vsys1_from_vlosfit, 38:double vsys2_from_vlosfit, 39:double perimeter
                        // set rings
                        MPI_CHAR, MPI_CHAR, MPI_INT,        // 40:double final_fit, 41:double fullFit, 42:double Nrings
                        MPI_INT, MPI_INT, MPI_INT,          // 43:double Nrings_intp, 44:double Nrings_to_semi_mx, 45:double N_reliable_rings
                        MPI_INT, MPI_DOUBLE, MPI_DOUBLE,    // 46:double tilted_ring, 47:double Npoints_in_tilted_ring, 48:double total_Npoints_allRings
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 49:double ring_radius, 50:double ring_intp, 51:double ring_s
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 52:double ring_e, 53:double ring_w, 54:double ring_s_for_einasto1Dfit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 56:double ring_e_for_einasto1Dfit, 57:double ring_w_for_einasto1Dfit, 58:double Npoints_in_tilted_ring_decim0
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 58:double npoints_inaring, 59:double npoints_inaring_decim0, 60:double rGalaxyPlane_pixel
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 61:double rGalaxyPlane_pixel_max, 62:double free_angle, 63:double wpow
                        MPI_DOUBLE,                         // 64:double rwpow
                        // xpos
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 65:char xpos_fix, 66:double xposF, 66:double xposF_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 68:double xpos0, 69:double xpos, 70:double xpos_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 71:double xpos_temp, 72:double xpos_mode, 73:double xpos_std
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 74:double xpos1, 75:double xpos2, 76:double xposF_EinastoFit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 77:double xposF_EinastoFit_e, 78:double xposF_EinastoFit_t, 79:double xposF_EinastoFit_e_t
                        // ypos
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 80:char ypos_fix, 81:double yposF, 82:double yposF_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 83:double ypos0, 84:double ypos, 85:double ypos_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 86:double ypos_temp, 87:double ypos_mode, 88:double ypos_std
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 89:double ypos1, 90:double ypos2, 91:double yposF_EinastoFit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 92:double yposF_EinastoFit_e, 93:double yposF_EinastoFit_t, 94:double yposF_EinastoFit_e_t
                        // vsys
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 95:char vsys_fix, 96:double vsysF, 97:double vsysF_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 98:double vsys0, 99:double vsys, 100:double vsys_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 101:double vsys_temp, 102:double vsys_mode, 103:double vsys_std
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 104:double vsys1, 105:double vsys2, 106:double vsysF_EinastoFit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 107:double vsysF_EinastoFit_e, 108:double vsysF_EinastoFit_t, 109:double vsysF_EinastoFit_e_t
                        // PA
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 110:char pa_fix, 96:double paF, 97:double paF_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 113:double pa0, 114:double pa, 115:double pa_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 116:double pa_temp, 117:double pa_temp_e, 118:double pa1
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 119:double pa2, 120:double pa1_for_TRfit, 121:double pa2_for_TRfit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 122:double _p1_tr, 123:double _p2_tr, 124:double PA_MAX_in_degree
                        MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR,   // 125:double paF_EinastoFit, 126:double paF_EinastoFit_e, 127:char pa_function
                        MPI_INT, MPI_INT, MPI_INT,          // 128:int n_coeffs_bspline_pa, 129:int pa_nbreak_bspline, 130:int pa_order_bspline
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 131:double _p_bs, 132:double _p_bs_e, 133:double _p_bs_t
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 134:double _p_bs_e_t, 135:double bspline1pa, 136:double bspline2pa
                        MPI_DOUBLE, MPI_DOUBLE,             // 137:double _p_bs_tr, 138:double _bspline_pa_hist_sigma
                        // INCL
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 139:char incl_fix, 140:double inclF, 141:double inclF_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 142:double incl0, 143:double incl, 144:double incl_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 145:double incl_temp, 147:double incl_temp_e, 147:double incl1
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 148:double incl2, 149:double incl1_for_TRfit, 150:double incl2_for_TRfit
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 151:double _i1_tr, 152:double _i2_tr, 153:double INCL_MAX_in_degree
                        MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR,   // 154:double inclF_EinastoFit, 155:double inclF_EinastoFit_e, 156:char incl_function
                        MPI_INT, MPI_INT, MPI_INT,          // 157:int n_coeffs_bspline_incl, 158:int incl_nbreak_bspline, 159:int incl_order_bspline
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 160:double _i_bs, 161:double _i_bs_e, 162:double _i_bs_t
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 163:double _i_bs_e_t, 164:double bspline1incl, 165:double bspline2incl
                        MPI_DOUBLE, MPI_DOUBLE,             // 166:double _i_bs_tr, 167:double _bspline_incl_hist_sigma
                        // VROT
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 168:char vrot_fix, 169:double vrotF, 170:double vrotF_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 171:double vrot1, 172:double vrot2, 173:double vrot0
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 174:double vrot, 175:double vrot_e, 176:double vrot_rec
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 177:double vrot_e_rec, 178:double vrot_app, 179:double vrot_e_app
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 180:double vrot_temp, 181:double vrot_temp_e, 182:double vrot_intp
                        MPI_DOUBLE,                         // 183:double vrot_e_intp
                        // VRAD
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 184:char vrad_fix, 185:double vrad1, 186:double vrad2
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 187:double vradF, 188:double vradF_e, 189:double vrad_max
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 190:double vrad, 191:double vrad_e, 192:double vrad_rec
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 193:double vrad_rec_e, 194:double vrad_app, 195:double vrad_app_e
                        MPI_CHAR, MPI_INT, MPI_INT,         // 196:char vrad_function, 197:int n_coeffs_bspline_vrad, 198:int vrad_nbreak_bspline
                        MPI_INT, MPI_DOUBLE, MPI_DOUBLE,    // 199:int vrad_order_bspline, 200:double bspline1vrad, 201:double bspline2vrad
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 202:double _vr_bs, 203:double _vr_bs_e, 204:double _vr_bs_t
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 205:double _vr_bs_e_t, 206:double _vr_bs_tr, 207:double vrad_temp
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 208:double vrad_temp_e, 209:double vrad_einastofit_bs, 210:double vrad_einastofit_bs_e
                        MPI_DOUBLE,                         // 211:double _bspline_vrad_hist_sigma
                        // parameters for Gaussian profile fit
                        MPI_INT, MPI_DOUBLE, MPI_DOUBLE,    // 212:int n_gauss, 213:double g_param, 214:double g01
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 215:double g02, 216:double gA1, 217:double gA2
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 218:double gS1, 219:double gS2, 220:double gX1
                        MPI_DOUBLE,                         // 221:double gX2
                        // line-of-sight velocities
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 222:double LOS_hist_Gfit_V0, 223:double LOS_hist_Gfit_sigma, 224:double LOS_vel_hist_rbm
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 225:double LOS_vel_hist_std, 226:double vlos_lower, 227:double vlos_upper
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 228:double vlos_lower_limit, 229:double vlos_upper_limit, 230:double vel_geo_app_side
                        MPI_DOUBLE, MPI_DOUBLE,             // 231:double vel_geo_rec_side, 232:double dispersion_VLOS
                        // sigma_factor
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 233:char sigma_factor_fix, 234:double sigma_factor, 235:double sigma_factor_e
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 236:double sigma_factor1, 237:double sigma_factor2, 238:double mean_mom4
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 239:double std_mom4, 240:double sigma_factor_mode, 241:double sigma_factor_std
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 242:double sigmaFactor_TR, 243:double scale_factor_const_vlose_w, 244:double scale_factor_var_vlose_w
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 245:double scale_factor_mom0weightedvf_e_to_mom2, 246:double hist_mean_vf_e_mom0_weighted_scaled, 247:double hist_std_vf_e_mom0_weighted_scaled
                        MPI_DOUBLE,                         // 248:double einastofit_BIC
                        // median box filter
                        MPI_INT, MPI_INT,                   // 249:int box_x, 250:int box_y
                        // Einasto profile
                        // n
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 251:char _n_fix, 252:double _n, 253:double _ne
                        MPI_DOUBLE, MPI_DOUBLE,             // 254:double _n1, 255:double _n2
                        // r_2
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 256:char r_2_fix, 257:double r_2, 258:double r_2e
                        MPI_DOUBLE, MPI_DOUBLE,             // 259:double r_21, 260:double r_22
                        // rho_2
                        MPI_CHAR, MPI_DOUBLE, MPI_DOUBLE,   // 261:char rho_2_fix, 262:double rho_2, 263:double rho_2e
                        MPI_DOUBLE, MPI_DOUBLE,             // 264:double rho_21, 265:double rho_22
                        // 
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 266:double _n_t, 267:double _r_2_t, 268:double _rho_2_t
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 269:double _ne_t, 270:double _r_2e_t, 271:double _rho_2e_t
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 272:double einasto1D_ne, 273:double einasto1D_r2e, 274:double einasto1D_rho2e
                        // Einasto profile : temporary
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 275:double _n1_t, 276:double _n2_t, 277:double r_21_t
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 278:double r_22_t, 279:double rho_21_t, 280:double rho_22_t
                        // Einasto halo params limits
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 281:double Ein_n_min, 282:double Ein_n_max, 283:double Ein_r_2_min
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 284:double Ein_r_2_max, 285:double Ein_rho_2_min, 286:double Ein_rho_2_max
                        // histogram
                        MPI_INT, MPI_DOUBLE, MPI_DOUBLE,    // 287:int n_hist_post, 288:double hist_x, 289:double hist_y
                        MPI_DOUBLE,                         // 290:double hist_ye
                        MPI_DOUBLE, MPI_DOUBLE, // 291:double decimX_einasto_halofit, 292:double decimY_einasto_halofit
                        MPI_INT, // 293: int use_allPixels
                        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 294:double Ein_r_2_max, 295:double Ein_rho_2_min, 296:double Ein_rho_2_max
                        MPI_DOUBLE, // 297:double _nu_studentT
                        MPI_DOUBLE, MPI_DOUBLE, // 298:double e_sigma_tr, 299:double e_sigma_e_tr
                        MPI_DOUBLE, MPI_DOUBLE, // 300:double vf_e_user, 301:double vdisp_user
                        MPI_INT // 302:int _nfilter
                        };

int blocklen[302] = {// Process + einasto output file + Parameters for multinest
                    1000, 1000, 1, 1, // 1:char wdir, 2:char txtfile_multinest_output, 3:int loop_check, 4:int n_freeParams
                    1, 1, 1, // 5:double e_sigma, 6:double e_sigma_fitted, 7:double e_sigma_fitted_t
                    999, 999, 999, // 8:double maxLogLike, 9:double logZ, 10:double logZerr
                    1, 1, 1, // 11:double maxLogLikeF, 12:double logZF, 13:double logZerrF
                    // VF fields
                    1, 1, 1, 1, // 14:int nax1, 15:int nax2, 16:int decim_x0, 17:int decim_y0
                    1, 1, 1, // 18:double decimX, 19:double decimY, 20:double decimX_einasto_halofit
                    1, 1, 1, // 21:double decimY_einasto_halofit, 22:double decimX_trfit, 23:double decimY_trfit
                    1, // 24:float pixelScale
                    // simple ellipse fit parametres
                    1, 1, 1, // 25:double N_all_pixels_in_a_ring, 25:double ellipse_xpos_boxfiltered, 26:double ellipse_ypos_boxfiltered
                    1, 1, 1, // 28:double ellipse_pa_boxfiltered, 29:double pa_EllipseFit_e, 30:double ellipse_incl_boxfiltered
                    1, 1, 1, // 31:double incl_EllipseFit_e, 32:double ellipse_semi_mx_boxfiltered, 33:double xpos1_from_ellipsefit
                    1, 1, 1, // 34:double xpos2_from_ellipsefit, 35:double ypos1_from_ellipsefit, 36:double ypos2_from_ellipsefit
                    1, 1, 1, // 37:double vsys1_from_vlosfit, 38:double vsys2_from_vlosfit, 39:double perimeter
                    // set rings
                    1, 1, 1,        // 40:double final_fit, 41:double fullFit, 42:double Nrings
                    1, 1, 1,          // 43:double Nrings_intp, 44:double Nrings_to_semi_mx, 45:double N_reliable_rings
                    4024*4024*2, 1, 1,    // 46:double tilted_ring, 47:double Npoints_in_tilted_ring, 48:double total_Npoints_allRings
                    999, 999, 1, // 49:double ring_radius, 50:double ring_intp, 51:double ring_s
                    1, 1, 1, // 52:double ring_e, 53:double ring_w, 54:double ring_s_for_einasto1Dfit
                    1, 1, 1, // 56:double ring_e_for_einasto1Dfit, 57:double ring_w_for_einasto1Dfit, 58:double Npoints_in_tilted_ring_decim0
                    999999, 999999, 1, // 58:double npoints_inaring, 59:double npoints_inaring_decim0, 60:double rGalaxyPlane_pixel
                    1, 1, 1, // 61:double rGalaxyPlane_pixel_max, 62:double free_angle, 63:double wpow
                    1,                         // 64:double rwpow
                    // xpos
                    1, 1, 1,   // 65:char xpos_fix, 66:double xposF, 66:double xposF_e
                    999, 999, 999, // 68:double xpos0, 69:double xpos, 70:double xpos_e
                    999, 1, 1, // 71:double xpos_temp, 72:double xpos_mode, 73:double xpos_std
                    1, 1, 1, // 74:double xpos1, 75:double xpos2, 76:double xposF_EinastoFit
                    1, 1, 1, // 77:double xposF_EinastoFit_e, 78:double xposF_EinastoFit_t, 79:double xposF_EinastoFit_e_t
                    // ypos
                    1, 1, 1,   // 80:char ypos_fix, 81:double yposF, 82:double yposF_e
                    999, 999, 999, // 83:double ypos0, 84:double ypos, 85:double ypos_e
                    999, 1, 1, // 86:double ypos_temp, 87:double ypos_mode, 88:double ypos_std
                    1, 1, 1, // 89:double ypos1, 90:double ypos2, 91:double yposF_EinastoFit
                    1, 1, 1, // 92:double yposF_EinastoFit_e, 93:double yposF_EinastoFit_t, 94:double yposF_EinastoFit_e_t
                    // vsys
                    1, 1, 1,   // 95:char vsys_fix, 96:double vsysF, 97:double vsysF_e
                    999, 999, 999, // 98:double vsys0, 99:double vsys, 100:double vsys_e
                    999, 1, 1, // 101:double vsys_temp, 102:double vsys_mode, 103:double vsys_std
                    1, 1, 1, // 104:double vsys1, 105:double vsys2, 106:double vsysF_EinastoFit
                    1, 1, 1, // 107:double vsysF_EinastoFit_e, 108:double vsysF_EinastoFit_t, 109:double vsysF_EinastoFit_e_t
                    // PA
                    1, 1, 1,   // 110:char pa_fix, 96:double paF, 97:double paF_e
                    999, 999, 999, // 113:double pa0, 114:double pa, 115:double pa_e
                    999, 999, 1, // 116:double pa_temp, 117:double pa_temp_e, 118:double pa1
                    1, 1, 1, // 119:double pa2, 120:double pa1_for_TRfit, 121:double pa2_for_TRfit
                    999, 999, 1, // 122:double _p1_tr, 123:double _p2_tr, 124:double PA_MAX_in_degree
                    1, 1, 13,   // 125:double paF_EinastoFit, 126:double paF_EinastoFit_e, 127:char pa_function
                    1, 1, 1,          // 128:int n_coeffs_bspline_pa, 129:int pa_nbreak_bspline, 130:int pa_order_bspline
                    99, 99, 99, // 131:double _p_bs, 132:double _p_bs_e, 133:double _p_bs_t
                    99, 99, 99, // 134:double _p_bs_e_t, 135:double bspline1pa, 136:double bspline2pa
                    99, 1,             // 137:double _p_bs_tr, 138:double _bspline_pa_hist_sigma
                    // INCL
                    1, 1, 1,   // 139:char incl_fix, 140:double inclF, 141:double inclF_e
                    999, 999, 999, // 142:double incl0, 143:double incl, 144:double incl_e
                    999, 999, 1, // 145:double incl_temp, 147:double incl_temp_e, 147:double incl1
                    1, 1, 1, // 148:double incl2, 149:double incl1_for_TRfit, 150:double incl2_for_TRfit
                    999, 999, 1, // 151:double _i1_tr, 152:double _i2_tr, 153:double INCL_MAX_in_degree
                    1, 1, 13,   // 154:double inclF_EinastoFit, 155:double inclF_EinastoFit_e, 156:char incl_function
                    1, 1, 1,          // 157:int n_coeffs_bspline_incl, 158:int incl_nbreak_bspline, 159:int incl_order_bspline
                    99, 99, 99, // 160:double _i_bs, 161:double _i_bs_e, 162:double _i_bs_t
                    99, 99, 99, // 163:double _i_bs_e_t, 164:double bspline1incl, 165:double bspline2incl
                    99, 1,             // 166:double _i_bs_tr, 167:double _bspline_incl_hist_sigma
                    // VROT
                    1, 1, 1,   // 168:char vrot_fix, 169:double vrotF, 170:double vrotF_e
                    1, 1, 999, // 171:double vrot1, 172:double vrot2, 173:double vrot0
                    999, 999, 999, // 174:double vrot, 175:double vrot_e, 176:double vrot_rec
                    999, 999, 999, // 177:double vrot_e_rec, 178:double vrot_app, 179:double vrot_e_app
                    999, 999, 999, // 180:double vrot_temp, 181:double vrot_temp_e, 182:double vrot_intp
                    999,                         // 183:double vrot_e_intp
                    // VRAD
                    1, 1, 1,   // 184:char vrad_fix, 185:double vrad1, 186:double vrad2
                    1, 1, 1, // 187:double vradF, 188:double vradF_e, 189:double vrad_max
                    999, 999, 999, // 190:double vrad, 191:double vrad_e, 192:double vrad_rec
                    999, 999, 999, // 193:double vrad_rec_e, 194:double vrad_app, 195:double vrad_app_e
                    13, 1, 1,         // 196:char vrad_function, 197:int n_coeffs_bspline_vrad, 198:int vrad_nbreak_bspline
                    1, 99, 99,    // 199:int vrad_order_bspline, 200:double bspline1vrad, 201:double bspline2vrad
                    99, 99, 99, // 202:double _vr_bs, 203:double _vr_bs_e, 204:double _vr_bs_t
                    99, 99, 999, // 205:double _vr_bs_e_t, 206:double _vr_bs_tr, 207:double vrad_temp
                    999, 999, 999, // 208:double vrad_temp_e, 209:double vrad_einastofit_bs, 210:double vrad_einastofit_bs_e 
                    1,                         // 211:double _bspline_vrad_hist_sigma
                    // parameters for Gaussian profile fit
                    1, 99, 1,    // 212:int n_gauss, 213:double g_param, 214:double g01
                    1, 1, 1, // 215:double g02, 216:double gA1, 217:double gA2
                    1, 1, 1, // 218:double gS1, 219:double gS2, 220:double gX1
                    1,                         // 221:double gX2
                    // line-of-sight velocities
                    1, 1, 1, // 222:double LOS_hist_Gfit_V0, 223:double LOS_hist_Gfit_sigma, 224:double LOS_vel_hist_rbm
                    1, 1, 1, // 225:double LOS_vel_hist_std, 226:double vlos_lower, 227:double vlos_upper
                    1, 1, 1, // 228:double vlos_lower_limit, 229:double vlos_upper_limit, 230:double vel_geo_app_side
                    1, 1,             // 231:double vel_geo_rec_side, 232:double dispersion_VLOS
                    // sigma_factor
                    1, 1, 1,   // 233:char sigma_factor_fix, 234:double sigma_factor, 235:double sigma_factor_e
                    1, 1, 1, // 236:double sigma_factor1, 237:double sigma_factor2, 238:double mean_mom4
                    1, 1, 1, // 239:double std_mom4, 240:double sigma_factor_mode, 241:double sigma_factor_std
                    999, 1, 1, // 242:double sigmaFactor_TR, 243:double scale_factor_const_vlose_w, 244:double scale_factor_var_vlose_w
                    1, 1, 1, // 245:double scale_factor_mom0weightedvf_e_to_mom2, 246:double hist_mean_vf_e_mom0_weighted_scaled, 247:double hist_std_vf_e_mom0_weighted_scaled
                    1,                         // 248:double einastofit_BIC
                    // median box filter
                    1, 1,                   // 249:int box_x, 250:int box_y
                    // Einasto profile
                    // n
                    1, 1, 1,   // 251:char _n_fix, 252:double _n, 253:double _ne
                    1, 1,             // 254:double _n1, 255:double _n2
                    // r_2
                    1, 1, 1,   // 256:char r_2_fix, 257:double r_2, 258:double r_2e
                    1, 1,             // 259:double r_21, 260:double r_22
                    // rho_2
                    1, 1, 1,   // 261:char rho_2_fix, 262:double rho_2, 263:double rho_2e
                    1, 1,             // 264:double rho_21, 265:double rho_22
                    // 
                    1, 1, 1, // 266:double _n_t, 267:double _r_2_t, 268:double _rho_2_t
                    1, 1, 1, // 269:double _ne_t, 270:double _r_2e_t, 271:double _rho_2e_t
                    1, 1, 1, // 272:double einasto1D_ne, 273:double einasto1D_r2e, 274:double einasto1D_rho2e
                    // Einasto profile : temporary
                    1, 1, 1, // 275:double _n1_t, 276:double _n2_t, 277:double r_21_t
                    1, 1, 1, // 278:double r_22_t, 279:double rho_21_t, 280:double rho_22_t
                    // Einasto halo params limits
                    1, 1, 1, // 281:double Ein_n_min, 282:double Ein_n_max, 283:double Ein_r_2_min
                    1, 1, 1, // 284:double Ein_r_2_max, 285:double Ein_rho_2_min, 286:double Ein_rho_2_max
                    // histogram
                    1, 9999, 9999,    // 287:int n_hist_post, 288:double hist_x, 289:double hist_y
                    9999,                        // 290:double hist_ye
                    1, 1, // 291:double decimX_einasto_halofit, 292:double decimY_einasto_halofit
                    1, // 293:int use_allPixels
                    999, 999, 999, // 294:double vrot_einasto_error, 295:double vrot_asymmetry_error, 296:double vrot_dispersion_error
                    1, // 297:double _nu_studentT
                    1, 1, // 298:double e_sigma_tr, 299:double e_sigma_e_tr
                    1, 1, // 300:double vf_e_user, 301:double vdisp_user
                    1 // 302:int _nfilter
                };

// --- End of line

#endif

