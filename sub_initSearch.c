#include "cap.h"
/*
  Function to search uniform moment tensors
       initSearchMT
 
  Function initSearchMT calls functions
       get_Rand_MT
       get_Grid_MT
       search_MT
 
  See uniformMT.m for examples to generate uniform moment tensors in matlab.
  See also "A uniform parameterization of moment tensors" (Tape & Tape, 2015)
 
  20160112 celso alvizuri - cralvizuri@alaska.edu 
 
 */

SOLN initSearchMT( int npar, // 3=mw; 2=iso; 1=clvd; 0=strike/dip/rake
        int nda,
        DATA *obs0,
        int nfm,
        FM *fm,
        float fm_thr,
        const int *max_shft,
        float tie,
        MTPAR *mt,
        GRID grid,
        int interp,
        int bootstrap,
        int search_type,
        int norm,
        SEARCHPAR * searchPar,
        ARRAYMT * arrayMT)
{
    int i, j, k, l, m, k1, kc, z0, z1, z2, ii, N, iso_len;
    float mw_ran; // half-range for magnitude search (previously int)
    int i_stk, i_dip, i_rak, i_iso;
    int perc, del_N, count_perc;
    float amp, rad[6], arad[4][3], x, x1, x2, y, y1, y2, cfg[NCP], s3d[9], temp[3], m_par, del_dip, del_iso;
    float *f_pt0, *f_pt1, *r_pt, *r_pt0, *r_pt1, *z_pt, *z_pt0, *z_pt1, *grd_err, *rnd_stk, *rnd_dip, *rnd_rak, *rnd_iso, *rnd_clvd, *iso;
    float dx, mtensor[3][3], *r_iso, *z_iso;
    DATA *obs;
    COMP *spt;
    SOLN sol, best_sol;

    switch(search_type)
    {
        case 1: 
            fprintf(stderr,"Preparing array of moment tensors using option %d: GRID\n", search_type);
            fprintf(stderr,"Solutions to prepare: %10d\n\n", searchPar->nsol);
            getGridMT(searchPar, arrayMT);
            break;
        case 2:
            fprintf(stderr,"Preparing array of moment tensors using option %d: RANDOM\n", search_type);
            fprintf(stderr,"Solutions to prepare: %10d\n\n", searchPar->nsol);
            getRandMT(searchPar, arrayMT);
            break;
        default:
            fprintf(stderr,"\nStop. Search type \"%d\" not available.\n", search_type);
            exit(-1);
            break;
    }

    // run FMT search
    best_sol = searchMT(npar, nda, obs0, nfm, fm, fm_thr, max_shft, tie, mt, grid, interp, bootstrap, search_type, norm, searchPar, arrayMT);
    return(best_sol);

} /* end function initSearchMT */

/* fill array with RANDOM moment tensors */
void getRandMT(SEARCHPAR * searchPar, ARRAYMT * arrayMT)
{
    int isol;
    float siso1, siso2;
    float h1, h2;

    /* temporary vectors for filling the moment tensor array */
    float * vec_v     = calloc(searchPar->nsol, sizeof(float));
    float * vec_gamma = calloc(searchPar->nsol, sizeof(float));
    float * vec_u     = calloc(searchPar->nsol, sizeof(float));
    float * vec_beta  = calloc(searchPar->nsol, sizeof(float));
    float * vec_delta = calloc(searchPar->nsol, sizeof(float));
    float * vec_kappa = calloc(searchPar->nsol, sizeof(float));
    float * vec_theta = calloc(searchPar->nsol, sizeof(float));
    float * vec_h     = calloc(searchPar->nsol, sizeof(float));
    float * vec_sigma = calloc(searchPar->nsol, sizeof(float));

    // convert w to u
    searchPar->u1 = w2u(searchPar->w1);
    searchPar->u2 = w2u(searchPar->w2);

    if(LUNE_GRID_INSTEAD_OF_UV == 1) {
        // This section runs only when compiling with flag LUNE_GRID_INSTEAD_OF_UV=1. As requested.
        // This aproach is to not have lune grid readily available. As requested.
        // The original version used flag -K3 from cap.pl. This option was deleted, as requested.

        fprintf(stderr,"\nWARNING. using lune grid.\n");
        fprintf(stderr,"gamma [ %10.5f, %10.5f]\n", searchPar->v1, searchPar->v2);
        fprintf(stderr,"delta [ %10.5f, %10.5f]\n", searchPar->w1, searchPar->w2);

        // gamma
        randvec(searchPar->v1 * d2r, searchPar->v2 * d2r, searchPar->nv, vec_gamma);

        // delta. siso = sin(iso)
        siso1 = sin(searchPar->w1 * d2r); siso2 = sin(searchPar->w2 * d2r);
        randvec(siso1, siso2, searchPar->nu, vec_beta);     // nu is calculated in cap.pl
                                                            // vec_beta is a dummy. no relation to beta here.
        siso2delta_vec(vec_beta, searchPar->nu, vec_delta);

        // strike
        randvec(searchPar->k1 * d2r, searchPar->k2 * d2r, searchPar->nk, vec_kappa);

        // dip
        h1 = cos(searchPar->h1 * d2r); h2 = cos(searchPar->h2 * d2r);
        randvec(h1, h2, searchPar->nh, vec_h);    // nh is calculated in cap.pl
        h2dip_vec(vec_h, searchPar->nh, vec_theta); 

        // rake
        randvec(searchPar->s1 * d2r, searchPar->s2 * d2r, searchPar->ns, vec_sigma);
    }
    else {
        // uniform grid. (default case)
        randvec(searchPar->v1, searchPar->v2, searchPar->nv, vec_v);    // gamma(v)
        randvec(searchPar->u1, searchPar->u2, searchPar->nu, vec_u);    // beta(u)
        randvec(searchPar->k1, searchPar->k2, searchPar->nk, vec_kappa);
        randvec(searchPar->h1, searchPar->h2, searchPar->nh, vec_h);    // dip(h)
        randvec(searchPar->s1, searchPar->s2, searchPar->ns, vec_sigma);
 
        u2beta_vec(vec_u, vec_beta, searchPar->nu);
        beta2delta_vec(vec_beta, vec_delta, searchPar->nu);
        v2gamma_vec(vec_v, searchPar->nv, vec_gamma);      
        h2dip_vec(vec_h, searchPar->nh, vec_theta);          
    }

    fprintf(stderr,"~~~~~~~~~~~input u1 %f \n", searchPar->u1);
    fprintf(stderr,"~~~~~~~~~~~DEBUG u[0] %f \n", vec_u[0]);
    fprintf(stderr,"~~~~~~~~~~~DEBUG beta[0] %f \n", vec_beta[0]);

    free(vec_v);
    free(vec_u);
    free(vec_beta);
    free(vec_h);

    /* fill the moment tensor array with moment tensors */
    fprintf(stderr,"\nFilling moment tensor table (nsol= %10d) ... ", searchPar->nsol);
    for(isol = 0; isol < searchPar->nsol; isol++) {
        arrayMT[isol].gamma = vec_gamma[isol];
        arrayMT[isol].delta = vec_delta[isol];
        arrayMT[isol].kappa = vec_kappa[isol];
        arrayMT[isol].theta = vec_theta[isol];
        arrayMT[isol].sigma = vec_sigma[isol];
//        fprintf(stdout,"index = %20d %f %f %f %f %f \n", isol, arrayMT[isol].gamma *r2d, arrayMT[isol].delta *r2d,  arrayMT[isol].kappa *r2d,  arrayMT[isol].theta *r2d,  arrayMT[isol].sigma *r2d);
    }
    fprintf(stderr,"done. nsol = %d \n", isol);
    if (isol == 1) {
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 
                0, arrayMT[0].gamma *r2d, arrayMT[0].delta *r2d,  arrayMT[0].kappa *r2d,  arrayMT[0].theta *r2d,  arrayMT[0].sigma *r2d); 
    } else if(isol > 1) {
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 
                0, arrayMT[0].gamma *r2d, arrayMT[0].delta *r2d,  arrayMT[0].kappa *r2d,  arrayMT[0].theta *r2d,  arrayMT[0].sigma *r2d); 
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n",
                isol-1, arrayMT[isol-1].gamma *r2d, arrayMT[isol-1].delta *r2d,  arrayMT[isol-1].kappa *r2d,  arrayMT[isol-1].theta *r2d,  arrayMT[isol-1].sigma *r2d);
    } 

    free(vec_gamma);
    free(vec_delta);
    free(vec_kappa);
    free(vec_theta);
    free(vec_sigma);
}   /* end function getRandMT */

/* fill array with a grid uniform moment tensors */
void getGridMT(SEARCHPAR * searchPar, ARRAYMT * arrayMT)
{
    int na, isol;
    int ig, id, ik, ih, is;

    float siso1, siso2;
    float h1, h2;

    /* temporary vectors for filling the moment tensor array */
    float * vec_v     = calloc(searchPar->nv, sizeof(float));
    float * vec_gamma = calloc(searchPar->nv, sizeof(float));
    float * vec_u     = calloc(searchPar->nu, sizeof(float));
    float * vec_beta  = calloc(searchPar->nu, sizeof(float));
    float * vec_delta = calloc(searchPar->nu, sizeof(float));
    float * vec_kappa = calloc(searchPar->nk, sizeof(float));       /* strike */
    float * vec_theta = calloc(searchPar->nh, sizeof(float));
    float * vec_h     = calloc(searchPar->nh, sizeof(float));
    float * vec_sigma = calloc(searchPar->ns, sizeof(float));       /* slip */

    // convert w to u
    searchPar->u1 = w2u(searchPar->w1);
    searchPar->u2 = w2u(searchPar->w2);

    if(LUNE_GRID_INSTEAD_OF_UV == 1) {
        // This section runs only when compiling with flag LUNE_GRID_INSTEAD_OF_UV=1. As requested.
        // This aproach is to not have lune grid readily available. As requested.
        // The original version used flag -K3 from cap.pl. This option was deleted, as requested.

        fprintf(stderr,"\nWARNING. using lune grid.\n");
        fprintf(stderr,"gamma [ %10.5f, %10.5f]\n", searchPar->v1, searchPar->v2);
        fprintf(stderr,"delta [ %10.5f, %10.5f]\n", searchPar->w1, searchPar->w2);

        // gamma
        gridvec(searchPar->v1 * d2r, searchPar->v2 * d2r, searchPar->nv, vec_gamma);

        // delta. siso = sin(iso)
        siso1 = sin(searchPar->w1 * d2r); siso2 = sin(searchPar->w2 * d2r);
        gridvec(siso1, siso2, searchPar->nu, vec_beta);     // nu is calculated in cap.pl
                                                            // vec_beta is a dummy. no relation to beta here.
        siso2delta_vec(vec_beta, searchPar->nu, vec_delta);

        // strike
        gridvec(searchPar->k1 * d2r, searchPar->k2 * d2r, searchPar->nk, vec_kappa);

        // dip
        h1 = cos(searchPar->h1 * d2r); h2 = cos(searchPar->h2 * d2r);
        gridvec(h1, h2, searchPar->nh, vec_h);    // nh is calculated in cap.pl
        h2dip_vec(vec_h, searchPar->nh, vec_theta); 

        // rake
        gridvec(searchPar->s1 * d2r, searchPar->s2 * d2r, searchPar->ns, vec_sigma);
    }
    else {
        // uniform grid. (default case)
        gridvec(searchPar->v1, searchPar->v2, searchPar->nv, vec_v);    // gamma(v)
        gridvec(searchPar->u1, searchPar->u2, searchPar->nu, vec_u);    // beta(u)
        gridvec(searchPar->k1, searchPar->k2, searchPar->nk, vec_kappa);
        gridvec(searchPar->h1, searchPar->h2, searchPar->nh, vec_h);    // dip(h)
        gridvec(searchPar->s1, searchPar->s2, searchPar->ns, vec_sigma);
 
        u2beta_vec(vec_u, vec_beta, searchPar->nu);
        beta2delta_vec(vec_beta, vec_delta, searchPar->nu);
        v2gamma_vec(vec_v, searchPar->nv, vec_gamma);      
        h2dip_vec(vec_h, searchPar->nh, vec_theta);          
    }

    /* fill the moment tensor array with moment tensors */
    fprintf(stderr,"\nFilling moment tensor table (nsol= %10d) ... ", searchPar->nsol);
    isol = 0;
    if(searchPar->nv == 0) {
        fprintf(stderr,"WARNING setting nv=1 to fill up arrayMT \n");
        searchPar->nv = 1;
    }
    if(searchPar->nu == 0) {
        fprintf(stderr,"WARNING setting nu=1 to fill up arrayMT \n");
        searchPar->nu = 1;
    }

    for(ig = 0; ig < searchPar->nv; ig++) {
    for(id = 0; id < searchPar->nu; id++) {
    for(ik = 0; ik < searchPar->nk; ik++) {
    for(ih = 0; ih < searchPar->nh; ih++) {
    for(is = 0; is < searchPar->ns; is++) {
                        arrayMT[isol].gamma = vec_gamma[ig];
                        arrayMT[isol].delta = vec_delta[id];
                        arrayMT[isol].kappa = vec_kappa[ik];
                        arrayMT[isol].theta = vec_theta[ih];
                        arrayMT[isol].sigma = vec_sigma[is];
//                        fprintf(stdout,"index= %20d %11.6f %11.6f %11.6f %11.6f %11.6f\n",
//                                isol, arrayMT[isol].gamma*r2d, arrayMT[isol].delta*r2d, arrayMT[isol].kappa*r2d, arrayMT[isol].theta*r2d, arrayMT[isol].sigma*r2d);
                        //fprintf(stdout,"index= %20d %6d %6d %6d %6d %6d\n", isol, ig, id, ik, ih, is);
                        //fprintf(stdout,"index= %20d %11.6f %11.6f %11.6f %11.6f %11.6f\n",
                        //        isol, vec_gamma[ig]*r2d, vec_delta[id]*r2d, vec_kappa[ik]*r2d, vec_theta[ih]*r2d, vec_sigma[is]*r2d);
                        isol++;
                    }
                }
            }
        }
    }
    fprintf(stderr,"done. nsol = %d \n", isol);
    if (isol == 1) {
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 
                0, arrayMT[0].gamma *r2d, arrayMT[0].delta *r2d,  arrayMT[0].kappa *r2d,  arrayMT[0].theta *r2d,  arrayMT[0].sigma *r2d); 
    } else if(isol > 1) {
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 
                0, arrayMT[0].gamma *r2d, arrayMT[0].delta *r2d,  arrayMT[0].kappa *r2d,  arrayMT[0].theta *r2d,  arrayMT[0].sigma *r2d); 
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 
                isol-1, arrayMT[isol-1].gamma *r2d, arrayMT[isol-1].delta *r2d,  arrayMT[isol-1].kappa *r2d,  arrayMT[isol-1].theta *r2d,  arrayMT[isol-1].sigma *r2d); 
    } 

    free(vec_v);
    free(vec_u);
    free(vec_beta);
    free(vec_h);

    free(vec_gamma);
    free(vec_delta);
    free(vec_kappa);
    free(vec_theta);
    free(vec_sigma);
}   /* end function getGridMT */

SOLN searchMT( int npar, // 3=mw; 2=iso; 1=clvd; 0=strike/dip/rake
        int nda,
        DATA *obs0,
        int nfm,
        FM *fm,
        float fm_thr,
        const int *max_shft,
        float tie,
        MTPAR *mt,
        GRID grid,
        int interp,
        int bootstrap,
        int search_type,
        int norm,
        SEARCHPAR * searchPar,
        ARRAYMT * arrayMT
        )
{
    float mtensor[3][3];
    float amp, rad[6], arad[4][3], x, x1, x2, y, y1, y2, cfg[NCP], s3d[9], temp[3], m_par, del_dip, del_iso;
    SOLN sol, best_sol;
    float *f_pt0, *f_pt1, *r_pt, *r_pt0, *r_pt1, *z_pt, *z_pt0, *z_pt1, *grd_err, *rnd_stk, *rnd_dip, *rnd_rak, *rnd_iso, *rnd_clvd, *iso;
    int misfit_fmp;

    int isol, nreject;
    int isol_best;
    int imag, nmag;
    float VR;
    int tid;

#ifdef WB
    // output files for postprocessing
    // mt = moment tensor elements
    // bb = gamma, delta, strike, dip, rake

    char outFileMT[255];
    char outFileBB[255];
    FILE *fidmt, *fidbb;
    OUTPUTBB outbb;
    OUTPUTMT outmt;

    if(search_type == 1) {
        sprintf(outFileMT, "capout_grid_mt.bin");
        sprintf(outFileBB, "capout_grid_bb.bin");
    } else if (search_type == 2) {
        sprintf(outFileMT, "capout_rand_mt.bin");
        sprintf(outFileBB, "capout_rand_bb.bin");
    } else {
        fprintf(stderr,"Abort. wrong search type.\n");
        exit(-1);
    }

    fidmt=fopen(outFileMT,"wb");
    fidbb=fopen(outFileBB,"wb");
#endif

    // count number of solutions computed
    // (this is ongoing, this variable may be deleted)
    // # 20151216 celso alvizuri - cralvizuri@alaska.edu
    int count=0;

    /* vars to track smallest misfit at each (gamma,delta) on the lune */
    FILE *fidmol;
    //    LUNE_MISFIT * bestmisfit;
    //    bestmisfit = (LUNE_MISFIT *) malloc(sizeof(LUNE_MISFIT));
    //    if(misfit_on_lune)
    //    {
    //        fidmol=fopen("out.misfit.wf_","w");
    //    }

    /* output of first motion polarity */
    FILE *fidfmp;
    if(only_first_motion)
    {
        fidfmp=fopen("out.misfit.fmp_","w");
    }

    /* TODO these will come from user input */
    temp[0] = 4.5; // ALASKA
    temp[0] = 2.8; // UTURUNCU
    temp[0] = 5.2; // ILLINOIS
    // temp[0] is set to starting magnitude (when no magnitude search).
    // another option is define point mag mw0 if point mag search
    temp[0] = searchPar->mw1;

    // prepare magnitude array
    // NOTE magnitude points include boundaries (unlike gridvec)
    float * vec_mag = calloc(searchPar->nmw, sizeof(float));
    magvec(searchPar->mw1, searchPar->mw2, searchPar->dmw, vec_mag);

    //amp = pow(10.,1.5 * +16.1-20);
    //amp = pow(10.,1.5*temp[0]+16.1-20);

    grd_err = grid.err;
    best_sol.err = FLT_MAX;

    fprintf(stderr,"\nRunning searchMT (nsol = %d) ... \n\n", searchPar->nsol);

    // LOOP MAGNITUDE
    nmag = searchPar->nmw; 
    for(imag = 0; imag < nmag; imag++) {
        fprintf(stderr,"Loop Mw[%d/%d] = %6.3f\n", imag+1, nmag, vec_mag[imag]);
        //        temp[0] = vec_mag[imag];
        amp = pow(10., 1.5 * vec_mag[imag] + 16.1 - 20);

        // LOOP THROUGH ALL SOLUTIONS
        // MAIN LOOP 
        nreject = 0;
#ifdef OMP
#pragma omp parallel for private(sol, mtensor)
#endif
        for(isol = 0; isol < searchPar->nsol; isol++) {
            /*
               TODO output coords (u,v) <--> (v,w)
               (v,w) should be straightforward from (u,v)
               w = 3pi/8 - u
               */

            //           temp[2] = arrayMT[isol].gamma * r2d;
            //           temp[1] = arrayMT[isol].delta * r2d;
            //           sol.meca.stk = arrayMT[isol].kappa * r2d;
            //           sol.meca.dip = arrayMT[isol].theta * r2d;
            //           sol.meca.rak = arrayMT[isol].sigma * r2d;

            //tt2cmt(temp[2], temp[1], 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);
            tt2cmt(arrayMT[isol].gamma * r2d, arrayMT[isol].delta * r2d, 1.0, arrayMT[isol].kappa * r2d, arrayMT[isol].theta * r2d, arrayMT[isol].sigma * r2d, mtensor);

            // KEY COMMAND. reject this solution if first motion polarities do not match.
            // NOTE the weight file needs to have polarity picks, otherwise this
            // command will not reject any solution.
            if (check_first_motion(mtensor,fm,nfm,fm_thr)<0) {
                nreject++;
                continue;
            }

            // compute misfit from first motion. data will be output to out.misfit.fmp_
            // misfit_fmp = misfit_first_motion(mtensor, nfm, fm, fidfmp, temp[2], temp[1], temp[0], sol.meca.stk, sol.meca.dip, sol.meca.rak);
            misfit_fmp = misfit_first_motion(mtensor, nfm, fm, fidfmp, arrayMT[isol].gamma * r2d, arrayMT[isol].delta * r2d, vec_mag[imag], sol.meca.stk, sol.meca.dip, sol.meca.rak);

            //--------------KEY COMMAND---call misfit function------
            sol=calerr(nda,obs0,max_shft,tie,norm,mtensor,amp,sol);

            // Ncomp = number of components.
            sol.err=sol.err/Ncomp;
            //*grd_err++ = sol.err; /*error for this solution*/
            //

            //        fprintf(stdout,"index= %10d\t%11.6f %11.6f %11.6f %11.6f %11.6f\n", isol, temp[2], temp[1], sol.meca.stk, sol.meca.dip, sol.meca.rak);
            if (best_sol.err > sol.err) {
                isol_best = isol;
#ifdef OMP
                tid = omp_get_thread_num();
#endif
                //    mt[0].par = temp[0];
                //     mt[1].par = temp[1];
                //      mt[2].par = temp[2];
                mt[0].par = vec_mag[imag];
                mt[2].par = arrayMT[isol].gamma * r2d;
                mt[1].par = arrayMT[isol].delta * r2d;
                sol.meca.stk = arrayMT[isol].kappa * r2d;
                sol.meca.dip = arrayMT[isol].theta * r2d;
                sol.meca.rak = arrayMT[isol].sigma * r2d;

                best_sol = sol; 

                /* output search status */
                VR = 100.*(1.-(sol.err/data2)*(sol.err/data2));
                fprintf(stderr,"(tid= %d) best sol isol=%9d (%3d%) mag=%5.2f %11.6f %11.6f %11.6f %11.6f %11.6f \t err=%13.6e VR=%6.1f%\n", 
                        tid,
                        isol, 100 * isol/searchPar->nsol,
                        temp[0], arrayMT[isol].gamma * r2d, arrayMT[isol].delta * r2d, sol.meca.stk, sol.meca.dip, sol.meca.rak,
                        sol.err, VR);

                /* output variables (only use for debug) */
                /*
                   fprintf(stderr,"best sol isol=%9d (%3d%) sol.err= %13.6e mag=%5.2f %11.6f %11.6f %11.6f %11.6f %11.6f\n", 
                   arrayMT[isol].gamma * r2d, arrayMT[isol].delta * r2d, 
                   arrayMT[isol].kappa * r2d, arrayMT[isol].theta * r2d, arrayMT[isol].sigma * r2d);
                   */
            }

            //  output binary data
#ifdef WB
//            outbb.g = arrayMT[isol].gamma * r2d;
//            outbb.d = arrayMT[isol].delta * r2d;
            outbb.v = gamma2v(arrayMT[isol].gamma);
            outbb.w = delta2w(arrayMT[isol].delta);
            outbb.kappa = arrayMT[isol].kappa * r2d; // output in degrees
            outbb.theta = arrayMT[isol].theta * r2d; // output in degrees
            outbb.sigma = arrayMT[isol].sigma * r2d; // output in degrees
            outbb.mag = vec_mag[imag];
            outbb.misfit_wf  = sol.err/data2;
            outbb.misfit_fmp = (float) misfit_fmp;

            outmt.mrr = mtensor[2][2];
            outmt.mtt = mtensor[0][0];
            outmt.mpp = mtensor[1][1];
            outmt.mrt = mtensor[0][2];
            outmt.mrp = -mtensor[1][2];
            outmt.mtp = -mtensor[0][1];
            outmt.mag = vec_mag[imag];
            outmt.misfit_wf = sol.err/data2;
            outmt.misfit_fmp = (float) misfit_fmp;

            fwrite(&outmt, sizeof outmt, 1, fidmt); 
            fwrite(&outbb, sizeof outbb, 1, fidbb);
#endif

        } /* end loop over solutions */
    } // end loop over magnitudes

    fprintf(stderr,"\nSearch completed.\n\n");

#ifdef WB
    // close binary output files
    fclose(fidmt);
    fclose(fidbb);
    fprintf(stderr,"binary data saved to file: %s \n", outFileMT);
    fprintf(stderr,"binary data saved to file: %s \n", outFileBB);
#endif

    free(vec_mag);

    if(nreject == searchPar->nsol) {
        fprintf(stderr,"----------------------------------------\n");
        fprintf(stderr,"Total solutions rejected nsol=  %10d (%6.2f%)\n", nreject, 100 * (float) nreject / searchPar->nsol);
        fprintf(stderr,"No solutions found. increase number of solutions or check first motion polarities.\n");
        fprintf(stderr,"----------------------------------------\n");
        exit(-1);
    } else {
        fprintf(stderr,"----------------------------------------\n");
        fprintf(stderr,"Best solution at index= %10d\n", isol_best);
        fprintf(stderr,"Total solutions processed nsol= %10d (%6.2f%)\n", isol, 100 *  (float) isol/searchPar->nsol);
        fprintf(stderr,"Total solutions rejected nsol=  %10d (%6.2f%)\n", nreject, 100 * (float) nreject / searchPar->nsol);
        fprintf(stderr,"----------------------------------------\n\n");
    }

    return(best_sol);
}

