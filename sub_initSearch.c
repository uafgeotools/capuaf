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
    FILE * fid;
    int isol;

    //TODO move this into LUT
    char interpfile[256] = "../values_u2beta_1e3.txt";
    int ninterp = 1e3;
    int na;

    /* prepare arrays to calculate beta(u) */
    fid = fopen(interpfile, "r");
    if(fid == NULL) {
        fprintf(stderr,"abort. unable to open datafile beta(u): %s\n", interpfile);
        exit(-1);
    }

    float * pxa = calloc(ninterp, sizeof(float));
    float * pya = calloc(ninterp, sizeof(float));

    na = 0;
    while(fscanf(fid, "%f %f", &pxa[na], &pya[na]) != EOF)
    {
        na++;
    }
    fclose(fid);

    /* temporary vectors for filling the moment tensor array */
    float * vec_beta  = calloc(searchPar->nsol, sizeof(float));
    float * vec_delta = calloc(searchPar->nsol, sizeof(float));
    float * vec_gamma = calloc(searchPar->nsol, sizeof(float));
    float * vec_dip   = calloc(searchPar->nsol, sizeof(float));
    float * vec_u     = calloc(searchPar->nsol, sizeof(float));
    float * vec_v     = calloc(searchPar->nsol, sizeof(float));
    float * vec_k     = calloc(searchPar->nsol, sizeof(float));
    float * vec_h     = calloc(searchPar->nsol, sizeof(float));
    float * vec_s     = calloc(searchPar->nsol, sizeof(float));

    // convert u to w
    // consider making a function
    searchPar->u0 = searchPar->w0 + (3.0 * PI / 8.0);   // 
    searchPar->uf = searchPar->wf + (3.0 * PI / 8.0);   // 
    if(searchPar->u0 < 0.0) {
        // some times  u0 == -0 then set to +0
        fprintf(stderr,"WARNING u0 = %e new value:\n");
        searchPar->u0 = 0.0;
        fprintf(stderr,"WARNING u0 = %e \n");
    }

    randvec(searchPar->v0, searchPar->vf, searchPar->nv, vec_v);
    randvec(searchPar->u0, searchPar->uf, searchPar->nu, vec_u);
    randvec(searchPar->k0, searchPar->kf, searchPar->nk, vec_k);
    randvec(searchPar->h0, searchPar->hf, searchPar->nh, vec_h);
    randvec(searchPar->s0, searchPar->sf, searchPar->ns, vec_s);

    u2beta_vec(pxa, pya, na, vec_u, vec_beta, searchPar->nu);     //NOTE needs pxa, pya, na from external file
    beta2delta_vec(vec_beta, vec_delta, searchPar->nu);
    v2gamma_vec(vec_v, searchPar->nv, vec_gamma);
    h2dip_vec(vec_h, searchPar->nh, vec_dip);

    fprintf(stderr,"~~~~~~~~~~~input u0 %f \n", searchPar->u0);
    fprintf(stderr,"~~~~~~~~~~~DEBUG u[0] %f \n", vec_u[0]);
    fprintf(stderr,"~~~~~~~~~~~DEBUG beta[0] %f \n", vec_beta[0]);

    free(vec_v);
    free(vec_u);
    free(vec_beta);
    free(vec_h);
    free(pxa);
    free(pya);

    /* fill the moment tensor array with moment tensors */
    fprintf(stderr,"\nFilling moment tensor table (nsol= %10d) ... ", searchPar->nsol);
    for(isol = 0; isol < searchPar->nsol; isol++) {
        arrayMT[isol].g = vec_gamma[isol];
        arrayMT[isol].d = vec_delta[isol];
        arrayMT[isol].k =     vec_k[isol];
        arrayMT[isol].t =   vec_dip[isol];
        arrayMT[isol].s =     vec_s[isol];
        fprintf(stdout,"index = %20d %f %f %f %f %f \n", isol, arrayMT[isol].g *r2d, arrayMT[isol].d *r2d,  arrayMT[isol].k *r2d,  arrayMT[isol].t *r2d,  arrayMT[isol].s *r2d);
    }
    fprintf(stderr,"done. nsol = %d \n", isol);
    if (isol == 1) {
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 0, arrayMT[0].g *r2d, arrayMT[0].d *r2d,  arrayMT[0].k *r2d,  arrayMT[0].t *r2d,  arrayMT[0].s *r2d); 
    } else if(isol > 1) {
    fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 0, arrayMT[0].g *r2d, arrayMT[0].d *r2d,  arrayMT[0].k *r2d,  arrayMT[0].t *r2d,  arrayMT[0].s *r2d); 
    fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", isol-1, arrayMT[isol-1].g *r2d, arrayMT[isol-1].d *r2d,  arrayMT[isol-1].k *r2d,  arrayMT[isol-1].t *r2d,  arrayMT[isol-1].s *r2d); 
    } 

    free(vec_gamma);
    free(vec_delta);
    free(vec_k);
    free(vec_dip);
    free(vec_s);
}   /* end function getRandMT */

/* fill array with a grid uniform moment tensors */
void getGridMT(SEARCHPAR * searchPar, ARRAYMT * arrayMT)
{
    int na, isol;
    int ig, id, ik, ih, is;

    FILE * fid;
    char interpfile[256] = "../values_u2beta_1e3.txt";
    int ninterp = 1e3;

    /* prepare arrays to calculate beta(u) */
    fid = fopen(interpfile, "r");
    if(fid == NULL)
    {
        fprintf(stderr,"abort. unable to open datafile beta(u): %s\n", interpfile);
        exit(-1);
    }

    float * pxa = calloc(ninterp, sizeof(float));
    float * pya = calloc(ninterp, sizeof(float));

    na = 0;
    while(fscanf(fid, "%f %f", &pxa[na], &pya[na]) != EOF)
    {
        na++;
    }
    fclose(fid);

    /* temporary vectors for filling the moment tensor array */
    float * vec_u     = calloc(searchPar->nu, sizeof(float));
    float * vec_beta  = calloc(searchPar->nu, sizeof(float));
    float * vec_delta = calloc(searchPar->nu, sizeof(float));
    float * vec_v     = calloc(searchPar->nv, sizeof(float));
    float * vec_gamma = calloc(searchPar->nv, sizeof(float));
    float * vec_k     = calloc(searchPar->nk, sizeof(float));       /* strike */
    float * vec_h     = calloc(searchPar->nh, sizeof(float));
    float * vec_dip   = calloc(searchPar->nh, sizeof(float));
    float * vec_s     = calloc(searchPar->ns, sizeof(float));       /* slip */

    // convert u to w
    // consider making a function
    searchPar->u0 = searchPar->w0 + (3.0 * PI / 8.0);   // 
    searchPar->uf = searchPar->wf + (3.0 * PI / 8.0);   // 
    if(searchPar->u0 < 0.0) {
        // some times  u0 == -0 then set to +0
        fprintf(stderr,"WARNING u0 = %e new value:\n");
        searchPar->u0 = 0.0;
        fprintf(stderr,"WARNING u0 = %e \n");
    }

    gridvec(searchPar->v0, searchPar->vf, searchPar->nv, vec_v);    /* gamma(v) */
    gridvec(searchPar->u0, searchPar->uf, searchPar->nu, vec_u);    /* beta(u) \\ delta = beta -90 */
    gridvec(searchPar->k0, searchPar->kf, searchPar->nk, vec_k);
    gridvec(searchPar->h0, searchPar->hf, searchPar->nh, vec_h);    /* dip(h) */
    gridvec(searchPar->s0, searchPar->sf, searchPar->ns, vec_s);

    u2beta_vec(pxa, pya, na, vec_u, vec_beta, searchPar->nu);     /* gets (pxa, pya, na) from external file */
    beta2delta_vec(vec_beta, vec_delta, searchPar->nu);                /* beta(u) --> delta, npts = nu */
    v2gamma_vec(vec_v, searchPar->nv, vec_gamma);                      /* gamma(v), npts = nv */
    h2dip_vec(vec_h, searchPar->nh, vec_dip);                          /* dip(h), npts = nh */

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
                        arrayMT[isol].g = vec_gamma[ig];
                        arrayMT[isol].d = vec_delta[id];
                        arrayMT[isol].k =     vec_k[ik];
                        arrayMT[isol].t =   vec_dip[ih];
                        arrayMT[isol].s =     vec_s[is];
                        fprintf(stdout,"index= %20d %11.6f %11.6f %11.6f %11.6f %11.6f\n",
                                isol, arrayMT[isol].g*r2d, arrayMT[isol].d*r2d, arrayMT[isol].k*r2d, arrayMT[isol].t*r2d, arrayMT[isol].s*r2d);
                        //fprintf(stdout,"index= %20d %6d %6d %6d %6d %6d\n", isol, ig, id, ik, ih, is);
                        //fprintf(stdout,"index= %20d %11.6f %11.6f %11.6f %11.6f %11.6f\n",
                        //        isol, vec_gamma[ig]*r2d, vec_delta[id]*r2d, vec_k[ik]*r2d, vec_dip[ih]*r2d, vec_s[is]*r2d);
                        isol++;
                    }
                }
            }
        }
    }
    fprintf(stderr,"done. nsol = %d \n", isol);
    if (isol == 1) {
        fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 0, arrayMT[0].g *r2d, arrayMT[0].d *r2d,  arrayMT[0].k *r2d,  arrayMT[0].t *r2d,  arrayMT[0].s *r2d); 
    } else if(isol > 1) {
    fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", 0, arrayMT[0].g *r2d, arrayMT[0].d *r2d,  arrayMT[0].k *r2d,  arrayMT[0].t *r2d,  arrayMT[0].s *r2d); 
    fprintf(stderr,"mt[%10d] = %11.6f %11.6f %11.6f %11.6f %11.6f \n", isol-1, arrayMT[isol-1].g *r2d, arrayMT[isol-1].d *r2d,  arrayMT[isol-1].k *r2d,  arrayMT[isol-1].t *r2d,  arrayMT[isol-1].s *r2d); 
    } 

    free(vec_v);
    free(vec_u);
    free(vec_beta);
    free(vec_h);
    free(pxa);
    free(pya);

    free(vec_gamma);
    free(vec_delta);
    free(vec_k);
    free(vec_dip);
    free(vec_s);
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

    int isol;
    int isol_best;
    float VR;

    OUTPUTGD outgd;
    OUTPUTMT outmt;

    // count number of solutions computed
    // (this is ongoing, this variable may be deleted)
    // # 20151216 celso alvizuri - cralvizuri@alaska.edu
    int count=0;

    /* vars to track smallest misfit at each (gamma,delta) on the lune */
    FILE *fidmol;
    LUNE_MISFIT * bestmisfit;
    bestmisfit = (LUNE_MISFIT *) malloc(sizeof(LUNE_MISFIT));
    if(misfit_on_lune)
    {
        fidmol=fopen("out.misfit.wf_","w");
    }

    /* output of first motion polarity */
    FILE *fidfmp;
    if(only_first_motion)
    {
        fidfmp=fopen("out.misfit.fmp_","w");
    }

    // output files for postprocessing
    // mt = moment tensor elements
    // gd = gamma, delta, strike, dip, rake
    FILE *fidmt, *fidgd;


    fprintf(stderr,"\nRunning searchMT (nsol = %d) ... \n\n", searchPar->nsol);
    //    fprintf(stderr,"************* check values . mag = %f %f %f \n", mt[0].min, mt[0].max, mt[0].par);

    /* TODO these will come from user input */
    temp[0] = 4.5; // ALASKA
    temp[0] = 2.8; // UTURUNCU
    temp[0] = 5.2; // ILLINOIS
    temp[0] = searchPar->mw0;
    //    searchPar->mw0 = 4.5;
    //    searchPar->mwf = 4.5;

    amp = pow(10.,1.5*temp[0]+16.1-20);
    grd_err = grid.err;
    best_sol.err = FLT_MAX;

    for(isol = 0; isol < searchPar->nsol; isol++) {
        /*
           TODO output coords (u,v) <--> (v,w)
           (v,w) should be straightforward from (u,v)
           w = 3pi/8 - u
           */

        /* output variables (only use for debug) */
        /* tt2cmt(arrayMT[isol].g * r2d, arrayMT[isol].d * r2d, 1.0, 
           arrayMT[isol].k * r2d, arrayMT[isol].t * r2d, arrayMT[isol].s * r2d, 
           mtensor);
           */

        temp[2] = arrayMT[isol].g * r2d;
        temp[1] = arrayMT[isol].d * r2d;
        sol.meca.stk = arrayMT[isol].k * r2d;
        sol.meca.dip = arrayMT[isol].t * r2d;
        sol.meca.rak = arrayMT[isol].s * r2d;

        tt2cmt(temp[2], temp[1], 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);

        // compute misfit from first motion. data will be output to out.misfit.fmp_
        misfit_fmp = misfit_first_motion(mtensor, nfm, fm, fidfmp, temp[2], temp[1], temp[0], sol.meca.stk, sol.meca.dip, sol.meca.rak);

        //--------------KEY COMMAND---call misfit function------
        sol=calerr(nda,obs0,max_shft,tie,norm,mtensor,amp,sol);

        // Nsta = number of components.
        sol.err=sol.err/Nsta;
        *grd_err++ = sol.err; /*error for this solution*/

        //        fprintf(stdout,"index= %10d\t%11.6f %11.6f %11.6f %11.6f %11.6f\n", isol, temp[2], temp[1], sol.meca.stk, sol.meca.dip, sol.meca.rak);
        if (best_sol.err > sol.err) {
            isol_best = isol;
            best_sol = sol; 
            mt[0].par = temp[0];
            mt[1].par = temp[1];
            mt[2].par = temp[2];

            /* output search status */
            VR = 100.*(1.-(sol.err/data2)*(sol.err/data2));
            fprintf(stderr,"best sol isol=%9d (%3d%) mag=%5.2f %11.6f %11.6f %11.6f %11.6f %11.6f \t VR=%6.1f%\n", 
                    isol, 100 * isol/searchPar->nsol,
                    temp[0], temp[2], temp[1], sol.meca.stk, sol.meca.dip, sol.meca.rak, VR);

            /* output variables (only use for debug) */
            /*
               fprintf(stderr,"best sol isol=%9d (%3d%) sol.err= %13.6e mag=%5.2f %11.6f %11.6f %11.6f %11.6f %11.6f\n", 
               arrayMT[isol].g * r2d, arrayMT[isol].d * r2d, 
               arrayMT[isol].k * r2d, arrayMT[isol].t * r2d, arrayMT[isol].s * r2d);
               */
        }
    } /* end search */

    // close output files

    fprintf(stderr,"\nTotal solutions processed nsol= %10d (%3d%)\n", isol, 100 *  isol/searchPar->nsol);
    fprintf(stderr,"Best solution at index= %10d\n", isol_best);
    fprintf(stderr,"Search completed.\n\n");

    return(best_sol);
}

