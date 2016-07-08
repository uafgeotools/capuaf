/****************************************************************
				cap.c
  Cut-and-Paste program. The code uses two windows (P and S) of
  3-component waveforms to determined the moment tensor
    M_ij = M0 * [ sqrt(1-iso*iso)*D_ij  + sqrt(2/3)*iso*I_ij ],
  where I_ij is the unit tensor and D_ij is a deviatoric tensor:
    D_ij = (1/sqrt(1-2*clvd+4*clvd^2)*[(1-clvd)*DC_ij + clvd*CLVD_ij],
  and
    DC_ij   = n_i v_j + n_j v_i,
    CLVD_ij = v_i v_j + n_i n_j - 2 N_i N_j,
    (n is the fault normal, v is the slip vector, N=nXv)
    iso = tr(M)/M0/sqrt(6),    -1<= iso <=1,
    clvd = -m2/(2*m1),         -0.5 <= clvd <=0.25,
  where m1 is the largest positive eigenvalue of the deviatoric tensor
  and m2 is the smallest in absolute value.

  The solution is given in terms of Mw, strike, dip,
  rake, iso, and clvd. The full moment tensor is also given in the
  output file in the formate of
    M0_in_dyncm Mxx Mxy Mxz Myy Myz Mzz
  where x=North, y=East, z=Down.

  For reference, see:
  Zhu and Helmberger, Advancements in source estimation
  techniques using broadband regional seismograms, BSSA, 86, 1634-1641, 1996.
  Zhao and Helmberger, 1994
  Ben-Zion and Zhu, 2012
  Tape and Tape, 2012, 2013

  requires:
	Green's function -- has P and S arrival time set (t1 and t2 in SAC header)
  optional:
	Data		 -- P pick (A in SAC header) --> align with t1
			 -- Pnl window (t1 and t2)
			 -- P-SV, SH window (t3 and t4)

  Modify history:
  June 19, 1998	Lupei Zhu	modified from srct.c
  July  9, 1998 Lupei Zhu	use different velocities for love and rayleigh 
  July 16, 1998 Lupei Zhu	improve m0 estimation using parabola near mimimum
  July 19, 1998 Lupei Zhu	allow inputs for time shifts
  July 26, 1998 Lupei Zhu	taper waveforms (w=0.4)
  Jan. 29, 1998 Lupei Zhu	add option of repeat inversion after discard bad comp.
  Nov.  9, 1999 Lupei Zhu	absorb shft_pnl into constant shift
				use rint() to convert float (t/dt) to int
  Dec   2, 1999 Lupei Zhu	taper waveform before conv() use w=0.3
  Dec  27, 1999 Lupei Zhu	compute windows using apparent Vp, Vs
  June 28, 2000 Lupei Zhu	switch to new greens function format
  Feb  15, 2001 Lupei Zhu	taper waveform after conv.
  June 27, 2001 Lupei Zhu	add fm_thr (firt-motion threshold)
  July 16, 2001 Lupei Zhu	add directivity option
  Jan. 02, 2002 Lupei Zhu	add an input for number of freedom per sec (nof_per_samp)
  Oct. 31, 2002 Lupei Zhu	use abolute time in the output sac files
  July 18, 2003 Lupei Zhu	use Butterworth filters to band-pass data
  July 30, 2003 Lupei Zhu	not absorb shft_pnl into constant shift
  Aug. 18, 2003 Lupei Zhu	normalize L2 of misfit by # of points
  Aug. 21, 2003 Lupei Zhu	tie SH and SV using variable tie (0-0.5)
  Sep. 28, 2003 Lupei Zhu	use P and S take-off angles in hd.user1/user2
  Oct. 12, 2003 Lupei Zhu	output other local minimums whose
				misfit-min is less than mltp*sigma.
  Apr. 06, 2004 Lupei Zhu	allow inputing a SAC source time function src.sac
  Sep. 06, 2004 Lupei Zhu	make cap to run above the event_dir level
  Jan. 28, 2007 Lupei Zhu	use C-wrapped buttworth filtering subroutines.
  Aug. 20, 2007 Lupei Zhu	use SAC's Butterworth filter routines.
  Jan. 28, 2008 Lupei Zhu	include teleseismic P and SH in inversion (by setting W_PnlR<0).
				The time windows of the z component and r/t components can be different
				for both the data and Greens functions.
  Mar. 11, 2010 Lupei Zhu	Correct a bug when computing m0 using interpolation.
				Change to no interpolation of FM untill the correct mw is found to avoid unstable interpolation in some cases.
  Mar. 12, 2010 Lupei Zhu       modified from cap.c by adding ISO.
  June 10, 2010 Lupei Zhu	Correct a bug introduced in Jan. 2008 which deleted
				distance compensation and Pnl weighting of the
				Greens functions.
  Feb. 13, 2012 Lupei Zhu       revise the decomposition of m_ij.
  Mar.  2, 2012 Lupei Zhu       correct a bug in using discard_bad_data().
  Mar. 25, 2012 Lupei Zhu       output misfit errors for bootstrapping.
  Sept 13, 2012 Lupei Zhu       correct a bug in iso interpolation.
  Oct. 29, 2012 Lupei Zhu	add CLVD and consolidate the searches for
				mw, iso, and clvd.
  Nov.  6, 2012 Lupei Zhu	correct a bug in CLVD_ij (in radiats.c).

  Known bugs:

****************************************************************/
#include "cap.h"

int total_n,loop=0,start=0,debug=0, Ncomp=0,Psamp[STN],Ssamp[STN],edep=-999;
float data2=0.0;

/* flags for computing uncertainty on the lune. 1==apply */
int only_first_motion=0;    // polarity misfit. runs ONLY polarity, no waveform misfit
int misfit_on_lune=0;       // waveform misfit. output misfit on the lune 
char filename_prefix[255];    // used for all output files

/* workaround for filter issues with small magnitude events (Uturuncu) */
// this has not been tested with DIRECTIVITY option
int FTC_data=1, FTC_green=0;// for original CAP set FTC_data=0, FTC_green=0

/* allows use of polarities even when weight=0.
 * Note CAP still needs at least 1 waveform for the inversion */
int skip_zero_weights=1;    // for original CAP set skip_zero_weights=1

// Flag to create regular grid as in Alvizuri & Tape (2016) and Silwal & Tape (2016).
// NOTE reproducibility may not be exact since grid spacing uses function gridvec.
// Function gridvec does not implement the discretization of the previous version of 
// cap.c which uses rules to account for special grid points.
// Function gridvec also avoids endpoints in all parameters.
// See NOTE flag LUNE_GRID_INSTEAD_OF_UV in function sub_inversion.c
int LUNE_GRID_INSTEAD_OF_UV = 0;    // default = 0 (ie do not run old grid mode)

int main (int argc, char **argv) {
  int 	i,j,k,k1,l,m,nda,npt,plot,kc,nfm,useDisp,dof,tele,indx,gindx,dis[STN],tsurf[STN],search_type,norm;
  int	n1,n2,ns, mltp, nup, up[3], n_shft, nqP, nqS,isurf=0,ibody=0,istat=0,Nsurf=0,Nbody=0,Nstat=0;
  int	mm[2],n[NCP],max_shft[NCP],npts[NRC];
  int	repeat;
  char	tmp[255],glib[128],dep[32],dst[16],eve[32],*c_pt;
  float	x,x1,x2,y,y1,amp,dt,rad[6],arad[4][3],fm_thr,tie,mtensor[3][3],rec2=0.,VR,evla,evlo,evdp;
  float	rms_cut[NCP], t0[NCP], tb[NRC], t1, t2, t3, t4, srcDelay;
  float	con_shft[STN], s_shft, shft0[STN][NCP],Pnl_win,ts, surf_win, P_pick[STN], P_win[STN], S_pick[STN], S_win[STN], S_shft[STN],maxamp_syn[200][NCP],maxamp_obs[200][NCP],kcc,lamp_thresh, ppick[200];
  float	tstarP, tstarS, attnP[NFFT], attnS[NFFT];
  float *data[NRC], *green[NGR];
  float	bs_body,bs_surf,bs[NCP],weight,nof_per_samp;
  float	w_pnl[NCP];
  float	distance,dmin=100.,vp,vs1,vs2,depSqr=25;
  float	*f_pt,*f_pt0,*f_pt1;
  float	*data_obs, *data_syn;   // save seismograms for plotting
  float *f_pt2;
  float *g_pt;  // for FTC_green
  int npt_data, offset_h=0;
  GRID	grid;
  MTPAR mt[3];  // DELETE
  COMP	*spt;
  DATA	*obs, *obs0;
  FM	*fm, *fm0;
  FM *fm_copy;  /*  copy of all first motions entered in weight file */

  fprintf(stderr,"\n----------------------------------------------------------\n");
  fprintf(stderr,"Initialize CAP routine\n\n");

  // start random number generator (randvec function). See also cap.h
  fprintf(stderr,"NOTE random seed = %d (randvec function)\n", RANDSEED);
  // option 1 seed is user defined
  srand(RANDSEED);
  // option 2 seed using current time
  //srand(time(NULL));

  // variables for uniform MT
  SEARCHPAR *searchPar = calloc(1, sizeof(SEARCHPAR));

  // variables to verify that Mw_best is not near search limits
  FILE * fid_warn;        // output file for warnings

  SOLN	sol;
  SACHEAD hd[NRC];
  FILE 	*f_out, *wt, *wt2, *wt3, *fid_srcfile ;
  float tau0, riseTime, *src;
  char type[2] = {'B','P'}, proto[2] = {'B','U'};
  double f1_pnl, f2_pnl, f1_sw, f2_sw;
  float pnl_sn[30], pnl_sd[30], sw_sn[30], sw_sd[30];
  long int order=4, nsects;
  void  principal_values(float *);

  fprintf(stderr,"\nReading input parameters for inversion\n");

#ifdef DIRECTIVITY
  int ns_pnl, ns_sw;
  float *src_pnl, *src_sw;
  float tau, faultStr, faultDip, rupDir, rupV, pVel, sVel, temp;
  scanf("%f%f%f%f%f",&pVel,&sVel,&riseTime,&tau0,&rupDir);
  rupDir *= DEG2RAD;
  rupV = 0.8*sVel;
#endif
  
  strcpy(eve,argv[1]);
  strcpy(dep,argv[2]);
 
  /* get station info and polarity */
  FILE *fidfmp;
  FMPDATA *fmpdata;
  fmpdata = (FMPDATA *) malloc(sizeof(FMPDATA));

  /****** input control parameters *************/
  char mod_dep[255]="-999";         /* for renaming .out file */
  char model[128];
  int depth=-999;
  scanf("%s %d",model, &depth); 
  edep=depth;
  
  // WRITE POLARITY AND STATION DATA
  // This section was used in previous CAP with flag only_first_motion=1
  // for generating  polarity misfit on the lune (Uturuncu FMT paper).
  // Now it's set to run for all inversions
  char filename_fmpdata[255];
  strcpy(fmpdata->evid, eve);
  strcpy(fmpdata->vmod, model);
  fmpdata->idep = depth;
  sprintf(mod_dep,"%s_%s_%03d",fmpdata->evid, fmpdata->vmod, fmpdata->idep );   
  if (stat("./OUTPUT_DIR") == -1) {mkdir("./OUTPUT_DIR", 0700);} // create directory is it doesn't exist
  sprintf(filename_prefix, "OUTPUT_DIR/%s_%s_%03d", fmpdata->evid, fmpdata->vmod, fmpdata->idep);
  sprintf(filename_fmpdata, "%s_fmpdata.txt", filename_prefix);
  fidfmp = fopen(filename_fmpdata, "w");
  // end

  scanf("%f%f%f%f%d%f%f",&x1,&y1,&x,&y,&repeat,&fm_thr,&tie);
  if (repeat) for(j=0;j<NCP;j++) scanf("%f",rms_cut+4-j);
  scanf("%f%f%f",&vp,&vs1,&vs2);
  scanf("%f%f%f%f",&bs_body,&bs_surf,&x2,&nof_per_samp);
  scanf("%d",&plot);
  scanf("%d%d",&useDisp,&mltp);
  scanf("%s",glib);
  scanf("%d",&search_type);
  scanf("%d",&norm);
  
  fprintf(stderr, "search type = %d, norm = %d\n", search_type, norm);

  /*** input source functions and filters for pnl and sw ***/
  scanf("%f",&dt);
  if (dt>0.) {
    scanf("%f%f",&tau0,&riseTime);
    if ((src = trap(tau0, riseTime, dt, &ns)) == NULL) {
      fprintf(stderr,"fail to make a trapzoid stf\n");
      return -1;
    }
    srcDelay = 0.;
  } else {
    scanf("%s",tmp); scanf("%f",&riseTime);
    if ((src = read_sac(tmp,hd)) == NULL) {
      fprintf(stderr,"fail to read in source time: %s\n",tmp);
      return -1;
    }
    dt = hd->delta;
    ns = hd->npts;
    srcDelay = -hd->b;
  }
  // write source function to a file
  fid_srcfile = fopen("./OUTPUT_DIR/srcfile","wb");
  for(i=0;i<=ns;i++) {
    fprintf(fid_srcfile,"%f \t %f\n",((double)i*dt),src[i]);
  }
  fclose(fid_srcfile);

  scanf("%lf%lf%lf%lf",&f1_pnl,&f2_pnl,&f1_sw,&f2_sw);
  if (f1_pnl>0.) design(order, type, proto, 1., 1., f1_pnl, f2_pnl, (double) dt, pnl_sn, pnl_sd, &nsects);
  if (f1_sw>0.)  design(order, type, proto, 1., 1., f1_sw, f2_sw, (double) dt, sw_sn, sw_sd, &nsects);

  /** max. window length, shift, and weight for Pnl portion **/
  mm[0]=rint(x1/dt);
  max_shft[3]=max_shft[4]=2*rint(x/dt);
  w_pnl[3]=w_pnl[4]=x2;
  /** max. window length, shift, and weight for P-SV, SH **/
  mm[1]=rint(y1/dt);
  max_shft[0]=max_shft[1]=max_shft[2]=2*rint(y/dt);
  w_pnl[0]=w_pnl[1]=w_pnl[2]=1;
  /** and tie of time shifts between SH and P-SV **/

  /** begin -- get range of search parameters **/
  fprintf(stderr, "\nInput parameter ranges \n");
  // NOTE magnitude parameters include dMw and number of points
  scanf("%f%f%d%f", &searchPar->mw1, &searchPar->mw2, &searchPar->nmw, &searchPar->dmw);
  // (v, w, k, h, s) format: (start, end, number of points and grid spacings (for regular grid))
  scanf("%f%f%d%d", &searchPar->v1, &searchPar->v2, &searchPar->nv, &searchPar->dv);
  scanf("%f%f%d%d", &searchPar->w1, &searchPar->w2, &searchPar->nw, &searchPar->dw);
  scanf("%f%f%d%d", &searchPar->k1, &searchPar->k2, &searchPar->nk, &searchPar->dk);
  scanf("%f%f%d%d", &searchPar->h1, &searchPar->h2, &searchPar->nh, &searchPar->dh);
  scanf("%f%f%d%d", &searchPar->s1, &searchPar->s2, &searchPar->ns, &searchPar->ds);
  // initialize u. It's used internally so is not passed from cap.pl
  searchPar->u1 = 0; searchPar->u2 = 0; searchPar->nu = 0; searchPar->du = 0;
  // total number of solutions
  scanf("%d", &searchPar->nsol);
  // u = [0, 3pi/4], w = [-3pi/8, +3pi/8]. u, w have the same number of points.
  searchPar->nu = searchPar->nw;

  // output values
  fprintf(stderr, "mw1= %6.3f mw2= %6.3f nmw= %d dmw= %6.3f\n", searchPar->mw1, searchPar->mw2, searchPar->nmw, searchPar->dmw);
  fprintf(stderr, "v1= %11.6f v2= %11.6f nv= %10d dv= %10d\n", searchPar->v1, searchPar->v2, searchPar->nv, searchPar->dv);
  fprintf(stderr, "w1= %11.6f w2= %11.6f nw= %10d dw= %10d\n", searchPar->w1, searchPar->w2, searchPar->nw, searchPar->dw);
  fprintf(stderr, "k1= %11.6f k2= %11.6f nk= %10d dk= %10d\n", searchPar->k1, searchPar->k2, searchPar->nk, searchPar->dk);
  fprintf(stderr, "h1= %11.6f h2= %11.6f nh= %10d dh= %10d\n", searchPar->h1, searchPar->h2, searchPar->nh, searchPar->dh);
  fprintf(stderr, "s1= %11.6f s2= %11.6f ns= %10d ds= %10d\n", searchPar->s1, searchPar->s2, searchPar->ns, searchPar->ds);
  fprintf(stderr, "\nNumber of solutions to prepare = %10d\n", searchPar->nsol);

  // allocate memory for (strike, dip, rake)
  // grid.err = (float *) malloc(grid.n[0]*grid.n[1]*grid.n[2]*sizeof(float));
  grid.err = (float *) calloc(searchPar->nsol, sizeof(float));
  if (grid.err == NULL ) {
    fprintf(stderr,"fail to allocate memory for storing misfit errors\n");
    return -1;
  }

    fprintf(stderr,"Allocating memory for moment tensors (nsol = %10d) ... ", searchPar->nsol);
    ARRAYMT * arrayMT = calloc(searchPar->nsol, sizeof(ARRAYMT));

    if (arrayMT == NULL) {
        fprintf(stderr,"Abort. unable to allocate.\n");
        return 0;
    } else {
        fprintf(stderr,"done.\n\n");
    }

  /** end -- get range of search parameters **/

#ifdef DIRECTIVITY
  faultStr = grid.x0[0]*DEG2RAD;
  faultDip = grid.x0[1]*DEG2RAD;
#endif

  /** input number of stations **/
  scanf("%d",&nda);

  if (nda > STN) {
    fprintf(stderr,"number of station, %d, exceeds max., some stations are discarded\n",nda);
    nda = STN;
  }
  obs = obs0 = (DATA *) malloc(nda*sizeof(DATA));
  fm = fm0 = (FM *) malloc(3*nda*sizeof(FM));
  
  /* used when not discarding stations with zero weight */
  if (skip_zero_weights==0){
      fm_copy = (FM *) malloc(nda*sizeof(FM));
  }

  if (obs == NULL || fm == NULL) {
    fprintf(stderr,"fail to allocate memory for data\n");
    return -1;
  }
  
  /**** loop over stations *****/
  total_n = 0;
  n_shft = 0;
  nfm = 0;

  fprintf(stderr,"Reading waveform data (nsta = %d) ... \n", nda);
 for(i=0;i<nda;i++) {

    fprintf(stderr,"%i ", i+1);

    /***** input station name and weighting factor ******/
    scanf("%s%s",tmp,dst);
    for(nup=0,j=0;j<NCP;j++) {
      scanf("%d",&obs->com[4-j].on_off);
      nup += obs->com[4-j].on_off;
    }
    scanf("%f%f%f%f%f",&x1,&Pnl_win,&ts,&surf_win,&s_shft);

    tsurf[i]=ts;
    tele = 0;
    bs[0] = bs[1] = bs[2] = bs_surf;
    bs[3] = bs[4] = bs_body;
    if (obs->com[3].on_off<0) {
      tele = 1;
      tstarS = obs->com[1].on_off;
      tstarP = obs->com[2].on_off;
      obs->com[1].on_off = obs->com[2].on_off = obs->com[3].on_off = 0;
      nup = obs->com[0].on_off + obs->com[4].on_off;
      bs[0] = bs[1] = bs[2] = bs_body;
      j = NFFT;
      if (tstarP>0.) fttq_(&dt, &tstarP, &j, &nqP, attnP);
      if (tstarS>0.) fttq_(&dt, &tstarS, &j, &nqS, attnS);
    }

    /*  original code: remove current station if all weights==0 */
    /*  updated code: don't remove station  */
    if (skip_zero_weights==1){
        if (nup==0) {	/* skip this station */
            nda--; i--;
            continue;
        }
    }

    /* up[i] unknown if not in weight file, so initialize*/
    up[0] = 0;
    up[1] = 0;
    up[2] = 0;
    nup = sscanf(tmp,"%[^/]/%d/%d/%d",obs->stn,&up[0],&up[1],&up[2]);
    if ( fm_thr > 1 ) nup = 1;

    /**************input waveforms************/
    strcat(strcat(strcat(strcpy(tmp,eve),"/"),obs->stn), ".t");
    c_pt = strrchr(tmp,(int) 't');
    for(j=0;j<NRC;j++){
      *c_pt = cm[j];
      if ((data[j] = read_sac(tmp,&hd[j])) == NULL) return -1;
      tb[j] = hd[j].b-hd[j].o;
      npts[j] = hd[j].npts;
    }
    obs->az = hd->az;
    obs->dist = distance = hd->dist;
    obs->tele = tele;
    if (x1<=0.) x1 = hd[2].a;
    if (ts<=0.) ts = hd[2].a;
    x1 -= hd[2].o;
    ts -= hd[2].o;
    if (tele && s_shft>0.) s_shft -= hd[0].o;
    t1 = hd[2].t1-hd[2].o;
    t2 = hd[2].t2-hd[2].o;
    t3 = hd[0].t3-hd[0].o;
    t4 = hd[0].t4-hd[0].o;
    if (dst[0]=='0' && dst[1]=='\0')
      snprintf(dst,10,"%1.0f", rint(obs->dist));     /*if 0 distance given use the distance from header files*/
    evla = hd->evla;
    evlo = hd->evlo;
    evdp = hd->evdp;
    /**************compute source time function***********/
#ifdef DIRECTIVITY
    temp = hd->az*DEG2RAD-faultStr;
    temp = rupV*cos(temp)*cos(rupDir)-sin(temp)*sin(rupDir)*cos(faultDip);
    tau = tau0*(1-temp/pVel);
    src_pnl = trap(tau, riseTime, dt, &ns_pnl);
    tau = tau0*(1-temp/sVel);
    src_sw  = trap(tau, riseTime, dt, &ns_sw);
    if (src_pnl == NULL || src_sw == NULL) {
      fprintf(stderr, "failed to make src for pnl or sw\n");
      return -1;
    }
    fprintf(stderr,"station %s %5.1f tauS %5.1f\n",obs->stn,hd->az,tau);
#endif
    
    /************input green's functions***********/
    strcat(strcat(strcat(strcat(strcpy(tmp,glib),dep),"/"),dst),".grn.0");
    c_pt = strrchr(tmp,(int) '0');
//    fprintf(stderr, "NOTE: convolving greens function with src time function (trapezoid) tau0=dura=%f riseTime=%f \n",
//            tau0, riseTime);

    // WRITE POLARITY AND STATION DATA
    // This section was used in previous CAP with flag only_first_motion=1
    // for generating  polarity misfit on the lune (Uturuncu FMT paper).
    // Now it's set to run for all inversions
    fmpdata->azim = hd->az;
    strcpy(fmpdata->stname, obs->stn);
    fmpdata->stlo = hd->stlo;
    fmpdata->stla = hd->stla;
    fmpdata->dist = hd->dist;
    // end

    for(j=0;j<NGR;j++) {
      *c_pt = grn_com[j];
      indx = 0; if (j>1) indx = 1; if (j>=kk[2]) indx=2;
      if ((green[j] = read_sac(tmp,&hd[indx])) == NULL) return -1;
      conv(src, ns, green[j], hd[indx].npts);
      if (tele) {
	if (tstarP>0. && j>=kk[2]) conv(attnP, nqP, green[j], hd[indx].npts);
	if (tstarS>0. && j< kk[2]) conv(attnS, nqS, green[j], hd[indx].npts);
      }
    }
    if (!tele) {hd[0].t2 = hd[2].t2; hd[0].user2 = hd[2].user2;}

    /* generate first-motion polarity data */
    if (nup>1 && (hd[2].user1<0. || hd[0].user2<0.)) {
      fprintf(stderr,"No P/S take-off angle in Greens' function %s\n",tmp);
    } else {
      obs->alpha = hd[2].user1;
      for(j=1;j<nup;j++) {
        /* type:  1=P; 2=SV; 3=SH; positive=up; negative=down */
        fm->type = up[j-1];
        fm->az = obs->az;
        if (abs(fm->type)==1)	fm->alpha = hd[2].user1;
	else 			fm->alpha = hd[0].user2;
        nfm++;
        fm++;
      }
    }

    /* make copy of station and polarity (if no polarity, set=0)  */
    /* note this is a workaround, only works with p-wave first motions */
    /* type:  1=P; 2=SV; 3=SH; positive=up; negative=down */
    if (skip_zero_weights==0){
        fm_copy->type = up[0];
        fm_copy->az = obs->az;
        fm_copy->alpha = hd[2].user1;
        fm_copy++;
    }

    // WRITE POLARITY AND STATION DATA
    // This section was used in previous CAP with flag only_first_motion=1
    // for generating  polarity misfit on the lune (Uturuncu FMT paper).
    // Now it's set to run for all inversions
    fmpdata->pol = up[0];
    fmpdata->toa = hd[2].user1;
    fmpdata->tp = hd[2].t1;
    fmpdata->ts = hd[2].t2;
    // end

    /*** calculate time shift needed to align data and syn approximately ****/
    /* positive shift means synthetic is earlier */
    con_shft[i] = -srcDelay;
    if ( x1 > 0.) {			/* if first-arrival is specified */
      con_shft[i] += x1 - hd[2].t1;	/* use it to align with greens' fn*/
    }
    if (tele && s_shft > x1 ) {
      s_shft -= hd[0].t2+con_shft[i];	/* align teleseismic S */
    }

    /** calculate time windows for Pnl and Surface wave portions **/

    /* for Pnl portion */
    if (t1 < 0 || t2 < 0 ) {	/* no time window in the data trace. use default time window in syn */
      if (!tele && vp>0.)
	t1 = sqrt(distance*distance+depSqr)/vp;	/* use vp to compute t1 */
      else
	t1 = hd[2].t1;					/* use tp as t1 */
      t1 = t1 - 0.4*mm[0]*dt + con_shft[i];
      t2 = hd[0].t2 + 0.2*mm[0]*dt + con_shft[i];	/* ts plus some delay */
      if (Pnl_win != 0)                                 /* for specific length of time window */
	t2 = t1 + Pnl_win;
    }

    /* do the same for the s/surface wave portion */
    if (ts<=0.)                                          /*if S wave arrical is not specified */
      ts= hd[0].t2;
    if (t3 < 0 || t4 < 0 ) {
      if (!tele && vs1>0. && vs2> 0.) {
	t3 = sqrt(distance*distance+depSqr)/vs1 - 0.3*mm[1]*dt;
	t4 = sqrt(distance*distance+depSqr)/vs2 + 0.7*mm[1]*dt;
      }
      else {
	t3 = ts - 0.3*mm[1]*dt;
	t4 = t3+mm[1]*dt;
      }
      if (ts > 0.){                                      /* if surface wave arrival time is given */
	t3 += s_shft;
        t4 += s_shft;
      }
      else{
	t3 += con_shft[i] + s_shft;              /* add con_shft only if surf arrival time is not specified*/
	t4 += con_shft[i] + s_shft;
	fprintf(stderr,"%f %f %f %f\n",t3,t4,hd[0].t2,con_shft[i]);
      }
      if (surf_win != 0)                                /* for specific length of time window */
	t4 = t3 + surf_win;
    }

    /*calculate the time windows */
    n1 = rint((t2 - t1)/dt);	/*Pnl*/
    n2 = rint((t4 - t3)/dt);	/*PSV/SH*/
    if (n1>mm[0]) n1=mm[0];
    if (n2>mm[1]) n2=mm[1];

    /* storing in array so that later it could saved in weight_cap.dat ouput file */
    P_pick[i] = t1;
    P_win[i] = n1*dt;
    S_pick[i] = t3;
    S_win[i] = n2*dt;
    S_shft[i] = s_shft;
    dis[i]=atoi(dst);
    Psamp[i] = n1;
    Ssamp[i] = n2;

    /***window data+Greens, do correlation and L2 norms **/
    t0[0]=t3;			/* love wave */
    t0[1]=t0[2]=t4-n2*dt;	/* rayleigh wave */
    t0[3]=t0[4]=t1;		/* Pnl */
    n[0]=n[1]=n[2]=n2;	n[3]=n[4]=n1;
    shft0[i][0] = shft0[i][1] = shft0[i][2] = s_shft;
    shft0[i][3] = shft0[i][4] = 0.;
    if (obs->com[0].on_off>0) n_shft++;
    if (obs->com[1].on_off>0 || obs->com[2].on_off>0) n_shft++;
    if (obs->com[3].on_off>0 || obs->com[3].on_off>0) n_shft++;
    isurf=0;
    ibody=0;
    istat=0;
    for(spt=obs->com,kc=2,j=0;j<NCP;j++,spt++,kc=NRF) {
        indx  = kd[j];
        gindx = kk[j];
        if (tele) {
            if (j==2) {indx=1; gindx=2;}		/* no vertical S, use the radial */
            if (j==3) {indx=2; gindx=kk[2];}	/* no radial P, use the vertical */
        }
        spt->npt = npt = n[j];
        spt->b = t0[j];
 
	weight = pow(distance/dmin,bs[j]);  // Caution: This weight scales the amplitude of waveforms. w_pnl = 1 ALWAYS (make changes in cap.pl)
	spt->on_off = (int)spt->on_off * w_pnl[j]; // multiply -Dflag to the weights
	if (spt->on_off) {total_n+=npt; Ncomp += spt->on_off;}  // Ncomp = number of all the components

	// count number of surface and body wave components
    if (j<3) {
        if (spt->on_off >= 1) isurf ++;
    }
    else {
        if (spt->on_off >= 1) ibody ++;
    }
        istat += spt->on_off;

        /* FILTER+CUT options for data */
        /* filter then cut */
        if(FTC_data) {
            npt_data = npts[indx]-offset_h;
            f_pt2 = cutTrace(data[indx], npts[indx], offset_h, npt_data);
            taper(f_pt2, npt_data);
            if (f_pt2 == NULL) {
                fprintf(stderr, "fail to window the data\n");
                return -1;
            }
        }
        /* cut then filter */
        else {
            f_pt = cutTrace(data[indx], npts[indx], (int) rint((t0[j]-tb[indx])/dt), npt);
            taper(f_pt, npt);
            if (f_pt == NULL) {
                fprintf(stderr, "fail to window the data\n");
                return -1;
            }
            spt->rec = f_pt; 
        }
        if (j<3) {
            if (f1_sw>0.) {
                /* filter then cut */
                if(FTC_data) {
                    apply(f_pt2,(long int) npt_data, 0,sw_sn,sw_sd,nsects);
                    f_pt = cutTrace(f_pt2, npt_data, (int) rint((t0[j]-tb[indx])/dt), npt);
                    taper(f_pt, npt);
                    spt->rec = f_pt; 
                }
                /* cut then filter */
                else {
                    apply(f_pt,(long int) npt,0,sw_sn,sw_sd,nsects); 
                }
            }
        }
        else {
            if (f1_pnl>0.) {
                /* filter then cut */
                if(FTC_data) {
                    apply(f_pt2,(long int) npt_data, 0,pnl_sn,pnl_sd,nsects);
                    f_pt = cutTrace(f_pt2, npt_data, (int) rint((t0[j]-tb[indx])/dt), npt);
                    taper(f_pt, npt);
                    spt->rec = f_pt; 
                }
                /* cut then filter */
                else {
                    apply(f_pt,(long int) npt,0,pnl_sn,pnl_sd,nsects);
                }
            }
        }
        if (useDisp==1) cumsum(f_pt, npt, dt); /*use displacement data*/
        for(x2=0.,l=0;l<npt;l++,f_pt++) {
            *f_pt *= weight;
            x2+=(*f_pt)*(*f_pt);
        }
        spt->rec2 = x2;
        if (norm==1) x2 = sqrt(x2);
        rec2 += spt->on_off*x2/(spt->npt);

        /* FILTER+CUT options for greens functions */
        for(m=0,k=0;k<kc;k++) {
            /* filter then cut */
            if(FTC_green){
                g_pt = cutTrace(green[gindx+k], hd[indx].npts, 0, hd[indx].npts);
                taper(g_pt, hd[indx].npts);
                if ( g_pt == NULL ) {
                    fprintf(stderr, "fail to window the Greens functions\n");
                    return -1;
                }
            }
            /* cut then filter */
            else{
                f_pt = cutTrace(green[gindx+k], hd[indx].npts, (int) rint((t0[j]-con_shft[i]-shft0[i][j]-hd[indx].b)/dt), npt);
                taper(f_pt, npt);
                if ( f_pt == NULL ) {
                    fprintf(stderr, "fail to window the Greens functions\n");
                    return -1;
                }
                spt->syn[k] = f_pt;
            }
            if (j<3) {
#ifdef DIRECTIVITY
                conv(src_sw, ns_sw, f_pt, npt);
#endif
                if (f1_sw>0.) {
                    /* filter then cut */
                    if(FTC_green){
                        apply(g_pt,(long int) hd[indx].npts, 0,sw_sn,sw_sd,nsects);
                        f_pt = cutTrace(g_pt, hd[indx].npts, (int) rint((t0[j]-con_shft[i]-shft0[i][j]-hd[indx].b)/dt), npt);
                        taper(f_pt, npt);
                        spt->syn[k] = f_pt;
                    }
                    /* cut then filter */
                    else {
                        apply(f_pt,(long int) npt,0,sw_sn,sw_sd,nsects);
                        taper(f_pt, npt);
                    }
                }
            } 
            else {
#ifdef DIRECTIVITY
                conv(src_pnl, ns_pnl, f_pt, npt);
#endif
                if (f1_pnl>0.) {
                    /* filter then cut */
                    if(FTC_green){
                        apply(g_pt,(long int) hd[indx].npts, 0,pnl_sn,pnl_sd,nsects);
                        f_pt = cutTrace(g_pt, hd[indx].npts, (int) rint((t0[j]-con_shft[i]-shft0[i][j]-hd[indx].b)/dt), npt);
                        taper(f_pt, npt);
                        spt->syn[k] = f_pt;
                    }
                    /* cut then filter */
                    else {
                        apply(f_pt,(long int) npt,0,pnl_sn,pnl_sd,nsects);
                        taper(f_pt, npt);
                    }
                }
            }
            if (useDisp) cumsum(f_pt, npt, dt);
            for(l=0;l<npt;l++) f_pt[l] *= weight;
            spt->crl[k] = crscrl(npt,spt->rec,f_pt,max_shft[j]);
            for(x=1.,k1=k;k1>=0;k1--,x=2.) {
                f_pt0=spt->syn[k];
                f_pt1=spt->syn[k1];
                for(x2=0.,l=0;l<npt;l++) {
                    x2+=(*f_pt0++)*(*f_pt1++);
                }
                spt->syn2[m++] = x*x2;
            }
        }
        //fprintf(stderr, "%s %e %e\n",obs->stn, spt->rec2, spt->syn2[j]);
        // fprintf(stderr, "%d %d %d %d \n",ibody,Nbody,isurf,Nsurf);
    } // end of loop over components
    Nsurf += isurf;
    Nbody += ibody;
    Nstat += (istat>0);

    obs++;
    for(j=0;j<NRC;j++) free(data[j]);
    for(j=0;j<NGR;j++) free(green[j]);

    // WRITE POLARITY AND STATION DATA
    // This section was used in previous CAP with flag only_first_motion=1
    // for generating  polarity misfit on the lune (Uturuncu FMT paper).
    // Now it's set to run for all inversions
    fmp_print_parameters(fidfmp, fmpdata);

  }	/*********end of loop over stations ********/
  fprintf(stderr,"\nDone reading waveforms.\n\n");

  // WRITE POLARITY AND STATION DATA
  // This section was used in previous CAP with flag only_first_motion=1
  // for generating  polarity misfit on the lune (Uturuncu FMT paper).
  // Now it's set to run for all inversions
  //fmp_print_parameters(fidfmp, fmpdata);
  fclose(fidfmp);
  free(fmpdata);
  // end

  fprintf(stderr,"Total number of stations Nda= %d\n", nda);
  fprintf(stderr,"Total number of components Ncomp= %d\n", Ncomp);
  fprintf(stderr,"Total body wave components Nbody= %d\n", Nbody);
  fprintf(stderr,"Total surf wave components Nsurf= %d\n", Nsurf);
  fprintf(stderr,"Total Nstat= %d\n",Nstat);

  data2=rec2/Ncomp;
  if (nda < 1) {
    fprintf(stderr,"No station available for inversion\n");
    return -1;
  }

 INVERSION:
  fprintf(stderr,"\n====================\n"); 
  fprintf(stderr,"BEGIN SEARCH\n"); 
  fprintf(stderr,"====================\n"); 

  // call to old error function for reference. can be deleted
  /* sol = error(3,nda,obs0,nfm,fm0,fm_thr,max_shft,tie,mt,grid,0,bootstrap,search_type,norm); */

  // call initSearchMT instead of "error". This call includes extra parameters (searchPar, arrayMT)
  sol = initSearchMT(nda,obs0,nfm,fm0,fm_thr,max_shft,tie,mt,grid,0,search_type,norm, searchPar, arrayMT);

  dof = nof_per_samp*total_n;
  x2 = sol.err/dof;		/* data variance */
  //fprintf(stderr,"\n=========total_n=%d \t dof=%d \t error=%f\t nof=%f===========\n",total_n,dof,sol.err, nof_per_samp);
  /* repeat grid search if needed */
  if ( repeat && discard_bad_data(nda,obs0,sol,x2,rms_cut) ) {
    repeat--;
    goto INVERSION;
  }

  //Compute variance reduction
  if (norm==1)
    VR = 100*(1.-(sol.err/data2)*(sol.err/data2));
  if (norm==2) 
    VR = 100*(1.-sol.err/data2);

/***** output waveforms for both data and synthetics ****/
  i = mm[1]; if(mm[0]>i) i=mm[0];
  data_obs = (float *) malloc(i*sizeof(float));
  data_syn = (float *) malloc(i*sizeof(float));

  if ((data_obs == NULL) || (data_syn == NULL)) {
    fprintf(stderr,"Failed to allocate memory for seismogram output.\n");
    return -1;
  }

  /**************output the results***********************/
  // START DELETE SECTION -- NOT APPLICABLE ANYMORE
// if (sol.flag) fprintf(stderr,"\nWarning: flag=%d => the minimum %5.1f/%4.1f/%5.1f is at boundary\n",sol.flag,sol.meca.stk,sol.meca.dip,sol.meca.rak);
// else principal_values(&(sol.dev[0]));
// for(i=0; i<3; i++) rad[i] = sqrt(2*x2/sol.dev[i]);
// if (sol.meca.dip>90.) {
//   fprintf(stderr,"Warning: dip corrected by %f\n",sol.meca.dip-90);
//   sol.meca.dip = 90.;
// }
// if (search_type) {
//   rad[0]=0.0;
//   rad[1]=0.0;
//   rad[2]=0.0;
//   sol.ms = 1;}
  // END DELETE SECTION -- NOT APPLICABLE ANYMORE

  // output warning if best magnitude = magnitude limit in magnitude search.
  if ( (sol.meca.mag <= searchPar->mw1 && searchPar->dmw != 0) || (sol.meca.mag >= searchPar->mw2 && searchPar->dmw != 0) ) {
      fid_warn = fopen("capout_error.txt","w");
      fprintf(stderr, "\n***********************************************************************\n");
      fprintf(stderr, "\tINVERSION STOPPED. See file capout_error.txt\n");
      fprintf(stderr, "***********************************************************************\n");
      fprintf(fid_warn, "***********************************************************************\n");
      fprintf(fid_warn, "INVERSION STOPPED. Best magnitude Mw = %4.1f is at a boundary.\n", sol.meca.mag);
      fprintf(fid_warn, "Boundaries [Mw1, Mw2] = [%4.1f, %4.1f]\n", searchPar->mw1, searchPar->mw2);
      fprintf(fid_warn, "Consider expanding the magnitude boundary.\n");
      fprintf(fid_warn, "***********************************************************************\n");
      fclose(fid_warn);
      return(-1);
  }

  fprintf(stderr,"Preparing out file ...\n");
  // sprintf(mod_dep,"%s_%s_%03d", eve, model, depth);                       // rename .out file
  // sprintf(mod_dep,"%s_%s_%03d",fmpdata->evid, fmpdata->vmod, fmpdata->idep);
  // strcat(strcat(strcat(strcpy(tmp,eve),"/"),mod_dep),".out");  
  sprintf(tmp, "%s.out", filename_prefix);
  f_out=fopen(tmp,"w");
  fprintf(f_out,"Event %s Model %s FM %4d %9.6f %4d Mw %4.2f rms %9.3e %5d ERR %3d %3d %3d CLVD %3.2f %3.2f ISO %10.6f %3.2f VR %3.1f data2 %9.3e\n",eve,mod_dep,
	  (int) rint(sol.meca.stk), sol.meca.dip, (int) rint(sol.meca.rak),
	  sol.meca.mag, sol.err, dof,
	  (int) rint(rad[0]), (int) rint(rad[1]), (int) rint(rad[2]),
      sol.meca.gamma, sqrt(sol.meca.gamma * x2),
	  sol.meca.delta, sqrt(sol.meca.delta * x2), VR , data2);
  fprintf(f_out,"# Hypocenter_sac_header elat %e elon %e edep %e\n",evla,evlo,evdp);

  // convert Mw to M0 using GCMT convention (also in Aki and Richards, 2002)
  // this is very close to Kanamori1977 (16.1 vs 16.1010)
  amp=pow(10., 1.5 * sol.meca.mag + 16.1 - 20);

  // compute moment tensor values again for best solution. this is used to
  // output the values to the cap out file.
  tt2cmt(sol.meca.gamma, sol.meca.delta, 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);

  // mtensor saved in output file shoudl be in M00, M11, M22, M01, M02, M12 order. FIX HERE and then perhaps also in the cap_plt! (FUTURE WORK)
  fprintf(f_out,"# tensor = %8.3e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",amp*1.0e20,mtensor[0][0],mtensor[0][1],mtensor[0][2],mtensor[1][1],mtensor[1][2],mtensor[2][2]);
  fprintf(f_out,"# norm L%d    # Pwin %g Swin %g    # N %d Np %d Ns %d\n",norm,mm[0]*dt,mm[1]*dt,y1,Nstat,Nbody,Nsurf);

  // START DELETE SECTION -- DOES NOT EXECUTE
// for(i=1;i<sol.ms;i++) {
//   j = sol.others[i];
//   if (grid.err[j]-grid.err[sol.others[0]]<mltp*x2) {
//     k = j/(grid.n[0]*grid.n[1]);
//     k1 = (j-k*grid.n[0]*grid.n[1])/grid.n[0];
//     fprintf(f_out,"# DELETE %3d %2d %3d %4.2f %9.3e %3.1f\n",
//         (int) rint(grid.x0[0]+(j-k1*grid.n[0]-k*grid.n[0]*grid.n[1])*grid.step[0]),
//         (int) rint(grid.x0[1]+k1*grid.step[1]),
//         (int) rint(grid.x0[2]+k*grid.step[2]),
//         sol.meca.mag,grid.err[j],(grid.err[j]-grid.err[sol.others[0]])/x2);
//   }
// } 
// // END DELETE SECTION -- DOES NOT EXECUTE

for(obs=obs0,i=0;i<nda;i++,obs++) {
    for(j=0;j<NCP;j++) {
      k = NCP - 1 - j;
      sol.shft[i][k]=sol.shft[i][k] - max_shft[k]/2;
    }
 }

if (plot==1) {
    fprintf(stderr,"Saving seismograms for stations ... \n");
    // generate synthetic and data waveforms for each 5 component (filtered and cut) 
    // these are deleted at the later stage (see cap.pl)
    for(obs=obs0,i=0;i<nda;i++,obs++) {
        fprintf(stderr,"%d ", i);
        mt_radiat(obs->az, mtensor, arad);
        rad[0]=arad[3][0];
        for(k=1;k<4;k++) rad[k]=arad[3-k][0];
        for(k=4;k<6;k++) rad[k]=arad[6-k][2];
        //strcat(strcat(strcat(strcat(strcat(strcpy(tmp,eve),"/"),dep), "_"),obs->stn),".0");
	strcat(strcat(strcat(strcat(strcat(strcpy(tmp,"OUTPUT_DIR"),"/"),dep), "_"),obs->stn),".0");
        c_pt = strrchr(tmp,(int) '0');
        for(kc=2,f_pt=rad+NRF,spt=obs->com,j=0;j<NCP;j++,spt++,kc=NRF,f_pt=rad) {
            npt=spt->npt;
            hd[0] = sachdr(dt, npt, spt->b);
            hd->dist = obs->dist; hd->az = obs->az; hd->user1 = obs->alpha;
            hd->a = hd->b;

            // find the max amplitude from all stations [i], from all components [j]
            // OBSERVED
            maxamp_obs[i][j]=0.;
            for(l=0;l<npt;l++) {
                data_obs[l] = spt->rec[l]; // amp is the scaled down M0
                if (fabs(data_obs[l]) > maxamp_obs[i][j]) {
                    maxamp_obs[i][j] = fabs(data_obs[l]); // OBSERVED maximum amplitude
                }
            }
            write_sac(tmp,hd[0],data_obs);   // write observed to file
            (*c_pt)++;

            // SYNTHETIC
            maxamp_syn[i][j]=0.;
            for(l=0;l<npt;l++) {
                for(x2=0.,k=0;k<kc;k++) {
                    x2 += f_pt[k] * spt->syn[k][l];
                }
                data_syn[l] = sol.scl[i][j] * x2 *amp;    // x2 is the green's funciton norm
                if (fabs(data_syn[l]) > maxamp_syn[i][j]) {
                    maxamp_syn[i][j] = fabs(data_syn[l]);  // SYNTHETIC maximum amplitude
                }
            }
            //kcc = log(maxamp_obs[i][j]/maxamp_syn[i][j]);
            //fprintf(stderr,"%.3f\t",kcc);
            hd->b -= (shft0[i][j]+con_shft[i]);
            hd->a = hd->b-sol.shft[i][j]*dt;
            write_sac(tmp,hd[0],data_syn);   // write synthetics to file
            (*c_pt)++;
        }
    }
    fprintf(stderr,"\tdone.\n");
} else {
    fprintf(stderr,"Option plot=%d. No plots to create for this inversion.\n", plot);
}

// set fm_copy to start at beginning of array
 if (skip_zero_weights==0){
   for(i=0;i<nda;i++, fm_copy--);
 }

 wt3 = fopen(strcat(strcat(strcpy(tmp,eve),"/"),"weight_run.dat"),"w");
 for(obs=obs0,i=0;i<nda;i++,obs++) {
   //             stname /  distance / shift (what)
   //              1     2     3
   fprintf(f_out,"%-9s %5.1f/%-5.2f",obs->stn, obs->dist, con_shft[i]);
   fprintf(wt3,"%s\t %d\t",obs->stn, dis[i]);
   for(j=0;j<NCP;j++) {
     k = NCP - 1 - j;
     kc = sol.cfg[i][k]; if (kc<0) kc = 0;
     //        on_off / station misfit % / kross-corr (what?) / t-shift / log(Aobs/Asyn) / Aobs / Asyn
     //               4    5    6    7     8     9     10 
     fprintf(f_out," %1d %6.2f %2d %5.2f %5.2f %8.2e %8.2e",
	     obs->com[k].on_off, sol.error[i][k]*100/(Ncomp*sol.err), kc, shft0[i][k]+dt*sol.shft[i][k], 
	     log(maxamp_obs[i][k]/maxamp_syn[i][k]), maxamp_obs[i][k], maxamp_syn[i][k]);
     if (k<3) lamp_thresh = 1.5;  // log(amplitude) threshold for surface waves
     else lamp_thresh = 2.5;   // log(amplitude) threshold for body waves
     if (log(maxamp_obs[i][k]/maxamp_syn[i][k]) > lamp_thresh){
       fprintf(wt3,"%d\t",0);}
     else{
       fprintf(wt3,"%d\t",obs->com[k].on_off);}  
   }
   fprintf(wt3,"%3.1f\t %3.1f\t %3.1f\t %3.1f\t %3.1f\n", ppick[i], 0., 0., 0., 0.);

   /* output observed polarity and predicted rad amplitude */
   if (skip_zero_weights==1){
     fprintf(f_out, " 0 0\n"); // note leading space
   } 
   else {
     /* if no polarity then output pol = 0 */
     if ( (fm_copy->type == 0) ) {
       fprintf(f_out," 0 %6.2f\n", radpmt(mtensor, fm_copy->alpha, fm_copy->az, 1)); 
     }    
     else {
       fprintf(f_out, " %2d ", fm_copy->type);
       fprintf(f_out, "%6.2f\n", radpmt(mtensor, fm_copy->alpha, fm_copy->az, 1)); 
     }    
     fm_copy++; 
   }
 }
 fclose(f_out);
 fclose(wt3);

 /**********ouput weight file **********/
 wt = fopen(strcat(strcat(strcpy(tmp,eve),"/"),"weight_capout.dat"),"w");
 wt2 = fopen(strcat(strcat(strcpy(tmp,eve),"/"),"weight_capin.dat"),"w");
 for(obs=obs0,i=0;i<nda;i++,obs++){
   fprintf(wt,"%s\t %d\t %d\t %d\t %d\t %d\t %d\t %3.1f\t %3.1f\t %3.1f\t %3.1f\t %3.1f\n",obs->stn, dis[i], obs->com[4].on_off, obs->com[3].on_off, obs->com[2].on_off, obs->com[1].on_off, obs->com[0].on_off, P_pick[i], P_win[i], S_pick[i], S_win[i], S_shft[i]);
   fprintf(wt2,"%s\t %d\t %d\t %d\t %d\t %d\t %d\t %3.1f\t %3.1f\t %3.1f\t %3.1f\t %3.1f\n",obs->stn, dis[i], obs->com[4].on_off, obs->com[3].on_off, obs->com[2].on_off, obs->com[1].on_off, obs->com[0].on_off, P_pick[i]+0.2*mm[0]*dt, P_win[i], S_pick[i]+0.3*mm[1]*dt, S_win[i], S_shft[i]);
 }
 fclose(wt);
 fclose(wt2);

 free(grid.err);
 free(arrayMT); 
 free(searchPar);
 free(obs0);
 free(fm0);
 free(fm_copy);
 free(data_syn);
 free(data_obs);

 fprintf(stderr,"\ncap.c: DONE\n");
 fprintf(stderr,"----------------------------------------\n\n");
 return 0;
}

