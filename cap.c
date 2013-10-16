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

int total_n,loop=0,start=0,debug=0, only_first_motion=0, misfit_on_lune=0, Npoints,Nsta=0,Psamp[STN],Ssamp[STN];
float data2=0.0,max_amp=0.0,synt2=0.0,synt,st2,err2,synt1,err1,st1,reco,synth;
int main (int argc, char **argv) {
  int 	i,j,k,k1,l,m,nda,npt,plot,kc,nfm,useDisp,dof,tele,indx,gindx,dis[STN],tsurf[STN],search;
  int	n1,n2,ns, mltp, nup, up[3], n_shft, nqP, nqS;
  int	mm[2],n[NCP],max_shft[NCP],npts[NRC];
  int	repeat, bootstrap;
  char	tmp[128],glib[128],dep[32],dst[16],eve[32],*c_pt;
  float	x,x1,x2,y,y1,amp,dt,rad[6],arad[4][3],fm_thr,tie,mtensor[3][3],rec2=0.;
  float	rms_cut[NCP], t0[NCP], tb[NRC], t1, t2, t3, t4, srcDelay;
  float	con_shft[STN], s_shft, shft0[STN][NCP],Pnl_win,ts, surf_win, P_pick[STN], P_win[STN], S_pick[STN], S_win[STN], S_shft[STN];
  float	tstarP, tstarS, attnP[NFFT], attnS[NFFT];
  float *data[NRC], *green[NGR];
  float	bs_body,bs_surf,bs[NCP],weight,nof_per_samp;
  float	w_pnl[NCP];
  float	distance,dmin=100.,vp,vs1,vs2,depSqr=25;
  float	*syn,*f_pt,*f_pt0,*f_pt1;
  GRID	grid;
  MTPAR mt[3];
  COMP	*spt;
  DATA	*obs, *obs0;
  FM	*fm, *fm0;
  SOLN	sol;
  SACHEAD hd[NRC];
  FILE 	*f_out, *log, *wt, *wt2 ;
  float tau0, riseTime, *src;
  char type[2] = {'B','P'}, proto[2] = {'B','U'};
  double f1_pnl, f2_pnl, f1_sw, f2_sw;
  float pnl_sn[30], pnl_sd[30], sw_sn[30], sw_sd[30];
  long int order=4, nsects;
  void  principal_values(float *);
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
  float stnpp;  // stn radius on beachball
  float stnaz;  // stn azim
  float toa;    // takeoff angle
  if(only_first_motion)
  {
      fidfmp = fopen("out.fmp_stndata_","w");
  }

  /****** input control parameters *************/
  char mod_dep[]="-999";         /* 20130102 calvizuri - for renaming .out file */
  char model[128];
  int depth=-999;
  scanf("%s %d",model, &depth);  /* 20130102 calvizuri - end */

  scanf("%f%f%f%f%d%d%f%f",&x1,&y1,&x,&y,&repeat,&bootstrap,&fm_thr,&tie);
  if (repeat) for(j=0;j<NCP;j++) scanf("%f",rms_cut+4-j);
  scanf("%f%f%f",&vp,&vs1,&vs2);
  scanf("%f%f%f%f",&bs_body,&bs_surf,&x2,&nof_per_samp);
  scanf("%d",&plot);
  scanf("%d%d",&useDisp,&mltp);
  scanf("%s",glib);
  scanf("%d",&search);
  
  fprintf(stderr,"=========search = %d==========\n",search);
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

  /** input grid-search range **/
  scanf("%f%f%f",&(mt[0].par),&(mt[0].dd),&(mt[1].dd)); mt[0].min =  1.;  mt[0].max = 10.;
  if(search) {      // use ranges for lune parameters
    scanf("%f%f",&(mt[1].min),&(mt[1].max));  // -90 to 90
    scanf("%f%f",&(mt[2].min),&(mt[2].max));  // -30 to 30
    if (search==1){
      if ((mt[1].min==0.) && (mt[1].max==0.) && (mt[2].min==0.) && (mt[2].max==0.) && (mt[1].dd != 0.)){
	fprintf(stderr,"Warning: Range not specified (-J flag); Doing full grid search\n");
	mt[1].min=-90.; mt[1].max=90.; mt[2].min=-30.; mt[2].max=30;}
    if (mt[1].dd==0. && (mt[1].min != mt[1].max) && (mt[2].min != mt[2].max)){ 
      mt[1].dd=10.; fprintf(stderr,"Warning: Increment not specified (-I flag); Setting it to 10\n");}
    }
    mt[2].dd=mt[1].dd;
    mt[1].par=mt[1].min;
    mt[2].par=mt[2].min;
  }
  else {
    scanf("%f%f",&(mt[1].par),&(mt[1].dd)); mt[1].min = -1.; mt[1].max = 1.;   // -1 to 1
    scanf("%f%f",&(mt[2].par),&(mt[2].dd)); mt[2].min = -0.5; mt[2].max = 0.25;  //-0.5 to 0.25
  }
  for(j=0;j<3;j++) {
    scanf("%f%f%f",&x1,&x2,&grid.step[j]);
    grid.n[j] = rint((x2-x1)/grid.step[j]) + 1;
    grid.x0[j] = x1;
  }
  grid.err = (float *) malloc(grid.n[0]*grid.n[1]*grid.n[2]*sizeof(float));
  if (grid.err == NULL ) {
    fprintf(stderr,"fail to allocate memory for storing misfit errors\n");
    return -1;
  }

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
  if (obs == NULL || fm == NULL) {
    fprintf(stderr,"fail to allocate memory for data\n");
    return -1;
  }
  
  /**** loop over stations *****/
  total_n = 0;
  n_shft = 0;
  nfm = 0;
  for(i=0;i<nda;i++) {

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
    if (nup==0) {	/* skip this station */
      nda--; i--;
      continue;
    }

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

    if(only_first_motion)
    {
        // get first motion parameters from each station
        stnaz = hd->az;
        fprintf(fidfmp,"%8s %10.6f %10.6f d= %10.6f az= %10.6f ",
                obs->stn, hd->stlo, hd->stla, hd->dist, stnaz);
    }
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
        fm->type = up[j-1];
        fm->az = obs->az;
        if (abs(fm->type)==1)	fm->alpha = hd[2].user1;
	else 			fm->alpha = hd[0].user2;
        nfm++;
        fm++;
      }
    }

    /* get data for first motion polarity */
    if(only_first_motion)
    {
        toa = hd->user1;
        if(hd->user1 > 90.0)
        {
            // project to lower hemisphere
            stnaz += 180.0;
            toa = 180.0-toa;
        }
        stnpp = sqrt(2) * sin(toa * PI/360.0);
        fprintf(fidfmp,"toa= %10.6f POL: %2d %10.6f %10.6f\n", 
                hd->user1, up[0], stnpp, stnaz);
    }

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
      t1 = t1 - 0.2*mm[0]*dt + con_shft[i];
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
    for(spt=obs->com,kc=2,j=0;j<NCP;j++,spt++,kc=NRF) {
      indx  = kd[j];
      gindx = kk[j];
      if (tele) {
	if (j==2) {indx=1; gindx=2;}		/* no vertical S, use the radial */
	if (j==3) {indx=2; gindx=kk[2];}	/* no radial P, use the vertical */
      }
      spt->npt = npt = n[j];
      spt->b = t0[j];
      if (spt->on_off) {total_n+=npt; Nsta += spt->on_off;}
      if (j<3) 
	weight = w_pnl[j]*pow(distance/dmin,bs[j])/sqrt(Ssamp[i]);
      else
	weight = w_pnl[j]*pow(distance/dmin,bs[j])/sqrt(Psamp[i]);
      f_pt = cutTrace(data[indx], npts[indx], rint((t0[j]-tb[indx])/dt), npt);
      if ( f_pt == NULL ) {
	fprintf(stderr, "fail to window the data\n");
	return -1;
      }
      spt->rec = f_pt;  
      if (j<3) {if (f1_sw>0.)  apply(f_pt,(long int) npt,0,sw_sn,sw_sd,nsects); }
      else     {if (f1_pnl>0.) apply(f_pt,(long int) npt,0,pnl_sn,pnl_sd,nsects);}
      if (useDisp==1) cumsum(f_pt, npt, dt); /*use displacement data*/
      for(x2=0.,l=0;l<npt;l++,f_pt++) {
	*f_pt *= weight;
	x2+=(*f_pt)*(*f_pt);
	if ((*f_pt) > max_amp)
	  max_amp = (*f_pt);
	//fprintf(stderr,"%f\n",max_amp);
      }
      spt->rec2 = x2;
      rec2 += spt->on_off*x2;
       for(m=0,k=0;k<kc;k++) {
	f_pt = cutTrace(green[gindx+k], hd[indx].npts, rint((t0[j]-con_shft[i]-shft0[i][j]-hd[indx].b)/dt), npt);
	if ( f_pt == NULL ) {
	  fprintf(stderr, "fail to window the Greens functions\n");
	  return -1;
	}
	spt->syn[k] = f_pt;
	if (j<3) {
#ifdef DIRECTIVITY
	  conv(src_sw, ns_sw, f_pt, npt);
#endif
	  if (f1_sw>0.)  apply(f_pt,(long int) npt,0,sw_sn,sw_sd,nsects);
	} else {
#ifdef DIRECTIVITY
	  conv(src_pnl, ns_pnl, f_pt, npt);
#endif
	  if (f1_pnl>0.) apply(f_pt,(long int) npt,0,pnl_sn,pnl_sd,nsects);
	}
	if (useDisp) cumsum(f_pt, npt, dt);
	for(l=0;l<npt;l++) f_pt[l] *= weight;
	spt->crl[k] = crscrl(npt,spt->rec,f_pt,max_shft[j]);
	for(x=1.,k1=k;k1>=0;k1--,x=2.) {
	  f_pt0=spt->syn[k];
	  f_pt1=spt->syn[k1];
	  for(x2=0.,l=0;l<npt;l++) x2+=(*f_pt0++)*(*f_pt1++);
	  spt->syn2[m++] = x*x2;
	}
      }
      //fprintf(stderr, "%s %e %e\n",obs->stn, spt->rec2, spt->syn2[j]);
    }

    obs++;
    for(j=0;j<NRC;j++) free(data[j]);
    for(j=0;j<NGR;j++) free(green[j]);

  }	/*********end of loop over stations ********/

  /*  first motion polarity - close file */
  if(only_first_motion)
  {
      fclose(fidfmp);
  }
  fprintf(stderr,"Nsta= %d",Nsta);
  data2=rec2/Nsta;
  if (nda < 1) {
    fprintf(stderr,"No station available for inversion\n");
    return -1;
  }

  /************grid-search for full moment tensor***********/
  log = fopen("log_diff","w");
  fclose(log);

  // maybe we label this section "cap messages"...
  if (search==0){
    if (mt[1].dd>=1 || mt[2].dd>=1)
      fprintf(stderr,"Warning: Possible error = Expecting grid-search\n(set search=1) or reduce the search increment (-J flag)\n");
    fprintf(stderr,"----------starting line-search-----------\n");}
  if (search==1){
    if ((mt[1].dd<=1&&mt[1].dd!=0) || (mt[2].dd<=1&&mt[2].dd!=0))
      fprintf(stderr,"Warning: Very fine grid-search = Expecting line-search\n(set search=0) or increase search increment (-J flag)\n");
    fprintf(stderr,"----------starting grid-search-----------\n");
    if(only_first_motion)
    {
        fprintf(stderr,"NOTE: Computing only first motion. results output to out.misfit.fmp_\n");
    }
    if(misfit_on_lune)
    {
        fprintf(stderr,"NOTE: misfit on lune output to file out.misfit.wf_\n");
    }
  }
  if (search==2)
     fprintf(stderr,"----------starting random-search-----------\n");

 INVERSION:
  sol = error(3,nda,obs0,nfm,fm0,fm_thr,max_shft,tie,mt,grid,0,bootstrap,search);

  dof = nof_per_samp*total_n;
  x2 = sol.err/dof;		/* data variance */
  //fprintf(stderr,"\n=========total_n=%d \t dof=%d \t error=%f\t nof=%f===========\n",total_n,dof,sol.err, nof_per_samp);
  /* repeat grid search if needed */
  if ( repeat && discard_bad_data(nda,obs0,sol,x2,rms_cut) ) {
    repeat--;
    goto INVERSION;
  }


  /**************output the results***********************/
  if (sol.flag) fprintf(stderr,"\nWarning: flag=%d => the minimum %5.1f/%4.1f/%5.1f is at boundary\n",sol.flag,sol.meca.stk,sol.meca.dip,sol.meca.rak);
  else principal_values(&(sol.dev[0]));
  for(i=0; i<3; i++) rad[i] = sqrt(2*x2/sol.dev[i]);
  if (sol.meca.dip>90.) {
    fprintf(stderr,"Warning: dip corrected by %f\n",sol.meca.dip-90);
    sol.meca.dip = 90.;
  }
  if (search) {
    rad[0]=0.0;
    rad[1]=0.0;
    rad[2]=0.0;
    sol.ms = 1;}
  //grid.n[0]=grid.n[1]=grid.n[2]=1;
  //grid.x0[0]=sol.meca.stk; grid.x0[1]=sol.meca.dip; grid.x0[2]=sol.meca.rak;
  //sol = error(nda,obs0,nfm,fm0,max_shft,m0,grid,fm_thr,tie);
  sprintf(mod_dep,"%s_%03d",model,depth);                     // 20120723 calvizuri - rename .out file
  strcat(strcat(strcat(strcpy(tmp,eve),"/"),mod_dep),".out"); // 20130102 calvizuri - end rename .out
  //  strcat(strcat(strcat(strcpy(tmp,eve),"/"),dep),".out");   // 20130102 calvizuri - original
  f_out=fopen(tmp,"w");
  fprintf(f_out,"Event %s Model %s FM %3d %2d %3d Mw %4.2f rms %9.3e %5d ERR %3d %3d %3d ISO %3.2f %3.2f CLVD %3.2f %3.2f\n",eve,mod_dep,
	  (int) rint(sol.meca.stk), (int) rint(sol.meca.dip), (int) rint(sol.meca.rak),
	  mt[0].par, sol.err, dof,
	  (int) rint(rad[0]), (int) rint(rad[1]), (int) rint(rad[2]),
	  mt[1].par, sqrt(mt[1].sigma*x2),mt[2].par, sqrt(mt[2].sigma*x2));
  fprintf(f_out,"# Variance reduction %4.1f\n",100*(1.-sol.err/data2));
  amp=pow(10.,1.5*mt[0].par+16.1-20);

  if(search) {
    tt2cmt(mt[2].par, mt[1].par, 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);
  }
  else {
    nmtensor(mt[1].par,mt[2].par,sol.meca.stk,sol.meca.dip,sol.meca.rak,mtensor); //original
  }
  fprintf(f_out,"# tensor = %8.3e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",amp*1.0e20,mtensor[0][0],mtensor[0][1],mtensor[0][2],mtensor[1][1],mtensor[1][2],mtensor[2][2]);
  for(i=1;i<sol.ms;i++) {
    j = sol.others[i];
    if (grid.err[j]-grid.err[sol.others[0]]<mltp*x2) {
      k = j/(grid.n[0]*grid.n[1]);
      k1 = (j-k*grid.n[0]*grid.n[1])/grid.n[0];
      fprintf(f_out,"# %3d %2d %3d %4.2f %9.3e %3.1f\n",
	      (int) rint(grid.x0[0]+(j-k1*grid.n[0]-k*grid.n[0]*grid.n[1])*grid.step[0]),
	      (int) rint(grid.x0[1]+k1*grid.step[1]),
	      (int) rint(grid.x0[2]+k*grid.step[2]),
	      mt[0].par,grid.err[j],(grid.err[j]-grid.err[sol.others[0]])/x2);
    }
  } 
  for(obs=obs0,i=0;i<nda;i++,obs++) {
    fprintf(f_out,"%-9s %5.1f/%-5.2f",obs->stn, obs->dist, con_shft[i]);
    for(j=0;j<NCP;j++) {
      k = NCP - 1 - j;
      sol.shft[i][k]=sol.shft[i][k] - max_shft[k]/2;
      kc = sol.cfg[i][k]; if (kc<0) kc = 0;
      fprintf(f_out," %1d %8.2e %2d %5.2f",obs->com[k].on_off,sol.error[i][k],kc,shft0[i][k]+dt*sol.shft[i][k]);
    }
    fprintf(f_out,"\n");
  }
  fclose(f_out);

  if ( ! plot ) return 0;

  /***** output waveforms for both data and synthetics ****/
  i = mm[1]; if(mm[0]>i) i=mm[0];
  syn = (float *) malloc(i*sizeof(float));
  if (syn == NULL) {
    fprintf(stderr,"fail to allocate memory for output\n");
    return -1;
  }
  for(obs=obs0,i=0;i<nda;i++,obs++){
    mt_radiat(obs->az, mtensor, arad);
    rad[0]=arad[3][0];
    for(k=1;k<4;k++) rad[k]=arad[3-k][0];
    for(k=4;k<6;k++) rad[k]=arad[6-k][2];
    strcat(strcat(strcat(strcat(strcat(strcpy(tmp,eve),"/"),dep), "_"),obs->stn),".0");
    c_pt = strrchr(tmp,(int) '0');
    for(kc=2,f_pt=rad+NRF,spt=obs->com,j=0;j<NCP;j++,spt++,kc=NRF,f_pt=rad) {
      npt=spt->npt;
      hd[0] = sachdr(dt, npt, spt->b);
      hd->dist = obs->dist; hd->az = obs->az; hd->user1 = obs->alpha;
      hd->a = hd->b;
      for(l=0;l<npt;l++) syn[l] = spt->rec[l]/amp;
      write_sac(tmp,hd[0],syn);
      (*c_pt)++;
      for(l=0;l<npt;l++) {
	for(x2=0.,k=0;k<kc;k++) x2 += f_pt[k]*spt->syn[k][l];
	syn[l] = sol.scl[i][j]*x2;
      }
      hd->b -= (shft0[i][j]+con_shft[i]);
      hd->a = hd->b-sol.shft[i][j]*dt;
      write_sac(tmp,hd[0],syn);
      (*c_pt)++;
    }
  }

  /**********ouput weight file -vipul**********/
  wt = fopen(strcat(strcat(strcpy(tmp,eve),"/"),"weight_capout.dat"),"w");
  wt2 = fopen(strcat(strcat(strcpy(tmp,eve),"/"),"weight_capin.dat"),"w");
  for(obs=obs0,i=0;i<nda;i++,obs++){
    fprintf(wt,"%s\t %d\t %d\t %d\t %d\t %d\t %d\t %3.1f\t %3.1f\t %3.1f\t %3.1f\t %3.1f\n",obs->stn, dis[i], obs->com[4].on_off, obs->com[3].on_off, obs->com[2].on_off, obs->com[1].on_off, obs->com[0].on_off, P_pick[i], P_win[i], S_pick[i], S_win[i], S_shft[i]);
    fprintf(wt2,"%s\t %d\t %d\t %d\t %d\t %d\t %d\t %3.1f\t %3.1f\t %3.1f\t %3.1f\t %3.1f\n",obs->stn, dis[i], obs->com[4].on_off, obs->com[3].on_off, obs->com[2].on_off, obs->com[1].on_off, obs->com[0].on_off, P_pick[i]+0.2*mm[0]*dt, P_win[i], S_pick[i]+0.3*mm[1]*dt, S_win[i], S_shft[i]);
  }
  fclose(wt);
  fclose(wt2);

  return 0;
}


// grid-search for the full moment tensor
SOLN	error(	int		npar,	// 3=mw; 2=iso; 1=clvd; 0=strike/dip/rake
		int		nda,
		DATA		*obs0,
		int		nfm,
		FM		*fm,
		float		fm_thr,
		const int	*max_shft,
		float		tie,
		MTPAR		*mt,
		GRID		grid,
		int		interp,
		int		bootstrap,
		int             search
		) {
  int	i, j, k, l, m, k1, kc, z0, z1, z2, ii, N, iso_len;
  float mw_ran; // 20130730 celso - half-range for magnitude search (previously int)
  int	i_stk, i_dip, i_rak, i_iso;
  float	amp, rad[6], arad[4][3], x, x1, x2, y, y1, y2, cfg[NCP], s3d[9], temp[3], m_par, del_dip, del_iso;
  float	*f_pt0, *f_pt1, *r_pt, *r_pt0, *r_pt1, *z_pt, *z_pt0, *z_pt1, *grd_err, *rnd_stk, *rnd_dip, *rnd_rak, *rnd_iso, *rnd_clvd, *iso;
  float dx, mtensor[3][3], *r_iso, *z_iso;
  DATA	*obs;
  COMP	*spt;
  SOLN	sol, sol1, sol2, best_sol;
  FILE *log,*st_err1, *st_err2;
  char logfile[16],st_err_body[20],st_err_surf[20];

  /* parameters that produce least misfit */
  float best_err, best_stk, best_dip, best_rak;
  FILE *fidmol;
  /* output of misfit on lune */
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
  
  if (debug) {
    sprintf(st_err_body,"err_file_body");
    st_err1=fopen(st_err_body,"w");
    fclose(st_err1);
    sprintf(st_err_surf,"err_file_surf");
    st_err2=fopen(st_err_surf,"w");
    fclose(st_err2);
    sprintf(logfile,"%s_%03d","log",loop);
    log = fopen(logfile,"w");
    fclose(log);
  }

  if (debug) fprintf(stderr, "loop=%d start=%d \n",loop,start);
  start++;

  if (search==2){

    //--------newly added section-------------
    fprintf(stderr,"Mw=%f\n",mt[0].par);
    mw_ran=1.0;
    mt[0].max = mt[0].par+mw_ran;
    mt[0].min = mt[0].par-mw_ran;;

    if (mt[0].dd==0){
      mt[0].max = mt[0].par;
      mt[0].min = mt[0].par;
      mt[0].dd=1.0;
    }

    N=10000;
    rnd_stk = (float*)malloc(sizeof(int) * N*sizeof(float));
    rnd_dip = (float*)malloc(sizeof(int) * N*sizeof(float));
    rnd_rak = (float*)malloc(sizeof(int) * N*sizeof(float));
    rnd_iso = (float*)malloc(sizeof(int) * N*sizeof(float));
    rnd_clvd = (float*)malloc(sizeof(int) * N*sizeof(float));
    
    if (rnd_stk==NULL || rnd_stk==NULL || rnd_stk==NULL){
      fprintf(stderr,"Cannot allocate space for random number generation");}

    for(i=0;i<N;i++){
      rnd_stk[i]=0.0+360.0*drand48();
      rnd_dip[i]=0.0+1.0*drand48();
      rnd_rak[i]=-90.0+180.0*drand48();
      if (mt[1].dd==0)
	rnd_iso[i]=sin(mt[1].par*(PI/180.0));
      else
	rnd_iso[i]=-1.0+2.0*drand48();
      if (mt[2].dd==0)
	rnd_clvd[i]=mt[2].par;
      else
	rnd_clvd[i]=-30.0+60.0*drand48();
    }

    best_sol.err = FLT_MAX;

     for(temp[0]=mt[0].min;temp[0]<=mt[0].max;temp[0]=temp[0]+mt[0].dd){ 

	  //--------newly added section ends here-------------

	  //==== the base case: grid-search for strike, dip, and rake =============
	  amp = pow(10.,1.5*temp[0]+16.1-20);
	  grd_err = grid.err;

	  //--------------random search loop-------------------------------------
	  for(ii=0; ii<N; ii++) {
	    temp[1]=asin(rnd_iso[ii])*(180.0/PI);
	    temp[2]=rnd_clvd[ii];
	    sol.meca.rak=rnd_rak[ii];
	    sol.meca.dip=acos(rnd_dip[ii])*(180.0/PI);
	    sol.meca.stk=rnd_stk[ii];

	    //nmtensor(mt[1].par,mt[2].par,sol.meca.stk,sol.meca.dip,sol.meca.rak,mtensor);
	    //nmtensor(temp[1],temp[2],sol.meca.stk,sol.meca.dip,sol.meca.rak,mtensor);
	    tt2cmt(temp[2], temp[1], 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);
	    if (check_first_motion(mtensor,fm,nfm,fm_thr)<0) {
	      *grd_err++ = sol.err = FLT_MAX;
	      continue;
	    }
	    synt=0.;
	    reco=0.;
	    if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %5.2f %5.2f %5.2f %5.1f %5.1f %5.1f\n", mt[0].par, mt[1].par, mt[2].par, sol.meca.stk, sol.meca.dip, sol.meca.rak);
	    for(obs=obs0,sol.err=0.,i=0;i<nda;i++,obs++){
	      synt1=0.;
	      synt2=0.;
	      st2=0.;
	      err2=0.;
	      st1=0.;
	      err1=0.;
	      mt_radiat(obs->az,mtensor,arad);
	      rad[0]=amp*arad[3][0];
	      for(k=1;k<4;k++) rad[k]=amp*arad[3-k][0];
	      for(k=4;k<6;k++) rad[k]=amp*arad[6-k][2];
	    
	      /*******find the time shift*************/
	      /**SH surface wave**/
	      spt = obs->com;
	      f_pt0 = spt->crl[0];
	      f_pt1 = spt->crl[1];
	      z0 = spt->on_off>0?1:0;
	      /**PSV surface wave**/
	      spt++;
	      r_pt0 = spt->crl[1];
	      r_pt1 = spt->crl[2];
	      r_pt  = spt->crl[3];
	      r_iso = spt->crl[0];
	      z1 = spt->on_off>0?1:0;
	      spt++;
	      z_pt0 = spt->crl[1];
	      z_pt1 = spt->crl[2];
	      z_pt  = spt->crl[3];
	      z_iso = spt->crl[0];
	      z2 = spt->on_off>0?1:0;
	      for(y1=y2=-FLT_MAX,l=0;l<=max_shft[1];l++) {
		x =rad[4]*(*f_pt0++)+rad[5]*(*f_pt1++);
		x1=rad[1]*(*r_pt0++)+rad[2]*(*r_pt1++)+rad[3]*(*r_pt++)+rad[0]*(*r_iso++);
		x2=rad[1]*(*z_pt0++)+rad[2]*(*z_pt1++)+rad[3]*(*z_pt++)+rad[0]*(*z_iso++);
		y = (1-tie)*z0*x + tie*(z1*x1 + z2*x2);
		if (y>y2) {y2=y;cfg[0]=x;sol.shft[i][0]=l;}
		y = tie*z0*x + (1-tie)*(z1*x1 + z2*x2);
		if (y>y1) {y1=y;cfg[1]=x1;cfg[2]=x2;m=l;}
	      }
	      sol.shft[i][1]=sol.shft[i][2]=m;
	      /**Pnl*/
	      spt++;
	      r_pt0 = spt->crl[1];
	      r_pt1 = spt->crl[2];
	      r_pt  = spt->crl[3];
	      r_iso = spt->crl[0];
	      z1 = spt->on_off>0?1:0;
	      spt++;
	      z_pt0 = spt->crl[1];
	      z_pt1 = spt->crl[2];
	      z_pt  = spt->crl[3];
	      z_iso = spt->crl[0];
	      z2 = spt->on_off>0?1:0;
	      for(y1=-FLT_MAX,l=0;l<=max_shft[3];l++) {
		x1=rad[1]*(*r_pt0++)+rad[2]*(*r_pt1++)+rad[3]*(*r_pt++)+rad[0]*(*r_iso++);
		x2=rad[1]*(*z_pt0++)+rad[2]*(*z_pt1++)+rad[3]*(*z_pt++)+rad[0]*(*z_iso++);
		y = z1*x1 + z2*x2;
		if (y>y1) {y1=y;cfg[3]=x1;cfg[4]=x2;m=l;}
	      }
	      sol.shft[i][3]=sol.shft[i][4]=m;
	      spt -= NCP - 1;
	    
	      /***error calculation*****/
	      for(kc=2,f_pt1=rad+NRF,j=0;j<NCP;j++,spt++,kc=NRF,f_pt1=rad) {

		/* compute the L2 norm of syn */
		for(x2=0.,f_pt0=spt->syn2,k=0;k<kc;k++)
		  for(k1=k;k1>=0;k1--)
		    x2+=f_pt1[k]*f_pt1[k1]*(*f_pt0++); 

		if (j<3) synt2+=spt->on_off*(x2);
		else synt1+=spt->on_off*x2;
		y1 = 1.;
		/* find out the scaling factor for teleseismic distances */
		if (obs->tele && spt->on_off) {
		  if (cfg[j]>0.) y1 = cfg[j]/x2;
		  else y1 = 0.;
		}
		sol.scl[i][j] = y1;

		x1 = (spt->rec2+x2*y1*y1-2.*cfg[j]*y1);
		sol.error[i][j] = x1;	/*L2 error for this com.*/
		sol.cfg[i][j] = 100*cfg[j]/sqrt(spt->rec2*x2);
		sol.err += spt->on_off*sol.error[i][j];
		if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %-10s %d %d %9.3e\n", obs->stn, j, spt->on_off, x1);
		if (j<3){
		  st2+=spt->on_off*spt->rec2;
		  err2+=spt->on_off*x1;
		}
		else{
		  st1+=spt->on_off*spt->rec2;
		  err1+=spt->on_off*x1;
		}	
	      }

	      if (0){
		st_err1 = fopen(st_err_body,"a");
		fprintf(st_err1, "%f %e %e %e\n",obs->dist,st1, synt1, err1);
		fclose(st_err1);
		st_err2 = fopen(st_err_surf,"a");
		fprintf(st_err2, "%f %e %e %e\n",obs->dist,st2, synt2, err2);
		fclose(st_err2);}
	      synt+=synt2+synt1;
	      reco+=st2+st1;
	    } /*-------------------------end of all stations*/
	    sol.err=sol.err/Nsta;
	    synt=synt/Nsta;
	    reco=reco/Nsta;
	    *grd_err++ = sol.err;		/*error for this solution*/
	    if (best_sol.err>sol.err) {best_sol = sol; synth=synt;
	      if (debug) fprintf(stderr,"misfit for best sol = %f; stk=%3.1f, dip=%3.1f, rak=%3.1f \n",best_sol.err,sol.meca.stk, sol.meca.dip, sol.meca.rak);
	      mt[0].par=temp[0];
	      mt[1].par=temp[1];
	      mt[2].par=temp[2];
	    }

	    if (debug) { 
	      log = fopen(logfile,"a");                 // output log file
	      fprintf(log,"%3.1f\t%3.1f\t%3.1f\t%e\t%2.2f\t%2.2f\t%2.2f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n",sol.meca.stk, sol.meca.dip, sol.meca.rak, sol.err/data2, temp[0], temp[1], temp[2], amp*1.0e20, mtensor[0][0], mtensor[0][1], mtensor[0][2], mtensor[1][1], mtensor[1][2], mtensor[2][2] );
	      fclose(log);
	    }

	    if (1){
	      loop++;
	      if (debug) {
		sprintf(logfile,"%s_000","log");  // changes the log file name for next sext search
		log = fopen(logfile,"a");
		fclose(log);
	      }
	      log = fopen("log_diff","a"); /*fprintf(stderr,"completed stk,dip,rake loop\n");        //summary log file*/
	      fprintf(log,"%d\t%d\t%3.1f\t%3.1f\t%3.1f\t%f\t%2.2f\t%2.2f\t%2.2f\n",loop,interp, best_sol.meca.stk, best_sol.meca.dip, best_sol.meca.rak, best_sol.err, mt[0].par, mt[1].par, mt[2].par);
	      fclose(log);
	    }
	    fprintf(stderr, "%d\t%3.2f\t%3.2f\t%3.2f\t%2.1f\t%2.1f\t%2.2f\n",ii+1,sol.meca.stk, sol.meca.dip,sol.meca.rak,temp[0],temp[1],temp[2]);
	  }
	  fprintf(stderr,"========================Minimum==================================\n");
	  fprintf(stderr, "%3.2f\t%3.2f\t%3.2f\t%2.1f\t%2.1f\t%2.2f\n",best_sol.meca.stk, best_sol.meca.dip,best_sol.meca.rak,temp[0],mt[1].par,mt[2].par);
     }
     // }
     //}
  

    if (debug) fprintf(stderr, "Mw=%5.2f  iso=%5.2f clvd=%5.2f misfit = %9.3e\n", mt[0].par, mt[1].par, mt[2].par, best_sol.err);
    if (interp=1) return(best_sol);
    /* do interpolation */
    best_sol.err = grid3d(grid.err,&(grid.n[0]),s3d,&(best_sol.flag),&(best_sol.ms),best_sol.others);
    if (debug) fprintf(stderr, " interpolation  misfit = %9.3e\n", best_sol.err);
    best_sol.meca.stk = grid.x0[0]+s3d[0]*grid.step[0];
    best_sol.meca.dip = grid.x0[1]+s3d[1]*grid.step[1];
    best_sol.meca.rak = grid.x0[2]+s3d[2]*grid.step[2];
    for(i=0;i<3;i++) best_sol.dev[i]  = s3d[3+i]/(grid.step[i]*grid.step[i]);
    best_sol.dev[3] = s3d[6]/(grid.step[0]*grid.step[1]);
    best_sol.dev[4] = s3d[7]/(grid.step[0]*grid.step[2]);
    best_sol.dev[5] = s3d[8]/(grid.step[1]*grid.step[2]);
    fprintf(stderr,"=======================");
    return(best_sol);

  }

  if (search==1){

    //--------newly added section-------------
    mw_ran = 0.5;
    mt[0].max = mt[0].par+mw_ran;
    mt[0].min = mt[0].par-mw_ran;

    if(only_first_motion) 
    {
        mt[0].min = mt[0].par;
        mt[0].max = mt[0].par;
    }
 
    for (ii=0; ii<3; ii++){
    if (mt[ii].dd==0){
      mt[ii].max = mt[ii].par;
      mt[ii].min = mt[ii].par;
      mt[ii].dd=1.0;
    }}
  
    iso_len = rint((mt[1].max - mt[1].min)/mt[1].dd) + 1;

    //Output search ranges
    fprintf(stderr,"=========GRID-SEARCH RANGE===========\n");
    for (ii=0; ii<3; ii++){
      if (ii==0)
	fprintf(stderr,"---------Mw--------\n");
      if (ii==1)
	fprintf(stderr,"---------sin(ISO)--------\n");
      if (ii==2)
	fprintf(stderr,"---------CLVD--------\n");      
      for(m_par = mt[ii].min; m_par<=mt[ii].max; m_par=m_par+mt[ii].dd){
	fprintf(stderr,"%f\n",m_par);
      }
    }

    for (ii=0; ii<3; ii++){
      if (ii==0)
	fprintf(stderr,"---------STK--------\n");
      if (ii==1)
	fprintf(stderr,"---------cos(DIP)--------\n");
      if (ii==2)
	fprintf(stderr,"---------RAK--------\n");      
      for(m_par = grid.x0[ii]; m_par<(grid.x0[ii]+(grid.n[ii])*grid.step[ii]); m_par=m_par+grid.step[ii]){
	fprintf(stderr,"%f\n",m_par);
      }
    }

    best_sol.err = FLT_MAX;

    for(temp[0]=mt[0].min;temp[0]<=mt[0].max;temp[0]=temp[0]+mt[0].dd){
      for(i_iso=0; i_iso<iso_len-1; i_iso++){
	if (iso_len==1)
	  del_iso=0.;
	else
	  del_iso=(sin(mt[1].max*PI/180.0)-sin(mt[1].min*PI/180.0))/(iso_len-1);
	temp[1]=asin(sin(mt[1].min*PI/180.0)+(i_iso*del_iso))*(180.0/PI);
	fprintf(stderr,"temp[1]=%f %d",temp[1],i_iso);
	if ((temp[1]==-90.) || (temp[1]==90.))
	  continue;
	fprintf(stderr,"-----------------------------------------------\n");
	for(temp[2]=mt[2].min;temp[2]<=mt[2].max;temp[2]=temp[2]+mt[2].dd)

	  //--------newly added section ends here-------------
       
       
	  {  // the base case: grid-search for strike, dip, and rake =============
        best_err = FLT_MAX; /* used when misfit_on_lune=1 */
	    amp = pow(10.,1.5*temp[0]+16.1-20);
	    grd_err = grid.err;
	    for(i_rak=0; i_rak<grid.n[2]; i_rak++) {
	      sol.meca.rak=grid.x0[2]+i_rak*grid.step[2];
	      for(i_dip=0; i_dip<grid.n[1]; i_dip++) {
		if (grid.n[1]==1)
		  del_dip=0.;
		else
		  del_dip=(cos(grid.x0[1]*PI/180.0)-cos((grid.x0[1]+(grid.n[1]-1)*grid.step[1])*PI/180.0))/(grid.n[1]-1);
		sol.meca.dip=acos(cos(grid.x0[1]*PI/180.0)-(i_dip*del_dip))*(180.0/PI);   //dip from -1 to 1
		if (sol.meca.dip==0.)
		  continue;
		//sol.meca.dip=grid.x0[1]+i_dip*grid.step[1];
		for(i_stk=0; i_stk<grid.n[0]; i_stk++) {
		  sol.meca.stk=grid.x0[0]+i_stk*grid.step[1];
		  if (sol.meca.stk == 360.)
		    continue;
		  tt2cmt(temp[2], temp[1], 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);

          // compute misfit from first motion. data will be output to out.misfit_fm.txt
          if(only_first_motion)
          {
              misfit_first_motion(mtensor, nfm, fm, fidfmp, temp[2], temp[1], temp[0], sol.meca.stk, sol.meca.dip, sol.meca.rak);
              continue;
          }

		  if (check_first_motion(mtensor,fm,nfm,fm_thr)<0) {
		    *grd_err++ = sol.err = FLT_MAX;
		    continue;
		  }
		  if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %5.2f %5.2f %5.2f %5.1f %5.1f %5.1f\n", mt[0].par, mt[1].par, mt[2].par, sol.meca.stk, sol.meca.dip, sol.meca.rak);
		  synt=0.;
		  reco=0.;
		  for(obs=obs0,sol.err=0.,i=0;i<nda;i++,obs++){
		    synt1=0.;
		    synt2=0.;
		    st2=0.;
		    err2=0.;
		    st1=0.;
		    err1=0.;
		    mt_radiat(obs->az,mtensor,arad);
		    rad[0]=amp*arad[3][0];
		    for(k=1;k<4;k++) rad[k]=amp*arad[3-k][0];
		    for(k=4;k<6;k++) rad[k]=amp*arad[6-k][2];
		 
		    /*******find the time shift*************/
		    /**SH surface wave**/
		    spt = obs->com;
		    f_pt0 = spt->crl[0];
		    f_pt1 = spt->crl[1];
		    z0 = spt->on_off>0?1:0;
		    /**PSV surface wave**/
		    spt++;
		    r_pt0 = spt->crl[1];
		    r_pt1 = spt->crl[2];
		    r_pt  = spt->crl[3];
		    r_iso = spt->crl[0];
		    z1 = spt->on_off>0?1:0;
		    spt++;
		    z_pt0 = spt->crl[1];
		    z_pt1 = spt->crl[2];
		    z_pt  = spt->crl[3];
		    z_iso = spt->crl[0];
		    z2 = spt->on_off>0?1:0;
		    for(y1=y2=-FLT_MAX,l=0;l<=max_shft[1];l++) {
		      x =rad[4]*(*f_pt0++)+rad[5]*(*f_pt1++);
		      x1=rad[1]*(*r_pt0++)+rad[2]*(*r_pt1++)+rad[3]*(*r_pt++)+rad[0]*(*r_iso++);
		      x2=rad[1]*(*z_pt0++)+rad[2]*(*z_pt1++)+rad[3]*(*z_pt++)+rad[0]*(*z_iso++);
		      y = (1-tie)*z0*x + tie*(z1*x1 + z2*x2);
		      if (y>y2) {y2=y;cfg[0]=x;sol.shft[i][0]=l;}
		      y = tie*z0*x + (1-tie)*(z1*x1 + z2*x2);
		      if (y>y1) {y1=y;cfg[1]=x1;cfg[2]=x2;m=l;}
		    }
		    sol.shft[i][1]=sol.shft[i][2]=m;
		    /**Pnl*/
		    spt++;
		    r_pt0 = spt->crl[1];
		    r_pt1 = spt->crl[2];
		    r_pt  = spt->crl[3];
		    r_iso = spt->crl[0];
		    z1 = spt->on_off>0?1:0;
		    spt++;
		    z_pt0 = spt->crl[1];
		    z_pt1 = spt->crl[2];
		    z_pt  = spt->crl[3];
		    z_iso = spt->crl[0];
		    z2 = spt->on_off>0?1:0;
		    for(y1=-FLT_MAX,l=0;l<=max_shft[3];l++) {
		      x1=rad[1]*(*r_pt0++)+rad[2]*(*r_pt1++)+rad[3]*(*r_pt++)+rad[0]*(*r_iso++);
		      x2=rad[1]*(*z_pt0++)+rad[2]*(*z_pt1++)+rad[3]*(*z_pt++)+rad[0]*(*z_iso++);
		      y = z1*x1 + z2*x2;
		      if (y>y1) {y1=y;cfg[3]=x1;cfg[4]=x2;m=l;}
		    }
		    sol.shft[i][3]=sol.shft[i][4]=m;
		    spt -= NCP - 1;
		 
		    /***error calculation*****/
		    for(kc=2,f_pt1=rad+NRF,j=0;j<NCP;j++,spt++,kc=NRF,f_pt1=rad) {

		      /* compute the L2 norm of syn */
		      for(x2=0.,f_pt0=spt->syn2,k=0;k<kc;k++)
			for(k1=k;k1>=0;k1--)
			  x2+=f_pt1[k]*f_pt1[k1]*(*f_pt0++);
		   
		      //fprintf(stderr, "%s %e\n",obs->stn,x2/Npoints);
		      if (j<3) synt2+=spt->on_off*(x2);
		      else synt1+=spt->on_off*x2;
		      y1 = 1.;
		      /* find out the scaling factor for teleseismic distances */
		      if (obs->tele && spt->on_off) {
			if (cfg[j]>0.) y1 = cfg[j]/x2;
			else y1 = 0.;
		      }
		      sol.scl[i][j] = y1;
		   
		      x1 = (spt->rec2+x2*y1*y1-2.*cfg[j]*y1);
		      //fprintf(stderr,"%e %e %e\n",spt->rec2,x2,cfg[j]);
		      sol.error[i][j] = x1;	/*L2 error for this com.*/
		      sol.cfg[i][j] = 100*cfg[j]/sqrt(spt->rec2*x2);
		      sol.err += spt->on_off*x1;
		      if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %-10s %d %d %9.3e\n", obs->stn, j, spt->on_off, x1);
		      //fprintf(stderr, "%s %e %e %d\n",obs->stn, x2/Npoints, x1,sol.cfg[i][j]);
		      if (j<3){
			st2+=spt->on_off*spt->rec2;  // syn surf wave norm, individual staion
			err2+=spt->on_off*x1;  // err surf wave norm, individual staion
		      }
		      else{
			st1+=spt->on_off*spt->rec2;   // syn body wave norm, individual staion
			err1+=spt->on_off*x1;   // err body wave norm, individual station
		      }	
		    }
		    
		    if (0){
		      st_err1 = fopen(st_err_body,"a");
		      fprintf(st_err1, "%f %e %e %e\n",obs->dist,st1, synt1, err1);
		      fclose(st_err1);
		      st_err2 = fopen(st_err_surf,"a");
		      fprintf(st_err2, "%f %e %e %e\n",obs->dist,st2, synt2, err2);
		      fclose(st_err2);}
		    synt+=synt2+synt1;  // syntheics norm (all stations)
		    reco+=st2+st1;      // data norm, same as data2
		  }
		   /*-------------------------end of all stations*/
		  //fprintf(stderr, "Nsta=%d\n",Nsta);
		  sol.err=sol.err/Nsta;
		  synt = synt/Nsta;
		  reco=reco/Nsta;
		  //fprintf(stderr, "data2=%e synt2=%e\n",data2,synt2);
		  *grd_err++ = sol.err;		/*error for this solution*/
          if(misfit_on_lune)
          {
              // remember smallest misfit for a given gamma/delta on lune
              if(best_err > sol.err)
              {
                  best_err = sol.err;
                  best_stk = sol.meca.stk;
                  best_dip = sol.meca.dip;
                  best_rak = sol.meca.rak;
              }
          }
		  if (best_sol.err>sol.err) {best_sol = sol; synth=synt; //fprintf(stderr, "synt=%e data2=%e err=%e\n",synth,data2,sol.err);
		    if (0) fprintf(stderr,"misfit for best sol = %f; stk=%3.1f, dip=%3.1f, rak=%3.1f \n",best_sol.err,sol.meca.stk, sol.meca.dip, sol.meca.rak);
		    mt[0].par=temp[0];
		    mt[1].par=temp[1];
		    mt[2].par=temp[2];
		    //fprintf(stderr, "data=%e syn=%e err=%e\n",data2, synt2, sol.err);
		  }
	       
		  if (debug) { 
		    log = fopen(logfile,"a");                 // output log file
		    fprintf(log,"%3.1f\t%3.1f\t%3.1f\t%e\t%2.2f\t%2.2f\t%2.2f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n",sol.meca.stk, sol.meca.dip, sol.meca.rak, sol.err/data2, temp[0], temp[1], temp[2], amp*1.0e20, mtensor[0][0], mtensor[0][1], mtensor[0][2], mtensor[1][1], mtensor[1][2], mtensor[2][2] );
		    //fprintf(stderr,"%3.1f\t%3.1f\t%3.1f\t%e\t%2.2f\t%2.2f\t%2.2f\n",sol.meca.stk, sol.meca.dip, sol.meca.rak, sol.err,  temp[0], temp[1], temp[2]);
		    fclose(log);	  
		  }
		} // stk loop
	      }   // dip loop
	    
	      if (sol.meca.stk==(grid.x0[0]+(grid.n[0]-1)*grid.step[0]) &&  sol.meca.dip==(grid.x0[1]+(grid.n[1]-1)*grid.step[1]) &&  sol.meca.rak==(grid.x0[2]+(grid.n[2]-1)*grid.step[2])){
		loop++;
		if (debug) {
		  sprintf(logfile,"%s_%03d","log",loop);  // changes the log file name for next sext search (for multiple log files - search over stk,dip and rake only)
		  log = fopen(logfile,"a");
		  fclose(log);
		}
		log = fopen("log_diff","a"); /*fprintf(stderr,"completed stk,dip,rake loop\n");        //summary log file*/
		fprintf(log,"%d\t%d\t%3.1f\t%3.1f\t%3.1f\t%f\t%2.2f\t%2.2f\t%2.2f\n",loop,interp, best_sol.meca.stk, best_sol.meca.dip, best_sol.meca.rak, best_sol.err, mt[0].par, mt[1].par, mt[2].par);
		fclose(log);
	      }
	    }  // rak loop
	    if (1) fprintf(stderr,"Mw=%2.1f \t iso=%2.2f \t clvd=%2.2f\n",temp[0],temp[1],temp[2]);

        if(misfit_on_lune)
        {
            fprintf(fidmol,"%5.1f %5.1f %5.1f %5.1f %5.1f %e %3.1f\n",
                    temp[1], temp[2], best_stk, best_dip, best_rak, best_err, temp[0]);
        }
	  }
      }
    }

    if(misfit_on_lune)
    {
        fclose(fidmol);
        fprintf(stdout, "\nFinished writing misfit to file out.misfit.wf_\n");
    }
 
    if(only_first_motion)
    {
        fclose(fidfmp);
        fprintf(stderr,"\nFinished writing to file out.misfit.fmp_\n");
        fprintf(stderr,"No figure should be created (no -P flag) in this mode.\n");
    }
 
    if (debug) fprintf(stderr, "Mw=%5.2f  iso=%5.2f clvd=%5.2f misfit = %9.3e\n", mt[0].par, mt[1].par, mt[2].par, best_sol.err);
    if (interp=1) return(best_sol);
    /* do interpolation */
    best_sol.err = grid3d(grid.err,&(grid.n[0]),s3d,&(best_sol.flag),&(best_sol.ms),best_sol.others);
    if (debug) fprintf(stderr, " interpolation  misfit = %9.3e\n", best_sol.err);
    best_sol.meca.stk = grid.x0[0]+s3d[0]*grid.step[0];
    best_sol.meca.dip = grid.x0[1]+s3d[1]*grid.step[1];
    best_sol.meca.rak = grid.x0[2]+s3d[2]*grid.step[2];
    for(i=0;i<3;i++) best_sol.dev[i]  = s3d[3+i]/(grid.step[i]*grid.step[i]);
    best_sol.dev[3] = s3d[6]/(grid.step[0]*grid.step[1]);
    best_sol.dev[4] = s3d[7]/(grid.step[0]*grid.step[2]);
    best_sol.dev[5] = s3d[8]/(grid.step[1]*grid.step[2]);
    fprintf(stderr,"=======================");
    return(best_sol);
 
  }
 
  else{
    if ( npar ) {	// line-search for mw, iso, and clvd ================================

      npar--;
      dx = mt[npar].dd;
      i = 1; if (dx>0.001) i = 0;
      sol = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,i,bootstrap,search);
      if (dx>0.001) {	/* do line search */
	mt[npar].par += dx;
	sol2 = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,0,bootstrap,search);
	if (sol2.err > sol.err) {	/* this is the wrong direction, turn around */
	  dx = -dx;
	  sol1 = sol2; sol2 = sol; sol  = sol1; /*swap sol, sol2 */
	  mt[npar].par += dx;
	}
	while(sol2.err < sol.err) {	/* keep going until passing by the mininum */
	  sol1 = sol;
	  sol = sol2;
	  mt[npar].par += dx;
	  if (mt[npar].par>mt[npar].max || mt[npar].par<mt[npar].min) sol2.err = sol1.err;
	  else sol2 = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,0,bootstrap,search);
	}
	mt[npar].sigma = 2*dx*dx/(sol2.err+sol1.err-2*sol.err);
	mt[npar].par -= dx+0.5*dx*(sol2.err-sol1.err)/(sol2.err+sol1.err-2*sol.err);
	sol = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,1,bootstrap,search);
      } else {
	mt[npar].sigma = 0.;
      }
      return(sol);
    } 
    else {	// the base case: grid-search for strike, dip, and rake =============
      amp = pow(10.,1.5*mt[0].par+16.1-20);
      best_sol.err = FLT_MAX;
      grd_err = grid.err;
      for(i_rak=0; i_rak<grid.n[2]; i_rak++) {
	sol.meca.rak=grid.x0[2]+i_rak*grid.step[2];
	for(i_dip=0; i_dip<grid.n[1]; i_dip++) {
	  sol.meca.dip=grid.x0[1]+i_dip*grid.step[1];
	  for(i_stk=0; i_stk<grid.n[0]; i_stk++) {
	    sol.meca.stk=grid.x0[0]+i_stk*grid.step[1];
	    nmtensor(mt[1].par,mt[2].par,sol.meca.stk,sol.meca.dip,sol.meca.rak,mtensor);
	    if (check_first_motion(mtensor,fm,nfm,fm_thr)<0) {
	      *grd_err++ = sol.err = FLT_MAX;
	      continue;
	    }
	    if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %5.2f %5.2f %5.2f %5.1f %5.1f %5.1f\n", mt[0].par, mt[1].par, mt[2].par, sol.meca.stk, sol.meca.dip, sol.meca.rak);
	    synt=0.;
	    reco=0.;
	    for(obs=obs0,sol.err=0.,i=0;i<nda;i++,obs++){
	      synt1=0.;
	      synt2=0.;
	      st2=0.;
	      err2=0.;
	      st1=0.;
	      err1=0.;
	      mt_radiat(obs->az,mtensor,arad);    // arad are green's function coefficient (see Wang & Hermann 1985)
	      rad[0]=amp*arad[3][0];
	      for(k=1;k<4;k++) rad[k]=amp*arad[3-k][0];
	      for(k=4;k<6;k++) rad[k]=amp*arad[6-k][2];
	    
	      /*******find the time shift*************/
	      /**SH surface wave**/
	      spt = obs->com;
	      f_pt0 = spt->crl[0];
	      f_pt1 = spt->crl[1];
	      z0 = spt->on_off>0?1:0;
	      /**PSV surface wave**/
	      spt++;
	      r_pt0 = spt->crl[1];
	      r_pt1 = spt->crl[2];
	      r_pt  = spt->crl[3];
	      r_iso = spt->crl[0];
	      z1 = spt->on_off>0?1:0;
	      spt++;
	      z_pt0 = spt->crl[1];
	      z_pt1 = spt->crl[2];
	      z_pt  = spt->crl[3];
	      z_iso = spt->crl[0];
	      z2 = spt->on_off>0?1:0;
	      for(y1=y2=-FLT_MAX,l=0;l<=max_shft[1];l++) {
		x =rad[4]*(*f_pt0++)+rad[5]*(*f_pt1++);
		x1=rad[1]*(*r_pt0++)+rad[2]*(*r_pt1++)+rad[3]*(*r_pt++)+rad[0]*(*r_iso++);
		x2=rad[1]*(*z_pt0++)+rad[2]*(*z_pt1++)+rad[3]*(*z_pt++)+rad[0]*(*z_iso++);
		y = (1-tie)*z0*x + tie*(z1*x1 + z2*x2);
		if (y>y2) {y2=y;cfg[0]=x;sol.shft[i][0]=l;}
		y = tie*z0*x + (1-tie)*(z1*x1 + z2*x2);
		if (y>y1) {y1=y;cfg[1]=x1;cfg[2]=x2;m=l;}
	      }
	      sol.shft[i][1]=sol.shft[i][2]=m;
	      /**Pnl*/
	      spt++;
	      r_pt0 = spt->crl[1];
	      r_pt1 = spt->crl[2];
	      r_pt  = spt->crl[3];
	      r_iso = spt->crl[0];
	      z1 = spt->on_off>0?1:0;
	      spt++;
	      z_pt0 = spt->crl[1];
	      z_pt1 = spt->crl[2];
	      z_pt  = spt->crl[3];
	      z_iso = spt->crl[0];
	      z2 = spt->on_off>0?1:0;
	      for(y1=-FLT_MAX,l=0;l<=max_shft[3];l++) {
		x1=rad[1]*(*r_pt0++)+rad[2]*(*r_pt1++)+rad[3]*(*r_pt++)+rad[0]*(*r_iso++);
		x2=rad[1]*(*z_pt0++)+rad[2]*(*z_pt1++)+rad[3]*(*z_pt++)+rad[0]*(*z_iso++);
		y = z1*x1 + z2*x2;
		if (y>y1) {y1=y;cfg[3]=x1;cfg[4]=x2;m=l;}
	      }
	      sol.shft[i][3]=sol.shft[i][4]=m;
	      spt -= NCP - 1;
	    
	      /***error calculation*****/
	      for(kc=2,f_pt1=rad+NRF,j=0;j<NCP;j++,spt++,kc=NRF,f_pt1=rad) {

		/* compute the L2 norm of syn */
		for(x2=0.,f_pt0=spt->syn2,k=0;k<kc;k++)
		  for(k1=k;k1>=0;k1--)
		    x2+=f_pt1[k]*f_pt1[k1]*(*f_pt0++);

		if (j<3) synt2+=spt->on_off*(x2);
		else synt1+=spt->on_off*x2;
		y1 = 1.;
		/* find out the scaling factor for teleseismic distances */
		if (obs->tele && spt->on_off) {
		  if (cfg[j]>0.) y1 = cfg[j]/x2;
		  else y1 = 0.;
		}
		sol.scl[i][j] = y1;

		x1 = (spt->rec2+x2*y1*y1-2.*cfg[j]*y1);
		sol.error[i][j] = x1;	/*L2 error for this com.*/
		sol.cfg[i][j] = 100*cfg[j]/sqrt(spt->rec2*x2);
		sol.err += spt->on_off*x1;
		if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %-10s %d %d %9.3e\n", obs->stn, j, spt->on_off, x1);
		if (j<3){
		  st2+=spt->on_off*spt->rec2;
		  err2+=spt->on_off*x1;
		}
		else{
		  st1+=spt->on_off*spt->rec2;
		  err1+=spt->on_off*x1;
		}
	      }
	      if (0){
		st_err1 = fopen(st_err_body,"a");
		fprintf(st_err1, "%f %e %e %e\n",obs->dist,st1, synt1, err1);
		fclose(st_err1);
		st_err2 = fopen(st_err_surf,"a");
		fprintf(st_err2, "%f %e %e %e\n",obs->dist,st2, synt2, err2);
		fclose(st_err2);}
	      synt+=synt2+synt1;  // syntheics norm (all stations)
	      reco+=st2+st1;      // data norm, same as data2
	    } /*-------------------------end of all stations*/
	    //fprintf(stderr, "Nsta=%d\n",Nsta);
	    sol.err=sol.err/Nsta;
	    synt = synt/Nsta;
	    reco = reco/Nsta;
	    *grd_err++ = sol.err;		/*error for this solution*/

	    if (best_sol.err>sol.err) {
	      best_sol=sol;
	      synth=synt;}

	    if (debug) { 
	      log = fopen(logfile,"a");
	      fprintf(log,"%3.1f\t%3.1f\t%3.1f\t%e\t%2.1f\t%2.2f\t%2.2f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n",sol.meca.stk, sol.meca.dip, sol.meca.rak, sol.err/data2, mt[0].par, mt[1].par, mt[2].par, amp*1.0e20, mtensor[0][0], mtensor[0][1], mtensor[0][2], mtensor[1][1], mtensor[1][2], mtensor[2][2] );
	      fclose(log);	  
	    }
	  } // loop for stk
	}   // loop for dip
	if (sol.meca.stk==360 &&  sol.meca.dip==90 &&  sol.meca.rak==90){
	  loop++;
	  if (debug) {
	    sprintf(logfile,"%s_%03d","log",loop);
	    log = fopen(logfile,"a");
	    fclose(log);
	  }
	  log = fopen("log_diff","a");
	  fprintf(log,"%d\t%d\t%3.1f\t%3.1f\t%3.1f\t%f\t%2.2f\t%2.2f\t%2.2f\n",loop,interp, best_sol.meca.stk, best_sol.meca.dip, best_sol.meca.rak, best_sol.err, mt[0].par, mt[1].par, mt[2].par);
	  fclose(log);
	}
      }
  
      if (debug) fprintf(stderr, "Mw=%5.2f  iso=%5.2f clvd=%5.2f misfit = %9.3e\n", mt[0].par, mt[1].par, mt[2].par, best_sol.err);
      if (interp == 0) return(best_sol);
      /* do interpolation */
      best_sol.err = grid3d(grid.err,&(grid.n[0]),s3d,&(best_sol.flag),&(best_sol.ms),best_sol.others);
      if (debug) fprintf(stderr, " interpolation  misfit = %9.3e\n", best_sol.err);
      best_sol.meca.stk = grid.x0[0]+s3d[0]*grid.step[0];
      best_sol.meca.dip = grid.x0[1]+s3d[1]*grid.step[1];
      best_sol.meca.rak = grid.x0[2]+s3d[2]*grid.step[2];
      for(i=0;i<3;i++) best_sol.dev[i]  = s3d[3+i]/(grid.step[i]*grid.step[i]);
      best_sol.dev[3] = s3d[6]/(grid.step[0]*grid.step[1]);
      best_sol.dev[4] = s3d[7]/(grid.step[0]*grid.step[2]);
      best_sol.dev[5] = s3d[8]/(grid.step[1]*grid.step[2]);
      return(best_sol);
    }
  }
}




