#include "cap.h"

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
		int             search,
		int             norm
		) {
  int	i, j, k, l, m, k1, kc, z0, z1, z2, ii, N, iso_len;
  float mw_ran; // half-range for magnitude search (previously int)
  int	i_stk, i_dip, i_rak, i_iso;
  int   perc, del_N, count_perc;
  float	amp, rad[6], arad[4][3], x, x1, x2, y, y1, y2, cfg[NCP], s3d[9], temp[3], m_par, del_dip, del_iso;
  float	*f_pt0, *f_pt1, *r_pt, *r_pt0, *r_pt1, *z_pt, *z_pt0, *z_pt1, *grd_err, *rnd_stk, *rnd_dip, *rnd_rak, *rnd_iso, *rnd_clvd, *iso;
  float dx, mtensor[3][3], *r_iso, *z_iso;
  DATA	*obs;
  COMP	*spt;
  SOLN	sol, sol1, sol2, best_sol;
  FILE *logf,*std_range;
  char logfile[16],range_file[20];

  OUTPUTGD outgd;
  OUTPUTMT outmt;

  int misfit_fmp;

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

  if (debug) {
    sprintf(logfile,"%s_%03d_%03d","log",edep,loop);
    logf = fopen(logfile,"w");
    fclose(logf);
  }

  if (debug) fprintf(stderr, "loop=%d start=%d \n",loop,start);
  start++;

  // search range file
  sprintf(range_file,"search_range");

  // RANDOM distribution of samples in MT space 
  // Distribution is homogeneous on lune 
  if (search==2){

    // output for postprocessing
    fidmt=fopen("capout_rand_mt.bin","wb");
    fidgd=fopen("capout_rand_gd.bin","wb");

    fprintf(stderr,"Mw=%f\n",mt[0].par);
    mw_ran=1.0; // Mw search range
    mt[0].max = mt[0].par+mw_ran;
    mt[0].min = mt[0].par-mw_ran;
 
    // if Mw search increment in -I flag is set to 0 then make Mw<min max> search range equal
    if (mt[0].dd==0){
      mt[0].max = mt[0].par;
      mt[0].min = mt[0].par;
      mt[0].dd=1.0;
    }

    N=100000; // Number of samples 
    del_N=(int)N/20; // To add 'percentage completed' marker; each del_N is 5%
    count_perc=1;  // counter to track percentage completed
    // preallocate arrays for MT paramters 
    rnd_stk = (float*)malloc(sizeof(int) * N*sizeof(float));   
    rnd_dip = (float*)malloc(sizeof(int) * N*sizeof(float));
    rnd_rak = (float*)malloc(sizeof(int) * N*sizeof(float));
    rnd_iso = (float*)malloc(sizeof(int) * N*sizeof(float));
    rnd_clvd = (float*)malloc(sizeof(int) * N*sizeof(float));
    
    // If space cannot be allocated 
    if (rnd_stk==NULL || rnd_stk==NULL || rnd_stk==NULL){
      fprintf(stderr,"Cannot allocate space for random number generation");}

    // Generate homogeneously distributed samples (Lune parameterization)
    // XXX same random samples are generated in every run (change drand (to srand) or add seed)
    // drand return values bteween 0 and 1
    for(i=0;i<N;i++){
      rnd_stk[i]=0.0+360.0*drand48();
      rnd_dip[i]=0.0+1.0*drand48();     // because cos(dip) is homogeneous 
      rnd_rak[i]=-90.0+180.0*drand48();
      if (mt[1].dd==0) // for direct search at a particular CLVD
	rnd_iso[i]=sin(mt[1].par*(PI/180.0));
      else
	rnd_iso[i]=-1.0+2.0*drand48();  // because sin(iso) is homogeneous 
      if (mt[2].dd==0) // for direct search at a particular ISO
	rnd_clvd[i]=mt[2].par;
      else
	rnd_clvd[i]=-30.0+60.0*drand48();
    }

    best_sol.err = FLT_MAX;
    
    /* mt[0] = Mw    (temp[0] searches for best Mw)
       mt[1] = ISO   (temp[1] searches for best ISO)
       mt[2] = CLVD  (temp[2] searches for best CLVD)  */ 
     for(temp[0]=mt[0].min;temp[0]<=mt[0].max;temp[0]=temp[0]+mt[0].dd){ 

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
	    tt2cmt(temp[2], temp[1], 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);   // get normalized mtensor from (CLVD,ISO,strike,dip, rake)
	    if (check_first_motion(mtensor,fm,nfm,fm_thr)<0) {
	      *grd_err++ = sol.err = FLT_MAX;
	      continue;
	    }
	    // This gets executed only when bootstrapping is performed
	    if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %5.2f %5.2f %5.2f %5.1f %5.1f %5.1f\n", mt[0].par, mt[1].par, mt[2].par, sol.meca.stk, sol.meca.dip, sol.meca.rak);
	   
	    //--------------KEY COMMAND---call misfit function------
	    sol=calerr(nda,obs0,max_shft,tie,norm,mtensor,amp,sol);

	    sol.err=sol.err/Nsta;               // normalize error by number of station
	    *grd_err++ = sol.err;		// error for this solution
	    
	    // save the sample if it has the misfit lower than previous
	    if (best_sol.err>sol.err) {best_sol = sol;
	      if (debug) fprintf(stderr,"misfit for best sol = %f; stk=%3.1f, dip=%3.1f, rak=%3.1f \n",best_sol.err,sol.meca.stk, sol.meca.dip, sol.meca.rak); // output on screen if in debug mode
	      mt[0].par=temp[0];
	      mt[1].par=temp[1];
	      mt[2].par=temp[2];
	    }

	    misfit_fmp =  misfit_first_motion(mtensor, nfm, fm, fidfmp, temp[2], temp[1], temp[0], sol.meca.stk, sol.meca.dip, sol.meca.rak);
	    //  output binary data
	    outgd.g = temp[2];
	    outgd.d = temp[1];
	    outgd.s = sol.meca.stk;
	    outgd.h = sol.meca.dip;
	    outgd.r = sol.meca.rak;
	    outgd.mag = temp[0];
	    outgd.misfit_wf  = sol.err/data2;
	    outgd.misfit_fmp = (float) misfit_fmp;

	    outmt.mrr = mtensor[2][2];
	    outmt.mtt = mtensor[0][0];
	    outmt.mpp = mtensor[1][1];
	    outmt.mrt = mtensor[0][2];
	    outmt.mrp = -mtensor[1][2];
	    outmt.mtp = -mtensor[0][1];
	    outmt.mag = temp[0];
	    outmt.misfit_wf = sol.err/data2;
	    outmt.misfit_fmp = (float) misfit_fmp;

	    //          fprintf(stderr,"DEBUG. %f %f %d \n", outgd.g, outgd.d, outmt.misfit_fmp);
	    fwrite(&outmt, sizeof outmt, 1, fidmt); 
	    fwrite(&outgd, sizeof outgd, 1, fidgd);

	    // in debug mode it will save all samples in a logfile
	    if (debug) { 
	      logf = fopen(logfile,"a");                 // output log file
	      fprintf(logf,"%3.1f\t%3.1f\t%3.1f\t%e\t%2.2f\t%2.2f\t%2.2f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n",sol.meca.stk, sol.meca.dip, sol.meca.rak, sol.err/data2, temp[0], temp[1], temp[2], amp*1.0e20, mtensor[0][0], mtensor[1][1], mtensor[2][2], mtensor[0][1], mtensor[0][2], mtensor[1][2] );
	      fclose(logf);
	    }
	    // Block to keep track of percentage of samples checked
	    if (ii==del_N*count_perc){
	      perc = 5*count_perc; // percentage is in multiples of 5
	      printf("%d%% done \n",perc);
	      count_perc = count_perc+1;
	    }
	  }  // number of samples (ii < N) loop ends

	  // change the log filename for every depth
	  if (debug) {
	    loop++;
	    sprintf(logfile,"%s_%03d_%03d","log",edep,loop);  // changes the log file name for next sext search
	    logf = fopen(logfile,"a");
	    fclose(logf);
	  }
	  
	  if (1){  // This saves only those samples which had misfit lower than previous sample
	    logf = fopen("log_diff","a"); /*fprintf(stderr,"completed stk,dip,rake loop\n");        //summary log file*/
	    fprintf(logf,"%d\t%d\t%3.1f\t%3.1f\t%3.1f\t%f\t%2.2f\t%2.2f\t%2.2f\n",loop,interp, best_sol.meca.stk, best_sol.meca.dip, best_sol.meca.rak, best_sol.err, mt[0].par, mt[1].par, mt[2].par);
	    fclose(logf);
	  }
	  fprintf(stderr, "%d\t%3.2f\t%3.2f\t%3.2f\t%2.1f\t%2.1f\t%2.2f\n",ii+1,sol.meca.stk, sol.meca.dip,sol.meca.rak,temp[0],temp[1],temp[2]);
	  
	  fprintf(stderr,"========================Minimum==================================\n"); // output the minimum misfit sample for every magnitude (Mw) searched
     fprintf(stderr, "%3.2f\t%3.2f\t%3.2f\t%2.1f\t%2.1f\t%2.2f\n",best_sol.meca.stk, best_sol.meca.dip,best_sol.meca.rak,temp[0],mt[1].par,mt[2].par);
     } // magnitude search (mt[0]) loop ends
  
    // end writing data for postprocessing
    fclose(fidmt);
    fclose(fidgd);

    if (debug) fprintf(stderr, "Mw=%5.2f  iso=%5.2f clvd=%5.2f misfit = %9.3e\n", mt[0].par, mt[1].par, mt[2].par, best_sol.err);
    if (interp==0) return(best_sol);
    /* do interpolation - But this is not implemented in random search yet - a gateway towards neighbourhood algorithm*/
    // Code does not reach this part - return(best_sol) is happening before
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

    // output for postprocessing
    fidmt=fopen("capout_grid_mt.bin","wb");
    fidgd=fopen("capout_grid_gd.bin","wb");

    //--------newly added section-------------
    mw_ran = 0.5;
    mt[0].max = mt[0].par+mw_ran;
    mt[0].min = mt[0].par-mw_ran;

    // reset magnitude range to so that cap runs only once in this mode
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

    std_range = fopen(range_file,"w");
    //Output search ranges
    fprintf(stderr,"=========GRID-SEARCH RANGE===========\n");
    for (ii=0; ii<3; ii=ii+2){
      if (ii==0)
	fprintf(std_range,"---------Mw--------\n");
      if (ii==2)
	fprintf(std_range,"---------CLVD--------\n");      
      for(m_par = mt[ii].min; m_par<=mt[ii].max; m_par=m_par+mt[ii].dd){
	fprintf(std_range,"%f\n",m_par);
      }
    }
    fprintf(std_range,"---------ISO--------\n");
    for(i_iso=0; i_iso<iso_len; i_iso++){
      if (iso_len==1)
	  del_iso=0.;
	else
	  del_iso=(sin(mt[1].max*PI/180.0)-sin(mt[1].min*PI/180.0))/(iso_len-1);
	temp[1]=asin(sin(mt[1].min*PI/180.0)+(i_iso*del_iso))*(180.0/PI);
	if (temp[1]==-90. || temp[1]==90. || temp[1] != temp[1])
	  continue;
	fprintf(std_range,"%f\t%f\t%f\n",(((float)i_iso)*mt[1].dd)+mt[1].min,sin(temp[1]*PI/180.),temp[1]);
    }

    for (ii=0; ii<3; ii=ii+2){
      if (ii==0)
	fprintf(std_range,"---------STK--------\n");
      if (ii==2)
	fprintf(std_range,"---------RAK--------\n");      
      for(m_par = grid.x0[ii]; m_par<(grid.x0[ii]+(grid.n[ii])*grid.step[ii]); m_par=m_par+grid.step[ii]){
	fprintf(std_range,"%f\n",m_par);
      }
    }
    fprintf(std_range,"---------DIP--------\n");
    for(i_dip=0; i_dip<grid.n[1]; i_dip++) {
      if (grid.n[1]==1)
	del_dip=0.;
      else
	del_dip=(cos(grid.x0[1]*PI/180.0)-cos((grid.x0[1]+(grid.n[1]-1)*grid.step[1])*PI/180.0))/(grid.n[1]-1); // equal increment in cosine(dip) space -> del_dip is in degrees, i.e. not equally spaced
      sol.meca.dip=acos(cos(grid.x0[1]*PI/180.0)-(i_dip*del_dip))*(180.0/PI);   //dip from -1 to 1     
      if (sol.meca.dip==0. || sol.meca.dip>90.)
      	continue;
      fprintf(std_range,"%f\t%f\t%f\n",((float)i_dip+1.)*grid.step[2],cos(sol.meca.dip*PI/180.),sol.meca.dip);
    }
    fclose(std_range);

    best_sol.err = FLT_MAX;

    for(temp[0]=mt[0].min;temp[0]<=mt[0].max;temp[0]=temp[0]+mt[0].dd){
      for(i_iso=0;i_iso<iso_len;i_iso++){
	if (iso_len==1)
	  del_iso=0.;
	else
	  del_iso=(sin(mt[1].max*PI/180.0)-sin(mt[1].min*PI/180.0))/(iso_len-1);
	temp[1]=asin(sin(mt[1].min*PI/180.0)+(i_iso*del_iso))*(180.0/PI);
	if (temp[1]==-90. || temp[1]==90. || temp[1] != temp[1])    // Do not include the limits, or if ISO is NaN; temp[1]!=temp[1] only if temp[1] is NaN (if sin(theta)>1)
	  continue;
	fprintf(stderr,"-----------------------------------------------\n");
	for(temp[2]=mt[2].min;temp[2]<=mt[2].max;temp[2]=temp[2]+mt[2].dd)

	  //--------newly added section ends here-------------
       
       	  {  // the base case: grid-search for strike, dip, and rake =============

        /*  variables to track smallest misfit at each point on the lune */
        /*  used only when misfit_on_lune=1 */
        bestmisfit->misfit = FLT_MAX;
        bestmisfit->gamma  = NAN;
        bestmisfit->delta  = NAN;
        bestmisfit->mrr = NAN;
        bestmisfit->mtt = NAN;
        bestmisfit->mpp = NAN;
        bestmisfit->mrt = NAN;
        bestmisfit->mrp = NAN;
        bestmisfit->mtp = NAN;
        bestmisfit->mag = NAN;
        bestmisfit->stk = NAN;
        bestmisfit->dip = NAN;
        bestmisfit->rak = NAN;

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
		if (sol.meca.dip==0. || sol.meca.dip>90.)
		  continue;
		//sol.meca.dip=grid.x0[1]+i_dip*grid.step[1];
		for(i_stk=0; i_stk<grid.n[0]; i_stk++) {
		  sol.meca.stk=grid.x0[0]+i_stk*grid.step[1];
		  if (sol.meca.stk == 360.)
		    continue;
		  tt2cmt(temp[2], temp[1], 1.0, sol.meca.stk, sol.meca.dip, sol.meca.rak, mtensor);

          // compute misfit from first motion. data will be output to out.misfit.fmp_
           misfit_fmp =  misfit_first_motion(mtensor, nfm, fm, fidfmp, temp[2], temp[1], temp[0], sol.meca.stk, sol.meca.dip, sol.meca.rak);

		  if (bootstrap && interp==0) fprintf(stderr,"BOOTSTRAPPING %5.2f %5.2f %5.2f %5.1f %5.1f %5.1f\n", mt[0].par, mt[1].par, mt[2].par, sol.meca.stk, sol.meca.dip, sol.meca.rak);

		  //--------------KEY COMMAND---call misfit function------
		  sol=calerr(nda,obs0,max_shft,tie,norm,mtensor,amp,sol);

		  //fprintf(stderr, "Nsta=%d\n",Nsta);
		  sol.err=sol.err/Nsta;
		  *grd_err++ = sol.err;		/*error for this solution*/

          /* track smallest misfit at each point on the lune */
          if(misfit_on_lune)
          {
              if(sol.err < bestmisfit->misfit)
              {
                  bestmisfit->gamma = temp[2];
                  bestmisfit->delta = temp[1];
                  bestmisfit->stk = sol.meca.stk;
                  bestmisfit->dip = sol.meca.dip;
                  bestmisfit->rak = sol.meca.rak;
                  bestmisfit->misfit = sol.err;
                  bestmisfit->mag = temp[0];

                  /* GCMT format */
                  bestmisfit->mrr = mtensor[2][2];
                  bestmisfit->mtt = mtensor[0][0];
                  bestmisfit->mpp = mtensor[1][1];
                  bestmisfit->mrt = mtensor[0][2];
                  bestmisfit->mrp = -mtensor[1][2];
                  bestmisfit->mtp = -mtensor[0][1];
              }
          }

		  if (best_sol.err>sol.err) {best_sol = sol; 
		    mt[0].par=temp[0];
		    mt[1].par=temp[1];
		    mt[2].par=temp[2];
		  }
	       
		  //  output binary data
          outgd.g = temp[2];
          outgd.d = temp[1];
          outgd.s = sol.meca.stk;
          outgd.h = sol.meca.dip;
          outgd.r = sol.meca.rak;
          outgd.mag = temp[0];
          outgd.misfit_wf  = sol.err/data2;
          outgd.misfit_fmp = (float) misfit_fmp;

          outmt.mrr = mtensor[2][2];
          outmt.mtt = mtensor[0][0];
          outmt.mpp = mtensor[1][1];
          outmt.mrt = mtensor[0][2];
          outmt.mrp = -mtensor[1][2];
          outmt.mtp = -mtensor[0][1];
          outmt.mag = temp[0];
          outmt.misfit_wf = sol.err/data2;
          outmt.misfit_fmp = (float) misfit_fmp;

//          fprintf(stderr,"DEBUG. %f %f %d \n", outgd.g, outgd.d, outmt.misfit_fmp);
          fwrite(&outmt, sizeof outmt, 1, fidmt); 
          fwrite(&outgd, sizeof outgd, 1, fidgd);

          // count number of solutions computed
          // (this is ongoing, this variable may be deleted)
          // # 20151216 celso alvizuri - cralvizuri@alaska.edu
          count++;

		}   /* end stk loop */
	      } /* end dip loop */
	    
          // this part checks floats (not reliable)... Any suggestions?
          if (sol.meca.stk==(grid.x0[0]+(grid.n[0]-1)*grid.step[0]) &&  sol.meca.dip==(grid.x0[1]+(grid.n[1]-1)*grid.step[1]) &&  sol.meca.rak==(grid.x0[2]+(grid.n[2]-1)*grid.step[2])){
              loop++;

	      // changes the log file name for next sext search (for multiple log files - search over stk,dip and rake only)
              if (debug) {
                  sprintf(logfile,"%s_%03d_%03d","log",edep,loop);
                  logf = fopen(logfile,"w");
                  fclose(logf);
              }
              
              logf = fopen("log_diff","a"); //fprintf(stderr,"completed stk,dip,rake loop\n");        //summary log file - Minimum after each Mw search (and all orientations within it)
              fprintf(logf,"%d\t%d\t%3.1f\t%3.1f\t%3.1f\t%f\t%2.2f\t%2.2f\t%2.2f\n",loop,interp, best_sol.meca.stk, best_sol.meca.dip, best_sol.meca.rak, best_sol.err, mt[0].par, mt[1].par, mt[2].par);
              fclose(logf);
          } // END XYZ TEST
          
	    }  /* end rake loop */

        // output values during grid search
	    fprintf(stderr,"Mw=%2.1f \t iso=%2.2f \t clvd=%2.2f\n",temp[0],temp[1],temp[2]);

        /* output smallest misfit at each (gamma,delta) on the lune */
        if(misfit_on_lune)
        {
            fprintf(fidmol,"%6.2f %6.2f %6.2f %6.2f %6.2f %9.6e %9.6e %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %4.1f\n",
                    bestmisfit->gamma, bestmisfit->delta, bestmisfit->stk, bestmisfit->dip, bestmisfit->rak,
                    bestmisfit->misfit, 100.0*(1.-(bestmisfit->misfit/data2)*(bestmisfit->misfit/data2)),
                    bestmisfit->mrr, bestmisfit->mtt, bestmisfit->mpp,
                    bestmisfit->mrt, bestmisfit->mrp, bestmisfit->mtp,
                    bestmisfit->mag);
        }
	  } /* end clvd loop */
      } /* end iso loop */
    }   /* end mag loop */

  // count number of solutions computed
  // (this is ongoing, this variable may be deleted)
  // # 20151216 celso alvizuri - cralvizuri@alaska.edu
    fprintf(stderr, "\ntotal solutions processed:\n%d\n", count);

    if(misfit_on_lune)
    {
        fclose(fidmol);
        free(bestmisfit);
        fprintf(stdout, "\nFinished writing waveform misfit to file out.misfit.wf_\n");
    }
 
    if(only_first_motion)
    {
        fclose(fidfmp);
        fprintf(stderr,"\nFinished writing to file out.misfit.fmp_\n");
        fprintf(stderr,"No figure should be created (no -P flag) in this mode.\n");
    }

    // end writing data for postprocessing
    fclose(fidmt);
    fclose(fidgd);
 
    if (debug) fprintf(stderr, "Mw=%5.2f  iso=%5.2f clvd=%5.2f misfit = %9.3e\n", mt[0].par, mt[1].par, mt[2].par, best_sol.err);
    if (interp==0) return(best_sol);
    /* do interpolation */
    // Code does not reach here - return(best_sol) is happening before
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
 
  }     // end loop for option: search=1
 
  else{
    if ( npar ) {	// line-search for mw, iso, and clvd ================================

      npar--;
      dx = mt[npar].dd;
      i = 1; if (dx>0.001) i = 0;
      sol = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,i,bootstrap,search,norm);
      if (dx>0.001) {	/* do line search */
	mt[npar].par += dx;
	sol2 = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,0,bootstrap,search,norm);
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
	  else sol2 = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,0,bootstrap,search,norm);
	}
	mt[npar].sigma = 2*dx*dx/(sol2.err+sol1.err-2*sol.err);
	mt[npar].par -= dx+0.5*dx*(sol2.err-sol1.err)/(sol2.err+sol1.err-2*sol.err);
	sol = error(npar,nda,obs0,nfm,fm,fm_thr,max_shft,tie,mt,grid,1,bootstrap,search,norm);
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

	    //--------------KEY COMMAND---call misfit function------
	    sol=calerr(nda,obs0,max_shft,tie,norm,mtensor,amp,sol);

	    //fprintf(stderr, "Nsta=%d\n",Nsta);
	    sol.err=sol.err/Nsta;
	    *grd_err++ = sol.err;		/*error for this solution*/

	    if (best_sol.err>sol.err)
	      best_sol=sol;
	    
	    if (debug) { 
	      logf = fopen(logfile,"a");
	      fprintf(logf,"%3.1f\t%3.1f\t%3.1f\t%e\t%2.1f\t%2.2f\t%2.2f\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n",sol.meca.stk, sol.meca.dip, sol.meca.rak, sol.err/data2, mt[0].par, mt[1].par, mt[2].par, amp*1.0e20, mtensor[0][0], mtensor[1][1], mtensor[2][2], mtensor[0][1], mtensor[0][2], mtensor[1][2]);
	      fclose(logf);	  
	    }
	  } // loop for stk
	}   // loop for dip
	if (sol.meca.stk==360. &&  sol.meca.dip==90. &&  sol.meca.rak==90.){
	  loop++;
	  if (debug) {
	    sprintf(logfile,"%s_%03d","log",loop);
	    logf = fopen(logfile,"a");
	    fclose(logf);
	  }
	  logf = fopen("log_diff","a");
	  fprintf(logf,"%d\t%d\t%3.1f\t%3.1f\t%3.1f\t%f\t%2.2f\t%2.2f\t%2.2f\n",loop,interp, best_sol.meca.stk, best_sol.meca.dip, best_sol.meca.rak, best_sol.err, mt[0].par, mt[1].par, mt[2].par);
	  fclose(logf);
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




