#include "cap.h"

SOLN get_tshift_corr_misfit(
        int   	nda,    // number of stations
        DATA       *obs0,  // pointer to data
        const int	*max_shft,
        float      tie,
        int        norm,// type of norm
        float      mtensor[3][3],
        float      amp,
        SOLN       sol
        ) {
  int i,j,k,l,m,k1,kc,z0,z1,z2;
  float x, x1, x2, y, y1, y2,cfg[NCP];
  DATA	*obs;
  // SOLN	sol;
  float rad[6], arad[4][3];
  COMP	*spt;
  float	*f_pt0, *f_pt1, *r_pt, *r_pt0, *r_pt1, *r_iso, *z_pt, *z_pt0, *z_pt1,*z_iso;

  for(obs=obs0,sol.wferr=0.,i=0;i<nda;i++,obs++){
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
    z0=1;
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
    if (z1==z2){z1=1;z2=1;}
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
    if (z1==z2){z1=1;z2=1;}
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
		   
      y1 = 1.;
      /* find out the scaling factor for teleseismic distances */
      if (obs->tele && spt->on_off) {
	if (cfg[j]>0.) y1 = cfg[j]/x2;
	else y1 = 0.;
      }
      sol.scl[i][j] = y1;
      
      x1 = spt->rec2+x2*y1*y1-2.*cfg[j]*y1; // MISFIT FUNCTION
       //fprintf(stderr,"%d\n",spt->npt);
      if (norm==1) x1 = sqrt(x1);
      x1 = x1/(spt->npt * spt->rew);
      sol.error[i][j] = spt->on_off*x1;	/*L2 error for this com.*/
      sol.cfg[i][j] = 100*cfg[j]/sqrt(spt->rec2*x2);
      sol.wferr += spt->on_off*x1;
    }
  }
  return(sol);
}
