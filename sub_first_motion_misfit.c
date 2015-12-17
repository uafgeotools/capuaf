/* =============================================================================
 *
 *       Filename:  sub_misfit_fm.c
 *
 *        Created:  08/08/2013 03:24:46 PM
 *         Author:  celso r. a.
 *
 *    Description:  this subroutine computes predicted P-first motion polarity
 *                  at each station
 *          Input:  mtensor, nfm, fm (struct), fid,
 *                  gamma, delta, strike, dip, rake
 *         Output:  (to file)
 *                  gamma, delta, strike, dip, rake, P1 P2 P3 ... PN
 *                  (N station polarities)
 *
 * =============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include "cap.h"

int misfit_first_motion(
            float mtensor[3][3],
              int nsta,
               FM *data,
             FILE *fid,
            float gamma,    // clvd. [-30,30] degrees
            float delta,    // isotropic. [0,180] degrees (note alternative delta [-90,90])
            float mw,       // not used but included here in case of future use
            float kappa,    // strike
            float theta,    // dip
            float sigma     // slip
        )
{
    float pol_theory;           // predicted polarity
    int ipol_theory, ipol_obs;  // integer polarities
    int misfit_fm=0;
    int debug_fm=0;             // 1: output predicted polarities. 0: don't output
    int i;

    // compute first motion misfit = sum(abs(obs_i-predicted_i))
    // output grid search parameters, polarity (in debug mode), and misfit
    // data written to fid (opened in cap.c)
    //fprintf(fid,"%5.1f %5.1f %5.1f %5.1f %5.1f  ", gamma, delta, kappa, theta, sigma);

//   if(debug_fm)
//   {
//       for(i=0; i<nsta; i++)
//       {
//           ipol_obs = (data+i)->type;
//           pol_theory = radpmt(mtensor, (data+i)->alpha, (data+i)->az, 1);
//           if(pol_theory>=0.0) ipol_theory = 1;
//           else ipol_theory = -1;
//           misfit_fm += abs(ipol_obs-ipol_theory);
//           fprintf(fid,"(%d %d) ", ipol_obs, ipol_theory);
//       }
//   }
//   else
//   {
        for(i=0; i<nsta; i++)
        {
            ipol_obs = (data+i)->type;
            pol_theory = radpmt(mtensor, (data+i)->alpha, (data+i)->az, 1);
            if(pol_theory>=0.0) ipol_theory = 1;
            else ipol_theory = -1;
            misfit_fm += abs(ipol_obs-ipol_theory);
        }
 //   }
    misfit_fm = misfit_fm/2;
//    fprintf(fid,"\t%d\n",misfit_fm);
    return misfit_fm;
}

