#include<math.h>
#include "cap.h"

void fmp_print_parameters(FILE *fid, FMPDATA *fmpdata)
{

    /* project to lower hemisphere */
    if(fmpdata->toa > 90.0)
    {
        fmpdata->azim += 180.0;
        if(fmpdata->azim > 360)
        {
            fmpdata->azim = fmpdata->azim - 360.0;
        }
        fmpdata->toa = 180.0-fmpdata->toa;
    }
    fmpdata->strad_lh = sqrt(2.0) * sin(fmpdata->toa * PI/360.0);

    /* output data to file */
    /* possibly add event lat/lon/dep */
    fprintf(fid, "%s %s %11.6f %11.6f %2d %s %02d %7.2f %7.2f %f %6.2f %6.2f %5.1f\n",
            fmpdata->evid, fmpdata->stname,
            fmpdata->stlo, fmpdata->stla,
            fmpdata->pol,
            fmpdata->vmod,
            fmpdata->idep,
            fmpdata->azim, 
            fmpdata->toa, fmpdata->strad_lh,
            fmpdata->tp, fmpdata->ts, fmpdata->dist);
}
