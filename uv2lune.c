/* 
   Functions to compute (gamma, delta) from (u,v)
   These are called from function initSearchMT

   For (u,v) parameterization see 
   "A uniform parameterization of moment tensors" (Tape & Tape, 2015)

   20160108 celso alvizuri - cralvizuri@alaska.edu 

*/

#include "uv2lune.h"

/* compute beta(u) in vector form */
void u2beta_vec(float *pxa, float *pya, int na, float *pArray_u, float *pArray_beta, int nsol)
{
    //interp1(xa, npts, ya, xi, nu2beta_out, yi);             // previous
    fprintf(stderr,"create vector u2beta...");
    interp_lin(pxa, pya, na, pArray_u, pArray_beta, nsol);   // best
    //    fprintf(stderr,"done\n");
}

float beta2delta(float beta)
{
    float delta;
    delta = beta - (PI / 2.0);     //delta = beta - 90.0;
//    fprintf(stdout,"CHECK DELTA. %f \n", delta * r2d);
    return delta;
}

/* compute delta=delta(beta) */
void beta2delta_vec(float *pArray_beta, float *pArray_delta, int nsol)
{
    int i;
    fprintf(stderr,"create vector beta2delta... ");
    for(i=0; i<nsol; i++) {
        // convert beta to delta
        pArray_delta[i] = beta2delta(pArray_beta[i]);
    }
    fprintf(stderr,"done. NPTS = %d\n", i);
}

/* compute gamma(v) */
float v2gamma(float v)
{
    float gamma;
    float k3 = 1./3;

    // sanity check range for v
    if((v < -k3) || (v > k3)) {
        fprintf(stderr,"WARNING. v= %5.2e is outside range +-1/3\n");
        return -1;
    }
    gamma = k3 * asinf(3. * v);
    // fprintf(stdout,"CHECK GAMMA %f \n", gamma * r2d);
    return gamma;
}

/* compute vector of gamma values */
void v2gamma_vec(float *pArray_v, int nsol, float *pArray_gamma)
{
    int i;
    fprintf(stderr,"create vector v2gamma ... ");
    for(i = 0; i < nsol; i++) {
        pArray_gamma[i] = v2gamma(pArray_v[i]);
    }
    fprintf(stderr,"done.\tNPTS = %d\n", i);
    //fprintf(stderr,"routine: %s line %d \n", __FILE__, __LINE__);
}

/* convert h to dip */
float h2dip(float h)
{
    return acosf(h);
}

/* get vector of dip values */
void h2dip_vec(float *pArray_h, int nsol, float *pArray_dip)
{
    int i;
    fprintf(stderr,"create vector h2dip... ");
    for(i=0; i<nsol; i++) {
        pArray_dip[i] = h2dip(pArray_h[i]);
    }
    fprintf(stderr,"done.\tNPTS = %d\n", i);
}

/* 
    Create vector of random points in a range.
    NOTE
    There may be more robust implementations.
    This function uses function RAND(). 
    Previous random implementations use DRAND48().
    The linux man pages say DRAND48(3) is obsolete and to use RAND() instead.
 */
void randvec(float xmin, float xmax, int nsol, float *pArray)
{
    int i;
    float width;
    float randval;
    width = (xmax - xmin);
    fprintf(stderr,"create RAND vector. xmin= %10.5f xmax= %10.5f NPTS = %d ... ", xmin, xmax, nsol);

    for(i = 0; i < nsol; i++) {
        randval = ((float) rand()) / RAND_MAX;
        randval = ((randval * width) + xmin);
        pArray[i] = randval;

        //check random values (for debugging)
        // fprintf(stdout,"CHECK RAND %10.6f %10.6f %22.18f \n", xmin, xmax, randval); 
    }
    fprintf(stderr,"\tdone. NSOL = %d\n", i);
}

/* create vector of grid points in a range */
/* 
    NOTE uniformMT.m applies offset to full range AND subrange
         shouldn't it apply only to full range?
*/
void gridvec(float xmin, float xmax, int npoints, float *pArray)
{
    int i;
    float dx, xoffset;
    if(npoints < 0) {
        fprintf(stderr,"STOP. npoints < 0\n");
        exit(-1);
    } else if (npoints == 1) {
        fprintf(stderr, "create GRID vector. xmin= %10.5f xmax= %10.5f dx= %12.5e ... ", xmin, xmax, dx);
        fprintf(stderr,"WARNING. this is a point solution.");
        fprintf(stderr,"\tdone. NPTS = %d\n", npoints);
        pArray[i] = xmin;   // or xmax. they should be identical. maybe add extra check
    } else {
        /* get dx based on original user input */
        dx = (xmax - xmin) / (float) (npoints - 1); // NOTE npts-1

        /* Offset from boundaries.
         * Factor 2.0 looks arbitrary. it could be made smaller.
         * Also it's not clear why dx should scale with the number of points. 
         * ie why not just make this a tiny constant. */
        xoffset = dx/2.0;

        fprintf(stderr, "create GRID vector. xmin= %10.5f xmax= %10.5f dx= %12.5e ... (old)\n", xmin, xmax, dx);
        /* apply offset at boundaries and recalculate dx based on new boundaries */
        xmin = xmin + xoffset;
        xmax = xmax - xoffset;
        dx = (xmax - xmin) / (float) (npoints - 1);

        fprintf(stderr, "create GRID vector. xmin= %10.5f xmax= %10.5f dx= %12.5e ... (new)", xmin, xmax, dx);
        for(i = 0; i < npoints; i++) {
            pArray[i] = xmin + dx * (float) i;
            //        fprintf(stdout,"CHECK GRIDVEC. %f \n", pArray[i] );
        }
        if((pArray[npoints-1] - xmax) > TOLERANCE) {
            fprintf(stderr,"WARNING. end point does not match expected end point!\n");
            fprintf(stderr,"xmax(actual) = %f. xmax(expected) = %f\n", pArray[npoints-1], xmax);
        }
        fprintf(stderr,"\tdone. NPTS = %d\n", i);
    }
}

void magvec(float xmin, float xmax, float dx, float *pArray)
{
    int i;
    int npoints = 0;
    float dmw;
    int count=0;

    if(dx < TOLERANCE) {
        npoints = 1;
    } else {
        // endpoints inclusive for magnitude vector
        npoints = (int) roundf(((xmax - xmin) / dx) + 1);
    }
    fprintf(stderr, "\ncreate Mw vector. xmin= %10.5f xmax= %10.5f dx= %10.5e npts = %d... ", xmin, xmax, dx, npoints);

    for(i = 0; i < npoints; i++) {
        pArray[i] = xmin + (dx * (float) i);
        count++;
    }
    fprintf(stderr,"\tdone. NPTS = %d\n", count);
}

/*
    This is brute force and slow. The gnu method is faster
    [timings]
    But the request is to not add more library dependencies (gnu) to CAP for now.
*/
int find_nearest_index(float value, float * pData, int len)
{
    float dist, newDist;
    int inearest;
    int i;

    inearest = -1;  // if negative then nearest index not found
    dist = DBL_MAX;

    for(i = 0; i < len; i++) {
        newDist = value - pData[i]; // negative dist means 'value' is to the left of pData[i]
        // check that distances are not negative (check when this becomes an issue)
        if ((newDist > 0) && (newDist < dist)) {
            dist = newDist;
            inearest = i;
        }
    }

    return inearest;
}

void interp_lin(float *x,     // x array  (input)
                float *y,     // y array  (input)
                int x_size,   // size     (input)
                float *xx,    // xx array (input)
                float *yy,    // yy array (output)
                int xx_size)  // size     (input)
{
    float dx, dy, *slope, *intercept;
    int i, index_nearest;

    slope = (float *) calloc(x_size, sizeof(float));
    intercept = (float *) calloc(x_size, sizeof(float));

    // Precompute slopes and intercepts for all elements in x and y
    // the last elements in x, y are made equal to the previous element
    for(i = 0; i < x_size; i++) {
        if(i < x_size-1) {
            dx = x[i + 1] - x[i];
            dy = y[i + 1] - y[i];
            slope[i] = dy / dx;
            intercept[i] = y[i] - x[i] * slope[i];
        } else {
            slope[i] = slope[i-1];
            intercept[i] = intercept[i-1];
        }
    }

    fprintf(stderr,"\tcalculating values beta(u) ... ");
    for (i = 0; i < xx_size; i++) {
        index_nearest = find_nearest_index(xx[i], x, x_size);

        // compute all interpolations yy for xx
        if (index_nearest != -1) {
            yy[i] = slope[index_nearest] * xx[i] + intercept[index_nearest];
//            fprintf(stderr,"CHECK beta(u). inearest= %d xi= %f yi= %f \n", index_nearest, x[i], y[i]); 
        }
        else
            yy[i] = DBL_MAX;
    }
    fprintf(stderr,"done. NPTS = %d\n", i);
    free(slope);
    free(intercept);
}

