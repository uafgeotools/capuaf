/* 
   This file contains functions used by routines initSearchMT, getRandMT, getGridMT,
   and searchMT in file sub_inversion.c

   Main references
   Tape, W., and C. Tape (2015), A uniform parametrization of moment tensors,
   Geophys. J. Int., 202, 2074-2081.
   Matlab routine uniformMT.m

   20160108 celso alvizuri - cralvizuri@alaska.edu 

*/

#include "uv2lune.h"

#include "generated_u2beta.c"   // data for interpolation of beta = beta(u)
                                // data generated using matlab script gen_u2beta.m

// compute beta(u) in vector form
void u2beta_vec(float *pArray_u, float *pArray_beta, int nsol)
{
    int i;
    float * pxa = calloc(NBETA, sizeof(float));
    float * pya = calloc(NBETA, sizeof(float));
    for(i = 0; i < NBETA; i++) {
        pxa[i] = generated_u2beta[i].u;
        pya[i] = generated_u2beta[i].beta;
    }

    fprintf(stderr,"create vector u2beta...");
    interp_lin(pxa, pya, NBETA, pArray_u, pArray_beta, nsol);

    free(pxa);
    free(pya);
}

// compute delta(beta)
float beta2delta(float beta)
{
    float delta;
    delta = beta - (PI / 2.0);     //delta = beta - 90.0;
//    fprintf(stdout,"CHECK DELTA. %f \n", delta * r2d);
    return delta;
}

// compute beta(delta)
float delta2beta(float delta)
{
    float beta;
    beta = delta + (PI / 2.0);     //beta = delta + 90;
    if((beta < 0) || (beta > PI)) {
        fprintf(stderr,"WARNING. beta is out of bounds [0, PI]. check input values.\n");
        fprintf(stderr,"beta = %f\n", beta);
    }
//    fprintf(stdout,"CHECK BETA. %f \n", delta * r2d);
    return beta;
}

// compute u(w)
float w2u(float w)
{
    float u;
    u = w + (3.0 * PI / 8.0);

    // if u=-0 set to u=+0
    // this is a numerical issue that can cause NANs when using var u
    if(u < 0.0) {
        fprintf(stderr,"WARNING u = %e new value:\n", u);
        u = 0.0;
        fprintf(stderr,"WARNING u = %e \n", u);
    }
    return u;
}

// compute w(u)
float u2w(float u)
{
    float w;

    // if u=-0 set to u=+0
    // this is a numerical issue that can cause NANs when using var u
    if(u < 0.0) {
        fprintf(stderr,"WARNING u = %e new value:\n", u);
        u = 0.0;
        fprintf(stderr,"WARNING u = %e \n", u);
    }

    w = u - (3.0 * PI / 8.0);
    return w;
}

// compute w(delta)
float delta2w(float delta)
{
    float u, w, beta;
    beta = delta2beta(delta);
    u = (3./4.) * beta - (1./2.) * sin(2.0 * beta) + (1./16.) * sin(4.0 * beta);
    w = u2w(u);
    return w;
}

// compute delta(beta)
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

// compute v(gamma)
float gamma2v(float gamma)
{
    float v;
    v = (1./3.) * sinf(3.0 * gamma);
    return v;
}

// compute gamma(v)
float v2gamma(float v)
{
    float gamma;
    float k3 = 1./3;

    // sanity check range for v
    if((v < -k3) || (v > k3)) {
        fprintf(stderr,"WARNING. v= %5.2e is outside range +-1/3\n", v);
        exit(-1);
    }
    gamma = k3 * asinf(3. * v);
    // fprintf(stdout,"CHECK GAMMA %f \n", gamma * r2d);
    return gamma;
}

// compute gamma(v) in vector form
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

// compute delta(sin(iso))
float siso2delta(float siso)
{
    float delta;
    delta = asinf(siso);
    // fprintf(stdout,"CHECK GAMMA %f \n", gamma * r2d);
    return delta;
}

// compute delta(sin(iso)) in vector form
void siso2delta_vec(float *pArray_s, int nsol, float *pArray_dip)
{
    int i;
    fprintf(stderr,"create vector siso2delta... ");
    for(i=0; i<nsol; i++) {
        pArray_dip[i] = siso2delta(pArray_s[i]);
    }
    fprintf(stderr,"done.\tNPTS = %d\n", i);
}

// compute dip(h)
float h2dip(float h)
{
    float dip;
    if((h < 0.0) || (h > 1.0)) {
        fprintf(stderr,"STOP. h= %f is out of bounds [0,1]\n", h);
        exit(-1);
    } else {
        dip = acosf(h);
    }
    return dip;
}

// compute dip(h) in vector form
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

    NOTE this function currently will skip the endpoints even when 
    searching a subset of the full parameter space. This is the same
    as in uniformMT.m
    It may make more sense to add special checks for each parameter
    so that subsets include endpoints.
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
            fprintf(stderr,"\nWARNING. end point does not match expected end point!\n");
            fprintf(stderr,"xmax(actual) = %f. xmax(expected) = %f\n", pArray[npoints-1], xmax);
        }
        fprintf(stderr,"\tdone. NPTS = %d\n", i);
    }
}

// create array of magnitudes. this option includes endpoints.
// NOTE dx here is decimal because changes in magnitude are usually <1
void magvec(float xmin, float xmax, float dx, float *pArray)
{
    int i;
    int npoints = 0;
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

