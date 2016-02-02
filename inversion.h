/**************************************************************************************
 * inversion.h:		head file for various inversion methods
 * revision history
 * 	10/23/1997	Lupei Zhu	initial coding
**************************************************************************************/
#ifndef __MY_INV__
  #define __MY_INV__

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MLTPLY 10.			// for the damping factor in Marquardt method
#define INV_NO_CHANGE	0.001		// stop iteration if no change in chi2

/*3D grid*/
typedef struct {
        int     n[3];           /* number of points */
        float   x0[3], step[3]; /* min and step */
        float   *err;
} GRID;

typedef struct
{
    float mw1;  float mw2;  float dmw; int nmw;

    float u1;   float u2;   int nu; int du;
    float v1;   float v2;   int nv; int dv;
    float w1;   float w2;   int nw; int dw;
    float k1;   float k2;   int nk; int dk;
    float h1;   float h2;   int nh; int dh;
    float s1;   float s2;   int ns; int ds;

    float gamma1;   float gamma2;   int dgamma;
    float delta1;   float delta2;   int ddelta;
    float dip1;     float dip2;     int ddip;

    int nsol;

} SEARCHPAR;

extern int	svdrs(float *, int, int, int, float *, int, float *);
extern float	iter(float *, float *, float *(*f)(float *), int, int, int, int, float, int);
float	marquardt(float *, float *, float *, float *(*)(float *), float *(*)(float *), int, int, int, float, float);
float	jump(float *, float *, float *, float *(*)(float *), float *(*)(float *), int, int, int, int, float);
float	ridge(float *, float *, float *, float *(*)(float *), float *(*)(float *), int, int, int, int, float);
float	grid2d(float *,int,int,float *,float *,float *,float *,float *,int *,int *);
float	grid3d(float *,int *,float *, int *, int *, int *);

#endif
