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
    float u0;
    float uf;
    float v0;
    float vf;
    float w0;
    float wf;
    float k0;
    float kf;
    float h0;
    float hf;
    float s0;
    float sf;
    float gamma0;
    float gammaf;
    float delta0;
    float deltaf;
    float dip0;
    float dipf;

    float mw0;
    float mwf;
    float nmw;

    int nu;
    int nv;
    int nw;
    int nk;
    int nh;
    int ns;
    int ngamma;
    int ndelta;
    int ndip;

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
