/***************************************************************
	cap.h		Head file for cap.c
***************************************************************/


#ifndef __CAP_HEAD__
  #define __CAP_HEAD__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "sac.h"
#include "Complex.h"
#define NRANSI
#include "inversion.h"
#include "radiats.h"

#include "sub_uv2lune.h"
// #include <time.h>   // when using OPTION 2 in the random seed generator

#include <omp.h>

/***********************Constants********************************/

#define STN	200		/* up to STN stations */

#define NRC	3		/* 3 components of records */
static char cm[NRC]={'t','r','z'};

#define NCP	5		/* 5 segments, 2*Pnl 2*PSV 1*SH */
static int kd[NCP]={0,1,2,1,2};	/* SH, SVr, SVz, Pnlr, Pnlz */
static int kk[NCP]={0,2,6,2,6};	/* index of 1st Greens' fn for each segment, see below */

#define NRF	4		/* number of fundamental sources: EX, SS, DS, DD. SH only has SS and DS*/
#define NGR	10		/* 8 com. of greens function for DC plus 2 for EX: SHx2, Rx4, Vx4 */
static char grn_com[NGR]={'8','5','b','7','4','1','a','6','3','0'};

#define NFFT	2048		/* for Q operator */

#define r2d 180.0 / PI
#define d2r PI / 180.0
#define RANDSEED 12345

#define TOLNMAG 1e-6    // check for point magnitudes
#define TOLGRID 1e-6    // check for grid endpoints
#define TOLBETA 1e-6    // check if close to the DC
#define TOLDELTA 1e-6   // check precision near 0

#define NBETA 1000  /*  number of points used for interpolation of u=u(beta)  */

/***********************global vars********************************/

extern int total_n,loop,start,debug, Ncomp,Psamp[STN],Ssamp[STN],edep;
extern float data2;

/* flags for computing uncertainty on the lune. 1==apply */
extern int only_first_motion;    // polarity misfit. runs ONLY polarity, no waveform misfit
extern int misfit_on_lune;       // waveform misfit. output misfit on the lune 
extern char filename_prefix[255];    // used for all output files

/* workaround for filter issues with small magnitude events (Uturuncu) */
// this has not been tested with DIRECTIVITY option
extern int FTC_data, FTC_green;// for original CAP set FTC_data=0, FTC_green=0

/* allows use of polarities even when weight=0.
 * Note CAP still needs at least 1 waveform for the inversion */
extern int skip_zero_weights;    // for original CAP set skip_zero_weights=1

// Flag to create regular grid as in Alvizuri & Tape (2016) and Silwal & Tape (2016).
// NOTE reproducibility may not be exact since grid spacing uses function gridvec.
// Function gridvec does not implement the discretization of the previous version of 
// cap.c which uses rules to account for special grid points.
// Function gridvec also avoids endpoints in all parameters.
// See NOTE flag LUNE_GRID_INSTEAD_OF_UV in function sub_inversion.c
extern int LUNE_GRID_INSTEAD_OF_UV;

/*********************** Data Structure***************************/

// focal mechanism updated to include the full moment tensor.
// (gamma, delta, strike, dip, rake, magnitude)
typedef struct {
	float	gamma;	/* CLVD */
	float	delta;	/* ISO */
	float	stk;	/* strike */
	float	dip;	/* dip */
	float	rak;	/* rake */
	float	mag;	/* Mw */
} MECA;

/* a portion of observed waveform and corresponding 3 components
of Green's functions, cross-correlations, and L2 norms, etc */
typedef struct {
	int	npt;		/* number of points */
	int	on_off;		/* on or off */
	float	*rec;		/* data */
	float	*syn[NRF];	/* Green's functions */
	float	b;		/* beginning time of the data */
/*	float	w; */
	float	rec2;		/* sum rec*rec */
	float	syn2[10];	/* c*sum syn[i]*syn[j], i=0..NRF-1, j=i..0, c=1 if i=j, 2 otherwise */
	float	*crl[NRF];	/* sum rec*syn[i], i=0..NRF-1 */
} COMP;

typedef struct {
	char	stn[10];
	float	*rec[NRC];
	float	*grn[NGR];
	float	az;
	float	dist;
	float	alpha;		/* take-off angle */
	int	tele;		/* 1 if the station is at teleseismic distances*/
	COMP	com[NCP];
} DATA;

typedef struct {
	MECA	meca;
	float	dev[6];			/* uncertainty ellipsis */
	float	err;			/* chi-square of waveform misfits for this solution */
	int	cfg[STN][NCP];		/* correlation for each comp. */
	int	shft[STN][NCP];		/* time shift for each comp. */
	float	error[STN][NCP];	/* chi-square of waveform misfits for each component */
	float   scl[STN][NCP];		/* amplifications to GF for each component */
	int	ms;			/* number of local minimums < 10 */
	int	others[10];		/* top 10 best solutions */
	int	flag;			/* =1 if the best soln is at boundary */
} SOLN;

typedef struct {
	float	par;	// can be mw, iso, or cvld
	float	sigma;	// variance
	float	min, max;	// range
	float	dd;	// search step
} MTPAR;

/* first-motion data */
typedef struct {
	float	az;	/* azimuth */
	float	alpha;	/* take-off angle */
	int	type;	/* 1=P; 2=SV; 3=SH; positive=up; negative=down */
} FM;

/* data for processing first-motion polarity */
typedef struct {
    char evid[256];
    char stname[256];
    char vmod[256]; // name of velocity model

    int idep;       // inversion depth
    float toa;      // takeoff angle

    float dist;     // source-stn distance

    float stlo;
    float stla;     
    float azim;
    float strad_lh; // station location on lower hemisphere (bb plot)
    int pol;        // observed polarity
    float tp, ts;   // arrival time for p and s waves
} FMPDATA;

/* for tracking misfit at each point on the lune */
typedef struct {
    /* beachball */
    float gamma;
    float delta;
    float stk;
    float dip;
    float rak;
    float misfit;

    /* moment tensor */
    float mrr;
    float mtt;
    float mpp;
    float mrt;
    float mrp;
    float mtp;

    float mag;
} LUNE_MISFIT;

typedef struct {
    float v;
    float w;
    float kappa;    // degrees
    float theta;    // degrees
    float sigma;    // degrees
    float mag;
    float misfit_wf;
    float misfit_fmp;
} OUTPUTBB;

typedef struct {
    float mrr;
    float mtt;
    float mpp;
    float mrt;
    float mrp;
    float mtp;
    float mw;
    float misfit_wf;
    float misfit_fmp;
} OUTPUTMT;

typedef struct
{
    float gamma;    // radians
    float delta;    // radians
    float kappa;    // radians
    float theta;    // radians
    float sigma;    // radians
    float mw; 
    float misfit_fmp;
    float misfit_wf;
    float VR;
    float v;    // radians
    float w;    // radians
} ARRAYMT; 

// for interpolation within u=u(beta)
typedef struct
{
    float u;
    float beta;
} GENU2BETA; 

/* function declaration */
  /* SOLN	error(int,int,DATA *,int,FM *,float,const int *,float,MTPAR *,GRID,int,int,int,int); */
SOLN initSearchMT(int,DATA *,int,FM *,float,const int *,float,MTPAR *,GRID,int,int,int, SEARCHPAR *, ARRAYMT *);
void    taper(float *aa, int n);
float	*trap(float, float, float, int *);
float	*cutTrace(float *, int, int, int);
int	discard_bad_data(int,DATA *,SOLN,float,float *);
int	check_first_motion(float mt[3][3], FM *fm, int n, float fm_thr);
void tt2cmt(float gamma, float delta, float m0, float kappa, float theta, float sigma, float mtensor[3][3]);
int misfit_first_motion(float mtensor[3][3], int nsta, FM *data, FILE *fid, float gamma, float delta, float mw, float kappa, float theta, float sigma);
void fmp_print_parameters(FILE *fid, FMPDATA *fmpdata);
SOLN get_tshift_corr_misfit(int,DATA *,const int *, float,int,float mtensor[3][3],float,SOLN);

/* functions for uniformMT */
void getRandMT(SEARCHPAR * searchPar, ARRAYMT * arrayMT);
void getGridMT(SEARCHPAR * searchPar, ARRAYMT * arrayMT);
SOLN searchMT(int,DATA *,int,FM *,float,const int *,float,MTPAR *,GRID,int,int,int, SEARCHPAR *, ARRAYMT *);

#endif
