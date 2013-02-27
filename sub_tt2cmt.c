/*
 * =============================================================================
 *    Description:  This function generates moment tensor elements given:
 *                  strike/dip/rake, gamma, delta, M0
 *                  (See Tape & Tape, 2012)
 *                  
 *                  This code based on function TT2CMT.m
 *
 *    20130123 --  Celso Alvizuri
 *  ----------------------------------------------------------------------------
 *
 *          Usage: tt2cmt(gamma, delta, m0, kappa, theta, sigma, tensor)
 *          input: gamma, delta, m0, kappa, theta, sigma
 *         output: tensor[3][3]
 *
 * =============================================================================
 */


#include<stdio.h>
#include<math.h>

void tt2cmt(float gamma,    // clvd. [-30,30] degrees
            float delta,    // isotropic. [0,180] degrees (note alternative delta [-90,90])
            float m0,       // not used but included here in case of future use
            float kappa,    // strike
            float theta,    // dip
            float sigma,    // slip
            float mtensor[3][3])    // output
{
    // define constants and common operations
    float k8R6, kR3, k2R6, k2R3,k4R6;          // constants(==k), R=square root. k8r6 = 8*sqrt(6);
    float Cb, Cg, Cs, Ct, Ck, C2k, C2s, C2t;   // definitions for trig operations
    float Sb, Sg, Ss, St, Sk, S2k, S2s, S2t;
    float beta;

    beta = 90-delta;

    float m0scale = 1.0;

    float pi  = 3.1415926535897932385;
    float deg2rad = pi/180.0;

//    fprintf(stderr, "debug. received values: %lf %f %f %f %f %f\n", gamma, beta, m0, kappa, theta, sigma);

    // define constants
    kR3  = sqrt(3);
    k2R6 = 2.0*sqrt(6);
    k2R3 = 2.0*sqrt(3);
    k4R6 = 4.0*sqrt(6);
    k8R6 = 8.0*sqrt(6);

    // define common trig operations
    Cb  = cos(beta*deg2rad);
    Cg  = cos(gamma*deg2rad);
    Cs  = cos(sigma*deg2rad);
    Ct  = cos(theta*deg2rad);
    Ck  = cos(kappa*deg2rad);
    C2k = cos(2.0*kappa*deg2rad);
    C2s = cos(2.0*sigma*deg2rad);
    C2t = cos(2.0*theta*deg2rad);

    Sb  = sin(beta*deg2rad);
    Sg  = sin(gamma*deg2rad);
    Ss  = sin(sigma*deg2rad);
    St  = sin(theta*deg2rad);
    Sk  = sin(kappa*deg2rad);
    S2k = sin(2.0*kappa*deg2rad);
    S2s = sin(2.0*sigma*deg2rad);
    S2t = sin(2.0*theta*deg2rad);

    /* 
     * ========================================================================
     * get moment tensor elements
     * ========================================================================
     */
    mtensor[0][0] = m0scale* (1./24.)*(k8R6*Cb + Sb*(-24.*Cg*(Cs*St*S2k + S2t*Sk*Sk*Ss) +
                kR3*Sg*(-(1. + 3.*C2k)*(-1. + 3.*C2s) + (12.*C2t*Cs*Cs*Sk*Sk) - (12.*Ct*S2k*S2s))));

    mtensor[1][1] = m0scale* (1./6.)*(k2R6*Cb + Sb*(kR3*Ct*Ct*Ck*Ck*(1. + 3.*C2s)*Sg - k2R3*Ck*Ck*Sg*St*St + 
                kR3*(1. - 3.*C2s)*Sg*Sk*Sk + 6.*Cg*Cs*St*S2k + 
                3.*Ct*(-4.*Cg*Ck*Ck*St*Ss + kR3*Sg*S2k*S2s)));

    mtensor[2][2] = m0scale* (1./12.)*(k4R6*Cb + Sb*(kR3*Sg*(-1. - 3.*C2t + 6.*C2s*St*St) + 12.*Cg*S2t*Ss));

    mtensor[0][1] = m0scale* (1./8.)*Sb*(4.*Cg*(2.*C2k*Cs*St + S2t*S2k*Ss) +
                kR3*Sg*((1. - 2.*C2t*Cs*Cs - 3.*C2s)*S2k + 4.*Ct*C2k*S2s));

    mtensor[0][2] = m0scale* (-1./2.)*Sb*(k2R3*Cs*Sg*St*(Ct*Cs*Sk - Ck*Ss) +
            2.*Cg*(Ct*Ck*Cs + C2t*Sk*Ss));

    mtensor[1][2] = m0scale* (1./2.)*Sb*(Ck*(kR3*Cs*Cs*Sg*S2t + 2.*Cg*C2t*Ss) +
            Sk*(-2.*Cg*Ct*Cs + kR3*Sg*St*S2s));

    mtensor[1][0] = mtensor[0][1];
    mtensor[2][0] = mtensor[0][2];
    mtensor[2][1] = mtensor[1][2];
}
