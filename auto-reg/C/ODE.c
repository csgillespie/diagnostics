#include "headers.h"

/* 
   Moment closure equations for the auto-reg network
*/

int func(double t, const double kappa[], double dkappa[], void *params)
{
 
  double *c = params;
  double k1, k2;
  k1 = c[12]; k2 = c[13];
  
  dkappa[11]=kappa[13]*c[4]-kappa[22]*kappa[9]*c[0]-kappa[11]*(c[1] + c[8] + kappa[20]*c[0]);
  dkappa[18]=-kappa[18]*(c[1] + c[11] + kappa[20]*c[0]) + c[7]*kappa[10]-kappa[9]*kappa[25]*c[0];
  dkappa[20]=-kappa[9]*c[1] + k1*c[3]-c[2]*kappa[23]-c[3]*kappa[5]- 
    kappa[20]*(c[10] + c[2]*kappa[5] + kappa[9]*c[0]) + k2*c[1] + c[5]*kappa[2]-c[0]*kappa[24];
  dkappa[23]=-kappa[23]*(c[10] + c[2]*kappa[5]-c[2] + kappa[9]*c[0] + c[3] +
  c[2]*kappa[20]) + k1*c[3] + c[2]*kappa[20]*kappa[5]-c[3]*kappa[5]
    -kappa[20]*c[0]*kappa[12]-c[1]*kappa[12]-c[3]*kappa[8] +
  kappa[7]*c[5]-c[2]*kappa[26]*kappa[5]-c[2]*kappa[8]*kappa[20];
  dkappa[4]=c[8]*kappa[2] + 2*kappa[11]*c[4]-2*c[8]*kappa[4] + kappa[9]*c[4];
  dkappa[10]=-kappa[21]*kappa[9]*c[0]-kappa[10]*(c[9] + c[1] + kappa[20]*c[0]) + c[6]*kappa[12];
  dkappa[13]=-kappa[9]*c[1]-2*kappa[13]*(kappa[20]*c[0] + c[1]) + k2*c[1] +
    c[0]*kappa[24] + kappa[9]*kappa[20]*c[0]-2*kappa[9]*c[0]*kappa[24];
  dkappa[12]=-kappa[12]*(c[3] + c[2]*kappa[20] + kappa[20]*c[0] +
                         c[1])-kappa[9]*kappa[23]*c[0]-c[2]*kappa[5]*kappa[24];
  dkappa[19]=2*kappa[15]*c[7] + kappa[0]*c[7]-2*c[11]*kappa[19] + kappa[14]*c[11];
  dkappa[24]=-kappa[12]*(c[3] + c[2]*kappa[20]) + kappa[11]*c[5] -
    kappa[9]*(c[1]+ kappa[26]*c[0])-
    kappa[24]*(c[10] + c[2]*kappa[5] + c[0]*(-1 + kappa[20] + kappa[9]) + c[1]) 
    + (kappa[9]-kappa[13])*kappa[20]*c[0] + (k2-kappa[13])*c[1];
  dkappa[9]=-kappa[9]*(c[1] + kappa[20]*c[0]) + k2*c[1]-c[0]*kappa[24];
  dkappa[7]=-c[2]*kappa[22]*kappa[5]-kappa[7]*(c[3] + c[2]*kappa[20] + c[8]) + kappa[12]*c[4];
  dkappa[6]=-kappa[6]*(c[9] + c[2]*kappa[20] + c[3])-kappa[21]*c[2]*kappa[5] + c[6]*kappa[8];
  dkappa[16]=-kappa[16]*(c[11] + c[8]) + c[7]*kappa[3] + kappa[18]*c[4];
  dkappa[3]=c[4]*kappa[10] + c[6]*kappa[7]-kappa[3]*(c[8] + c[9]);
  dkappa[15]=-kappa[15]*(c[9] + c[11]) + kappa[17]*c[6] + kappa[1]*c[7];
  dkappa[5]=k1*c[3]-c[2]*kappa[23]-kappa[5]*(c[3] + c[2]*kappa[20]);
  dkappa[17]=-kappa[17]*(c[2]*kappa[20] + c[3] + c[11])-c[2]*kappa[25]*kappa[5] + kappa[6]*c[7];
  dkappa[0]=c[6]*kappa[5]-c[9]*kappa[0];
  dkappa[25]=-kappa[25]*(c[10] + c[2]*kappa[5] + c[11] +
  kappa[9]*c[0])-kappa[18]*(c[1] + kappa[20]*c[0])-kappa[17]*(c[2]*kappa[20] + c[3]) + kappa[21]*c[7] + c[5]*kappa[16];
  dkappa[22]=-kappa[22]*(c[2]*kappa[5] + c[10] + kappa[9]*c[0] + c[8]) + c[4]*kappa[24] + kappa[4]*c[5]-kappa[11]*(c[1] + kappa[20]*c[0])-kappa[7]*(c[3] + c[2]*kappa[20]);
  dkappa[26]=-2*kappa[26]*(c[0]*kappa[9] + c[10] + c[2]*kappa[5])-kappa[9]*c[1] + k1*c[3] + kappa[23]*(c[2]-2*c[3]-2*c[2]*kappa[20])-c[3]*kappa[5] + kappa[20]*(c[10] + kappa[9]*c[0] + c[2]*kappa[5]) + k2*c[1] + c[5]*kappa[2] + kappa[24]*(c[0]-2*kappa[20]*c[0]-2*c[1]) + 2*kappa[22]*c[5];
  dkappa[2]=-c[8]*kappa[2] + kappa[9]*c[4];
  dkappa[8]=k1*c[3] + kappa[23]*(c[2]-2*c[2]*kappa[5])-kappa[5]*(c[3]-c[2]*kappa[20])-2*kappa[8]*(c[3] + c[2]*kappa[20]);
  dkappa[14]=kappa[0]*c[7]-kappa[14]*c[11];
  dkappa[21]=-kappa[21]*(kappa[9]*c[0] + c[9] + c[2]*kappa[5] + c[10]) + c[6]*kappa[23]-kappa[20]*(c[0]*kappa[10] + kappa[6]*c[2])-c[3]*kappa[6]-c[1]*kappa[10] + kappa[3]*c[5];

  dkappa[1]=-2*c[9]*kappa[1] + c[6]*(kappa[5] + 2*kappa[6]) + c[9]*kappa[0];

  return GSL_SUCCESS;
}

void ode(double  *params, double maxtime, double * sps)
{
    double *y;
    y = sps;
    const gsl_odeiv_step_type *T = gsl_odeiv_step_gear2;
    gsl_odeiv_step *stp = gsl_odeiv_step_alloc(T, 27);
    gsl_odeiv_control *c = gsl_odeiv_control_y_new(1e-4, 0.0);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(27);
    gsl_odeiv_system sys = {func, 0, 27, params};

    double t = 0.0, t1 = maxtime;
    double h = 1e-6;


    while (t < t1) {
        gsl_odeiv_evolve_apply (e, c, stp, &sys, &t, t1, &h, y);
    }
  
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (stp);
}
