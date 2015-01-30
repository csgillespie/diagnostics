#include <math.h>
#include <gsl/gsl_randist.h>

/* This is a particularly nasty piece of code. I can claim no credit in creating
   this evil thing

   Direct method of the auto-reg model
*/

void gillespie(double *sps, double *pars, double time_max, gsl_rng *r) {

  double rg,ri,g,i,G,I,t=0.0,ran,total;
  double k1, k2;
  k1 = pars[12]; k2 = pars[13];

  /* k1,k2 are global and set at the start */

  /* set species */
  rg=sps[0]; ri = sps[1]; g = sps[2];
  i = sps[3]; G = sps[4]; I = sps[5];

  while(t<time_max) {
    ran=gsl_ran_flat(r,0,1);
    total= I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + 
      i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8] +
      rg*pars[9] + I*pars[10] + G*pars[11];

    t += gsl_ran_exponential(r,1/(total));
    if(t>time_max){
      break;
    }
    if(ran<((I*i*pars[0])/total)) {
      I--;
      i--;
    } else if(ran >((I*i*pars[0])/total) && ran<((I*i*pars[0] + (k2-i)*pars[1])/total)) { 
      I++;
      i++;
    } else if(ran>((I*i*pars[0] + (k2-i)*pars[1])/total) && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2])/total){
      I--;
      g--;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3])/total) {
      I++;
      g++;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4])/total) {
      ri++;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5])/total) {
      I++;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6])/total) {
      rg++;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7])/total) {
      G++;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8])/total) {
      ri--;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8] + rg*pars[9])/total) {
      rg--;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8] + rg*pars[9])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8] + rg*pars[9] + I*pars[10])/total) {
      I--;
    } else if(ran>(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8] + rg*pars[9] + I*pars[10])/total && ran<(I*i*pars[0] + (k2-i)*pars[1] + I*g*pars[2] + (k1-g)*pars[3] + i*pars[4] + ri*pars[5] + g*pars[6] + rg*pars[7] + ri*pars[8] + rg*pars[9] + I*pars[10] + G*pars[11])/total) {
      G--;
    } else {
      printf("end of options reached in simulate_t\n ran=%f\n total=%f",ran,total);
      exit(EXIT_FAILURE);
    }
  }
    
       
    sps[0] = rg; sps[1] = ri; sps[2] = g;
    sps[3] = i; sps[4] = G; sps[5] = I;

}
