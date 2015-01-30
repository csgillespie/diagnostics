#include "headers.h"
#include "ODE.h"
#include "gillespie.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define MSET gsl_matrix_set
#define MGET gsl_matrix_get

/* 
   C code for reading in a csv file
 */
gsl_matrix *readMatrix(char *filename)
{
  int line_length = 500;    
  int nrows, ncols;
       
  FILE* f;  
  char *pch;
  char line[line_length];
  gsl_matrix *particles;
    
  f = fopen(filename, "r");
  if(NULL==f) {
    fprintf(stderr, "Cannot open file %s\n", filename);
    exit(1);
  }
  nrows = 0; ncols = 0;
  /*Scan once to get the dimensions
    there doesn't seem to be a realloc matrix function
  */

  while(fgets(line, line_length, f) != NULL){
    pch = strtok(line,",");
    while(nrows == 0 && pch != NULL ) {
      ncols++;
      pch = strtok(NULL,",");
    }
    nrows++;
  }
  
  fclose(f);
        
  /*Create matrix and fill up*/
  particles = gsl_matrix_alloc(nrows, ncols);
  nrows = 0; ncols = 0;
  f=fopen(filename, "r");
      
  while(fgets(line, line_length, f) != NULL){
    pch = strtok(line,",");
    while(pch != NULL ) {
      MSET(particles, nrows, ncols, atof(pch));
      ncols++;
      pch = strtok(NULL,",");
    }
    ncols = 0;
    nrows++;
  }
  fclose(f);
    
  return(particles);    
}  





int main(int argc, char *argv[])
{
  /*no_mc == #MC ODEs; tstep = time step */
  int i, j, no_mc = 27;
  double tstep=1;
  int par_index, sim_index;
  gsl_matrix *pars_mat;

  int prior = 0;
  if(prior==1) {
    pars_mat = readMatrix("../data/prior.csv");
  } else {
    pars_mat = readMatrix("../data/post.csv");
  }
  gsl_matrix *sim_mat = readMatrix("../data/sim.csv");
  
  double *par_wts, *sim_wts;
  par_wts = calloc(sizeof(double), pars_mat->size1);
  sim_wts = calloc(sizeof(double), sim_mat->size1);

  double *sps, *pars, *sps_gil;
  sps = malloc(no_mc*sizeof(double));
  sps_gil = malloc(6*sizeof(double));

  /*Add in (fixed) k1 and k2 to pars*/
  pars = malloc((pars_mat->size2+1)*sizeof(double));
  /*total g & i*/
  pars[pars_mat->size2-1] = 10;  pars[pars_mat->size2] = 2;
  gsl_ran_discrete_t *sim_sample, *par_sample;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);  
  gsl_rng_set (r, 1);   

  par_sample = gsl_ran_discrete_preproc(pars_mat->size1, (const double *) par_wts);
  sim_sample = gsl_ran_discrete_preproc(sim_mat->size1, (const double *) sim_wts);
  
  /* 1. Select parameters from post (select row)
     2. Select time point of interest (select column)
     3. Simulate from MC and Gillespie one time unit
     4. Compare
     5. Profit
  */
  
  for(i=0; i<10000; i++) {
    if(prior == 1) {
      par_index = i;
    } else {
      par_index = gsl_ran_discrete(r, par_sample);
    }
    sim_index = gsl_ran_discrete(r, sim_sample);

    /*Reset MC */
    for(j=0; j<no_mc;j++) {
      sps[j] = 0;
    }
    /* Sps order: rg, ri, g, i, G, I 
       column zero iteration number
    */
    sps[0] = MGET(sim_mat, sim_index, 1);
    sps[2] = MGET(sim_mat, sim_index, 2);
    sps[5] = MGET(sim_mat, sim_index, 3);
    sps[9] = MGET(sim_mat, sim_index, 4);
    sps[14] = MGET(sim_mat, sim_index, 5);
    sps[20] = MGET(sim_mat, sim_index, 6);

    for(j=0; j<(pars_mat->size2-1); j++) {
      pars[j] = MGET(pars_mat, par_index, j+1);
    }

    sps_gil[0] = sps[0]; sps_gil[1] = sps[2]; sps_gil[2] = sps[5];
    sps_gil[3] = sps[9]; sps_gil[4] = sps[14]; sps_gil[5] = sps[20];

    /* simulate from MC */
    ode(pars, tstep, sps);
    gillespie(sps_gil, pars, tstep, r);

    
    printf("%f, %f, %f, %f, %f, %f, %d, %d\n", 
           (sps[0] - sps_gil[0])/(pow(sps[1], 0.5)),
           (sps[2] - sps_gil[1])/(pow(sps[4], 0.5)),
           (sps[5] - sps_gil[2])/(pow(sps[8], 0.5)),
           (sps[9] - sps_gil[3])/(pow(sps[13], 0.5)),
           (sps[14] - sps_gil[4])/(pow(sps[19], 0.5)),
           (sps[20] - sps_gil[5])/(pow(sps[26], 0.5)),
           par_index, sim_index
           );
    
  }


  return(EXIT_SUCCESS);
}

