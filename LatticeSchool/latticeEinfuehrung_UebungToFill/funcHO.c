/**********************************************************************************
 *
 * File: funcHO.c
 *
 *
 * It contains: all the needed EXTERNAL FUNCTIONS (except the initialization ones)
 *
 * Author: Karl Jansen  
 * Adopted from a HMC programm by Jenifer Gonzalez Lopez, Andreas Nube and Dru Renner
 *
 * 
 * Berlin...Summer_2015
 * ********************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "init.h"
#include "funcHO.h"
#include "global.h"


/*************************************************
 * MATHEMATICAL functions needed for the action
 ************************************************/



/* return random number in the interval [0,1] */

double random_number(){

  return ((double)rand())/((double)RAND_MAX);

}

/* generates one gaussian distributed random number (taken from the Swinger code of Andreas) */

double gauss(){

  return (sqrt(-2.*log(random_number()))*cos(2.*M_PI*random_number()));
    
}



/* calculates the sum of the square == square norm of a vector ==> returns a scalar*/

double square_norm (double *path){

  int i;
  double sum=0;


  for(i=0;i<N;i++){

    sum += path[i]*path[i];

  }
  return sum;
}


/* scalar product of two vectors ==> returns a scalar */

double scalar_product (double *path_1,double *path_2){

  int i;
  double sum=0;

  for(i=0;i<N;i++){

    sum += path_1[i]*path_2[i];

  }

  return sum;
}



/* Calculates the product: (scalar)*(vector) ==> returns a vector */

void scalar_times_vector ( double *vec_out, double scal,double *vec_in){

  int i;

  for(i=0;i<N;i++){

    vec_out[i] = scal*vec_in[i];

  }

}



/* Adds two vectors ==> returns a vector */

void add_vectors (double *vec_out,double *vec_in2,double *vec_in1){

  int i;

  for(i=0;i<N;i++){

    vec_out[i] = vec_in1[i] + vec_in2[i];

  }
}


/* calculates: vector_out = vector_in1 + scalar * vector_in2 ==> returns a vector */
void v_eq_v_p_s_t_v(double *A,double *B,double c,double *D){

  int i;

  for(i=0;i<N;i++){

    A[i]= B[i] + c * D[i];

  }
}


/* Goes onces through whole the path g_X and changes it 
   according to harmonic oscillator action by performing
   Metropolis accept-reject steps,
   returns the mean acceptance of this Metropolis step */

double S_ho (){

  int i;
  double Sold;
  double Snew;
  double xnew;
  int leftneighbor;
  int rightneighbor;
  int itest;

  double acceptance = 0.0;

  for(i=0;i<N;i++){

    ///////* TO DO *//////////
    /* specify neighbouring points and boundary conditions */
    if(i == 0){
      leftneighbor  =  N-1;
    }
    else{
      leftneighbor  =  i-1;
    }
    if(i == N-1){
      rightneighbor =  0;
    }
    else{
      rightneighbor =  i+1;
    }
    
    ///////* TO DO *//////////
    /* Determine old action */ 
    Sold =  ;
    
    ///////* TO DO *//////////
    /* Compute new proposal for coordinate */ 
    xnew =  ;
    
    ///////* TO DO *//////////
    /* Determine new action */
    Snew =  ;
    
    ///////* TO DO *//////////  
    /* perform Metropolis test: modify itest_metropolis() function below */
    itest = itest_metropolis(Sold,Snew);

    
    if(itest == 1) {
      g_X[i]     = xnew;
      acceptance = acceptance + 1.0;
    }

    /** for debugging **/
    /*   printf(" Sold = %e Snew= %e \n",Sold,Snew); */  
    /*   printf(" itest = %i \n",itest); */  

  }

  /** for debugging **/
  /* printf(" Acceptance = %e \n", *acceptance); */  
  return acceptance / (double)N;
}



/* Metropolis test: 
     Accept/Reject new variable depending on calculated old and new action
     accept: return 1
     reject: return 0
 */

int itest_metropolis(double Sold, double Snew) {

  double R;
  double r;
  int itest;

   /** Accept/reject step (Metropolis) **/

   R = exp(Sold-Snew);

   ///////* TO DO *//////////
   /* Perform Metropolis test */
 


   /** for debugging **/
   /*printf(" R = %e \n",R);*/
   
   return itest;
   
}

