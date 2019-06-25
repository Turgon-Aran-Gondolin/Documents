/**************************************************************************
 *
 * File: harm_osc.c
 *
 *
 * It calculates: generates trajectories and calculates observables on the fly
 * 
 * Author Karl Jansen
 * adopted from a HMC code by Gonzalez Lopez, Jenifer, Nube, Andreas and Renner, Dru
 *
 *
 *  Berlin Summer ... 2015
 * **************************************************************************/


#define MAINPROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#include "init.h"
#include "funcHO.h"
#include "global.h"

int main( int argc, char* argv[] ){

  /* "argc"   is the number of arguments you pass to the programm: integer number */
  /* "argv[]" are the arguments you pass to the programm: string */
  /* argv[0] is always the name of the programm so if you pass "n" arguments then "argc = n+1" */


  int n_MCsteps;
  int n_MCmeas=0;

  int seed; /* used for generating the random number */

  double x2_CF;
  double x_N,R;
  double acceptance=0.0;
  double total_x=0.0;
  double total_sdx=0.0;
  double total_x2=0.0;
  double total_sdx2=0.0;
  double total_acc=0.0;
  double total_sdacc=0.0;
  double error=0.0;
  double x1=0.0;
  double x2=0.0;
  FILE *cfout;
  FILE *Rout;
  int t;


  /* the next if checks the number of arguments that have been passed to the program when calling it */
  /* if the number of arguments is different from the one expected then it finishes and prints the error message */
  /* for calling the program we use the script file: "run.sh" */

  if( argc != 9 )
    {
      printf( "%s <N> <seed> <g_num_trajectories> <g_init_method> <g_a> <g_M_0> <g_mu2> \n", argv[0] );
      return 1;
    }

  /** define arguments whose value will be asigned in the call of the programm in the "run.sh" script file **/

  N = atoi(argv[1]);
  seed = atoi(argv[2]);
  g_num_MCsteps = atoi(argv[3]);
  g_init_method = atoi(argv[4]); /* 0 = cold  ;  1 = hot */
  
  /* (global) parameters for HO (Creutz) action */
  g_a = atof(argv[5]);
  g_M_0= atof(argv[6]);
  g_mu2= atof(argv[7]);
  
  /* number of skipped trajetories */
  g_num_MCskip = atoi(argv[8]);

  printf( "N = %i\n", N );
  printf( "seed = %i\n", seed );
  printf( "Numb_skip = %i\n", g_num_MCskip );
  printf( "Numb_MCsteps = %i\n", g_num_MCsteps );
  printf( "g_init_method = %i\n", g_init_method );
  printf( "g_a = %f\n", g_a );
  printf( "g_M_0 = %f\n", g_M_0 );
  printf( "g_mu2 = %f\n", g_mu2 );


  /** initialize random number generator **/

  srand(seed);


  /************************
   * allocation of memory
   ***********************/

  /** allocation of memory for PATHS                   **/
  /** g_X is global                                    **/
  /** g_X[i] with i in 0 ... (N-1) are entries of g_X  **/

  g_X=(double*)malloc(sizeof(double)*N);

  if(g_X==NULL) {
    fprintf(stderr,"allocation of memory failed");
    exit(-1);
  }

  /*****************************************/
  /* Compute theoretical value of harmonic */
  /* oscillator from Creutz-Freedman       */
  /*****************************************/

   x_N=(double) N;

   R=1. + 0.50*g_a*g_a*g_mu2 - g_a*sqrt(g_mu2) * sqrt(1.0 + 0.25 *g_a*g_a*g_mu2);

   x2_CF=(1.0 + pow(R,x_N))/(1-pow(R,x_N)); 

   x2_CF=x2_CF/(2.0*g_M_0*sqrt(g_mu2)*sqrt(1.0 + 0.25*g_a*g_a*g_mu2));

   cfout=fopen("cfoutput","w");
   fprintf(cfout," X^2= %14.7e\n",x2_CF);
   fclose(cfout); 

   
   /***********************/
   /* initialize path g_X */
   /***********************/

  if(g_init_method==0){

    init_cold();  /* for all i: g_X[i] = 1 */

  }


  else if(g_init_method==1){

    init_hot(1);  /* all g_X[i] are set randomly */

  }



  /********************************/
  /********************************/
  /* Main Loop over MC steps (MD) */
  /********************************/
  /********************************/

  Rout=fopen("Routput","w");
  fprintf(Rout,"No. of path, <X>, <X^2>, Acceptance\n");

  for(n_MCsteps=1; n_MCsteps <= g_num_MCsteps; n_MCsteps++){

    ///////* TO DO */////////
    /* change path g_X according to the harmonic oscillator action:
       modify function S_ho() in file funcH0.c */
    acceptance = S_ho();


    /**********************************************************/
    /**  perform measurements, calculate observables on g_X  **/
    /**********************************************************/ 

    ///////* TO DO *//////////
    /* Initialization */
    x1 =   ;
    x2 =   ;

    for( t=0; t<N; t++ )
    {
      ///////* TO DO *//////////
      /* Compute here average x and average x^2 for the given trajectory */


    }


    /** for debugging **/
    /* printf( "x^1 %i %14.7e\n", ntraj, x1 ); */
    /* printf( "x^2 %i %14.7e\n", ntraj, x2 ); */
    /* printf( "acceptance %i %14.7e\n", ntraj, acceptance ); */

    if (n_MCsteps > g_num_MCskip ) 
    {
      /* print observables to file */
      fprintf(Rout," %i %14.7e %14.7e %14.7e\n", n_MCsteps, x1, x2, acceptance);

      /* sum observables over all measuring MC steps */
      total_x     = total_x     + x1;
      total_sdx   = total_sdx   + x1*x1;
      total_x2    = total_x2    + x2;
      total_sdx2  = total_sdx2  + x2*x2;
      total_acc   = total_acc   + acceptance;
      total_sdacc = total_sdacc + acceptance*acceptance;
      n_MCmeas    = n_MCmeas    + 1;
    }
    
    

  } /* end of loop over trajectories */


  fclose(Rout); 
  printf( "total number of measurements %i \n", n_MCmeas );
  


  /****************************************************/
  /* Calculate and print final results of observables */
  /****************************************************/

  ///////* TO DO *//////////
  /* Compute total average of x,x^2 and acceptance including errors
     and print values */
  


  printf( "<x>= %14.7e d(<x>)= %14.7e \n", total_x, error );



  printf( "<x^2>= %14.7e d(<x^2>)= %14.7e theor= %14.7e \n", total_x2, error, x2_CF );



  printf( "<Acc>= %14.7e d(<Acc>)= %14.7e \n", total_acc, error );


  /***************/
  /* free memory */
  /***************/ 

  /** free memory for path **/
  free(g_X);

  return 0;

}




