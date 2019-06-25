/**************************************************************************
 *
 * File: init.c
 *
 * It calculates: INITIALIZATION of paths and momenta
 *
 *
 * The externally accessible functions are:
 *
 *    random_number ()    defined in funcHO.c
 *    gauss ()            defined in funcHO.c
 *
 *
 *
 * Author: Karl Jansen 
 * from a HMC code by Jenifer Gonzalez Lopez, Andreas Nube
 *
 * 
 * Berlin...Summer_2015
 * **************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "init.h"
#include "funcHO.h"


/** Paths initialization **/

/* initialize paths with zeroes */
void init_cold(){

  int i;

  for(i=0;i<N;i++){

    g_X[i] = 1.;

  }


}


/* initialize paths with random values from -a .. a */
void init_hot(double a){

  int i;

  for(i=0;i<N;i++){

    g_X[i] = a*(2.*random_number()-1.);

  }

}

