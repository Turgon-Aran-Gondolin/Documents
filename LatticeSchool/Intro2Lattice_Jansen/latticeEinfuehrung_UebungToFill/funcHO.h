/*********************************************************************************
 *
 * File: funcHO.h
 *
 *
 * It is the HEADER file for funcHO.c, where the external functions are defined
 *
 * Karl Jansen 
 * adopted from a HMC code by Jenifer Gonzalez Lopez, Andreas Nube, Dru Renner
 *
 * 
 * Berlin...Summer_2015
 * *******************************************************************************/


#ifndef FUNCHO_H
#define FUNCHO_H



/*************************************************
 * MATHEMATICAL functions needed for the action
 ************************************************/

double random_number();
double gauss();


void scalar_times_vector ( double *vec_out, double scal,double *vec_in);
double square_norm (double *path);
double scalar_product (double *path_1,double *path_2);

void add_vectors (double *vec_out,double *vec_in2,double *vec_in1);
void v_eq_v_p_s_t_v(double *A,double *B,double c,double *D);

/*************************************
 * Specific functions for HO action
 ************************************/

double S_ho();
void metropolis(double Snew, double Sold, int itest);

void momentum_update(double dtau);
void path_update(double dtau);

int itest_metropolis(double Sold, double Snew);


#endif /*FUNCHO_H*/


