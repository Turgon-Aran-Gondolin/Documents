/**************************************************************************
 *
 * File: init.h
 *
 * This is a: header file for init.c
 *
 *
 * Author/s: Jenifer Gonzalez Lopez 
 *           Andreas Nube
 *
 * 
 * Berlin...Summer_2008
 * **************************************************************************/



#ifndef INIT_H
#define INIT_H

/* initialize paths cold */
void init_cold();

/* initialize paths hot */
void init_hot(double a);

/* initizlize momentum gaussian */
void init_momentum(double *P);

#endif /*INIT_H*/
