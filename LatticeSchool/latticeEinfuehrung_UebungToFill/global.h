
#ifdef MAINPROGRAM
#define EXTERN
#else
#define EXTERN extern
#endif


/** global parameters **/

/* Time Extend */

EXTERN int N;


/* Parameters for HO (Creutz) Action */

EXTERN double g_M_0;
EXTERN double g_a;
EXTERN double g_mu2;

/* algorithmic parameters */

EXTERN int g_init_method; /* 0= cold  ;  1=hot */
EXTERN int g_num_MCsteps;
EXTERN int g_num_MCskip;


/* Number of Paths */


/* #define NOPATHS 2 */


/* Paths */

EXTERN double *g_X;
/* EXTERN double *g_X_copy; */


/* /\* conjugated momenta (paths) *\/ */

/* EXTERN double *g_P; */
/* EXTERN double *g_P_copy; */


/* /\* derivatives of paths *\/ */

/* EXTERN double *g_d2X;   /\* second order derivative*\/ */







