/*
CC -o Nintegrate Nintegrate.cpp 
*/

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>
double mass=1.0;
double omega=1.0;
double x00=0.0;
double x01=0.0;
double xu=3.0;
double xl=-3.0;

inline double Z(double x0,double x1,double t)
{
	//Code in the propagator here;
}


double g (std::vector<double> &k,double dt)
{
	//Code in the integrand here;
}

double multi_dim_integration(int idim,std::vector<double> &x,double dx,double (*pfunc)(std::vector<double> &k,double dt),double dt)
{
	double value=0;
	for(double x_tmp=xl+0.5*dx;x_tmp<xu;x_tmp+=dx)
	{
		x[idim]=x_tmp;
		if(idim==x.size()-1)
		{
			// do the calculation here;
		}
		else
		{
			// call itself with idim=idim +1;
		}
	}
	return value;
}

double multi_dim_integration(int N,double dx,double (*pfunc)(std::vector<double> &k,double dt),double dt)
{
	std::vector<double> x(N);
	return multi_dim_integration(0,x,dx,pfunc,dt);
}

int main (int argc, char* argv[])
{

        size_t calls=0;
        if(argc==2) calls=atoi(argv[1]);
        if(calls<1000&&calls>0) calls=1000;
        if(argc>2) {printf("Usage: ./Nintegrate number_of_points\n");return 1;}
        double value;

        double x_01=0.0;
        double x_00=0.0;
        double tau=4.0;

        printf("Z(0.0,0.0,4.0) = %14.6f\n",Z(x_01,x_00,tau));

        std::vector<double> x(10);
        for(int i=0;i<10;i++) x[i]=0.1*i;
        double dt=0.1;
        printf("Integrand = %14.6f\n",g(x,dt));

	if(calls>0)
	{

//-----------
		for(int T=4;T<10;T+=2)
		for(int dimen=1;dimen<10;dimen+=3)	    
		{
			printf("T=%6d, N=%6d\n",T,dimen);
			printf("Value (Analytical ) = %14.6f\n",Z(x00,x01,1.0*T));		
// set dimension of integration;	
			double dt=T*1.0/(dimen+1);
			int calls_per_dimen=(int)pow(calls,1.0/dimen)+1;
			double dx=(xu-xl)/calls_per_dimen;

// Sandard numerical integration;
			value=multi_dim_integration(dimen,dx,g,dt);
			printf("Value (Standard   ) = %14.6f\n",value);
		}	
	}
    return 0;
}
