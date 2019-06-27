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

struct stat_obj
{
	double aver;
	double error;
};

struct stat_obj get_aver(std::vector<double> &res)
{
	stat_obj tmp;

	//calculate the average and standard deviation here
	
	return tmp;
}

double PI_1D(std::vector<double> &k)
{
	if(k.size()!=1) return 0;
	return sqrt(1-k[0]*k[0]);
	//Define the integrand of 1D Monte Carlo here
}

double PI_2D(std::vector<double> &k)
{
	//Define the integrand of 2D Monte Carlo here
}

void MonteCarlo(std::vector<double> &res,int dim,double (*pfunc)(std::vector<double> &k))
{
	if(dim<1||dim>20){printf("the dimension of the integration is incorrect\n"); return;}
	std::vector<double> x(dim);
	for(int i=0;i<res.size();i++)
	{
		//set the value of each sample here
	}
}

int main (int argc, char* argv[])
{
        size_t calls=0;
        if(argc!=2) {printf("Usage: ./Nintegrate number_of_points\n");return 1;}
        if(argc==2) calls=atoi(argv[1]);
        if(calls<100) {printf("Number of the points should be at least 100\n");calls=100;}
        if(calls>1000000) {printf("Very large number is not recommended, will use 1000000\n");calls=1000000;}
	
	srand((unsigned)time(NULL));
	std::vector<double> res(calls);
	MonteCarlo(res,1,PI_1D);
	stat_obj tmp=get_aver(res);
	printf("aver=%13.5f(%8.5f)\n",tmp.aver*4,tmp.error*4);

        MonteCarlo(res,2,PI_2D);
        stat_obj tmp2=get_aver(res);
        printf("aver=%13.5f(%8.5f)\n",tmp2.aver*4,tmp2.error*4);	

	return 0;
}
