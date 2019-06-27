/****************************
Author: Liuming Liu
Data: 2019.01.06
Description:
Using metropolis algorithm to generator configurations for 1D harmornic oscillator.
Weight factor exp(-S/hbar)
S/hbar = \sum_n mass*(x(n+1) - x(n))^2/(2.0*a*hbar) + a*mass*omega^2*(x(n+1)^2 + x(n)^2)/(4.0*hbar)
here a is the lattice spacing in time direction.
Introduce a dimensionless variable u = sqrt(mass*omega/hbar) x, rewrite the action
S/hbar = \sum_n (u(n+1) - u(n))^2/(2.0*a*omega) + a*omega*(u(n+1)^2 + u(n)^2)/4.0 

The code read in 5 parameters from an input file:
a*omega
Nt: number of lattice site
Nconf: number of configurations
delta: step size for update
output_path: the directory of the output file 
*****************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <string.h>

void read_int(char const * name, int * value){
        char check[60];
        scanf("%s %i", check, value);

        if(strcmp(check,name)!=0){
        printf("error input %s \n", check);
        exit(0);
        }
}

void read_double(char const * name, double * value){
        char check[60];
        scanf("%s %lf", check, value);

        if(strcmp(check,name)!=0){
        printf("error input %s \n", check);
        exit(0);
        }
}

void read_string(char const * name, char * value){
        char check[60];
        scanf("%s %s", check, value);

        if(strcmp(check,name)!=0){
        printf("error input %s %s %s\n", check, name, value);
        exit(0);
        }
}



// compute the weight factor, update one site x1, the change of the action only involve the two neighboring sites x0(backward) and x2 (forward)
double weight(double x0, double x1, double x2, double aomega){
		return exp(-(pow(x0 - x1,2)/(2.0*aomega) + aomega*(x1*x1 + x0*x0)/4.0 ))*exp(-(pow(x2 - x1,2)/(2.0*aomega) + aomega*(x1*x1 + x2*x2)/4.0 ));
} 

int main(){
	int Nt,Nconf; 
	double aomega; // a(lattice spacing)*omega(angular frequency)
	double delta;  //update step size
	char output_path[200]; // output file path

// read in input
	read_double("aomega", &aomega);
	read_int("Nt", &Nt);
	read_int("Nconf", &Nconf);
	read_double("delta", &delta);
	read_string("output_path", output_path);
	//printf("%d %f %d",Nt,aomega,Nconf);
	printf("read complete\n");
	
	double* x1 = new double[Nt];
	double* x2 = new double[Nt];
	double* x = new double[Nt*Nconf];
	std::ofstream out;

        std::default_random_engine seed;
        seed.seed(time(NULL));
        std::uniform_real_distribution<double> ran1(-0.1,0.1);
        std::uniform_real_distribution<double> ran2(-1.0,1.0);
        std::uniform_real_distribution<double> ran3(0,1.0);
	

// initialize the configuration to be a random number between [-0.1, 0.1]	
	for(int i=0; i<Nt; i++)
//		*(x1+i) = *(x2+i) = ran1(seed);
		*(x1+i) = *(x2+i) = 0.0;  // cold start

// update the configuration with metropolis algorithm
	std::vector<int> t_iter;
	for(int i=0; i<Nt; i++) t_iter.push_back(i);
	for(int n=0; n<Nconf; n++){ 
		std::random_shuffle (t_iter.begin(), t_iter.end());  //shuffle the order of lattice sites for update
		#pragma omp parallel for
		for(int t=0; t<Nt; t++){
			int i=t_iter[t];
			*(x1+i)=*(x2+i)+ran2(seed)*delta;
			int i0=(i-1+Nt)%Nt;
			int i2=(i+1+Nt)%Nt;
			//printf("%f\n", (weight(*(x2+i0),*(x2+i),*(x2+i2),aomega)/weight(*(x1+i0),*(x1+i),*(x1+i2),aomega)));
			*(x2+i)=(ran3(seed)<(weight(*(x1+i0),*(x1+i),*(x1+i2),aomega)/weight(*(x2+i0),*(x2+i),*(x2+i2),aomega)))
				?(*(x1+i)):(*(x2+i));
			*(x1+i)=*(x2+i);
// choose a candidate value *(x2+i)

// if condition to accept or reject the candidate 

		}
	
		for(int i=0; i<Nt; i++)
			*(x + n*Nt + i) = *(x2 + i);
	}
	printf("build complete\n");

// output
	FILE *fp = NULL;
	char outfile[300];	
	sprintf(outfile, "%s/conf_aomega%.4f_Nt%i_delta%.4f", output_path, aomega, Nt, delta);
        if((fp = fopen(outfile, "wb")) == NULL){
                std::cout<<"fail to open output file."<<std::endl;
                exit(0);
        }
        fwrite((double*) x, sizeof(double), Nt*Nconf, fp);
        fclose(fp);

}


