#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
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

int main(){
	int Nt, Nconf;
	int Nskip, Ndump;
	double aomega, delta;
	char input_path[200], output_path[200];

// read in input
	read_double("aomega", &aomega);
	read_int("Nt", &Nt); // number of lattice sites 
	read_int("Nconf", &Nconf);
	read_int("Nskip", &Nskip);
	read_int("Ndump", &Ndump);
	read_double("delta", &delta);
	read_string("input_path", input_path);
	read_string("output_path", output_path);

	int nsample = (Nconf - Ndump)/Nskip;
	double* corr = new double[Nt*nsample];
	double* x = new double[Nt*Nconf];

// read in configurations
	char infile[300];
	sprintf(infile, "%s/conf_aomega%.4f_Nt%i_delta%.4f", input_path, aomega, Nt, delta);
	FILE *fp = NULL;
	if((fp = fopen(infile, "rb")) == NULL){
		std::cout << "failed to open input file: " << infile << "\n" << std::endl;
        	exit(0);
	}

	fread(x, sizeof(double), Nt*Nconf, fp);
	for(int i=0; i<Nt; i++)
		std::cout<<*(x+i)<<std::endl;

	char outfile[300];
	std::ofstream out;
	sprintf(outfile, "%s/corr_aomega%.4f_Nt%i_delta%.4f", output_path, aomega, Nt, delta);
	out.open(outfile);
	if(out.is_open()) {
		out<<Nsample<<" "<<Nt<<std::endl;
	}
	else{std::cout<<"error opening output file"<<std::endl; exit(0);}
	
	for( int n=Ndump; n<Nconf; n++){
		if((n-Ndump)%Nskip == 0){		 
			for(int i=0; i<Nt; i++){
				*(corr + i) = 0.0;
				for(int t=0; t<Nt; t++)
					*(corr+i) += *(x+n*Nt+i)*(*(x+n*Nt+(i+t)%Nt));
				*(corr+i)/=(double)Nt;
			}

			for(int t=0; t<Nt; t++)
				out<<std::setprecision(16)<<t<<" "<<*(corr+t)<<std::endl;
		}	
	}
	out.close();
}


