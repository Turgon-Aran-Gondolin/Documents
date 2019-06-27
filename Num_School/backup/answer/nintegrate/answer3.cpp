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
        return sqrt(mass*omega/(2*M_PI*sinh(omega*t)))*exp(-mass*omega/(2*sinh(omega*t))*((x0*x0+x1*x1)*cosh(omega*t)-2*x0*x1));
}

inline double rand_gaussian(double val,double sig)//Box-Muller, P(x)=1/sqrt(2M_PI)*exp(-x^2/2)
{
	double u=rand()*1.0/RAND_MAX,v=rand()*1.0/RAND_MAX;
	return sqrt(-2*log(u))*cos(2*M_PI*v)*sig+val;
}


inline double S(double x0, double x1,double dt)
{
        return + 0.5*log(2*M_PI*dt)
                 + 0.5*mass*pow(x0-x1,2)/dt
                 + 0.25*mass*pow(omega,2)*(x0*x0+x1*x1)*dt;
}

double g (std::vector<double> &k,double dt)
{
// This part BELOW can be charged
        int dim=k.size();
        double value=-S(x00,k[0],dt);
        for(int i=0;i<dim-1;i++)
                value-=S(k[i],k[i+1],dt);
        return exp(value-S(k[dim-1],x01,dt));
// This part ABOVE can be charged
}

struct stat_obj
{
        double aver;
        double error;
};

void jackknife(std::vector<double> &res,std::vector<double> &dest)
{
        int size=res.size();
        dest.resize(size);
        double sum=0;
        for(int i=0;i<size;i++)
                sum+=res[i];

        for(int i=0;i<size;i++)
                dest[i]=(sum-res[i])/(size-1);
}

struct stat_obj get_aver(std::vector<double> &res,int type)
//type=0 for the standard samples, 1 for jackknife samples
{
        int size=res.size();
        double aver=0;
        for(int i=0;i<size;i++)
                aver+=res[i]/size;
        double sig=0;
        for(int i=0;i<size;i++)
                sig+=pow(res[i]-aver,2);

        if(type==1)sig*=pow(size-1,2);
        stat_obj tmp;
        tmp.aver=aver;
        tmp.error=sqrt(sig/(size*(size-1)));
        return tmp;
}

int main (int argc, char* argv[])
{
	size_t calls=0;
	if(argc!=2) {printf("Usage: ./Nintegrate number_of_points\n");return 1;}
	if(argc==2) calls=atoi(argv[1]);
	if(calls<1000) calls=1000;
	double value;

//-----------
// set dimension of integration here;	
        int Nsample=1000;
	int NpSample=calls/Nsample;
	calls=NpSample*Nsample;
        std::vector<double> data[9];
	
	for(int T=1;T<10;T++)
	{
        	printf("T=%6d\n",T);
//        	int dimen=10;
//		double dt=T*1.0/(dimen+1);

                double dt=0.2;
                int dimen=(int)(T/dt)-1;
                if(dimen>10)
                {   
                    dimen=10;
                    dt=T*1.0/(dimen+1);
                }
                
		int calls_per_dimen=(int)pow(calls,1.0/dimen)+1;
        	int size=calls/Nsample;
		double dx=(xu-xl)/calls_per_dimen;
        	double dx_all=pow(xu-xl,dimen)/size;
		printf("Value (Analytical ) = %14.6f\n",Z(x00,x01,1.0*T));

//MC numerical integration;
        	srand((unsigned)time(NULL));
		data[T-1].resize(Nsample);
        	std::vector<double> x(dimen);
        	for(int i=0;i<calls;i++)
        	{
        		for(int id=0;id<dimen;id++)
        			x[id]=rand()*(xu-xl)/RAND_MAX+xl;
        		data[T-1][i%Nsample]+=g(x,dt)*dx_all;
        	}
        	stat_obj res1=get_aver(data[T-1],0);
        	printf("Value (Monte Carlo) = %14.6f(%9.6f)\n",res1.aver,res1.error);

//Vegas_simple;
        	double sigma=sqrt(dt);
        	for(int i=0;i<Nsample;i++)data[T-1][i]=0.0;
        	for(int i=0;i<calls;i++)
        	{
        		for(int id=0;id<dimen;id++)
        			x[id]=rand_gaussian(0.0,sigma);
        		double P=g(x,dt);
        		for(int id=0;id<dimen;id++)
        			P*=(sqrt(2*M_PI)*sigma)*exp(x[id]*x[id]/(2*sigma*sigma));
        		data[T-1][i%Nsample]+=P/size;
        	}
        	stat_obj res2=get_aver(data[T-1],0);
        	printf("Value (VegasSimple) = %14.6f(%9.6f)\n",res2.aver,res2.error);
	}

//effective mass;
	for(int T=0;T<9;T++)
	{
		std::vector<double> tmp;
		jackknife(data[T],tmp);
		for(int i=0;i<Nsample;i++)
			data[T][i]=tmp[i];
		if(T>0)
		{
			for(int i=0;i<Nsample;i++)	
				tmp[i]=log(data[T-1][i]/data[T][i]);
			stat_obj res=get_aver(tmp,1);
			printf("Mass[%3d] = %14.6f(%9.6f), real=%14.6f\n",T,res.aver,res.error,log(Z(x00,x01,T)/Z(x00,x01,T+1)));
		}
	}
	return 0;
}
