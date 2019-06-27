/*
CC -o Markov Markov.cpp 
*/

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <complex>

double sig_delta=0.000001;
double x00=0.0;
double x01=0.0;
double mass=1.0;
double omega=1.0;
double dt=1.0;

inline double S(double x0, double x1)
{
        return + 0.5*log(2*M_PI*dt)
                 + 0.5*mass*pow(x0-x1,2)/dt
                 + 0.25*mass*pow(omega,2)*(x0*x0+x1*x1)*dt;
}

double g (double *k, size_t dim,void *params)
{
// This part BELOW can be charged
//        double value=-S(x00,k[0]);
        double value=0.0;
        for(int i=0;i<dim-1;i++)
                value-=S(k[i],k[i+1]);
//    return exp(value-S(k[dim-1],x01));
    return exp(value-S(k[dim-1],k[0]));;
// This part ABOVE can be charged
}

inline double Z(double x0,double x1,double t)
{
        return sqrt(mass*omega/(2*M_PI*sinh(omega*t)))*exp(-mass*omega/(2*sinh(omega*t))*((x0*x0+x1*x1)*cosh(omega*t)-2*x0*x1));
}

double accept(std::vector<double> &x,int pos,double dx)
{
        int nt=x.size();
        double x0=x[pos]+rand()*(2*dx)/RAND_MAX-dx;
        double dS=0;
        if(pos<nt-1) dS+=S(x0,x[pos+1])-S(x[pos],x[pos+1]);
        else dS+=S(x0,x[0])-S(x[pos],x[0]);
        if(pos>0) dS+=S(x[pos-1],x0)-S(x[pos-1],x[pos]);
        else dS+=S(x[nt-1],x0)-S(x[nt-1],x[pos]);
        double P=(dS>0)?exp(-dS):1.0;
        double judge=rand()*1.0/RAND_MAX;
        if(P>judge)
        {
            x[pos]=x0;
            return 1.0;
        }
        return 0.0;
}

double iteration(std::vector<double> &x,double dx)
{
        int nt=x.size();
        std::vector<int> order;
        for(int i=0;i<nt;i++) order.push_back(i);
        std::random_shuffle(order.begin(), order.end());
        double rate=0;
        for(int i=0;i<nt;i++)
                rate+=accept(x,order[i],dx);
        order.clear();
        return rate/nt;
}

double iteration2(std::vector<double> &x,double dx)
{
        int nt=x.size();
        double dS=0;
        std::vector<double> tmp(nt);
    for(int pos=0;pos<nt;pos++)
    {
        double x0=x[pos]+rand()*(2*dx)/RAND_MAX-dx;
        tmp[pos]=x0;
        if(pos<nt-1) dS+=S(x0,x[pos+1])-S(x[pos],x[pos+1]);
        else dS+=S(x0,x[0])-S(x[pos],x[0]);
        if(pos>0) dS+=S(x[pos-1],x0)-S(x[pos-1],x[pos]);
        else dS+=S(x[nt-1],x0)-S(x[nt-1],x[pos]);
    }
        double P=(dS>0)?exp(-dS):1.0;
        double judge=rand()*1.0/RAND_MAX;
        if(P>judge)
        {
            for(int pos=0;pos<nt;pos++)
                 x[pos]=tmp[pos];
            return 1.0;
        }
        return 0.0;
}


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

void bootstrap(std::vector<double> &res,std::vector<double> &dest,int nboot)
{
     int size=res.size();
     dest.resize(nboot);

     for(int ib=0;ib<nboot;ib++)
     {
         double sum=0;
         for(int i=0;i<size;i++)
              sum+=res[rand()%size];
         dest[ib]=sum/size;
     }
}

struct stat_obj
{
     double aver;
     double error;
};

struct stat_obj get_aver(std::vector<double> &res,int type)
//type=0 for the standard samples, 1 for jackknife samples, 2 for bootstrap samples
{
     int size=res.size();
     double aver=0;
     for(int i=0;i<size;i++)
          aver+=res[i]/size;
     double sig=0;
     for(int i=0;i<size;i++)
          sig+=pow(res[i]-aver,2);

     if(type==1)sig*=pow(size-1,2);
     if(type==2)sig*=size;
     stat_obj tmp;
     tmp.aver=aver;
     tmp.error=sqrt(sig/(size*(size-1)));
     return tmp;
}


int main (int argc, char* argv[])
{
//        if(argc==1) {printf("HarmonicOscillator points nt\n");return 0;}
//        if(argc>3) {printf("The argument should be two: points, nt\n");return 1;}
	size_t calls = atoi(argv[1]);
	if(calls<10) calls=10;
	double D=atof(argv[2]);
	int sep=25;
        
//-----------
// set dimension of integration here;	
	std::vector<double> xl,xu,dx;
	int nt=atoi(argv[3]);
	int bin_size=10;
//	int dimen=2000;
//	dt=nt*1.0/dimen;
	dt=0.05;
	int dimen=(int)(nt/dt);
	int Dt=(int)(1/dt);
//	int Dt=0;
	printf("%4d\n",Dt);
        std::vector<double> x;x.resize(dimen);
	printf("Value (Analytical ) = %14.7e\n",Z(x00,x01,(dimen+1)*dt));

//Markov Chain
        for(int id=0;id<dimen;id++)
             x[id]=0.0;
        //warn up;
        int iter=0;
        double acp=0.0;
        double x20=0.0;
        do
//        for(int i=0;i<30000;i++)
        {
            acp=0.0;
            x20=0.0;
            for(int i=0;i<100;i++)
            {
                acp+=iteration(x,D)/100;
//                acp+=iteration2(x,D)/100;
                for(int id=Dt;id<dimen-Dt;id++)x20+=mass*pow(omega,2)*x[id]*x[id]/((dimen-2*Dt))/100;
            }
            if(iter%10==0) printf("warn up %6d:acp=%8.2f, x2=%8.2f\n",iter,acp,x20);
            iter++;
        }
        while(x20<0.55&&iter<10000);
        printf("warn up %6d:acp=%8.2f, x2=%8.2f\n",iter,acp,x20);
        std::vector<double> x2(calls/bin_size),x_corr(calls/bin_size),x_corr_p1(calls/bin_size);
        for(int i=0;i<calls;i++)
        {
            for(int iter=0;iter<25/acp;iter++) acp=iteration(x,D);
//            for(int iter=0;iter<25/acp;iter++) acp=iteration2(x,D);
//             double S=g(x.data(),dimen,NULL);
               double x2_tmp=0.0;
               for(int id=Dt;id<dimen-Dt;id++)x2_tmp+=mass*pow(omega,2)*x[id]*x[id]/((dimen-2*Dt));
               if(i%100==0) printf("running %6d:acp=%8.2f, x2=%8.2f\n",i,acp,x2_tmp);
               x2[i/bin_size]+=x2_tmp/bin_size;
               double ga1=0,ga2=0;
               for(int id=Dt;id<dimen-2*Dt-1;id++)
               {
                  x_corr[i/bin_size]+=x[id]*x[id+Dt]/bin_size;
                  x_corr_p1[i/bin_size]+=x[id]*x[id+Dt+1]/bin_size;
               }
        }
        std::vector<double> x2_j,x_corr_j,x_corr_p1_j;
        jackknife(x2,x2_j);
        stat_obj e0=get_aver(x2_j,1);
        jackknife(x_corr,x_corr_j);
        jackknife(x_corr_p1,x_corr_p1_j);
        for(int i=0;i<calls/bin_size;i++)
            x_corr_j[i]=log(x_corr_j[i]/x_corr_p1_j[i])/dt+x2_j[i];
        stat_obj e1=get_aver(x_corr_j,1);
        printf("e0=%14.5f(%8.5f) and e1=%14.5f(%8.5f)\n",e0.aver,e0.error,e1.aver,e1.error);

	return 0;
}
