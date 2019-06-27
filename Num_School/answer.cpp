#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <complex>

using namespace std;

#define cmplx complex<double>
#define vec vector<double>
#define mat vector<vector<double> >
#define tensor vector<vector<vector<double> > >

double a=0.1;

struct stat_obj
{
        double aver;
        double error;
};

void jackknife(vec &res,vec &dest)
{
        int size=res.size();
        dest.resize(size);
        double sum=0;
        for(int i=0;i<size;i++)
                sum+=res[i];

        for(int i=0;i<size;i++)
                dest[i]=(sum-res[i])/(size-1);
}

struct stat_obj get_aver(vec &res,int type)
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

double part_of_S(mat &lat, int T, int X, double lam, double msq)
//the part of action which depends on lat[T][X]
{
	int nt = lat.size();
	int nx = (lat[0]).size();
        int Xb=(X==0)?(nx-1):(X-1);
        int Xf=(X==nx-1)?0:(X+1);
        int Tb=(T==0)?(nt-1):(T-1);
        int Tf=(T==nt-1)?0:(T+1);
        return  pow(lat[T][X]-lat[Tf][X],2)+pow(lat[Tb][X]-lat[T][X],2)
	       +pow(lat[T][X]-lat[T][Xf],2)+pow(lat[T][Xb]-lat[T][X],2)
	       +msq*pow(lat[T][X],2)+lam*pow(lat[T][X],4);
}

int sample(mat &lat, mat &lat_tmp, double lam, double msq, int T, int X, double D)
//Determine whether the update will be accepted
{
	lat_tmp[T][X]=lat[T][X]+rand()*D*2.0/RAND_MAX-D;

	//calculate the change of the action, and calculation the acceptation rate;	
	double dS=part_of_S(lat_tmp,T,X,lam,msq)-part_of_S(lat,T,X,lam,msq);

        //deside whether the update can be accepted.
        double judge = rand()*1.0/RAND_MAX;     
        if (judge < exp(-dS))
        {
        	//update the configuration
                lat[T][X] = lat_tmp[T][X];
                return 1;
        }
        else
        {
        	//discard the value of the temporary configuration
        	lat_tmp[T][X]=lat[T][X];
        	return 0;
	}
}

double iteration(mat &lat, double lam, double msq, double D)
{
	//shuffle the order of the sites on the lattice
	int nt = lat.size();
	int nx = (lat[0]).size();
	vector<int> ordert(nt);
	vector<int> orderx(nx);
	for(int it=0;it<nt;it++) ordert[it]=it;
	for(int ix=0;ix<nx;ix++) orderx[ix]=ix;
	random_shuffle(ordert.begin(), ordert.end());
	random_shuffle(orderx.begin(), orderx.end());

	//the temporary configuration
	mat lat_tmp(nt, vec(nx, 0.0));
	for(int it=0; it < nt; it++)
	for(int ix=0; ix < nx; ix++)
		lat_tmp[it][ix]=lat[it][ix];
	
	//update the configuration
	double acp = 0.;
	for(int it=0; it < nt; it++)
	for(int ix=0; ix < nx; ix++)
	        acp += sample(lat, lat_tmp, lam, msq, ordert[it], orderx[ix], D);
	orderx.clear();
	ordert.clear();
	
	//return the averaged acceptation rate
	return acp *1.0 / (nt * nx);
}

void measure(int nt, int nx, double a, double m, double lambda, double D)
{
	int Ns = 1000; //the number of the samples 
	int swpl = 100; //sample interval 
  	int nk = nx/8+1; // the maxinum momentum, p_i=i*M_PI*2.0/nx, i=0,1,...nk-1;
  	
  	int nt_max = (int)(1.0/a);

	double msq = pow(m*a,2);
	double lam = lambda *a*a; 
	
	mat lat(nt, vec(nx, 0.0));
	
	srand((unsigned)time(NULL));

  	//Warm up
  	double acp = 1.;
	acp = iteration(lat, lam, msq, D);
  	for (int i = 0; i <= 100/acp; i++)//if the acceptation rate is low, do a longer warm up.
  	{
  		acp = iteration(lat, lam, msq, D);
  		if(i%20==0)
  		{
  			double phi2=0.0;
  			for(int it=0;it<nt;it++)
  			for(int ix=0;ix<nx;ix++)
  				phi2+=pow(lat[it][ix],2);
  			printf("acp=%13.2f, phi2=%13.2f\n",acp,phi2/(nt*nx));
		}
	}
  	
  	tensor c2(nt_max,mat(nk,vec(Ns,0.0)));
  	vector<vector<cmplx> > oper(nt, vector<cmplx>(nk, 0.0));
	//sampling
	for(int count=0;count<Ns;count++)
	{
		for(int iter=0;iter<swpl;iter++) iteration(lat, lam, msq, D);
		for(int it = 0; it < nt; it++)
		for(int ik = 0; ik < nk; ik++)
		{
			//Calculate the operator O(t,k);
		        oper[it][ik]=0.0;
        		for(int ix = 0; ix < nx; ix++)
			{
			        double eta=ik*ix*M_PI*2.0/nx;
			        oper[it][ik]+=lat[it][ix]*cmplx(cos(eta),sin(eta));
                        }
                }
                for(int it2 = 0; it2 < nt_max; it2++)
                for(int it = 0; it < nt; it++)
                for(int ik = 0; ik < nk; ik++)
			c2[it2][ik][count]+=cmplx(conj(oper[it][ik])*oper[(it+it2)%nt][ik]).real();
                        
                //Calculate the correlation function;
		if(count%100==0)printf("%6d%13.5f%13.5f\n",count,c2[0][0][count]/(nt*nx),c2[1][0][count]/(nt*nx));fflush(stdout);
	}
	
	//Calculate the averaged correlation function;
	printf("Correlation function: t S(t,0), S(t,2Pi/nx)...\n");
	vec tmp(Ns);
	for(int it2 = 0; it2 < nt_max; it2++)
	for(int ik = 0; ik < nk; ik++)
	{
                if(ik==0)printf("%4d",it2);
	        jackknife(c2[it2][ik],tmp);
	        c2[it2][ik]=tmp;
	        stat_obj corr=get_aver(tmp,1);
	        printf("%13.3f(%7.3f)",corr.aver/(nx*nt),corr.error/(nx*nt));
	        if(ik==nk-1)printf("\n");
	        
        }
        printf("\n");
	printf("Effective mass: t m(t,0), m(t,2Pi/nx)...\n");

	//effective mass;
        stat_obj mass_0;
        for(int it2 = 0; it2 < nt_max-1; it2++)
        for(int ik = 0; ik < nk; ik++)
        {
                if(ik==0)printf("%4d",it2);
                for(int is = 0; is < Ns; is++)
                        tmp[is] = log(c2[it2][ik][is]/c2[it2+1][ik][is]);
                stat_obj mass=get_aver(tmp,1);
                printf("%13.3f(%7.3f)",mass.aver,mass.error);
                if(ik==0&&it2==nt_max/2) mass_0=mass;
                if(ik==nk-1)printf("\n");
	}
	
	
	printf("\n");
	printf("Expected energy from the mass: t E(0), E(2Pi/nx)...\n");
	printf("%4d",-1);
	for(int ik = 0; ik < nk; ik++)
	{
		double p=ik*M_PI*2.0/nx;
		double E=sqrt(mass_0.aver*mass_0.aver+p*p);
		double dE=mass_0.error*mass_0.aver/E;
		printf("%13.3f(%7.3f)",E,dE);
	}
	printf("\n");
}

int main(int argc, char * argv[])
{
	a = 0.1;
	if(argc>1)a=atof(argv[1]);
	int nx = 25;
	if(argc>2)nx=atoi(argv[2]);
	int nt = 10/a;
	double m = 0.1/a;
	double lambda = 0.0;
	if(argc>3)lambda=atof(argv[3]);
	double D = 1.5;
	measure(nt, nx, a, m, lambda, D);
	
	return 0;
}
