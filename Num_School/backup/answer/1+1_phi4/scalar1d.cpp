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

void printvector(mat &v)
{
	int Nt = v.size();
	int Nx = (v[0]).size();
	printf("Scalar Field is:\n");
	for(int it = 0; it < Nt; it++)
	{
		for(int ix = 0; ix < Nx; ix++)
		        printf("%13.4f",v[it][ix]);
                printf("\n");
	}
	printf("over\n");
}

double S(mat &lat, int T, int X, double lam, double msq, double latnew)
{
	int nt = lat.size();
	int nx = (lat[0]).size();
        int Xb=(X==0)?(nx-1):(X-1);
        int Xf=(X==nx-1)?0:(X+1);
        int Tb=(T==0)?(nt-1):(T-1);
        int Tf=(T==nt-1)?0:(T+1);
        return -latnew*(lat[X][Tb]+lat[X][Tf]+lat[Xb][T]+lat[Xf][T])/2+(2-msq/2)*pow(latnew,2)-lam*pow(latnew,4);
}

int sample(mat &lat, double lam, double msq, int T, int X, double D)
{
        double smp = lat[T][X] + rand()*(2*D)/RAND_MAX-D;
        double dS=S(lat,T,X,lam,msq,smp)-S(lat,T,X,lam,msq,lat[T][X]);

        double P = min(1., exp(-dS));

        double judge = rand()*1.0/RAND_MAX;     
        if (judge < P)
        {
                lat[T][X] = smp;
                return 1;
        }
        return 0;
}

double sweep(mat &lat, double lam, double msq, double D)
{
	//shuffle the order of the item in the lattice in one single sweep
	int nt = lat.size();
	int nx = (lat[0]).size();
	vector<int> ordert(nt);
	vector<int> orderx(nx);
	for(int it=0;it<nt;it++) ordert.push_back(it);
	for(int ix=0;ix<nx;ix++) orderx.push_back(ix);
	std::random_shuffle(ordert.begin(), ordert.end());
	std:random_shuffle(orderx.begin(), orderx.end());

	double acp = 0.;
	for(int it=0; it < nt; it++)
	{
		for(int ix=0; ix < nx; ix++)
		        acp += sample(lat, lam, msq, ordert[it], orderx[ix], D);
	}
	orderx.clear();
	ordert.clear();
	return acp *1.0 / (nt * nx);
}

void qft(int nt, int nx, double a, double m, double g, double D)
{
	int Ns = 1000; //the amount of samples 
	int swpl = 50; //sample interval 
  	int nk = 5;
  	int nt_max = (int)(2.0/a);

	double msq = pow(m*a,2);
	double lam = g *a*a; 
	
	mat lat(nt, vec(nx, 0.0));
	
//	srand((unsigned)time(NULL));

  	double acp = 1.;
  	//warm up
  	for (int i = 0; i < 1000; i++)
  	{
  		acp = sweep(lat, lam, msq, D);
  		if(i%100==0)
  		{
  			double phi2=0.0;
  			for(int it=0;it<nt;it++)
  			for(int ix=0;ix<nx;ix++)
  				phi2+=pow(lat[it][ix],2);
  			printf("acp=%13.2f, phi2=%13.2f\n",acp,phi2/(nt*nx));
		}
		  //cout << acp << endl;
	}
  	
  	
  	tensor c2(nt_max,mat(nk,vec(Ns,0.0)));
  	vector<vector<cmplx> > oper(nt, vector<cmplx>(nk, 0.0));
	//sampling
	int count = 0;
	int Nswp = 0;
	while(count < Ns)
	{
//		acp = sweep(lat, lam, msq, D);
		if(Nswp % swpl == 0)
		{
       			for(int it = 0; it < nt; it++)
		        for(int ik = 0; ik < nk; ik++)
			{
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
			if(count%100==0)printf("%6d%13.5f%13.5f\n",count,c2[0][0][count]/(nt*nx),c2[1][0][count]/(nt*nx));fflush(stdout);
			count++;
		}
		Nswp++;
	}
	
	//effective mass;
	vec tmp(Ns);
	for(int it2 = 0; it2 < nt_max; it2++)
	for(int ik = 0; ik < nk; ik++)
	{
                if(ik==0)printf("%4d",it2);
	        jackknife(c2[it2][ik],tmp);
	        c2[it2][ik]=tmp;
	        stat_obj corr=get_aver(tmp,1);
	        printf("%13.4f(%7.4f)",corr.aver,corr.error);
	        if(ik==nk-1)printf("\n");
	        
        }
        printf("\n");
        for(int it2 = 0; it2 < nt_max-1; it2++)
        for(int ik = 0; ik < nk; ik++)
        {
                if(ik==0)printf("%4d",it2);
                for(int is = 0; is < Ns; is++)
                        tmp[is] = log(c2[it2][ik][is]/c2[it2+1][ik][is]);
                stat_obj mass=get_aver(tmp,1);
                printf("%13.4f(%7.4f)",mass.aver,mass.error);
                if(ik==nk-1)printf("\n");
	}
}

int main()
{
	double a = 0.1;
	int nt = 100;
	int nx = 100;
	double m = 1.;
//	double g = 1e-3;
	double g = 1;
	double D = 1.3;
	qft(nt, nx, a, m, g, D);
	
	return 0;
}
