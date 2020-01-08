#include "imd.h"
// --------------------------------------------
// ZWECK: 1d-minimierer. Die GSL-alternative ist nicht gut, denn hier muss man unbedingt einen Schätzwert mitliefern, 
//	  dessen Bestimmung oftmals nicht einfach ist. Passt der Schätzwert nicht, meckert gsl und bricht ab.
// QELLE: fminbnd.m von octave
double fminbnd(double a, double b,double (* f)(double,double,double), double tolx,double rho,double eng)
{
	tolx=1e-6;
	double c=0.5*(3-sqrt(5));
	double v=a+c*(b-a);
	double w=v;
	double x=v;
	double e=0;
	double fv=(*f) (x,rho,eng);
	double fw=fv;
	double fval=fv;
	double fu=0.0;
	
	int i;
	int imax=500;
	double xm=0.0;
	double tol=0.0;
	double sqrteps=2.220446049250313e-16; 
	int dogs=0;
	double r,q,p,d,u;
	for(i=0;i<imax+1;i++)
	{
		xm=0.5*(a+b);
		tol=2.0*sqrteps*ABS(x)+ tolx/3.0;
//printf("i:%d, x:%.4e,a:%.4e,b:%.4e\n", i,x,a,b);			
		if(ABS(x-xm)<= (2*tol-0.5*(b-a)) || i==imax)
		{
if(i==imax)
printf("myid:%d, iter>iter_max in fminbnd for rho=%.4e. Return x=%.4e \n",myid,rho,x);
			return x;
		}		

		if(ABS(e)>tol)
		{
			dogs=0;
			r=(x-w)*(fval-fv);
			q=(x-v)*(fval-fw);
			p=(x-v)*q-(x-w)*r;
			q=2*(q-r);
			p*=-SIGNUM(q);
			q=ABS(q);
			r=e;
			e=d;
			if(ABS(p)<ABS(0.5*q*r) && p > q*(a-x) && p < q*(b-x))
			{
				d=p/q;
				u=x+d;
			
				if(MIN(u-a,b-u) < 2*tol)
				{
					d=tol*(SIGNUM(xm-x) + (xm==x));
				}
			}
			else
			{
				dogs=1;
			}		
		}
		else
		{
			dogs=1;
		}
		if(dogs==1)
		{
			e= (x >= xm ? a-x : b-x);
			d=c*e;
		}
		u=x+MAX(ABS(d),tol)*(SIGNUM(d)+(d==0));
		fu=(*f)(u,rho,eng);
		
		if(fu<fval)
		{
			if(u<x)
				b=x;
			else
				a=x;
			v=w;fv=fw;
			w=x;fw=fval;
			x=u;fval=fu;			
		}
		else
		{
			if(u<x)
				a=u;
			else
				b=u;
			if(fu<=fw  || w==x)
			{
				v=w;fv=fw;
				w=u;fw=fu;
			}
			else if(fu <= fv || v==x || v==w)
			{
				v=u;
				fv=fu;
			}
		}
		
	}
}
