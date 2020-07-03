#include "imd.h"
#include <sys/time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


// #define USEFLOAT  // hauptsächlich in der funktion genexptint. Profiling zeigte, dass
                  // hier die meiste zeit verbraucht wird -> float verdoppelt performance

#ifdef USEFLOAT
typedef float Real;
#define REALTYPE MPI_FLOAT
#else
typedef double Real;
#define REALTYPE MPI_DOUBLE
#endif

#ifdef USEFLOAT
  #define EXPR     expf //exp zu floaten ist eine ganz mieeese idee
  #define SQRTR    sqrtf
  #define POWR     powf
  #define LOGR     logf
#else
  #define EXPR     exp
  #define SQRTR    sqrt
  #define POWR     pow
#define LOGR       log
#endif


// *********************************************************
// PHYSICAL CONSTANTS
// *********************************************************
// const double eV2J=1.6021766E-19;
const Real eV2H=0.03674932; //eV to Hartree
const Real colrad_reltol=1e-5;
const Real colrad_abstol=10.0;

// const Real J2eV=6.2415091E18;
const Real planck=6.62607004E-34; // J/s
const Real bohr_radius=0.52917721067E-10; // m
const Real bohr_radius_sq=2.800285202924816e-21;
const Real hbar_cub=1.172812163789953e-102; //hbar^3 
const Real double_emass_pow_3_2 = 2.459112949719466e-45; // (2*emass)^3/2
const int MAXLINE = 255;
const Real  pi=3.141592653589793;
const Real  pi_sq=9.869604401089358;
const Real E_ion_H=13.6; // eV
const Real E_ion_H_J=2.178960176000000e-18; // J
const Real E_ion_H_sq_J=4.747867448593952e-36;

const Real colrad_tequi=1e-12;//TEST// 1e-12; //bei initial equi ohne Temperatur-variation erst einmal 
                                //die Saha-besetzungsdichten equilibrieren

//const double  LIGHTSPEED=2.997925458e8; // m/s
Real  LASERFREQ;




int colrad_ydot(double t, N_Vector u, N_Vector udot, void *user_data);
void do_Saha(Real Te,Real totalc,Real ne,N_Vector y);
int colrad_GetCoeffs(N_Vector y,Real It, void * user_data);

// Die Zwei müssen nach Prototypes.h
// void do_colrad(double dt);
// void colrad_init(void);

void colrad_read_states(void);
void colrad_Saha_init(int i,int j,int k);




// ******************************************************************************
// *                CROSS SECTION INTEGRATION STUFF
// ******************************************************************************

gsl_integration_workspace * winteg_inner=NULL;
gsl_integration_workspace * winteg_outer=NULL;
gsl_integration_workspace * winteg_fermi=NULL;
gsl_integration_workspace * winteg_exc=NULL; //excitation

gsl_integration_romberg_workspace * winteg_rb_inner=NULL;
gsl_integration_romberg_workspace * winteg_rb_outer=NULL;


struct my_f_params { Real ne; Real T;Real mu; Real E;Real DeltaE; int allowed;};
// struct my_f_params fparams_inner; //For inner integrand
// struct my_f_params fparams_outer; //outer integrand
// struct my_f_params fparams_fermi; 
// struct my_f_params fparams_exc; 

double inner_integrand_ionization(double x, void *p); // integrate along E'
double outer_integrand_ionization(double x,void *p);  // integrate along E


Real double_integral_ionization(Real ne,Real T, Real mu, Real DeltaE); //evaluates double integral


double inner_integrand_recombination(double x, void *p);
double outer_integrand_recombination(double x,void *p);
Real double_integral_recombination(Real ne,Real T, Real mu, Real DeltaE);

double integrand_excitation(double x,void *p);
Real eval_excitation_integral(Real ne,Real T,Real mu, Real DeltaE, int allowed); 

Real eval_dexcitation_integral(Real ne,Real T,Real mu, Real DeltaE, int allowed);
double integrand_deexcitation(double x,void *p);


double fermi_integrand(double x, void *p);
Real eval_fermi_integrand(Real ne,Real T, Real mu);


double integrand_excitation_debug(double x,void *p);


double outer_integrand_ionization2(double x,struct my_f_params* p);
Real double_integral_ionization2(Real ne,Real T, Real mu, Real DeltaE); //evaluates double integral
double inner_integrand_ionization2(double x, struct my_f_params* p);

// **********************************************************************************************
// *          PAR INTEGRAL STUFF
// **********************************************************************************************
int terminate_gkq;
int terminate_gkq_outer;
int terminate_gkq_inner;
int terminate_serial;
int gkq_iter_serial; // nr of iterations

const double gkq_alpha=0.816496580927726;
const double gkq_beta=0.447213595499958;          
static const double xgkq[12] =   
{
  0.0,
  -0.942882415695480,
  -0.816496580927726,
  -0.641853342345781,
  -0.447213595499958,
  -0.236383199662150,
  0.0,
  0.236383199662150,
  0.447213595499958,
  0.641853342345781,
  0.816496580927726,
  0.942882415695480
};


Real integral_simpson(Real (*f)(Real, void*), Real a, Real b,int n,void* p);

int simpson_error;
const Real tolmax=1e-20;
const Real simpson_itermax=120;
#define INITIAL_STACK_SIZE 128   /* initial size of new stacks */




/* the stack structure */
struct stack_s{               
  int el_count;            /* count of elements on stack */    
  int el_size;             /* size of an element */
  int mem_reserve;         /* allocated memory for stack */
  void* elements;          /* pointer to begin of stack */
};

typedef struct _work_t{
  double a;
  double b;
  double tol;
  double S;
  double fa;
  double fb;
  double fm;
  double rec;
  int iter;
  struct my_f_params * p; //pointer auf params
} work_t;

typedef struct _work_t_gkq{
  double a;
  double b;
  double toler;
  double I_13;
  double I_prev;
  double fa;
  double fb;
  struct my_f_params * p; //pointer auf params
  shortint iter;
} work_gkq;


typedef struct stack_s* stack_t;
double integral_simpson_par(double (*f)(double, struct my_f_params*), stack_t stack);

double gkq_adapt_OMP(double (*f)(double, struct my_f_params*), stack_t stack);
double gkq_OMP(double (*f)(double, struct my_f_params*), double a, double b, double TOL, struct my_f_params* p,stack_t stack);

double gkq_serial(double (*f)(double, struct my_f_params*), double a, double b, double TOL, struct my_f_params* p);
double gkq_adapt_serial(double (*f)(double, struct my_f_params*), double a, double b, 
                        double fa,double fb, double toler,double I_13, struct my_f_params* p);

// void create_stack(stack_t* stack, int element_size);
// int  empty_stack(stack_t stack);
// void push_stack(stack_t stack, void* element);
// void pop_stack(stack_t stack, void* element);





/******************************************
 * create new stack
 ******************************************/
void create_stack(
             stack_t* stack,     /* stack to create */
             int element_size)   /* size of a stack element */
{
  int initial_size = INITIAL_STACK_SIZE;

  /* allocate memory for new stack struct */
  (*stack) = (stack_t) malloc(sizeof(struct stack_s));
  if (!(*stack)){
    char errstr[255];
    sprintf(errstr, "error: could not allocate memory for stack.. Abort.\n"); 
    error(errstr);
    // exit(1);
  }    

  /* allocate memory for stack elements */
  (*stack)->elements = (void*) malloc(element_size * initial_size);
  (*stack)->mem_reserve = initial_size; 
  if (!(*stack)->elements){
    char errstr[255];
    sprintf(errstr, "error: could not allocate memory for stack.. Abort.\n");
    error(errstr);
  }

  (*stack)->el_size = element_size;
  (*stack)->el_count = 0;

}

/*****************************************
 * check if the stack is empty 
 *****************************************/
int empty_stack(stack_t stack)
{
  return stack->el_count <= 0;
}


/*****************************************
 * push a element on stack
 *****************************************/
void push_stack(stack_t stack,    /* target stack */
           void* element)    /* element to push */
{
  int i, new_reserve;
  int log2_count;

  /* check if we need more memory for stack */    
  if (stack->el_count >= stack->mem_reserve)
  {

      /* calculate new size for the stack
         it should be a power of two */
      for (i = stack->el_count, log2_count = 0; 
           i > 0; 
           i>>1, log2_count++);
      new_reserve = 1 << log2_count;
 
      /* reallocate memory for phase thread tables 
         and nullify new values */
      stack->elements = (void *) realloc(stack->elements, 
                      stack->el_size * new_reserve);
      if (!stack->elements){
        char errstr [255];
        sprintf(errstr, "error: can't reallocate stack.. Aborting\n");
        error(errstr);
        // exit(1);
      }

      stack->mem_reserve = new_reserve;
  }
  
  /* now push the element on top of the stack */
  memcpy((char*)stack->elements + stack->el_count*stack->el_size, 
            element, stack->el_size);
  stack->el_count++;

}


/*****************************************
 * pop an element from stack
 *****************************************/
void pop_stack(
          stack_t stack,    /* target stack */
          void* element)    /* where poped el. should be stored */
{
  if (stack->el_count <= 0){
    char errstr[255];
    sprintf(errstr, "error: trying to pop from empty stack.\n");
    error(errstr);
    // exit(2);
  }

  stack->el_count--;
  memcpy(element, 
          (char*)stack->elements + stack->el_count*stack->el_size, 
          stack->el_size);

}



// ***************************************************************************
// *      Gauss-kronard quadrature, parallel
// ***************************************************************************
double gkq_OMP(double (*f)(double, struct my_f_params*), double a, double b, double TOL, struct my_f_params* p,stack_t stack)
{
  //1st integration
  double result=0.0;
// *********************************************                

  double m=0.5*(a+b);
  double h=0.5*(b-a);
  double y[13];
  double fa=y[0]=f(a,p);
  double fb=y[12]=f(b,p);
  int i;
  for(i=1;i<12;i++)
    y[i]=f(m+xgkq[i]*h,p);

  double I_4= (h/6.0)*(y[0]+y[12]+5.0*(y[4]+y[8]));  // 4-point gauss-lobatto
  double I_7= (h/1470.0)*(77.0*(y[0]+y[12])+432.0*(y[2]+y[10])+  // 7-point kronrod
               625.0*(y[4]+y[8])+672.0*y[6]);

  double I_13= h*(0.0158271919734802*(y[0]+y[12])+0.0942738402188500*(y[1]+y[11])+0.155071987336585*(y[2]+y[10])+
                  0.188821573960182*(y[3]+y[9])+0.199773405226859*(y[4]+y[8])+0.224926465333340*(y[5]+y[7])+
                  0.242611071901408*y[6]);   //13-point Kronrod

  
  double Err1=fabs(I_7-I_13);
  double Err2=fabs(I_4-I_13);

  double r=(Err2 != 0.0) ? Err1/Err2 : 1.0;
  double toler=(r > 0.0 && r < 1.0) ? TOL/r : TOL; 

  if(I_13 == 0)
    I_13=b-a;
  I_13=fabs(I_13);

  
  //Prepare work and push onto stack
  work_gkq work;

  work.a = a;
  work.b = b;
  work.toler = toler;
  work.I_13=I_13;
  work.fa=fa;
  work.fb=fb;
  work.p=p;
  work.I_prev=I_7;

  //ANTI-FOLGENDES:
  //OUT OF TOLERANCE !!!, mll:3.0162e-18, a:3.0162e-18, b:3.0162e-18, mrr:3.0162e-18,I_7-I_4:0.0000e+00, tol:1.6002e-315,I_13:7.0585e-313
  if(I_13 < 1e-150) 
    return 0;

  push_stack(stack, &work); 
  result=gkq_adapt(f,stack);
  return result;
}

double gkq_serial(double (*f)(double, struct my_f_params*), double a, double b, double TOL, struct my_f_params* p)
{
  //1st integration

  double result=0.0;
  gkq_iter_serial=0;
// *********************************************                

  double m=0.5*(a+b);
  double h=0.5*(b-a);
  double y[13];
  double fa=y[0]=f(a,p);
  double fb=y[12]=f(b,p);
  int i;
  for(i=1;i<12;i++)
    y[i]=f(m+xgkq[i]*h,p);

  double I_4= (h/6.0)*(y[0]+y[12]+5.0*(y[4]+y[8]));  // 4-point gauss-lobatto
  double I_7= (h/1470.0)*(77.0*(y[0]+y[12])+432.0*(y[2]+y[10])+  // 7-point kronrod
               625.0*(y[4]+y[8])+672.0*y[6]);

  double I_13= h*(0.0158271919734802*(y[0]+y[12])+0.0942738402188500*(y[1]+y[11])+0.155071987336585*(y[2]+y[10])+
                  0.188821573960182*(y[3]+y[9])+0.199773405226859*(y[4]+y[8])+0.224926465333340*(y[5]+y[7])+
                  0.242611071901408*y[6]);   //13-point Kronrod

  
  double Err1=fabs(I_7-I_13);
  double Err2=fabs(I_4-I_13);

  double r=(Err2 != 0.0) ? Err1/Err2 : 1.0;
  double toler=(r > 0.0 && r < 1.0) ? TOL/r : TOL; 

  if(I_13 == 0)
    I_13=b-a;
  I_13=fabs(I_13);


  result=gkq_adapt_serial(f,a,b,fa,fb,toler,I_13, p);  

  return result;
}


// ***********************************************
// * RECURSIVE ADAPTION ROUTINE FOR PARALLEL-GK-QUADRATURE 
// **********************************************
double gkq_adapt_OMP(double (*f)(double, struct my_f_params*), stack_t stack)
{
  work_gkq work;
  work.iter=0;
  int ready, idle, busy;
  double integral_result = 0.0;  

  busy = 0;
  terminate_gkq=0;

#pragma omp parallel default(none) \
    shared(stack, integral_result,f,busy,terminate_gkq,myid) \
    private(work, idle, ready)
  {  
// printf("me:%d, err:%d\n",omp_get_thread_num(),simpson_error);    

    ready = 0;
    idle = 1;

    while(!ready) // && !terminate_gkq)//  && !simpson_error) //<-- so NICHT!
    {
      #pragma omp critical (stack)
      {
        if (!empty_stack(stack))
        {
          /* we have new work */ 
          pop_stack(stack, &work);
          if (idle)
          {
            /* say others i'm busy */
            busy += 1;
            idle = 0;
          }
        }
        else
        {
          /* no new work on stack */
          if (!idle){
            busy -= 1;
            idle = 1;
          }

          /* nobody has anything to do; let us leave the loop */
          if (busy == 0)
          {
            ready = 1;        
          }
        }
      } /* end critical(stack) */

      if (idle)
        continue; //if ready==1 --> leave loop
double I_prev=work.I_prev;

      double a = work.a;
      double b = work.b;      
      double toler = work.toler;    
      double I_13=work.I_13; 
      double fa=work.fa;
      double fb=work.fb;
      int iter=work.iter;
      // double *y= work.y; // brauch ich nicht!
      struct my_f_params * p = work.p;
      
      double m = (a+b)/2;
      double h = (b -a)/2;
      double mll=m-gkq_alpha*h;
      double ml=m-gkq_beta*h;
      double mr=m+gkq_beta*h;
      double mrr=m+gkq_alpha*h;

      double fmll=f(mll,p);
      double fml=f(ml,p);
      double fm=f(m,p);
      double fmr=f(mr,p);
      double fmrr=f(mrr,p);
      double I_4=h/6.0*(fa+fb+5.0*(fml+fmr));   // 4-point Gauss-Lobatto formula.
      double I_7=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
      
// if(myid==1)
//   printf("I_7:%.4e, I_13:%.4e,I_4:%.4e, minus:%.4e, to:%.4e\n",I_7,I_13,I_4,I_7-I_4, toler*I_13);
int maxiter=50; //max. subdivisions
double abstol=1e-30;
work.I_prev=I_7; // für abstolcheck in nächster recursion

      if (fabs(I_7-I_4) <= toler*I_13 || mll <= a || b <= mrr || iter > maxiter || fabs(I_7-I_prev) < abstol ) 
      {
        if ((mll <= a || b <= mrr)) //Error
         {
           // out_of_tolerance=true; // Interval contains no more machine numbers
           // printf("OUT OF TOLERANCE !!!, mll:%.4e, a:%.4e, b:%.4e, mrr:%.4e,I_7-I_4:%.4e, tol:%.4e,I_13:%.4e\n", 
           //        mll,b,b,mrr,I_7-I_4, toler*I_13,I_13);
           terminate_gkq=1;  
         }
        #pragma omp critical (integral_result)            
        {
          integral_result += I_7;      //Terminate recursion.  
        }        
// printf("me ok:%d, a:%f,b:%f, tler:%.5e,I_4:%f,I_7:%f,ubteg;%.4e\n", omp_get_thread_num(), a,b,toler,I_4,I_7,integral_result);        

      }
      else  //subdivide interval and push new work on stack
      {
        #pragma omp critical (stack)
        {

          // printf("me NOOOO:%d, a:%f,b:%f, tler:%.5e,I_4:%f,I_7:%f\n", omp_get_thread_num(), a,b,toler,I_4,I_7);

          work.iter=iter+1;

          work.a=a;
          work.b=mll;
          work.fa=fa;
          work.fb=fmll;
          push_stack(stack, &work);   

          work.a=mll;
          work.b=ml;
          work.fa=fmll;
          work.fb=fml;
          push_stack(stack, &work);   

          work.a=ml;
          work.b=m;
          work.fa=fml;
          work.fb=fm;
          push_stack(stack, &work);             

          work.a=m;
          work.b=mr;
          work.fa=fm;
          work.fb=fmr;
          push_stack(stack, &work);             

          work.a=mr;
          work.b=mrr;
          work.fa=fmr;
          work.fb=fmrr;
          push_stack(stack, &work);             

          work.a=mrr;
          work.b=b;
          work.fa=fmrr;
          work.fb=fb;
          push_stack(stack, &work);                       

        } // pragma critical stack
      }   // else ..non-acceptable error
    } // while
  } /* end omp parallel */
  return integral_result;    
}


double gkq_adapt_serial(double (*f)(double, struct my_f_params*), double a, double b, double fa,
                        double fb, double toler,double I_13, struct my_f_params* p)
{
  double m = (a+b)/2;
  double h = (b -a)/2;
  double mll=m-gkq_alpha*h;
  double ml=m-gkq_beta*h;
  double mr=m+gkq_beta*h;
  double mrr=m+gkq_alpha*h;

  double fmll=f(mll,p);
  double fml=f(ml,p);
  double fm=f(m,p);
  double fmr=f(mr,p);
  double fmrr=f(mrr,p);
  double I_4=h/6.0*(fa+fb+5.0*(fml+fmr));   // 4-point Gauss-Lobatto formula.
  double I_7=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
  
  gkq_iter_serial++;


  if ( (fabs(I_7-I_4) <= toler*I_13 || mll <= a || b <= mrr) && gkq_iter_serial) 
  {
    if ((mll <= a || b <= mrr) && !terminate_serial) //Error
     {
      // out_of_tolerance=true; // Interval contains no more machine numbers
      printf("OUT OF TOLERANCE !!!, mll:%.4e, a:%.4e, b:%.4e, mrr:%.4e\n", mll,b,b,mrr);
      terminate_serial=1;  
     }

// printf("me ok:%d, a:%f,b:%f, tler:%.5e,I_4:%f,I_7:%f\n", omp_get_thread_num(), a,b,toler,I_4,I_7);      
    return I_7;
  }
  else
  {

    // printf("me NOOOO:%d, a:%f,b:%f, tler:%.5e,I_4:%f,I_7:%f\n", omp_get_thread_num(), a,b,toler,I_4,I_7);

    return  gkq_adapt_serial(f, a,mll,fa,fmll,toler,I_13,p) +
            gkq_adapt_serial(f, mll,ml,fmll,fml,toler,I_13,p) +
            gkq_adapt_serial(f, ml,m,fml,fm,toler,I_13,p) +
            gkq_adapt_serial(f, m,mr,fm,fmr,toler,I_13,p) +
            gkq_adapt_serial(f, mr,mrr,fmr,fmrr,toler,I_13,p) +
            gkq_adapt_serial(f, mrr,b,fmrr,fb,toler,I_13,p);
  }


}