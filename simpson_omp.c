#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <sys/time.h>

// ***************************************************************************************
// * Pre-Port Entwicklung und Testen eines OMP-parallisierten, adaptiven Simpson-integrators für 
// * die fiesen doppelintgrale in imd_colrad.c
// ***************************************************************************************
int num_threads;
int simpson_error;
int funcevals;
int iter;

const double tolmax=1e-30;
#define INITIAL_STACK_SIZE 128; //128   /* initial size of new stacks */


const double alpha=0.816496580927726;
const double beta=0.447213595499958;          
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






/* the stack structure */
struct stack_s{               
  int el_count;            /* count of elements on stack */    
  int el_size;             /* size of an element */
  int mem_reserve;         /* allocated memory for stack */
  void* elements;          /* pointer to begin of stack */
};


struct my_f_params { double ne; double T;double mu; double E;double DeltaE; int allowed;};



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


typedef struct stack_s* stack_t;

typedef struct _work_t_gkq{
  double a;
  double b;
  double toler;
  double I_13;
  double fa;
  double fb;
  struct my_f_params * p; //pointer auf params
  double (*f)(double, struct my_f_params*);
  stack_t stack_inner;
  short int is_parent;
  short int subtask_nr; //inner task nr
  short int task_nr; //outer task nr
  short int subtasks_left; //erst pop'n wenn aller inner tasks (intergals) vollständig!
  double inner_integrals[5];
} work_gkq;





// **************************************** FUNCTION DEFS *************************************************
void create_stack(stack_t* stack, int element_size);
int empty_stack(stack_t stack);
void push_stack(stack_t stack, void* element);
void pop_stack(stack_t stack, void* element);


double  
integral2(     double (*f)(double, struct my_f_params*), /* function to integrate */
     stack_t stack);


double gkq_adapt_double(stack_t stack);
double gkq_adapt_single(stack_t stack); //for single integral, first integral

double gkq_double(double (*f)(double, struct my_f_params*), double (*finner)(double, struct my_f_params*), 
                  double a, double b, 
                  double TOL, struct my_f_params* p,stack_t stack);

double gkq_adapt_serial(double (*f)(double, struct my_f_params*), double a, double b, double fa,double fb, double toler,double I_13, struct my_f_params* p);

work_gkq gkq_create_inner_task(double (*f)(double, struct my_f_params*), double a, double b, double TOL, struct my_f_params* p);
double gkq_single(double (*f)(double, struct my_f_params*), double a, double b, 
                       double TOL, struct my_f_params* p,stack_t stack);

int terminate_serial;
int terminate_gkq;

// static double myfun(double x,struct my_f_params* p)
// {  
//   double T=p->T;  
//   double fermi=1/(1+exp(-x/T));
//   double fermi2=1- 1/(1+exp(-x/T));
//   double sigma=x/T*log(x/T*2)/pow(1/x,2.0); //einfach nur ärger machen

//   // return  exp(-x*x)/p->T+log(x*pow(p->T,2.0))/p->ne;
//   return fermi*fermi2*sigma;
// }

// static double myfun2(double x,void* pv)
// { 
//   struct my_f_params *p= (struct my_f_params*) pv;  

//   double T=p->T;  
//   double fermi=1/(1+exp(-x/T));
//   double fermi2=1- 1/(1+exp(-x/T));
//   double sigma=x/T*log(x/T*2)/pow(1/x,2.0); //einfach nur ärger machen

  
//   // return  exp(-x*x)/p->T+log(x*pow(p->T,2.0))/p->ne;
//   return fermi*fermi2*sigma;
// }


static double inner_integrand(double x,void *pv)
{
  struct my_f_params *p= (struct my_f_params*) pv;  
  double T=p->T;  
  double fermi=1/(1+exp(-x/T));
  double fermi2=1- 1/(1+exp(-x/T));    
  return exp(-x*x); //fermi*fermi2;
}
static double inner_integrand2(double x,struct my_f_params* p)
{
  double T=p->T;  
  double fermi=1/(1+exp(-x/T));
  double fermi2=1- 1/(1+exp(-x/T));
  return fermi*fermi2;
}
void get_integ_bounds_inner(double *integ_bnd, double x, struct my_f_params *p)
{
  integ_bnd[0]=-4;
  integ_bnd[1]=4;
}

gsl_integration_workspace * gsinner;
gsl_integration_workspace * gsinner2;

static double myfun(double x,struct my_f_params* p)
{  
  double T=p->T;  
  double fermi=1/(1+exp(-x/T));
  double fermi2=1- 1/(1+exp(-x/T));
  double sigma=x/T*log(x/T*2)/pow(1/x,2.0); //einfach nur ärger machen

  // double result, err;
  // gsl_function gslfun;
  // gslfun.function=&inner_integrand;
  // gslfun.params=p;    

  // gsl_integration_qag(&gslfun, x, 300, 0.0, 1e-3, 1000,1,
  //                      gsinner2, &result, &err);  

  // gsl_integration_workspace_free(gsinner2);
  // stack_t stack2;  
  // create_stack(&stack2, sizeof(work_gkq));

  //double result=gkq(inner_integrand2, x, 600, 1e-3, p);//, stack);
  // free(stack2->elements);


  return sigma;
  //return  exp(-x*x);
}; 


static double myfun_gsl(double x,void* pv)
{ 
  struct my_f_params *p= (struct my_f_params*) pv;  

  double T=p->T;  
  double fermi=1/(1+exp(-x/T));
  double fermi2=1- 1/(1+exp(-x/T));
  double sigma=x/T*log(x/T*2)/pow(1/x,2.0); //einfach nur ärger machen


  gsl_function gslfun;
  gslfun.function=&inner_integrand;
  gslfun.params=pv;    
  double result, err;
  gsl_integration_qag(&gslfun, x, 600, 0.0, 1e-3, 1000,1,
                       gsinner, &result, &err);  

  return sigma*result;
}


// static double myfun2(double x,void* pv)
// { 
//   struct my_f_params *p= (struct my_f_params*) pv;  
// static float sinfc(float x) { return sinf(x); } 
// static double frand48c(double x) { funcevals++; return (double) drand48(); } 

// *************************************************************************************************************
int main(int argc,char** argv)
{
  iter=0;
  num_threads=omp_get_num_threads();
  funcevals=0;
  simpson_error=0;

  struct timeval start, end;
    
  double gsltime, simpsontime;
  // ************************************************
  gsinner=   gsl_integration_workspace_alloc (1000);   
  
  

  // double gslinteg=0.0;
  // double integ_err;

  // gsl_function gslfun;
  // gslfun.function=&myfun;


  // **********************************************

  double pi=3.141592653589793;
  double xmin, xmax;
  double answer = 0.0;

  int i;

  xmin=2;
  xmax=100; 
  
  // ***********************************************
  double T0=1;
  double ne0=1e5;

  struct my_f_params ptest;
  ptest.T=T0;
  ptest.ne=ne0;

  // *************************************************
  int method=GSL_INTEG_GAUSS61;
  gsl_function gslfun;
  gslfun.function=&myfun_gsl;
  gslfun.params=&ptest;  
  gsl_integration_workspace * gswork=NULL;
  gswork=   gsl_integration_workspace_alloc (1000); 
  double gslinteg=0.0;
  double integ_err;

  gettimeofday(&start, NULL); 
  gsl_integration_qag(&gslfun, xmin, xmax, 0.0, 1e-8, 1000,1,
                       gswork, &gslinteg, &integ_err);  
  gettimeofday(&end, NULL);
  gsltime = ((end.tv_sec  - start.tv_sec) * 1000000u + 
              end.tv_usec - start.tv_usec) /1.e6;

  printf("gslres:%.12e, gsltime:%.4e\n", gslinteg, gsltime);
  // *************************************************

  stack_t stack;  //global stack
  create_stack(&stack, sizeof(work_gkq));

  gettimeofday(&start, NULL);     
  double integ=gkq_double(myfun,inner_integrand2, xmin, xmax, 1e-2, &ptest,stack);
  
  gettimeofday(&end, NULL);
  free(stack->elements);
  free(stack);

  double partime = ((end.tv_sec  - start.tv_sec) * 1000000u + 
                     end.tv_usec - start.tv_usec) /1.e6;

  printf("my:%.12e, partime:%.4e\n", integ,partime);

  return 0;
}



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
    fprintf(stderr, "error: could not allocate memory for stack.. Abort.\n"); 
    exit(1);
  }    

  /* allocate memory for stack elements */
  (*stack)->elements = (void*) malloc(element_size * initial_size);
  (*stack)->mem_reserve = initial_size; 
  if (!(*stack)->elements){
    fprintf(stderr, "error: could not allocate memory for stack.. Abort.\n");
    exit(1);
  }

  (*stack)->el_size = element_size;
  (*stack)->el_count = 0;

}

/*****************************************
 * check if the stack is empty 
 *****************************************/
int empty_stack(stack_t stack)
{
  //MYMOD
  if(stack==NULL)
    return 1;
  //ENDOF MYMOD
  return stack->el_count <= 0;
}


/*****************************************
 * push a element on stack
 *****************************************/
void push_stack(
           stack_t stack,    /* target stack */
           void* element)    /* element to push */
{
  int i, new_reserve;
  int log2_count;

  /* check if we need more memory for stack */      
  if (stack->el_count >= stack->mem_reserve)
  {

      /* calculate new size for the stack
         it should be a power of two */
      // for (i = stack->el_count, log2_count = 0; i > 0;  i>>1, log2_count++);
      //   new_reserve = 1 << log2_count;

      log2_count=0;
      for(i=stack->el_count;i>0; i=i*0.5) //i>>1)
      {
        log2_count++;
      }
      new_reserve= 1 << log2_count;
 
      /* reallocate memory for phase thread tables 
         and nullify new values */
      stack->elements = (void *) realloc(stack->elements, 
                      stack->el_size * new_reserve);
      if (!stack->elements){
        fprintf(stderr, "error: can't reallocate stack.. Aborting\n");
        exit(1);
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
    fprintf(stderr, "error: trying to pop from empty stack.\n");
    exit(2);
  }

  stack->el_count--;
  memcpy(element, 
         (char*)stack->elements + stack->el_count*stack->el_size, 
         stack->el_size);

  
}





// ***************************************************************************
// *      Gauss-kronard quadrature
// ***************************************************************************
work_gkq gkq_create_inner_task(double (*f)(double, struct my_f_params*), double a, double b, double TOL, struct my_f_params* p)
{
  work_gkq work;

  work.f=f;
  struct my_f_params* pinner;
  pinner=(struct my_f_params*) malloc(sizeof(struct my_f_params)); //ACHTUNG: Muss wieder gefreed werden
  pinner->mu= p->mu;
  pinner->T=p->T;
  pinner->ne=p->ne;
  work.p=pinner;
  work.a=a;
  work.b=b;
  work.subtasks_left=0;//beliebig!

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

// printf("create inner task: a:%f, b:%f, I4:%f,I7:%f,I13:%f\n",a,b,I_4,I_7,I_13);
  work.toler = toler;
  work.I_13=I_13;
  work.fa=fa;
  work.fb=fb;
  
  work.is_parent=0;
  for(i=0;i<5;i++)
    work.inner_integrals[i]=0.0;
  return work;

}



double gkq_double(double (*f)(double, struct my_f_params*), double (*finner)(double, struct my_f_params*), 
                  double a, double b, 
                  double TOL, struct my_f_params* p,stack_t stack)
{
  //1st integration
  //create parent task

  double result=0.0;
// *********************************************                

  double m=0.5*(a+b);
  double h=0.5*(b-a);
  int i;


  double y[13];

  stack_t st;
  create_stack(&st, sizeof(work_gkq));

  double integ_bnd[2];

  get_integ_bounds_inner(integ_bnd, a, p);
  y[0]=f(a,p)* gkq_single(finner,integ_bnd[0],integ_bnd[1],1e-3, p, st);

  for(i=1;i<12;i++)
  {
    get_integ_bounds_inner(integ_bnd, m+xgkq[i]*h, p);
    y[i]=f(m+xgkq[i]*h,p)*gkq_single(finner,integ_bnd[0],integ_bnd[1],1e-3,p,st);
  }

  get_integ_bounds_inner(integ_bnd, b, p);
  y[12]=f(b,p)*gkq_single(finner,integ_bnd[0],integ_bnd[1],1e-3, p, st);


  free(st->elements);
  free(st);


  
  double fa=y[0];
  double fb=y[12];


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

  printf("First approx I_13:%.4e\n",I_13); 
  
  //Prepare work and push onto stack
  work_gkq work;
  work.a = a;
  work.b = b;
  work.toler = toler;
  work.I_13=I_13;
  work.fa=fa;
  work.fb=fb;

  work.p=p;
  work.f=f;
  work.is_parent=1;
  work.task_nr=0; //fange bei 0 das zählen an
  

  for(i=0;i<5;i++)
    work.inner_integrals[i]=0.0;
  
  // ALLOC INNER WS FOR INTEGR.
  // gsl_integration_workspace *gsinner2=  gsl_integration_workspace_alloc (1000); 
  
  stack_t stack_inner; 
  create_stack(&stack_inner, sizeof(work_gkq));
  work.stack_inner=stack_inner;

  //push work on inner stack
  //for comp of 5 inner intgrals
  work_gkq winner;
  double m_inner, h_inner, mll_inner,ml_inner, mr_inner, mrr_inner;

  m_inner=  (work.a+work.b)/2;
  h_inner=  (work.a+work.b)/2;
  mll_inner=(m_inner-alpha*h_inner);
  ml_inner= (m_inner-beta*h_inner);
  mr_inner= (m_inner+beta*h_inner);
  mrr_inner=(m_inner+alpha*h_inner);            


  //1st subtask
  get_integ_bounds_inner(integ_bnd, mll_inner, work.p);            
  winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0], integ_bnd[1],1e-3, work.p);
  winner.subtask_nr=0;  
  winner.task_nr=work.task_nr;

  push_stack(work.stack_inner,&winner);

  //2nd subtask
  get_integ_bounds_inner(integ_bnd, ml_inner, work.p);            
  winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0], integ_bnd[1],1e-3, work.p);
  winner.subtask_nr=1;
  winner.task_nr=work.task_nr;
  push_stack(work.stack_inner,&winner); 

  //3rd subtask
  get_integ_bounds_inner(integ_bnd, m_inner, work.p);            
  winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0], integ_bnd[1],1e-3, work.p);
  winner.subtask_nr=2;
  winner.task_nr=work.task_nr;
  push_stack(work.stack_inner,&winner);     

  //4th subtask
  get_integ_bounds_inner(integ_bnd, mr_inner, work.p);            
  winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0], integ_bnd[1],1e-3, work.p);
  winner.subtask_nr=3;
  winner.task_nr=work.task_nr;
  push_stack(work.stack_inner,&winner);     

  //5th subtask
  get_integ_bounds_inner(integ_bnd, mrr_inner, work.p);            
  winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0], integ_bnd[1],1e-3, work.p);
  winner.subtask_nr=4;
  winner.task_nr=work.task_nr;
  push_stack(work.stack_inner,&winner);

  //push outer integtand work on outer stack
  work.subtasks_left=5;
  push_stack(stack, &work);                   


/*
  int elemsinner=work.stack_inner->el_count;
  int elout=stack->el_count;


  work_gkq* top=(work_gkq*) stack->elements +(stack->el_count-1); //get top element
  
  int inner2=0;
  inner2=top->stack_inner->el_count;
*/
  // printf("elemsinner:%d,ou:%d,in2:%d\n",elemsinner,elout, inner2);
  

  // int in2=  top->stack_inner->el_count;


  result=gkq_adapt_double(stack);//,stack);



  // result=gkq_adapt_serial(f,a,b,fa,fb,toler,I_13, p);  
  // gsl_integration_workspace_free(gsinner2);
  

  // free(stack->elements);

  return result;
}



double gkq_single(double (*f)(double, struct my_f_params*), double a, double b, 
                       double TOL, struct my_f_params* p,stack_t stack)
{
  double result=0.0;
// *********************************************                

  double m=0.5*(a+b);
  double h=0.5*(b-a);
  int i;


  double y[13];
  double fa=y[0]=f(a,p);
  double fb=y[12]=f(b,p);
  
  for(i=1;i<12;i++)           //hier müssen direkt 13 integrale berechnet werden!!!
    y[i]=f(m+xgkq[i]*h,p);    //auch das sollte parallel erfolgen (pragma parfor?)


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
  work.a=a;
  work.f=f;


  push_stack(stack, &work); 
  result=gkq_adapt_single(stack);//,stack);

  return result;
}

double gkq_adapt_serial(double (*f)(double, struct my_f_params*), double a, double b, 
                        double fa,double fb, double toler,double I_13, struct my_f_params* p)
{
  double m = (a+b)/2;
  double h = (b-a)/2;
  double mll=m-alpha*h;
  double ml=m-beta*h;
  double mr=m+beta*h;
  double mrr=m+alpha*h;

  double fmll=f(mll,p);
  double fml=f(ml,p);
  double fm=f(m,p);
  double fmr=f(mr,p);
  double fmrr=f(mrr,p);
  double I_4=h/6.0*(fa+fb+5.0*(fml+fmr));   // 4-point Gauss-Lobatto formula.
  double I_7=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
  

  if (fabs(I_7-I_4) <= toler*I_13 || mll <= a || b <= mrr) 
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


double gkq_adapt_single(stack_t stack)
{
  work_gkq work;
  // work.iter=0;
  int ready, idle, busy;
  double integral_result = 0.0;  

  busy = 0;
  terminate_gkq=0;

#pragma omp parallel default(none) \
    shared(stack, integral_result,busy) \
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
// double I_prev=work.I_prev;

double (*f)(double, struct my_f_params*)=work.f; 

      double a = work.a;
      double b = work.b;      
      double toler = work.toler;    
      double I_13=work.I_13; 
      double fa=work.fa;
      double fb=work.fb;
      // int iter=work.iter;
      // double *y= work.y; // brauch ich nicht!
      struct my_f_params * p = work.p;
      
      double m = (a+b)/2;
      double h = (b -a)/2;
      double mll=m-alpha*h;
      double ml=m-beta*h;
      double mr=m+beta*h;
      double mrr=m+alpha*h;

      double fmll=f(mll,p);
      double fml=f(ml,p);
      double fm=f(m,p);
      double fmr=f(mr,p);
      double fmrr=f(mrr,p);
      double I_4=h/6.0*(fa+fb+5.0*(fml+fmr));   // 4-point Gauss-Lobatto formula.
      double I_7=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
      
// if(myid==1)
//   printf("I_7:%.4e, I_13:%.4e,I_4:%.4e, minus:%.4e, to:%.4e\n",I_7,I_13,I_4,I_7-I_4, toler*I_13);
// int maxiter=50; //max. subdivisions
// double abstol=1e-30;
// work.I_prev=I_7; // für abstolcheck in nächster recursion

      if (fabs(I_7-I_4) <= toler*I_13 || mll <= a || b <= mrr) // || iter > maxiter || fabs(I_7-I_prev) < abstol ) 
      {
        if ((mll <= a || b <= mrr)) //Error
         {
           // out_of_tolerance=true; // Interval contains no more machine numbers
           // printf("OUT OF TOLERANCE !!!, mll:%.4e, a:%.4e, b:%.4e, mrr:%.4e,I_7-I_4:%.4e, tol:%.4e,I_13:%.4e\n", 
           //        mll,b,b,mrr,I_7-I_4, toler*I_13,I_13);

         }
        #pragma omp critical (integral_result)            
        {
          integral_result += I_7;      //Terminate recursion.  
        }        
      }
      else  //subdivide interval and push new work on stack
      {
        #pragma omp critical (stack)
        {          
          // work.iter=iter+1;

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

double gkq_adapt_double(stack_t stack)
{
  work_gkq work;
  work_gkq* pwork_outer;

  int ready, idle, busy;
  double integral_result = 0.0;    
  busy = 0;
  int myid;
  int elcnt; //nr of outer task elements on stack

#pragma omp parallel default(none)  \
    shared(stack, integral_result,busy,iter,elcnt) \
    private(work, idle, ready,myid,pwork_outer)
  {      

    myid=omp_get_thread_num();

    ready = 0;
    idle = 1;

    while(!ready) // && !terminate_gkq)//  && !simpson_error) //<-- so NICHT!
    {


// printf("\n NEXT, curapprox:%.4e\n",integral_result);
      #pragma omp critical (stack)
      {
        //pointer to outer work element
        //work_gkq* pwork=(work_gkq*) stack->elements + stack->el_count*stack->el_size;

        elcnt=stack->el_count;
     

// if(elcnt>1)
//  getchar();
        
        stack_t stack_inner;
        if(elcnt>0)
        {
          pwork_outer=(work_gkq*) stack->elements +(elcnt-1); //get top element

          //IST WOHL SUBOPTIMAL
            while(pwork_outer->task_nr > 0 && pwork_outer->subtasks_left==0) //this task is complete, goto next outer task
            {
              printf("myid:%d, outer task nr:%d, complete:subs:%d, decrement\n",myid,pwork_outer->task_nr,pwork_outer->subtasks_left);
              pwork_outer--;
            }
            stack_inner=pwork_outer->stack_inner;
        }
        else
        {
          stack_inner=NULL;
          pwork_outer=NULL;
        }

        printf("myid:%d,elcnt:%d,sinner:%p,pwout:%p\n",myid,elcnt,stack_inner,pwork_outer);   

        { //stackinner gehört zu work, und work ist private!          
          if(!empty_stack(stack_inner))
          {           
            printf("myid:%d,iter:%d, inner stack not empty,pop..\n",myid,iter);              
            pop_stack(stack_inner, &work);
//if letztes inner-work elem, blockiere outer work elem bis ich das integral hab

            iter++;

            if (idle)
            {          
              busy += 1;
              idle = 0;
            }              
          }
          else //inner stack is empty, pop from outer
          {
            printf("myid:%d,iter:%d, inner stack is empty?\n",myid,iter);
            if (!empty_stack(stack)) //work elem ist blockiert solange inner integs nicht vollständig
            {              
              // printf("myid:%d,iter:%d, inner stack empty. work on outer with cnt=%d\n",myid,iter,pwork_outer->subtasks_left);
              if(pwork_outer->subtasks_left==0)
              {
                printf("myid:%d,iter:%d, outer stack not empty\n",myid,iter);
                pop_stack(stack, &work);
                iter++;
                if (idle)
                {            
                  busy += 1;
                  idle = 0;
                }              
              }
              else
              {
                printf("myid:%d,iter:%d, inner integs not complete:%d\n",myid,iter,pwork_outer->subtasks_left);
                if (!idle)
                {
                  busy -= 1;
                  idle = 1;
                }


              }
            }
            else //auch outer stack ist leer
            {           
              printf("myid:%d,iter:%d, outer stack empty too\n",myid,iter);   
              if (!idle)
              {
                busy -= 1;
                idle = 1;
              }
              if (busy == 0)
              {
                ready = 1;        
              }
            }
          }            
        } // critical inner stack
      } //critical outer stack


      if (idle)
      {
        printf("myid:%d,iter:%d, noth8ing to do\n",myid,iter);
        continue; //if ready==1 --> leave loop
      }


      //work on inner tasks first,  if available
      if(work.is_parent==0)
      {
        printf("myid:%d,iter:%d,tasknr:%d,innertask:%d, left subs:%d\n",
                myid,iter,work.task_nr,work.subtask_nr,pwork_outer->subtasks_left); //subtasks nr nicht immer aktuell!

// getchar();

        double (*f)(double, struct my_f_params*) = work.f;
        double a = work.a;
        double b = work.b;      
        double toler = work.toler;    
        double I_13=work.I_13; 
        double fa=work.fa;
        double fb=work.fb;
        // double *y= work.y; // brauch ich nicht!
        struct my_f_params * p = work.p;
        
        double m = (a+b)/2;
        double h = (b -a)/2;
        double mll=m-alpha*h;
        double ml=m-beta*h;
        double mr=m+beta*h;
        double mrr=m+alpha*h;

        double fmll=f(mll,p);
        double fml=f(ml,p);
        double fm=f(m,p);
        double fmr=f(mr,p);
        double fmrr=f(mrr,p);
        double I_4=h/6.0*(fa+fb+5.0*(fml+fmr));   // 4-point Gauss-Lobatto formula.
        double I_7=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
      
// printf("myid:%d,iter:%d non parent checkpoint 1 passed\n",myid,iter);

        if (fabs(I_7-I_4) <= toler*I_13 || mll <= a || b <= mrr) 
        {
          if ((mll <= a || b <= mrr))// && !terminate_gkq) //Error
           {             
             printf("OUT OF TOLERANCE !!!, mll:%.4e, a:%.4e, b:%.4e, mrr:%.4e\n", mll,a,b,mrr);
           }
            
// printf("myid:%d,iter:%d non parent checkpoint 2.1 passed\n",myid,iter);

          int tasknr=work.subtask_nr;

          printf("myid:%d,iter:%d, winner error acceptable: I_7:%.4e, I_4:%.4e, I_13:%.4e,toler*I_13:%.4e,toler:%.4e\n",
                  myid,iter,I_7,I_4, I_13,toler*I_13,toler);


          #pragma omp critical (innerinteg) // workers greifen auf selbes array zu.nur an anderer Stelle..evtl. besser 6 versch.doubles statt array?
          {
            double *inner_integrals=pwork_outer->inner_integrals;
            inner_integrals[tasknr]+=I_7;
            pwork_outer->subtasks_left=pwork_outer->subtasks_left-1; //hier sollte subtasks_left aktuell sein

            printf("myid:%d,iter:%d, inner_integral[%d]:%.4e added:%.4e,subtask left:%d\n",
                    myid,iter,tasknr,inner_integrals[tasknr],I_7,pwork_outer->subtasks_left);
          }            
        }
        else  //subdivide interval and push new work on stack
        {
          
          stack_t stack_inner;         
          #pragma omp critical //(stack) //(stack_inner)
          { 
            stack_inner =  pwork_outer->stack_inner;
            pwork_outer->subtasks_left=pwork_outer->subtasks_left-1;//weil dieser stack durch 6 andere ersetzt wird


            work.task_nr=pwork_outer->task_nr; //inner task soll auch wissen wer sein outer task is (bisher nur für debug)
            //subtasknr bleibt dieselbe (bezogen auf parent task)
            work.a=a;
            work.b=mll;
            work.fa=fa;
            work.fb=fmll;
            
            push_stack(stack_inner, &work);   

            work.a=mll;
            work.b=ml;
            work.fa=fmll;
            work.fb=fml;
            
            push_stack(stack_inner, &work);   

            work.a=ml;
            work.b=m;
            work.fa=fml;
            work.fb=fm;
            
            push_stack(stack_inner, &work);             

            work.a=m;
            work.b=mr;
            work.fa=fm;
            work.fb=fmr;
            
            push_stack(stack_inner, &work);             

            work.a=mr;
            work.b=mrr;
            work.fa=fmr;
            work.fb=fmrr;
            
            push_stack(stack_inner, &work);             

            work.a=mrr;
            work.b=b;
            work.fa=fmrr;
            work.fb=fb;
            
            push_stack(stack_inner, &work);                       

            pwork_outer->subtasks_left=pwork_outer->subtasks_left+6;

            printf("myid:%d,iter:%d, inner_integral approx not accurate..subdivided: I_4:%f,I_7:%f,I_13:%f,left subs:%d\n",
                    myid,iter,I_4,I_7,I_13,pwork_outer->subtasks_left);

          } // pragma critical stack_inner
        }   // else ..non-acceptable error
      }// !isparent
      else  //parent task without any inner tasks
      {

printf("myid:%d,iter:%d, work is parent,subs left:%d\n",myid,iter,work.subtasks_left);

        // stack_t stack_inner = work.stack_inner; //brauch ich hier nicht
        // free(stack_inner->elements);
        // free(stack_inner);

        double (*f)(double, struct my_f_params*) = work.f;
        double a = work.a;
        double b = work.b;      
        double toler = work.toler;    
        double I_13=work.I_13; 
        double fa=work.fa;
        double fb=work.fb;
        // double *y= work.y; // brauch ich nicht!
        struct my_f_params * p = work.p;
        
        double m = (a+b)/2;
        double h = (b -a)/2;
        double mll=m-alpha*h;
        double ml=m-beta*h;
        double mr=m+beta*h;
        double mrr=m+alpha*h;

        //the inner intergrals are pre-calculated and saved in an array
        //the function only calculates a prefactor...
        double fmll=f(mll,p)*work.inner_integrals[0];
        double fml=  f(ml,p)*work.inner_integrals[1];
        double fm=   f(m,p)*work.inner_integrals[2];
        double fmr=  f(mr,p)*work.inner_integrals[3];
        double fmrr= f(mrr,p)*work.inner_integrals[4];

        double I_4=h/6.0*(fa+fb+5.0*(fml+fmr));   // 4-point Gauss-Lobatto formula.
        double I_7=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);        

        if (fabs(I_7-I_4) <= toler*I_13 || mll <= a || b <= mrr) 
        {
          if ((mll <= a || b <= mrr))// && !terminate_gkq) //Error
           {
             // out_of_tolerance=true; // Interval contains no more machine numbers
             printf("outer task OUT OF TOLERANCE !!!, mll:%.4e, a:%.4e, b:%.4e, mrr:%.4e\n", mll,b,b,mrr);
           }
           //#pragma omp critical(integal_result)
           {
              #pragma omp atomic
              integral_result += I_7;
              printf("myid:%d,iter:%d, I7 added, integs: 0=%f,1=%f,2=%f,3=%f,4=%f\n", myid,iter, 
                    work.inner_integrals[0],work.inner_integrals[1],work.inner_integrals[2],
                    work.inner_integrals[3],work.inner_integrals[4]);
           }
        }
        else  //subdivide interval and push new work on stack
        {

printf("myid:%d,iter:%d, iouter_integral approx not accurate: I_7=%.4e,I_4=%.4e,I_13=%.4e, toler*I_13:%.4e\n"
        "integs: 0=%f,1=%f,2=%f,3=%f,4=%f..subdivide\n",
        myid,iter,I_7,I_4,I_13,I_13*toler,work.inner_integrals[0],work.inner_integrals[1],work.inner_integrals[2],
                    work.inner_integrals[3],work.inner_integrals[4] );

          stack_t stack_inner = work.stack_inner;
          #pragma omp critical (stack)          
          { 
            //create sub-stack for inner-integrals
            
            work.subtasks_left=5; //for all new outer work elements

            create_stack(&stack_inner, sizeof(work_gkq));            
            work.is_parent=1;

            work_gkq winner;         

            //outer task elem
            work.a=a;
            work.b=mll;
            work.fa=fa;
            work.fb=fmll;          
            
            double m_inner=  (work.a+work.b)/2;
            double h_inner=  (work.a+work.b)/2;
            double mll_inner=(m_inner-alpha*h_inner);
            double ml_inner= (m_inner-beta*h_inner);
            double mr_inner= (m_inner+beta*h_inner);
            double mrr_inner=(m_inner+alpha*h_inner);

            //Jedes subinterval bekommt 5 inner tasks für die inneren integrale!

            int curcnt=stack->el_count;
            work.task_nr=curcnt; //weil ich bei 0 anfange zu zählen und push_stack gleich el_count incrementiert

            double integ_bnd[2];            

            //1st subtask
            get_integ_bounds_inner(integ_bnd, mll_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.task_nr=work.task_nr;
            winner.subtask_nr=0;
            push_stack(work.stack_inner,&winner);

            //2nd subtask
            get_integ_bounds_inner(integ_bnd, ml_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=1;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);

            //3rd subtask
            get_integ_bounds_inner(integ_bnd, m_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=2;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);     

            //4th subtask
            get_integ_bounds_inner(integ_bnd, mr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=3;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);     

            //5th subtask
            get_integ_bounds_inner(integ_bnd, mrr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=4;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
              
            push_stack(stack, &work);


            // *************************
            // *   2nd subinterval
            // **************************
            work.a=mll;
            work.b=ml;
            work.fa=fmll;
            work.fb=fml;
            work.task_nr++;

            m_inner=  (work.a+work.b)/2;
            h_inner=  (work.a+work.b)/2;
            mll_inner=(m_inner-alpha*h_inner);
            ml_inner= (m_inner-beta*h_inner);
            mr_inner= (m_inner+beta*h_inner);
            mrr_inner=(m_inner+alpha*h_inner);


            //1st subtask
            get_integ_bounds_inner(integ_bnd, mll_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=0;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
            

            //2nd subtask
            get_integ_bounds_inner(integ_bnd, ml_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=1;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
              

            //3rd subtask
            get_integ_bounds_inner(integ_bnd, m_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=2;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //4th subtask
            get_integ_bounds_inner(integ_bnd, mr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=3;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //5th subtask
            get_integ_bounds_inner(integ_bnd, mrr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=4;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);

            // #pragma omp critical (stack)
            // {
              push_stack(stack, &work);
            // }                 

            // *************************
            // *   3rd subinterval
            // **************************
            work.a=ml;
            work.b=m;
            work.fa=fml;
            work.fb=fm;
            work.task_nr++;

            m_inner=  (work.a+work.b)/2;
            h_inner=  (work.a+work.b)/2;
            mll_inner=(m_inner-alpha*h_inner);
            ml_inner= (m_inner-beta*h_inner);
            mr_inner= (m_inner+beta*h_inner);
            mrr_inner=(m_inner+alpha*h_inner);            

            //1st subtask
            get_integ_bounds_inner(integ_bnd, mll_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=0;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
            

            //2nd subtask
            get_integ_bounds_inner(integ_bnd, ml_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=1;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
              

            //3rd subtask
            get_integ_bounds_inner(integ_bnd, m_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=2;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //4th subtask
            get_integ_bounds_inner(integ_bnd, mr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=3;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //5th subtask
            get_integ_bounds_inner(integ_bnd, mrr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2, integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=4;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);

            // #pragma omp critical (stack)
            // {
              push_stack(stack, &work);
            // }
            



            // *************************
            // *   4th subinterval
            // **************************
            work.a=m;
            work.b=mr;
            work.fa=fm;
            work.fb=fmr;
            work.task_nr++;

            m_inner=  (work.a+work.b)/2;
            h_inner=  (work.a+work.b)/2;
            mll_inner=(m_inner-alpha*h_inner);
            ml_inner= (m_inner-beta*h_inner);
            mr_inner= (m_inner+beta*h_inner);
            mrr_inner=(m_inner+alpha*h_inner);            

            //1st subtask
            get_integ_bounds_inner(integ_bnd, mll_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=0;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
            

            //2nd subtask
            get_integ_bounds_inner(integ_bnd, ml_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=1;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
              

            //3rd subtask
            get_integ_bounds_inner(integ_bnd, m_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=2;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //4th subtask
            get_integ_bounds_inner(integ_bnd, mr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=3;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //5th subtask
            get_integ_bounds_inner(integ_bnd, mrr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=4;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);

            // #pragma omp critical (stack)
            {
              push_stack(stack, &work);
            }                              

            // *************************
            // *   5th subinterval
            // **************************
            work.a=mr;
            work.b=mrr;
            work.fa=fmr;
            work.fb=fmrr;
            work.task_nr++;

            m_inner=  (work.a+work.b)/2;
            h_inner=  (work.a+work.b)/2;
            mll_inner=(m_inner-alpha*h_inner);
            ml_inner= (m_inner-beta*h_inner);
            mr_inner= (m_inner+beta*h_inner);
            mrr_inner=(m_inner+alpha*h_inner);            

            //1st subtask
            get_integ_bounds_inner(integ_bnd, mll_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=0;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
            

            //2nd subtask
            get_integ_bounds_inner(integ_bnd, ml_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=1;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
              

            //3rd subtask
            get_integ_bounds_inner(integ_bnd, m_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=2;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //4th subtask
            get_integ_bounds_inner(integ_bnd, mr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=3;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //5th subtask
            get_integ_bounds_inner(integ_bnd, mrr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=4;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);

            // #pragma omp critical (stack)
            {
              push_stack(stack, &work);
            }

            // *************************
            // *   6th subinterval
            // **************************
            work.a=mrr;
            work.b=b;
            work.fa=fmrr;
            work.fb=fb;
            work.task_nr++;

            m_inner=  (work.a+work.b)/2;
            h_inner=  (work.a+work.b)/2;
            mll_inner=(m_inner-alpha*h_inner);
            ml_inner= (m_inner-beta*h_inner);
            mr_inner= (m_inner+beta*h_inner);
            mrr_inner=(m_inner+alpha*h_inner);            

            //1st subtask
            get_integ_bounds_inner(integ_bnd, mll_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=0;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
            

            //2nd subtask
            get_integ_bounds_inner(integ_bnd, ml_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=1;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
              

            //3rd subtask
            get_integ_bounds_inner(integ_bnd, m_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=2;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //4th subtask
            get_integ_bounds_inner(integ_bnd, mr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=3;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);
                          

            //5th subtask
            get_integ_bounds_inner(integ_bnd, mrr_inner, work.p);            
            winner=gkq_create_inner_task(inner_integrand2,integ_bnd[0],integ_bnd[1],1e-3, work.p);
            winner.subtask_nr=4;
            winner.task_nr=work.task_nr;
            push_stack(work.stack_inner,&winner);

            // #pragma omp critical (stack)
            {
              push_stack(stack, &work);
            }           

            printf("myid:%d,iter:%d, 6 subinvtervals with 5 subtasks each pushed to stack outer\n",myid,iter);


          } // pragma critical
        }   // else ..non-acceptable error


      }


    } // while
  } /* end omp parallel */



  return integral_result;    
}









//SIMPSON

  // *****************************************************************************************************************************************
// double integral2(double (*f)(double, struct my_f_params*), stack_t stack)
// {
//   work_t work;
//   int ready, idle, busy;
//   double integral_result = 0.0;

//   busy = 0;
  
// #pragma omp parallel default(none) \
//     shared(stack, integral_result,f,busy,simpson_error) \
//     private(work, idle, ready)
//   {

// // printf("me:%d, err:%d\n",omp_get_thread_num(),simpson_error);    

//     ready = 0;
//     idle = 1;

//     while(!ready  && !simpson_error) //<-- so NICHT!
//     {
// #pragma omp critical (stack)
//       {
//         if (!empty_stack(stack))
//         {
//           /* we have new work */ 
//           pop_stack(stack, &work);
//           if (idle)
//           {
//             /* say others i'm busy */
//             busy += 1;
//             idle = 0;
//           }
//         }
//         else{
//           /* no new work on stack */
//           if (!idle){
//             busy -= 1;
//             idle = 1;
//           }

//           /* nobody has anything to do; let us leave the loop */
//           if (busy == 0)
//           {
//             ready = 1;        
//           }
//         }
//       } /* end critical(stack) */

//       if (idle)
//         continue; //if ready==1 --> leave loop

//       double b = work.b;
//       double a = work.a;
//       double tol = work.tol;
//       double S=work.S; // previous TOTAL integral
//       double fa=work.fa;
//       double fb=work.fb;
//       double fm=work.fm;
//       int    rec=work.rec;

//       double h = (b - a)/2;
//       double mid = (a+b)/2;
//       double lm=(a+mid)/2;
//       double rm=(mid+b)/2;      

//       double flm=f(lm,work.p);
//       double frm=f(rm,work.p);
      
//       double Sl=h/6*(fa+4*flm+fm);
//       double Sr=h/6*(fm+4*frm+fb);

//       double delta=Sl+Sr-S;

//       // serious numerical trouble: it won't converge
//       if ((tol/2 == tol) || fabs(a-lm) <= tol) // || tol < tolmax) 
//       { 
//              simpson_error = 1; 
// #pragma omp critical (integral_result)            
//              integral_result = S;      
//       }
//       //need something against spurious convergence
//       //for now: iter > 5 (c.f. numerical recipes)
//       if( work.iter > 5 &&  (rec <= 0 || fabs(delta) <= 15*tol))  //error acceptable
//       {
// #pragma omp critical (integral_result)
//         integral_result += Sl+Sr+delta/15;
//       }
//       else // error not acceptable
//       {          
//         //push new subintervals to stack
//         work.a = a;
//         work.b = mid;
//         work.tol = tol/2;
//         work.S = Sl;
//         work.fa=fa;
//         work.fb=fm;
//         work.fm=flm;
//         work.rec=rec-1;
//         work.iter=work.iter+1;

//   #pragma omp critical (stack)
//         {
//           //LEFT
//           push_stack(stack, &work);      
//           //prepare RIGHT side and push to stack
//           work.a = mid;
//           work.b = b;
//           work.tol = tol/2;
//           work.S=Sr;
//           work.fa=fm;
//           work.fb=fb;
//           work.fm=frm;
//           work.rec=rec-1;
//           push_stack(stack, &work);

//         }        
//       }
//     } /* while */
//   } /* end omp parallel */


//   return integral_result;
// }