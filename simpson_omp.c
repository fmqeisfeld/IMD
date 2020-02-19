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
const double tolmax=1e-30;
#define INITIAL_STACK_SIZE 128   /* initial size of new stacks */


static const double xgk[8] =    /* abscissae of the 15-point kronrod rule */
{
  0.991455371120812639206854697526329,
  0.949107912342758524526189684047851,
  0.864864423359769072789712788640926,
  0.741531185599394439863864773280788,
  0.586087235467691130294144838258730,
  0.405845151377397166906606412076961,
  0.207784955007898467600689403773245,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule */

static const double wg[4] =     /* weights of the 7-point gauss rule */
{
  0.129484966168869693270611432679082,
  0.279705391489276667901467771423780,
  0.381830050505118944950369775488975,
  0.417959183673469387755102040816327
};

static const double wgk[8] =    /* weights of the 15-point kronrod rule */
{
  0.022935322010529224963732008058970,
  0.063092092629978553290700663189204,
  0.104790010322250183839876322541518,
  0.140653259715525918745189590510238,
  0.169004726639267902826583426598550,
  0.190350578064785409913256402421014,
  0.204432940075298892414161999234649,
  0.209482141084727828012999174891714
};

void 
qkrule (const int n, 
                    const double xgk[], const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function * f, double a, double b,
                    double *result, double *abserr,
                    double *resabs, double *resasc);




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


void create_stack(stack_t* stack, int element_size);
int empty_stack(stack_t stack);
void push_stack(stack_t stack, void* element);
void pop_stack(stack_t stack, void* element);


double  
integral2(
     double (*f)(double, struct my_f_params*), /* function to integrate */
     stack_t stack);


static double myfun(double x,struct my_f_params* p)
{  
  double T=p->T;  
  double fermi=1/(1+exp(-x/T));
  double fermi2=1- 1/(1+exp(-x/T));
  double sigma=x/T*log(x/T*2)/pow(1/x,2.0); //einfach nur ärger machen

  // return  exp(-x*x)/p->T+log(x*pow(p->T,2.0))/p->ne;
  return fermi*fermi2*sigma;
}

static double myfun2(double x,void* pv)
{ 
  struct my_f_params *p= (struct my_f_params*) pv;  

  double T=p->T;  
  double fermi=1/(1+exp(-x/T));
  double fermi2=1- 1/(1+exp(-x/T));
  double sigma=x/T*log(x/T*2)/pow(1/x,2.0); //einfach nur ärger machen

  
  // return  exp(-x*x)/p->T+log(x*pow(p->T,2.0))/p->ne;
  return fermi*fermi2*sigma;
}
// static double myfun(double x,struct my_f_params* p){  funcevals++; return  exp(-x*x);}; ///p->T+log(x*pow(p->T,2.0))/p->ne;}
// static float sinfc(float x) { return sinf(x); } 
// static double frand48c(double x) { funcevals++; return (double) drand48(); } 

// *************************************************************************************************************
int main(int argc,char** argv)
{
  num_threads=omp_get_num_threads();
  funcevals=0;
  simpson_error=0;

  struct timeval start, end;
    
  double gsltime, simpsontime;
  // ************************************************
  gsl_integration_workspace * gswork=NULL;
  gswork=   gsl_integration_workspace_alloc (1000); 


  double gslinteg=0.0;
  double integ_err;

  gsl_function gslfun;
  gslfun.function=&myfun2;


  // **********************************************

  double pi=3.141592653589793;
  double xmin, xmax;
  double answer = 0.0;

  int i;

  // xmin = -4.0;
  // xmax = 4.0;

  xmin=2;
  xmax=100; 

  /* prepare stack */
  stack_t stack;
  work_t work;  
  create_stack(&stack, sizeof(work_t));
  
  // ***********************************************
double T0=1;
double ne0=1e5;
for(i=0;i<10;i++)  
{  
  answer=0.0;
  struct my_f_params ptest;
  ptest.T=T0+i;
  ptest.ne=pow(ne0,((double) i)/2.0);


// *****************************************************
  gslfun.params=&ptest;  
  int method=GSL_INTEG_GAUSS61;

  gettimeofday(&start, NULL); 
  gsl_integration_qag(&gslfun, xmin, xmax, 0.0, 1e-6, 1000,method,
                       gswork, &gslinteg, &integ_err);  
  gettimeofday(&end, NULL);
  gsltime = ((end.tv_sec  - start.tv_sec) * 1000000u + 
              end.tv_usec - start.tv_usec) /1.e6;

// *********************************************                

  double fa=myfun(xmin,&ptest);
  double fb=myfun(xmax,&ptest);
  double fm=myfun((xmin+xmax)/2, &ptest);
  double h=xmax-xmin;
  double Sinit=h/6*(fa+4*fm+fb);
  // printf("Sinit:%.4e\n",Sinit);

  work.a = xmin;
  work.b = xmax;
  work.tol = 1e-6;
  work.S=Sinit;
  work.fa=fa;
  work.fb=fb;
  work.rec=100;
  work.fm=fm;
  work.p=&ptest;
  work.iter=0;


  push_stack(stack, &work);

  //#pragma omp parallel
  {
    gettimeofday(&start, NULL); 

    answer = integral2(myfun, stack);

    gettimeofday(&end, NULL);
    simpsontime = ((end.tv_sec  - start.tv_sec) * 1000000u + 
             end.tv_usec - start.tv_usec) /1.e6;

  } /* omp parallel */

  if(simpson_error)
  {
    // answer=0.0;
    printf("ERROR: Integral wont converge..cleanup\n");        

    while(!empty_stack(stack))
    {
      pop_stack(stack, &work);
    }    
    simpson_error=0;
  }  

  // 466 enum
  // 467   {
  // 468     GSL_INTEG_GAUSS15 = 1,      /* 15 point Gauss-Kronrod rule */
  // 469     GSL_INTEG_GAUSS21 = 2,      /* 21 point Gauss-Kronrod rule */
  // 470     GSL_INTEG_GAUSS31 = 3,       31 point Gauss-Kronrod rule 
  // 471     GSL_INTEG_GAUSS41 = 4,      /* 41 point Gauss-Kronrod rule */
  // 472     GSL_INTEG_GAUSS51 = 5,      /* 51 point Gauss-Kronrod rule */
  // 473     GSL_INTEG_GAUSS61 = 6       /* 61 point Gauss-Kronrod rule */
  // 474   };



  printf("simpson: %.12f,dt:%.4e,  gsl:%.12f, dt:%.4e\n",answer,simpsontime, gslinteg, gsltime);  

}   // for i
free(stack);
  




//   stack->el_count--;
//   memcpy(element, 
//          (char*)stack->elements + stack->el_count*stack->el_size, 
//          stack->el_size);





  return 0;
}


double integral2(double (*f)(double, struct my_f_params*), stack_t stack)
{
  work_t work;
  int ready, idle, busy;
  double integral_result = 0.0;

  busy = 0;
  
#pragma omp parallel default(none) \
    shared(stack, integral_result,f,busy,simpson_error) \
    private(work, idle, ready)
  {

// printf("me:%d, err:%d\n",omp_get_thread_num(),simpson_error);    

    ready = 0;
    idle = 1;

    while(!ready  && !simpson_error) //<-- so NICHT!
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
        else{
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

      double b = work.b;
      double a = work.a;
      double tol = work.tol;
      double S=work.S; // previous TOTAL integral
      double fa=work.fa;
      double fb=work.fb;
      double fm=work.fm;
      int    rec=work.rec;

      double h = (b - a)/2;
      double mid = (a+b)/2;
      double lm=(a+mid)/2;
      double rm=(mid+b)/2;      

      double flm=f(lm,work.p);
      double frm=f(rm,work.p);
      
      double Sl=h/6*(fa+4*flm+fm);
      double Sr=h/6*(fm+4*frm+fb);

      double delta=Sl+Sr-S;

      // serious numerical trouble: it won't converge
      if ((tol/2 == tol) || fabs(a-lm) <= tol) // || tol < tolmax) 
      { 
             simpson_error = 1; 
#pragma omp critical (integral_result)            
             integral_result = S;      
      }
      //need something against spurious convergence
      //for now: iter > 5 (c.f. numerical recipes)
      if( work.iter > 5 &&  (rec <= 0 || fabs(delta) <= 15*tol))  //error acceptable
      {
#pragma omp critical (integral_result)
        integral_result += Sl+Sr+delta/15;
      }
      else // error not acceptable
      {          
        //push new subintervals to stack
        work.a = a;
        work.b = mid;
        work.tol = tol/2;
        work.S = Sl;
        work.fa=fa;
        work.fb=fm;
        work.fm=flm;
        work.rec=rec-1;
        work.iter=work.iter+1;

  #pragma omp critical (stack)
        {
          //LEFT
          push_stack(stack, &work);      
          //prepare RIGHT side and push to stack
          work.a = mid;
          work.b = b;
          work.tol = tol/2;
          work.S=Sr;
          work.fa=fm;
          work.fb=fb;
          work.fm=frm;
          work.rec=rec-1;
          push_stack(stack, &work);

        }        
      }
    } /* while */
  } /* end omp parallel */


  return integral_result;
}

/******************************************
 * create new stack
 ******************************************/
void 
create_stack(
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
int 
empty_stack
(stack_t stack)
{
  return stack->el_count <= 0;
}


/*****************************************
 * push a element on stack
 *****************************************/
void 
push_stack(
           stack_t stack,    /* target stack */
           void* element)    /* element to push */
{
  int i, new_reserve;
  int log2_count;

  /* check if we need more memory for stack */    
  if (stack->el_count >= stack->mem_reserve){

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





void 
qkrule (const int n, 
                    const double xgk[], const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function * f, double a, double b,
                    double *result, double *abserr,
                    double *resabs, double *resasc)
{

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs (half_length);
  const double f_center = GSL_FN_EVAL (f, center);

  double result_gauss = 0;
  double result_kronrod = f_center * wgk[n - 1];

  double result_abs = fabs (result_kronrod);
  double result_asc = 0;
  double mean = 0, err = 0;

  int j;

  if (n % 2 == 0)
    {
      result_gauss = f_center * wg[n / 2 - 1];
    }

  for (j = 0; j < (n - 1) / 2; j++)
    {
      const int jtw = j * 2 + 1;  /* in original fortran j=1,2,3 jtw=2,4,6 */
      const double abscissa = half_length * xgk[jtw];
      const double fval1 = GSL_FN_EVAL (f, center - abscissa);
      const double fval2 = GSL_FN_EVAL (f, center + abscissa);
      const double fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss += wg[j] * fsum;
      result_kronrod += wgk[jtw] * fsum;
      result_abs += wgk[jtw] * (fabs (fval1) + fabs (fval2));
    }

  for (j = 0; j < n / 2; j++)
    {
      int jtwm1 = j * 2;
      const double abscissa = half_length * xgk[jtwm1];
      const double fval1 = GSL_FN_EVAL (f, center - abscissa);
      const double fval2 = GSL_FN_EVAL (f, center + abscissa);
      fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk[jtwm1] * (fval1 + fval2);
      result_abs += wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
    };

  mean = result_kronrod * 0.5;

  result_asc = wgk[n - 1] * fabs (f_center - mean);

  for (j = 0; j < n - 1; j++)
    {
      result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
    }

  /* scale by the width of the integration region */

  err = (result_kronrod - result_gauss) * half_length;

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  *result = result_kronrod;
  *resabs = result_abs;
  *resasc = result_asc;
  *abserr = rescale_error (err, result_abs, result_asc);

}




static int
qag (const gsl_function * f,
     const double a, const double b,
     const double epsabs, const double epsrel,
     const size_t limit,
     gsl_integration_workspace * workspace,
     double *result, double *abserr)
{
  double area, errsum;
  double result0, abserr0, resabs0, resasc0;
  double tolerance;
  size_t iteration = 0;
  int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

  double round_off;     

  /* Initialize results */
  // initialise (workspace, a, b);

  *result = 0;
  *abserr = 0;

  /* perform the first integration */
  q (f, a, b, &result0, &abserr0, &resabs0, &resasc0);
  set_initial_result (workspace, result0, abserr0);
  *result = result0;
  *abserr = abserr0;

  area = result0;
  errsum = abserr0;

  iteration = 1;

  do
    {
      double a1, b1, a2, b2;
      double a_i, b_i, r_i, e_i;
      double area1 = 0, area2 = 0, area12 = 0;
      double error1 = 0, error2 = 0, error12 = 0;
      double resasc1, resasc2;
      double resabs1, resabs2;

      /* Bisect the subinterval with the largest error estimate */

      retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

      a1 = a_i; 
      b1 = 0.5 * (a_i + b_i);
      a2 = b1;
      b2 = b_i;

      q (f, a1, b1, &area1, &error1, &resabs1, &resasc1);
      q (f, a2, b2, &area2, &error2, &resabs2, &resasc2);

      area12 = area1 + area2;
      error12 = error1 + error2;

      errsum += (error12 - e_i);
      area += area12 - r_i;

      if (resasc1 != error1 && resasc2 != error2)
        {
          double delta = r_i - area12;

          if (fabs (delta) <= 1.0e-5 * fabs (area12) && error12 >= 0.99 * e_i)
            {
              roundoff_type1++;
            }
          if (iteration >= 10 && error12 > e_i)
            {
              roundoff_type2++;
            }
        }

      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));

      update (workspace, a1, b1, area1, error1, a2, b2, area2, error2);
      retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

      iteration++;

    }
  while (iteration < limit && !error_type && errsum > tolerance);

  *result = sum_results (workspace);
  *abserr = errsum;
}


static double
rescale_error (double err, const double result_abs, const double result_asc)
{
  err = fabs(err) ;

  if (result_asc != 0 && err != 0)
      {
        double scale = pow((200 * err / result_asc), 1.5) ;
        
        if (scale < 1)
          {
            err = result_asc * scale ;
          }
        else 
          {
            err = result_asc ;
          }
      }
  if (result_abs > GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))
    {
      double min_err = 50 * GSL_DBL_EPSILON * result_abs ;

      if (min_err > err) 
        {
          err = min_err ;
        }
    }
  
  return err ;
}



static inline
void update (gsl_integration_workspace * workspace,
             double a1, double b1, double area1, double error1,
             double a2, double b2, double area2, double error2)
{
  double * alist = workspace->alist ;
  double * blist = workspace->blist ;
  double * rlist = workspace->rlist ;
  double * elist = workspace->elist ;
  size_t * level = workspace->level ;

  const size_t i_max = workspace->i ;
  const size_t i_new = workspace->size ;

  const size_t new_level = workspace->level[i_max] + 1;

  /* append the newly-created intervals to the list */
  
  if (error2 > error1)
    {
      alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;
      
      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
    }
  else
    {
      blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;
      
      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
    }
  
  workspace->size++;

  if (new_level > workspace->maximum_level)
    {
      workspace->maximum_level = new_level;
    }

  qpsrt (workspace) ;
}



static inline void
retrieve (const gsl_integration_workspace * workspace, 
          double * a, double * b, double * r, double * e)
{
  const size_t i = workspace->i;
  double * alist = workspace->alist;
  double * blist = workspace->blist;
  double * rlist = workspace->rlist;
  double * elist = workspace->elist;

  *a = alist[i] ;
  *b = blist[i] ;
  *r = rlist[i] ;
  *e = elist[i] ;
}


static inline double
sum_results (const gsl_integration_workspace * workspace)
{
  const double * const rlist = workspace->rlist ;
  const size_t n = workspace->size;

  size_t k;
  double result_sum = 0;

  for (k = 0; k < n; k++)
    {
      result_sum += rlist[k];
    }
  
  return result_sum;
}