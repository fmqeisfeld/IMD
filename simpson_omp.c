#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>



// ***************************************************************************************
// * Pre-Port Entwicklung und Testen eines OMP-parallisierten, adaptiven Simpson-integrators für 
// * die fiesen doppelintgrale in imd_colrad.c
// ***************************************************************************************
int num_threads;
int simpson_error;
int funcevals;
const double tolmax=1e-30;
#define INITIAL_STACK_SIZE 128   /* initial size of new stacks */




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


static double myfun(double x,struct my_f_params* p){  funcevals++; return  exp(-x*x)/p->T+log(x*pow(p->T,2.0))/p->ne;}
//static double myfun(double x,struct my_f_params* p){  funcevals++; return  exp(-x*x);}; ///p->T+log(x*pow(p->T,2.0))/p->ne;}
// static float sinfc(float x) { return sinf(x); } 
// static double frand48c(double x) { funcevals++; return (double) drand48(); } 

int main(int argc,char** argv)
{
  num_threads=omp_get_num_threads();
  funcevals=0;
  simpson_error=0;

  double pi=3.141592653589793;
  double xmin, xmax;
  double answer = 0.0;

  int i;

  xmin = -4.0;
  xmax = 4.0;

  xmin=60;
  xmax=200; 

  /* prepare stack */
  stack_t stack;
  work_t work;  
  create_stack(&stack, sizeof(work_t));
  
  // ***********************************************
double T0=10;
double ne0=1e10;
for(i=0;i<100;i++)  
{  
  struct my_f_params ptest;
  ptest.T=T0+10*i;
  ptest.ne=pow(ne0,i/2);

  double fa=myfun(xmin,&ptest);
  double fb=myfun(xmax,&ptest);
  double fm=myfun((xmin+xmax)/2, &ptest);
  double h=xmax-xmin;
  double Sinit=h/6*(fa+4*fm+fb);
  // printf("Sinit:%.4e\n",Sinit);

  work.a = xmin;
  work.b = xmax;
  work.tol = 1e-16;
  work.S=Sinit;
  work.fa=fa;
  work.fb=fb;
  work.rec=1000;
  work.fm=fm;
  work.p=&ptest;


  push_stack(stack, &work);

  //#pragma omp parallel
  {
    answer = integral2(myfun, stack);
  } /* omp parallel */
  printf("The integral is approximately %.12f,sqrt(pi):%.12f\n",answer,sqrt(pi));  

  if(simpson_error)
  {
    printf("ERROR: Integral wont converge..cleanup\n");        
    while(!empty_stack(stack))
    {
      pop_stack(stack, &work);
    }    
    simpson_error=0;
  }  
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

printf("me:%d, err:%d\n",omp_get_thread_num(),simpson_error);    

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
            int i=0;            
            for(i=0;i<stack->el_count;i++)
              free(stack->elements+i);            
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

      if(rec <= 0 || fabs(delta) <= 15*tol) //error acceptable
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



