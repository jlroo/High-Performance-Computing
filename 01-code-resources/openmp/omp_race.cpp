#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <my_timer.h>
#include <my_omp.h>

int global_counter = 0; // Shared (or global) variable.
double global_sum = 0;

int doSomeThing (int i)
{
   double x = sqrt(i*2);
   return int(x);
}

int test_critical(const int n)
{
   myTimer_t t0 = getTimeStamp();
   global_counter = 0;
   global_sum = 0;

   #pragma omp parallel for default(shared) 
   for (int i = 0; i < n; ++i)
   {
      double myValue = doSomeThing(i);
      #pragma omp critical
      {
         global_sum += myValue;
         if (myValue > 10)
            global_counter++;
      }
   }

   myTimer_t t1 = getTimeStamp();
   printf("critical: global_sum = %e %d %f\n", global_sum, global_counter, 1000*getElapsedTime(t0,t1));

   return 0;
}

int test_manual_reduction(const int n)
{
   myTimer_t t0 = getTimeStamp();
   global_counter = 0;
   global_sum = 0;

   #pragma omp parallel default(shared)
   {

   int my_counter = 0; // Private to each thread.
   double my_sum = 0;

   #pragma omp for nowait
   for (int i = 0; i < n; ++i)
   {
      double myValue = doSomeThing(i);
      my_sum += myValue;
      if (myValue > 10)
         my_counter++;
   }

   #pragma omp critical
   {
      global_counter += my_counter;
      global_sum += my_sum;
   }

   } // end parallel

   myTimer_t t1 = getTimeStamp();
   printf("manual: global_sum = %e %d %f\n", global_sum, global_counter, 1000*getElapsedTime(t0,t1));

   return 0;
}

int test_reduction(const int n)
{
   myTimer_t t0 = getTimeStamp();
   global_counter = 0;
   global_sum = 0;

   #pragma omp parallel for default(shared) reduction(+:global_sum,global_counter)
   for (int i = 0; i < n; ++i)
   {
      double myValue = doSomeThing(i);
      global_sum += myValue;
      if (myValue > 10)
         global_counter++;
   }

   myTimer_t t1 = getTimeStamp();
   printf("reduction: global_sum = %e %d %f\n", global_sum, global_counter, 1000*getElapsedTime(t0,t1));

   return 0;
}

#ifdef _MY_OPENMP_HAS_ATOMICS
int test_atomic(const int n)
{
   myTimer_t t0 = getTimeStamp();
   global_counter = 0;
   global_sum = 0;

   #pragma omp parallel for default(shared)
   for (int i = 0; i < n; ++i)
   {
      double myValue = doSomeThing(i);
      #pragma omp atomic
      global_sum += myValue;
      if (myValue > 10)
      {
         #pragma omp atomic
         global_counter++;
      }
   }

   myTimer_t t1 = getTimeStamp();
   printf("atomic: global_sum = %e %d %f\n", global_sum, global_counter, 1000*getElapsedTime(t0,t1));

   return 0;
}
#endif

int main (int argc, char* argv[])
{
#ifdef _OPENMP
   printf("OpenMP specification c. %s %d.%d\n", TOSTRING(_OPENMP), _MY_OPENMP_MAJOR, _MY_OPENMP_MINOR);
   #ifdef _MY_OPENMP_HAS_ATOMICS
   printf("OpenMP has atomics!\n");
   #endif

   int num_threads = 1;

   #pragma omp parallel default(shared)
   {
      #pragma omp master
      num_threads = omp_get_num_threads();
   }

   printf("num_threads = %d\n", num_threads);
#endif

   int n = 100000;

   for (int i = 1; i < argc; ++i)
   {
#define check_index(i,str) \
   if ((i) >= argc) \
      { fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }

      if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
      {
         fprintf(stderr,"omp_list --help|-h --length|-n\n");
         return 1;
      }
      else if (strcmp(argv[i],"--length") == 0 || strcmp(argv[i],"-n") == 0)
      {
         check_index(i+1,"--length|-n");
         i++;
         if (isdigit(*argv[i]))
            n = atoi( argv[i] );
      }
   }

   printf("length: %d\n", n);

   for (int iter = 0; iter < 10; iter++)
   {
      test_critical(n);

#ifdef _MY_OPENMP_HAS_ATOMICS
      test_atomic(n);
#endif

      test_manual_reduction(n);
      test_reduction(n);
   }

   return 0;
}
