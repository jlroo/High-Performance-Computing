#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <my_timer.h>
#include <my_omp.h>
#include <aligned_allocator.h>

double drand (void)
{
   return drand48();
}

// Estimate PI by integration: Int{ dx / (1 + x^2) } ~= pi/4
template <typename ValueType>
void calcPiInt (const int n, const int numThreads)
{
   myTimer_t t0 = getTimeStamp();

   ValueType step = 1.0 / ValueType(n);
   ValueType sum = 0;

   #pragma omp parallel for default(shared) reduction(+:sum)
   for (int i = 0; i < n; ++i)
   {
      ValueType x = (i + 0.5) * step;
      sum += 1.0 / (1.0 + x*x);
   }

   ValueType pi = sum * (4 * step);
   myTimer_t t1 = getTimeStamp();

   printf("Integration: pi = %e %f(%%) %d %f (ms) %d\n", pi, 100*fabs(M_PI-pi)/M_PI, n, 1000.*getElapsedTime(t0, t1), numThreads);

   return;
}

template <typename ValueType>
void calcPiInt_FS (const int n, const int numThreads)
{
   myTimer_t t0 = getTimeStamp();

   ValueType step = 1.0 / ValueType(n);

   ValueType threadSums[32][4];
#ifdef _OPENMP
   if (numThreads > 32)
   {
      fprintf(stderr,"omp_get_max_threadS() > 32 here %d %d %s\n", numThreads, __LINE__, __FILE__);
      return;
   }
#endif

   //for (int k = 0; k < numThreads; ++k)
   //   threadSums[k][0] = 0;

   #pragma omp parallel default(shared) shared(threadSums)
   {
      int thread_id = 0;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#endif

      threadSums[thread_id][0] = 0;

      #pragma omp for nowait
      for (int i = 0; i < n; ++i)
      {
         ValueType x = (i + 0.5) * step;
         threadSums[thread_id][0] += 1.0 / (1.0 + x*x);
      }
   }

   // Sum up the private accumulators.
   ValueType sum = 0;
   for (int k = 0; k < numThreads; ++k)
      sum += threadSums[k][0];

   ValueType pi = sum * (4 * step);
   myTimer_t t1 = getTimeStamp();

   printf("IntegrationFS: pi = %e %f(%%) %d %f (ms) %d\n", pi, 100*fabs(M_PI-pi)/M_PI, n, 1000.*getElapsedTime(t0, t1), numThreads);

   return;
}

// Estimate PI by Monte Carlo method.
void calcPiMC (const int n, const int numThreads)
{
   myTimer_t t0 = getTimeStamp();
   int counter = 0;

   #pragma omp parallel for default(shared) reduction(+:counter)
   for (int i = 0; i < n; ++i)
   {
      double x = drand();
      double y = drand();
      if ((x*x + y*y) <= 1.0)
         counter ++;
   }

   double pi = 4.0 * (double(counter) / n);

   myTimer_t t1 = getTimeStamp();

   printf("Monte-Carlo: pi = %e %f(%%) %d %f (ms) %d\n", pi, 100*fabs(M_PI-pi)/M_PI, n, 1000.*getElapsedTime(t0, t1), numThreads);

   return;
}

int main(int argc, char *argv[])
{
   int len = 1024; // Length of the array.
   int method = 0; // Which method?

   for (int i = 1; i < argc; ++i)
   {
#define check_index(i,str) \
   if ((i) >= argc) \
      { fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }

      if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
      {
         fprintf(stderr,"omp_false_sharing --help|-h --length|-n --method|-m\n");
         return 1;
      }
      else if (strcmp(argv[i],"--length") == 0 || strcmp(argv[i],"-n") == 0)
      {
         check_index(i+1,"--length|-n");
         i++;
         if (isdigit(*argv[i]))
            len = atoi( argv[i] );
      }
      else if (strcmp(argv[i],"--method") == 0 || strcmp(argv[i],"-m") == 0)
      {
         check_index(i+1,"--method|-m");
         i++;
         if (isdigit(*argv[i]))
            method = atoi( argv[i] );
      }
   }

   int max_threads = 1;
#ifdef _OPENMP
   max_threads = omp_get_max_threads();
#endif

   printf("OpenMP Parallel False-sharing example: length= %d, num_threads= %d, method= %d\n", len, max_threads, method);

   if (method == 1)
      calcPiInt<double>(len, max_threads);
   else if (method == 2)
      calcPiInt_FS<double>(len, max_threads);
   else
      calcPiMC(len, max_threads);

   return 0;
}
