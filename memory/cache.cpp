#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm> // std::min
#include <cmath> // std::log, etc.

#define HighResTimer
#include <my_timer.h>

#ifndef MMAX
#define MMAX (13)
#endif
#ifndef NMAX
#define NMAX (20480)
#endif
#ifndef NPAD
#define NPAD (0)
#endif

//#define INDEX(i,j,m,n) ((j) + (i)*(n+NPAD))
#define INDEX(i,j,m,n) ((j) + (i)*(n+npad))

template <typename ValueType, int npad>
void assoc_test (const int miters)
{
   const double nticks_per_sec = getTicksPerSecond();

   ValueType *matrix = (ValueType *) malloc(sizeof(ValueType) * MMAX * (NMAX+npad));

   // Run a series of tests touching part of the matrix to stress the cache.
   for (int m = 2; m < MMAX; ++m)
   {
      for (int n = 128; n < NMAX; n = 2*n)
      {
         int niters = miters;
         //while (niters <= miters)
         {
            // Flush the matrix.
            for (int i = 0; i < m; ++i)
               for (int j = 0; j < n; ++j)
                  matrix[INDEX(i,j,m,n)] = (double) INDEX(i,j,m,n);

            // Touch the matrix columns.
            uint64_t ticks_start = getClockTicks();
            //myTimer_t t0 = getTimeStamp();
            for (int iter = 0; iter < niters; ++iter)
            {
               for (int j = 0; j < n; ++j)
                  for (int i = 1; i < m; ++i)
                      matrix[INDEX(0,j,m,n)] += matrix[INDEX(i,j,m,n)];
            }

            uint64_t ticks_stop = getClockTicks();
            //double time = (ticks_stop - ticks_start) / nticks_per_sec;
            //myTimer_t t1 = getTimeStamp();
            //double time = getElapsedTime(t0,t1);
            //if (time > minTime or niters == miters)
            {
               printf("m = %d n = %d pad = %d ticks = %g kbytes= %.3f\n", m, n, npad, (double(ticks_stop-ticks_start)/(1.*niters))/(1.*n)/(1.*m), sizeof(ValueType)*(n+npad)*m/1024.);
               //printf("m = %d n = %d pad = %d ticks = %g niters= %d\n", m, n, npad, (double(ticks_stop-ticks_start)/(1.*niters))/(1.*n)/(1.*m), niters);
               //printf("m = %d n = %d pad = %d ticks = %g niters= %d %e\n", m, n, npad, (time/(1.*niters))/(1.*n)/(1.*m), niters, time);
               //break;
            }
            //else
            //   niters *= 2;
         }
      }
   }

   free(matrix);
}

template <typename ValueType>
void stride_test (const int niters)
{
   const double nticks_per_sec = getTicksPerSecond();
   const int MSTRIDE = 64;
   const int L1SIZE = 32 * 1024; // bytes
   const int L1WORDS = L1SIZE / sizeof(ValueType); // # of words

   const size_t melems = 2 * L1WORDS * MSTRIDE;
   ValueType *array = (ValueType *) malloc(sizeof(ValueType) * melems);

   const ValueType alpha = 1.1;
   const ValueType beta = 2.2;

   // Run a series of tests touching part of the matrix to stress the cache.
   for (int stride = 1; stride <= MSTRIDE; ++stride)
   //for (int stride = 1; stride <= MSTRIDE; stride *= 2)
   {
      // Flush the array and clear the cache.
      for (int i = 0; i < melems; ++i)
         array[i] = ValueType(0);

      // Touch the matrix columns.
      uint64_t ticks_start = getClockTicks();
      for (int iter = 0; iter < niters; ++iter)
      {
         for (int i = 0, n = 0; i < L1WORDS; ++i, n += stride)
            array[n] = alpha * array[n] + beta;
      }

      uint64_t ticks_stop = getClockTicks();
      printf("stride = %d ticks = %g\n", stride, double(ticks_stop-ticks_start)/(1.*niters)/(1.*L1WORDS));
   }

   free(array);
}

int main (int argc, char * argv[])
{
   int niters = 1000; // Number of samples for each test.

   if (argc > 1)
      if (isdigit(*argv[1]))
         niters = atoi(argv[1]);

   typedef double ValueType;

   // Check the timer accuracy.
   double tDelta = 1e50;
   {
      const int size = 10000;
      double *y = new double [size];

      for (int i = 0; i < size; ++i)
         y[i] = 0.0;

      myTimer_t t0 = getTimeStamp();
      for (int i = 0; i < size; ++i)
      {
         myTimer_t t1 = getTimeStamp();
         y[i] = getElapsedTime(t0,t1);
         t0 = t1;
      }

      for (int i = 0; i < size-1; ++i)
         if (std::abs(y[i+1] - y[i]) > 0.0)
            tDelta = std::min(tDelta, std::abs(y[i+1]-y[i]));

      printf("Smallest detectable time = %e (ms)\n", tDelta*1000);

      printf("getTicksPerSecond = %e\n", getTicksPerSecond());
   }

   //assoc_test<float >( niters );
   assoc_test<double,0>( niters );
   assoc_test<double,1>( niters );

   //stride_test<float >( niters );
   stride_test<double>( niters );

   return 0;
}
