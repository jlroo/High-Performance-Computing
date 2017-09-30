#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm> // std::min
#include <cmath> // std::log, etc.
#include <string> // std::log, etc.

#define HighResTimer
#include <my_timer.h>

#include <aligned_allocator.h>

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
         {
            // Flush the matrix.
            for (int i = 0; i < m; ++i)
               for (int j = 0; j < n; ++j)
                  matrix[INDEX(i,j,m,n)] = (double) INDEX(i,j,m,n);

            // Touch the matrix columns.
            uint64_t ticks_start = getClockTicks();

            for (int iter = 0; iter < niters; ++iter)
            {
               for (int j = 0; j < n; ++j)
                  for (int i = 1; i < m; ++i)
                      matrix[INDEX(0,j,m,n)] += matrix[INDEX(i,j,m,n)];
            }

            uint64_t ticks_stop = getClockTicks();

            printf("m = %d n = %d pad = %d ticks = %g kbytes= %.3f\n", m, n, npad, (double(ticks_stop-ticks_start)/(1.*niters))/(1.*n)/(1.*m), sizeof(ValueType)*(n+npad)*m/1024.);
         }
      }
   }

   free(matrix);
}

template <typename ValueType, int npad>
void stride_test (const int miters, const int stride, const int L1SIZE = 32)
{
   const double nticks_per_sec = getTicksPerSecond();
   const int L1WORDS = L1SIZE / sizeof(ValueType); // # of words

   const int melems = L1WORDS * stride * 2;
   ValueType *array = (ValueType *) aligned_alloc<ValueType>(melems+npad, 64 /*bytes*/);

   const ValueType alpha = 1.1;
   const ValueType beta = 2.2;

   // Flush the array and clear the cache.
   for (int i = 0; i < melems; ++i)
      array[i+npad] = ValueType(0);

   int niters = miters;
   while (1)
   {
      ValueType sum(0);

      // Touch the matrix columns.
      uint64_t ticks_start = getClockTicks();
      //myTimer_t t0 = getTimeStamp();
      for (int iter = 0; iter < niters; ++iter)
      {
         for (int i = 0; i < L1WORDS; ++i)
            array[i*stride+npad] = alpha * array[i*stride+npad] + beta;

         sum += array[iter%L1WORDS+npad];
      }

      //myTimer_t t1 = getTimeStamp();
      uint64_t ticks_stop = getClockTicks();
      double ticks = ticks_stop-ticks_start;
      double tCalc = ticks/nticks_per_sec;
      if (tCalc > 0.01)
      {
         //printf("stride = %d ticks = %g\n", stride, 1000.*tCalc/(1.*niters)/(1.*L1WORDS));
         printf("stride = %d ticks = %g %g %d\n", stride, ticks/(niters*L1WORDS), tCalc/(niters*L1WORDS), niters);
         break;
      }
      else
         niters *= 2;
   }

   free(array);
}

void show_usage( const char* prog )
{
   printf("Usage for %s\n", prog);
   printf("\t--niters | -n <int value> : Number of iterations. (1000)\n");
   printf("\t--assoc  | -a             : Run associativity test. (default)\n");
   printf("\t--stride | -s             : Run stride test.\n");
}

int main (int argc, char * argv[])
{
   int niters = 1000; // Number of samples for each test.
   int test = 0;

   // Get user inputs.
   {
      #define check_index(_i) { if ((_i) >= argc){ fprintf(stderr,"Missing value for argument %s\n", for
      for (int i = 1; i < argc; i++)
      {
         std::string arg = argv[i];
         if (arg == "--niters" || arg == "-n")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for --niters\n"); show_usage(argv[0]); return 1; }
            niters = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "--assoc" || arg == "-a")
         {
            test = 0;
         }
         else if (arg == "--stride" || arg == "-s")
         {
            test = 1;
         }
         else if (arg == "--help" || arg == "-h")
         {
            show_usage(argv[0]); return 0;
         }
      }
   }

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

   typedef double ValueType;

   if (test == 0)
   {
      assoc_test<ValueType,/*npad*/0>( niters );
      assoc_test<ValueType,/*npad*/8>( niters );
   }
   else
   {
      for (int i = 1; i <= 129; ++i)
         stride_test<ValueType,0>( niters, i, 20480 );
   }

   return 0;
}
