#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm> // std::min
#include <cmath> // std::log, etc.

#include <my_timer.h>

int main (int argc, char * argv[])
{
   const int nmin = 30;
   int nsizes = 34; // Number of tests between [nmin,nmax]
   int nmax = 16000000;
   int niters = 50; // Number of samples for each test.
   const int npad = 5; // Array padding.

   // Allow the user to change the defaults.
   if (argc > 1)
      if (isdigit(*argv[1]))
         nmax = atoi(argv[1]);

   if (argc > 2)
      if (isdigit(*argv[2]))
         nsizes = atoi(argv[2]);

   if (argc > 3)
      if (isdigit(*argv[3]))
         niters = atoi(argv[3]);

   typedef double ValueType;

   // Allocate arrays with a little padding.
   ValueType *a = new ValueType[nmax+5];
   ValueType *b = new ValueType[nmax+5];

   // Timing and statistics.
   double (*time)[4] = new double[niters][4];

   // Check the timer accuracy.
   double tDelta = 1e50;
   {
      for (int i = 0; i < std::min(10000,nmax); ++i)
         a[i] = 0.0;

      myTimer_t t0 = getTimeStamp();
      for (int i = 0; i < std::min(10000,nmax); ++i)
      {
         myTimer_t t1 = getTimeStamp();
         a[i] = getElapsedTime(t0,t1);
         //a[i] = getElapsedTime();
         t0 = t1;
      }

      for (int i = 0; i < std::min(10000,nmax)-1; ++i)
         if (std::abs(a[i+1] - a[i]) > 0.0)
            tDelta = std::min(tDelta, std::abs(a[i+1]-a[i]));

      printf("Smallest detectable time = %e (ms)\n", tDelta*1000);

      printf("getTicksPerSecond = %e\n", getTicksPerSecond());
   }

   printf("Size, Fill, Copy, AXPY, DOT, Uncertainty\n");

   // Loop over the problem sizes.
   for (int size = 0; size < nsizes; size++)
   {
      // Determine the size for this test.
      const double expon = log(double(nmin)) +
                     (double(size)/(nsizes-1)) * (log(double(nmax)) - log(double(nmin)));

      const int m = rint(exp(expon));
      //printf("%d %f %d %f\n", m, exp(expon), size, log(double(nmax)) - log(double(nmin)));

      // Zero the arrays.
      for (int i = 0; i < m; ++i)
         a[i] = b[i] = 0.0;

      // Run this test many times depending on the size.
      for (int k = 0; k < niters; ++k)
      {
         const int nloop = nmax / m;

         // FILL (a[] := s): set to scalar value.
         {
            myTimer_t t0 = getTimeStamp();
            for (int loop = 0; loop < nloop; ++loop)
            {
               const ValueType s(k+loop);
               for (int i = 0; i < m; ++i)
                  a[i] = s;
            }
            double t1 = getElapsedTime( t0, getTimeStamp());
            time[k][0] = t1/double(nloop);
         }

         // COPY (a[] := b[]): copy arrays.
         {
            myTimer_t t0 = getTimeStamp();
            for (int loop = 0; loop < nloop; ++loop)
            {
               a[loop] = ValueType(k+loop);
               for (int i = 0; i < m; ++i)
                  a[i] = b[i];
            }
            double t1 = getElapsedTime( t0, getTimeStamp());
            time[k][1] = t1/double(nloop);
         }

         // AXPY (b[] := b[] + s*a[]): Vector + Vector*scalar.
         {
            myTimer_t t0 = getTimeStamp();
            for (int loop = 0; loop < nloop; ++loop)
            {
               const ValueType s(k+loop);
               for (int i = 0; i < m; ++i)
                  b[i] = b[i] + s*a[i];
            }
            double t1 = getElapsedTime( t0, getTimeStamp());
            time[k][2] = t1/double(nloop);
         }

         // DOT (s = Sum_i (a[]*b[]): Inner product of two vectors.
         {
            ValueType sumAll(0);
            myTimer_t t0 = getTimeStamp();
            for (int loop = 0; loop < nloop; ++loop)
            {
               b[loop] = ValueType(k+loop);
               ValueType sum(0);
               for (int i = 0; i < m; ++i)
                  sum += a[i]*b[i];
               sumAll += sum;
            }
            double t1 = getElapsedTime( t0, getTimeStamp());
            time[k][3] = t1/double(nloop);
            srand((int)sumAll); // Force the compiler to keep sumAll.
         }
      }

      // Pick the best times for each test.
      double bestTime[4], aveTime[4], testTime[4], rate[4];
      double numBytes[] = {8,16,24,16};
      for (int test = 0; test < 4; test++)
      {
         bestTime[test] = 1e50;
         aveTime[test] = 0.0;
         for (int k = 0; k < niters; ++k)
         {
            bestTime[test] = std::min(bestTime[test], time[k][test]);
            aveTime[test] += time[k][test];
            //printf("%d %d %f\n", test, k, time[k][test]*1000);
         }

         aveTime[test] /= double(niters);
         testTime[test] = aveTime[test];
         rate[test] = double(m) * numBytes[test] / testTime[test] / (1024.*1024.);
      }
      printf("%d, %f, %f, %f, %f, %.1f%%\n", m, rate[0], rate[1], rate[2], rate[3], 100*tDelta / testTime[0]);
   }

   return 0;
}
