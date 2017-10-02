#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm> // std::min
#include <cmath> // std::log, etc.
#include <string>

#include <my_timer.h>

void show_usage( const char* prog )
{
   printf("Usage for %s\n", prog);
   printf("\t-min    <int value> : Minimum array size to start. (30B)\n");
   printf("\t-max    <int value> : Maximum array size to finish. (16MB)\n");
   printf("\t-ntests <int value> : # of samples between min/max.   (34)\n");
}
int main (int argc, char * argv[])
{
   int nmin = 30;
   int nsizes = 34; // Number of tests between [nmin,nmax]
   int nmax = 16000000;
   int niters = 50; // Number of samples for each test.
   const int npad = 5; // Array padding.

   // Allow the user to change the defaults.
   {
      #define check_index(_i) { if ((_i) >= argc){ fprintf(stderr,"Missing value for argument %s\n", for
      for (int i = 1; i < argc; i++)
      {
         std::string arg = argv[i];
         if (arg == "-min")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for -min\n"); show_usage(argv[0]); return 1; }
            nmin = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "-max")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for -max\n"); show_usage(argv[0]); return 1; }
            nmax = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "-ntests")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for -ntests\n"); show_usage(argv[0]); return 1; }
            nsizes = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "--help" || arg == "-h")
         {
            show_usage(argv[0]); return 0;
         }
         else
         {
            fprintf(stderr,"Unknown option %s\n", arg.c_str());
            show_usage(argv[0]);
            return 1;
         }
      }
   }

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

   printf("Size (iterations), Fill (MB/s), Copy, AXPY, DOT, Uncertainty (%%)\n");

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

         bool allPassed = true;
         for (int test = 0; test < 4; test++)
            allPassed = allPassed && (time[k][test] > 1.0);

         if (allPassed)
            break;
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
