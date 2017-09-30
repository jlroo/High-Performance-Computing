#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void help()
{
   fprintf(stderr,"omp_for --help|-h --length|-n\n");

   return;
}

int main(int argc, char* argv[])
{
   int len = 1024; // Length of the array.

   for (int i = 1; i < argc; ++i)
   {
#define check_index(i,str) \
   if ((i) >= argc) \
      { fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }

      if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
      {
         help();
         return 1;
      }
      else if (strcmp(argv[i],"--length") == 0 || strcmp(argv[i],"-n") == 0)
      {
         check_index(i+1,"--length|-n");
         i++;
         if (isdigit(*argv[i]))
            len = atoi( argv[i] );
      }
   }

   int max_threads = 1;
#ifdef _OPENMP
   max_threads = omp_get_max_threads();
#endif

   printf("OpenMP Parallel For example: length= %d, num_threads= %d\n", len, max_threads);

   float *x = new float [len];
   float *y = new float [len];
   float *z = new float [len];

   float xsum = 0;

   #pragma omp parallel \
           default(none) shared(x,y,z,xsum,len)
           /* Control the data sharing mode:
                default(none) = no implicit definitions.
                shared(...) = explicitly declare this as shared. */
   {

      #pragma omp for
      for (int i = 0; i < len; ++i)
      {
         x[i] = y[i] = i;
      }
      // <-- implicit barrier here.

      #pragma omp for reduction(+:xsum)
      for (int i = 0; i < len; ++i)
      {
         z[i] = x[i] + 2*y[i];
         xsum += x[i];
      }

      #pragma omp single
      {
         int thread_id = 0;
#ifdef _OPENMP
         thread_id = omp_get_thread_num();
#endif
         printf("xsum = %g %g %d\n", xsum, (len-1.0)*len/2.0, thread_id);
      }

   } // <-- implicit barrier here.

   printf("All done\n");

   return 0;
}
