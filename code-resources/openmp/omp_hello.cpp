#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
   int num_threads = 1; // When serial, there's still one thread.

   printf("Hello from master\n");
   #pragma omp parallel \
           default(none) shared(num_threads)
           /* Control the data sharing mode:
                default(none) = no implicit definitions.
                shared(num_threads) = explicitly declare this as shared. */
   {
      // Only one thread inquires about the team state.
      #pragma omp single
      {
#ifdef _OPENMP
         // This is an API function so guard in case OpenMP not enabled.
         num_threads = omp_get_num_threads();
#endif
      } // <-- implicit barrier here.

      printf("Hello from worker %d\n", omp_get_thread_num());

   } // <-- implicit barrier here.

   printf("Goodbye from master: total team size was %d\n", num_threads);

   return 0;
}
