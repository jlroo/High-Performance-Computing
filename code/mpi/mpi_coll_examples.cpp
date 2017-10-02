#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <my_mpi_header.h>
 
/************************************************************
 * This is a simple collective program in MPI
 * ************************************************************/

void demo_broadcast (int count, MPI_Comm comm)
{
   int numProcs, myRank; // Withint this communicator
   MPI_Comm_size(comm,&numProcs);
   MPI_Comm_rank(comm,&myRank);

   int *array = new int[count];
   for (int i = 0; i < count; ++i)
      array[i] = -(1 + myRank);

   const int root = 0;

   MPI_Bcast (array, count, MPI_INT, root, comm);

   for (int i = 0; i < count; ++i)
      printf("I am %d of %d in comm: array[%d] = %d\n", myRank, numProcs, i, array[i]);

   delete [] array;
}
void demo_reduce (int count, MPI_Comm comm)
{
   int numProcs, myRank; // Withint this communicator
   MPI_Comm_size(comm,&numProcs);
   MPI_Comm_rank(comm,&myRank);

   int *arrayIn = new int[count];
   int *arrayOut = new int[count];
   for (int i = 0; i < count; ++i)
      arrayIn[i] = arrayOut[i] = -(1 + myRank);

   const int root = 0;

   MPI_Reduce (arrayIn, arrayOut, count, MPI_INT, MPI_SUM, root, comm);

   for (int i = 0; i < count; ++i)
      printf("I am %d of %d in comm: arrayOut[%d] = %d\n", myRank, numProcs, i, arrayOut[i]);

   delete [] arrayIn;
   delete [] arrayOut;
}
void demo_allreduce (int count, MPI_Comm comm)
{
   int numProcs, myRank; // Withint this communicator
   MPI_Comm_size(comm,&numProcs);
   MPI_Comm_rank(comm,&myRank);

   int *arrayIn = new int[count];
   int *arrayOut = new int[count];
   for (int i = 0; i < count; ++i)
      arrayIn[i] = arrayOut[i] = -(1 + myRank);

   MPI_Allreduce (arrayIn, arrayOut, count, MPI_INT, MPI_SUM, comm);

   for (int i = 0; i < count; ++i)
      printf("I am %d of %d in comm: arrayOut[%d] = %d\n", myRank, numProcs, i, arrayOut[i]);

   delete [] arrayIn;
   delete [] arrayOut;
}
void demo_scan (int count, MPI_Comm comm)
{
   int numProcs, myRank; // Withint this communicator
   MPI_Comm_size(comm,&numProcs);
   MPI_Comm_rank(comm,&myRank);

   int *arrayIn = new int[count];
   int *arrayOut = new int[count];
   for (int i = 0; i < count; ++i)
      arrayIn[i] = arrayOut[i] = -(1 + myRank);

   MPI_Scan (arrayIn, arrayOut, count, MPI_INT, MPI_SUM, comm);

   for (int i = 0; i < count; ++i)
      printf("I am %d of %d in comm: arrayOut[%d] = %d\n", myRank, numProcs, i, arrayOut[i]);

   delete [] arrayIn;
   delete [] arrayOut;
}
void demo_gather (int count, MPI_Comm comm)
{
   int numProcs, myRank; // Withint this communicator
   MPI_Comm_size(comm,&numProcs);
   MPI_Comm_rank(comm,&myRank);

   const int root = 0;

   int *arrayIn = new int[count];
   int *arrayOut = new int[count* numProcs];

   for (int i = 0; i < count; ++i)
      arrayIn[i] = i + myRank * count;

   for (int n = 0; n < numProcs; ++n)
   {
      for (int i = 0; i < count; ++i)
         arrayOut[i + n*count] = -1;

      if (n == myRank)
      {
         printf("%d: ", myRank);
         for (int i = 0; i < count; ++i)
            printf("%d ", arrayIn[i]);
         printf("\n");
      }
      MPI_Barrier (comm);
   }

   MPI_Gather (arrayIn , count, MPI_INT,
               arrayOut, count, MPI_INT,
               root, comm);

   for (int n = 0; n < numProcs; ++n)
   {
      if (n == myRank)
      {
         printf("ArrayOut[%d]: ", myRank);
         for (int i = 0; i < numProcs*count; ++i)
            printf("%d ", arrayOut[i]);
         printf("\n");
      }
      MPI_Barrier (comm);
   }

   delete [] arrayIn;
   delete [] arrayOut;
}
void demo_scatter (int count, MPI_Comm comm)
{
   int numProcs, myRank; // Withint this communicator
   MPI_Comm_size(comm,&numProcs);
   MPI_Comm_rank(comm,&myRank);

   const int root = 0;

   int *arrayIn = new int[count* numProcs];
   int *arrayOut = new int[count];

   if (myRank == 0)
   {
      for (int n = 0; n < numProcs; ++n)
         for (int i = 0; i < count; ++i)
            arrayIn[i + n*count] = i + n*count;

      printf("ArrayIn[%d]: ", myRank);
      for (int i = 0; i < numProcs * count; ++i)
         printf("%d ", arrayIn[i]);
      printf("\n");
   }
   else
   {
      for (int n = 0; n < numProcs; ++n)
         for (int i = 0; i < count; ++i)
            arrayIn[i + n*count] = -1;
   }

   MPI_Scatter (arrayIn , count, MPI_INT,
                arrayOut, count, MPI_INT,
                root, comm);

   for (int n = 0; n < numProcs; ++n)
   {
      if (n == myRank)
      {
         printf("ArrayOut[%d]: ", myRank);
         for (int i = 0; i < count; ++i)
            printf("%d ", arrayOut[i]);
         printf("\n");
      }
      MPI_Barrier (comm);
   }

   delete [] arrayIn;
   delete [] arrayOut;
}

int main (int argc, char* argv[])
{
   callMPI( MPI_Init(&argc,&argv) );

   int myRank, numProcs;

   callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
   callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );

   const char *testname[] = {
      "demo_broadcast",
      "demo_reduce",
      "demo_allreduce",
      "demo_scan",
      "demo_gather",
      "demo_scatter" };

   typedef void (*func_ptr_t)(int count, MPI_Comm comm);

   func_ptr_t testfptr[] = {
      demo_broadcast,
      demo_reduce,
      demo_allreduce,
      demo_scan,
      demo_gather,
      demo_scatter };

   const int numTests = sizeof(testfptr)/sizeof(func_ptr_t);

   int len = 1; // Length of the message to send/recv.
   int test = 2; // Which test?

   for (int i = 1; i < argc; ++i)
   {
#define check_index(i,str) \
   if ((i) >= argc) \
      { fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }

      if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
      {
         if (myRank == 0)
         {
            fprintf(stderr,"mpi_coll_examples --help|-h --length|-n --test|-t\n");
            fprintf(stderr,"available tests: \n");
            for (int k = 0; k < numTests; ++k)
               fprintf(stderr,"\t%-2d: %s\n", k, testname[k]);
         }

         callMPI( MPI_Finalize() );

         return 0;
      }
      else if (strcmp(argv[i],"--length") == 0 || strcmp(argv[i],"-n") == 0)
      {
         check_index(i+1,"--length|-n");
         i++;
         if (isdigit(*argv[i]))
            len = atoi( argv[i] );
      }
      else if (strcmp(argv[i],"--test") == 0 || strcmp(argv[i],"-t") == 0)
      {
         check_index(i+1,"--test|-m");
         i++;
         if (isdigit(*argv[i]))
            test = atoi( argv[i] );
      }
   }

   if (test >= 0 && test < numTests)
   {
      MPI_Comm comm = MPI_COMM_NULL;
      MPI_Comm_dup (MPI_COMM_WORLD, &comm);

      if (myRank == 0)
         printf("Calling %s with %d processes\n", testname[test], numProcs);

      testfptr[test]( len, comm );
   }
   else
   {
      if (myRank == 0)
         fprintf(stderr,"Test (%d) is out of range [0,%d)\n", test, numTests);
   }

   callMPI( MPI_Finalize() );

   return 0;
}
