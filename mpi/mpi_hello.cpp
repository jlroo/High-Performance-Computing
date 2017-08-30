#include <stdio.h>
#include <stdlib.h>

#include <my_mpi_header.h>

/************************************************************
 * This is a simple hello world program. Each processor prints out 
 * it's rank and the size of the current MPI run (Total number of
 * processors).
 * ************************************************************/
int main (int argc, char* argv[])
{
   // Initialize the MPI system before anything else.
   callMPI( MPI_Init(&argc,&argv)                   );

   // Get the rank and size of the MPI world.
   int myRank, numProcs;
   callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
   callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );

   // Get the name of the local host.
   char mpi_name[MPI_MAX_PROCESSOR_NAME];
   int mpi_name_len;
   callMPI( MPI_Get_processor_name( mpi_name, &mpi_name_len ) );

   // Loop over all of the ranks and print the NUMA hardware.
   for (int rank = 0; rank < numProcs; rank++)
   {
      if (myRank == rank)
      {
         /* print out my rank and this run's PE size*/
         printf("Hello from %d of %d on %s\n", myRank, numProcs, mpi_name);

         int iret = system("numactl --show | grep physcpu");
      }

      // Block the rest of the ranks for progressing so the stdout is managable.
      callMPI( MPI_Barrier(MPI_COMM_WORLD) );
   }

   // Close the MPI library.
   callMPI( MPI_Finalize() );

   return 0;
}
