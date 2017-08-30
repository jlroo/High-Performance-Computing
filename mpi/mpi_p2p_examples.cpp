#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <my_mpi_header.h>
 
/************************************************************
 * This is a simple send/receive program in MPI
 * ************************************************************/

void always_deadlocks (int count, int myRank, int numProcs)
{
   // This will always deadlock since the Recv's are blocking and the
   // Send's will never be posted.
   if (myRank <= 1)
   {
      if (myRank == 0)
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -1;
         int *toRecv = new int[count];

         MPI_Status status;
         MPI_Recv (toRecv, count, MPI_INT, 1, 2, MPI_COMM_WORLD, &status);
         printf("process %d recv %d\n", myRank, toRecv[0]);

         MPI_Send (toSend, count, MPI_INT, 1, 1, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         delete [] toSend;
         delete [] toRecv;
      }
      else
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -2;
         int *toRecv = new int[count];

         MPI_Status status;
         MPI_Recv (toRecv, count, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
         printf("process %d got %d\n", myRank, toRecv[0]);

         MPI_Send (toSend, count, MPI_INT, 0, 2, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         delete [] toSend;
         delete [] toRecv;
      }
   }
}

void sometimes_deadlocks (int count, int myRank, int numProcs)
{
   // This sometimes will deadlock if the Send too larger
   // to fit into a system buffer. Otherwise, small messages will
   // be buffered and Send will return.
   if (myRank <= 1)
   {
      if (myRank == 0)
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -1;
         int *toRecv = new int[count];

         MPI_Send (toSend, count, MPI_INT, 1, 1, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         MPI_Status status;
         MPI_Recv (toRecv, count, MPI_INT, 1, 2, MPI_COMM_WORLD, &status);
         printf("process %d recv %d\n", myRank, toRecv[0]);

         delete [] toSend;
         delete [] toRecv;
      }
      else
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -2;
         int *toRecv = new int[count];

         MPI_Send (toSend, count, MPI_INT, 0, 2, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         MPI_Status status;
         MPI_Recv (toRecv, count, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
         printf("process %d got %d\n", myRank, toRecv[0]);

         delete [] toSend;
         delete [] toRecv;
      }
   }
}

void never_deadlocks (int count, int myRank, int numProcs)
{
   // this never deadlocks since the Send on 0 and the Recv on 1
   // are matched (ordered). But, this is hard to do when multiple
   // message must be sent to different neighbors.
   if (myRank <= 1)
   {
      if (myRank == 0)
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -1;
         int *toRecv = new int[count];

         MPI_Send (toSend, count, MPI_INT, 1, 1, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         MPI_Status status;
         MPI_Recv (toRecv, count, MPI_INT, 1, 2, MPI_COMM_WORLD, &status);
         printf("process %d recv %d\n", myRank, toRecv[0]);

         delete [] toSend;
         delete [] toRecv;
      }
      else
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -2;
         int *toRecv = new int[count];

         MPI_Status status;
         MPI_Recv (toRecv, count, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
         printf("process %d got %d\n", myRank, toRecv[0]);

         MPI_Send (toSend, count, MPI_INT, 0, 2, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         delete [] toSend;
         delete [] toRecv;
      }
   }
}
void sendrecv_never_deadlocks (int count, int myRank, int numProcs)
{
   // This never deadlocks since the SendRecv function handles this
   // internally. Good choice if you have send and recv messages
   // regularly from a few processes.
   if (myRank <= 1)
   {
      if (myRank == 0)
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -1;
         int *toRecv = new int[count];

         MPI_Status status;
         MPI_Sendrecv (toSend, count, MPI_INT, 1, 1,
                       toRecv, count, MPI_INT, 1, 2, MPI_COMM_WORLD, &status);

         printf("process %d sent %d, recv %d\n", myRank, toSend[0], toRecv[0]);

         delete [] toSend;
         delete [] toRecv;
      }
      else
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -2;
         int *toRecv = new int[count];

         MPI_Status status;
         MPI_Sendrecv (toSend, count, MPI_INT, 0, 2,
                       toRecv, count, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

         printf("process %d sent %d, recv %d\n", myRank, toSend[0], toRecv[0]);

         delete [] toSend;
         delete [] toRecv;
      }
   }
}
void irecv_never_deadlocks (int count, int myRank, int numProcs)
{
   // This will never deadlock since it uses non-blocking Recv. The
   // Irecv function returns immediately so the Send can proceed. But,
   // we have to test when the Irecv is finished with a Wait.
   if (myRank <= 1)
   {
      if (myRank == 0)
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -1;
         int *toRecv = new int[count];

         MPI_Request request;
         MPI_Irecv (toRecv, count, MPI_INT, 1, 2, MPI_COMM_WORLD, &request);
         printf("process %d posted recv\n", myRank);

         MPI_Send (toSend, count, MPI_INT, 1, 1, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         MPI_Status status;
         MPI_Wait (&request, &status);
         printf("process %d waited and finished recv %d\n", myRank, toRecv[0]);

         delete [] toSend;
         delete [] toRecv;
      }
      else
      {
         int *toSend = new int[count];
         for (int i = 0; i < count; ++i)
            toSend[i] = -2;
         int *toRecv = new int[count];

         MPI_Request request;
         MPI_Irecv (toRecv, count, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
         printf("process %d posted recv\n", myRank);

         MPI_Send (toSend, count, MPI_INT, 0, 2, MPI_COMM_WORLD);
         printf("process %d sent %d\n", myRank, toSend[0]);

         MPI_Status status;
         MPI_Wait (&request, &status);
         printf("process %d waited and finished recv %d\n", myRank, toRecv[0]);

         delete [] toSend;
         delete [] toRecv;
      }
   }
}
void nonblocking_never_deadlocks (int count, int myRank, int numProcs)
{
   // Same as above but a simplier code. This is more common in practice.
   if (myRank <= 1)
   {
      int tag = 1000;

      int *toSend = new int[count];
      for (int i = 0; i < count; ++i)
         toSend[i] = -(1 + myRank);
      int *toRecv = new int[count];

      int dest = 1 - myRank;

      MPI_Request request;
      MPI_Irecv (toRecv, count, MPI_INT, dest, tag, MPI_COMM_WORLD, &request);
      printf("process %d posted recv from %d\n", myRank, dest);

      MPI_Send (toSend, count, MPI_INT, dest, tag, MPI_COMM_WORLD);
      printf("process %d sent %d to %d\n", myRank, toSend[0], dest);

      MPI_Status status;
      MPI_Wait (&request, &status);
      printf("process %d waited and finished recv %d\n", myRank, toRecv[0]);

      delete [] toSend;
      delete [] toRecv;
   }
}

void nonblocking2_never_deadlocks (int count, int myRank, int numProcs)
{
   // More advanced non-blocking. Uses non-blocking Send and non-blocking Recv
   // and Wait's for both. I staggered the delete's to show when it is safe
   // to reuse (or delete) the message data.
   if (myRank <= 1)
   {
      int tag = 1000;

      int *toSend = new int[count];
      for (int i = 0; i < count; ++i)
         toSend[i] = -(1 + myRank);
      int *toRecv = new int[count];

      int dest = 1 - myRank;

      MPI_Request recv_req;
      MPI_Irecv (toRecv, count, MPI_INT, dest, tag, MPI_COMM_WORLD, &recv_req);
      printf("process %d posted recv from %d\n", myRank, dest);

      MPI_Request send_req;
      MPI_Isend (toSend, count, MPI_INT, dest, tag, MPI_COMM_WORLD, &send_req);
      printf("process %d posted send %d to %d\n", myRank, toSend[0], dest);

      MPI_Status status;
      MPI_Wait (&recv_req, &status);
      printf("process %d waited and finished recv %d\n", myRank, toRecv[0]);

      // Recv is done so it's safe to delete this.
      delete [] toRecv;

      MPI_Wait (&send_req, &status);
      printf("process %d waited and finished send\n", myRank);

      // Send is done so it's safe to delete this.
      delete [] toSend;
   }
}

int main (int argc, char* argv[])
{
   callMPI( MPI_Init(&argc,&argv) );

   int myRank, numProcs;

   callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
   callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );

   const char *testname[] = { "always_deadlocks",
                              "sometimes_deadlocks",
                              "never_deadlocks",
                              "sendrecv_never_deadlocks",
                              "irecv_never_deadlocks",
                              "nonblocking_never_deadlocks",
                              "nonblocking2_never_deadlocks" };

   typedef void (*func_ptr_t)(int count, int myRank, int numProcs);

   func_ptr_t testfptr[] = {
                        always_deadlocks,
                        sometimes_deadlocks,
                        never_deadlocks,
                        sendrecv_never_deadlocks,
                        irecv_never_deadlocks,
                        nonblocking_never_deadlocks,
                        nonblocking2_never_deadlocks };

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
            fprintf(stderr,"mpi_p2p_examples --help|-h --length|-n --test|-t\n");
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
      printf("I am %d of %d processes; len = %d; test = %s.\n", myRank, numProcs, len, testname[test]);
      testfptr[test]( len, myRank, numProcs );
   }
   else
   {
      if (myRank == 0)
         fprintf(stderr,"Test (%d) is out of range [0,%d)\n", test, numTests);
   }

   callMPI( MPI_Finalize() );

   return 0;
}
