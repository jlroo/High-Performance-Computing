#include <mpi.h>

#include <stdio.h>
#include <string.h>
#include <string>

const int TAG_TASK = 1;
const int TAG_RESULT = 2;
const int TAG_TERMINATE = 3;

/* This is the manager */
int manager (MPI_Comm manager_comm, MPI_Comm comm)
{
   int  rank,size;
   char buf[256];

   MPI_Comm_rank ( comm, &rank );
   MPI_Comm_size ( comm, &size );
   printf("Hello from manager in split comm %d %d\n", rank, size);

   MPI_Comm_rank ( manager_comm, &rank );
   MPI_Comm_size ( manager_comm, &size );
   printf("Hello from manager in manager_comm %d %d\n", rank, size);

   int ntasks = 5*size;
   int ntasks_remaining = ntasks;
   int nworkers_active = size-1;

   while (ntasks_remaining > 0 or nworkers_active > 0)
   {
      MPI_Status status;

      // Listen if a message is pending.
      MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, manager_comm, &status);

      int worker_rank = status.MPI_SOURCE;
      int tag         = status.MPI_TAG;

      // What did the work send?
      if (tag == TAG_TASK)
      {
         // Worker is requesting a task. Accept the message but ignore it.
         MPI_Recv (buf, 256, MPI_CHAR, worker_rank, TAG_TASK, manager_comm, &status );
         printf("Worker %d sent task request %s\n", worker_rank, buf);

         // Send the next task or the terminate signal.
         if (ntasks_remaining > 0)
         {
            MPI_Send (&ntasks_remaining, 1, MPI_INT, worker_rank /*dest*/, TAG_TASK, manager_comm);
            ntasks_remaining--;
         }
         else
         {
            MPI_Send (&ntasks_remaining, 1, MPI_INT, worker_rank /*dest*/, TAG_TERMINATE, manager_comm);
            nworkers_active--;
         }
      }
      else if (tag == TAG_RESULT)
      {
         // Worker is sending a result.
         // 1. Figure out how much data they are sending me.
         int count = 0;
         MPI_Get_count ( &status, MPI_CHAR, &count );

         MPI_Recv (buf, count, MPI_CHAR, worker_rank, TAG_RESULT, manager_comm, &status );
         printf("Worker %d sent results count=%d: %s\n", worker_rank, count, buf);
      }
   }

   return 0;
}

std::string toBinary (int decimal)
{
   std::string binary;
   int i = decimal;
   while (i > 0)
   {
      int digit = i % 2;
      binary = ((digit == 0) ? "0" : "1") + binary;
      i /= 2;
   }

   return binary;
}

std::string toOctal (int decimal)
{
   std::string binary;
   int i = decimal;
   while (i > 0)
   {
      int digit = i % 8;
      char odig[2];
      sprintf(odig,"%d", digit);
      binary = odig[0] + binary;
      i /= 8;
   }

   return binary;
}

std::string toHexademical (int decimal)
{
   std::string binary;
   int i = decimal;
   while (i > 0)
   {
      int digit = i % 16;
      if (digit < 10)
      {
         char odig[2];
         sprintf(odig,"%d", digit);
         binary = odig[0] + binary;
      }
      else if (digit == 10) binary = "A" + binary;
      else if (digit == 11) binary = "B" + binary;
      else if (digit == 12) binary = "C" + binary;
      else if (digit == 13) binary = "D" + binary;
      else if (digit == 14) binary = "E" + binary;
      else if (digit == 15) binary = "F" + binary;

      i /= 16;
   }

   return binary;
}

/* This is the worker */
int worker (MPI_Comm manager_comm, MPI_Comm worker_comm)
{
   char buf[256];
   int  worker_rank;

   MPI_Comm_rank( worker_comm, &worker_rank );

   printf("Hello from worker %d.\n", worker_rank);

   int ntasks_processed = 0;

   while(true)
   {
      // Tell the manager I'm available to work.
      sprintf( buf, "Worker %d is ready ...", worker_rank );
      MPI_Send( buf, strlen(buf) + 1, MPI_CHAR, 0, TAG_TASK, manager_comm );

      // Receive a task command from the manager.
      MPI_Status status;
      int task;
      //printf("Worker %d is waiting for command ...\n", worker_rank);
      MPI_Recv ( &task, 1, MPI_INT, 0, MPI_ANY_TAG, manager_comm, &status );

      int tag = status.MPI_TAG;

      if (tag == TAG_TERMINATE)
      {
         break;
      }
      else
      {
         std::string result = toBinary(task);
         //std::string result = toOctal(task);
         //std::string result = toHexademical(task);
         ntasks_processed++;
         sprintf( buf, "Worker %d processed task %d %s %d", worker_rank, task, result.c_str(), ntasks_processed );
         MPI_Send( buf, strlen(buf) + 1, MPI_CHAR, 0, TAG_RESULT, manager_comm );
      }
   }

   printf("Worker %d completed %d tasks and ended.\n", worker_rank, ntasks_processed );

   return 0;
}

int main (int argc, char* argv[])
{
   int rank, size;
   MPI_Comm new_comm;

   MPI_Init ( &argc, &argv );
   MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
   MPI_Comm_size ( MPI_COMM_WORLD, &size );

   if (size > 1)
   {
      // Create a new communicator that splits the manager and the workers into unique groups.
      MPI_Comm_split ( MPI_COMM_WORLD, (rank == 0) /*color*/, rank/*rank key*/, &new_comm );

      if (rank == 0)
         manager ( MPI_COMM_WORLD, new_comm );
      else
         worker ( MPI_COMM_WORLD, new_comm );
   }
   else
   {
      if (rank == 0)
         fprintf(stderr,"Must use > 1 processes\n");
   }

   MPI_Finalize( );
   return 0;
}
