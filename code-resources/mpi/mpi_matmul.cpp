#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <cmath>
#include <algorithm>
#include <vector> // std::vector container.

#include <my_timer.h>
#include <my_mpi_header.h>

#ifndef RESTRICT
#define RESTRICT
#endif

typedef double ValueType;
const MPI_Datatype MPI_ValueType = MPI_DOUBLE;

// BLAS library call.
extern "C" void dgemm_ (char* transa, char* transb, int* _m, int* _n, int* _k, double* _alpha, double A[], int* _lda, double B[], int* _ldb, double* _beta, double C[], int* _ldc);
extern "C" void sgemm_ (char* transa, char* transb, int* _m, int* _n, int* _k, float* _alpha, float A[], int* _lda, float B[], int* _ldb, float* _beta, float C[], int* _ldc);

void matmul_blas (const int m, const int n, const int k, const double alpha, double A[], const int lda, double B[], const int ldb, const double beta, double C[], const int ldc)
{
   char transa[] = {'N'};
   char transb[] = {'N'};

   dgemm_ (transa, transb, const_cast<int*>(&m), const_cast<int*>(&n), const_cast<int*>(&k), const_cast<double*>(&alpha), A, const_cast<int*>(&lda), B, const_cast<int*>(&ldb), const_cast<double*>(&beta), C, const_cast<int*>(&ldc));
}

typedef void (*matmul_ptr) (const int m, const int n, const int k, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc);

void matmul_vec (const int m, const int n, const int k, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc)
{
   // Restructed to improved vectorization.
   for (int j = 0; j < n; ++j)
   {
      if (beta == 0.0)
         for (int i = 0; i < m; ++i) C[i + j*ldc] = 0.0;
      else if (beta != 1.0)
         for (int i = 0; i < m; ++i) C[i + j*ldc] *= beta;

      for (int l = 0; l < k; ++l)
      {
         //if (B[l + j*ldb] != 0.0)
         {
            ValueType temp = alpha * B[l + j*ldb];
            for (int i = 0; i < m; ++i)
               C[i + j*ldc] += temp * A[i + l*lda];
         }
      }
   }
}
void matmul_block (const int m, const int n, const int k, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc)
{
#ifndef blockSize
#define blockSize (32)
#warning 'Setting blockSize = 32'
#endif

   for (int i = 0; i < m; i += blockSize)
      for (int j = 0; j < n; j += blockSize)
      {
         const int m_blk = std::min( blockSize, m-i);
         const int n_blk = std::min( blockSize, n-j);

         ValueType *C_blk = C + i + j*ldc;

         if (beta == 0.0)
            for (int jj = 0; jj < n_blk; ++jj)
               for (int ii = 0; ii < m_blk; ++ii)
                  C_blk[ii + jj*ldc] = 0.0;
         else if (beta != 1.0)
            for (int jj = 0; jj < n_blk; ++jj)
               for (int ii = 0; ii < m_blk; ++ii)
                  C_blk[ii + jj*ldc] *= beta;

         for (int l = 0; l < k; l += blockSize)
         {
            const int k_blk = std::min( blockSize, k-l);

            ValueType *A_blk = A + i + l*lda;
            ValueType *B_blk = B + l + j*ldb;

            const ValueType one(1.0);
            matmul_vec (m_blk, n_blk, k_blk, alpha, A_blk, lda, B_blk, ldb, one, C_blk, ldc);
         }
      }
}

int partition_range (const int global_start, const int global_end,
                     const int num_partitions, const int rank,
                     int& local_start, int& local_end)
{
   // Total length of the iteration space.
   const int global_length = global_end - global_start;

   // Simple per-partition size ignoring remainder.
   const int chunk_size = global_length / num_partitions;

   // And now the remainder.
   const int remainder = global_length - chunk_size * num_partitions;

   // We want to spreader the remainder around evening to the 1st few ranks.
   // ... add one to the simple chunk size for all ranks < remainder.
   if (rank < remainder)
   {
      local_start = global_start + rank * chunk_size + rank;
      local_end   = local_start + chunk_size + 1;
   }
   else
   {
      local_start = global_start + rank * chunk_size + remainder;
      local_end   = local_start + chunk_size;
   }

   return 0;
}

struct mpiMatMulType
{
   // [Begin,End) indices for the block partitions.
   std::vector<int> blockBegin, blockEnd;

   // # of row/column partitions (assuming equal).
   int numBlocks;

   // Local row/column block rank.
   int rowRank, colRank;

   // Global communicator.
   MPI_Comm allComm;
};

int matmul_block_mpi (const int M, const int N, const int K, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc, const mpiMatMulType &mpiData)
{
   int myRank;

   callMPI( MPI_Comm_rank( mpiData.allComm, &myRank)   );

   const int myRow = mpiData.rowRank;
   const int myCol = mpiData.colRank;

   const int myRowBegin = mpiData.blockBegin[myRow];
   const int myRowEnd   = mpiData.blockEnd[myRow];
   const int myColBegin = mpiData.blockBegin[myCol];
   const int myColEnd   = mpiData.blockEnd[myCol];

   const int myRowSize  = myRowEnd - myRowBegin;
   const int myColSize  = myColEnd - myColBegin;

   // We only communicate with other processes on the same block row and column.
   // Create a row and column communicator.
   MPI_Comm rowComm = MPI_COMM_NULL;
   MPI_Comm colComm = MPI_COMM_NULL;

   callMPI( MPI_Comm_split (mpiData.allComm, myRow, myRank, &rowComm) );
   callMPI( MPI_Comm_split (mpiData.allComm, myCol, myRank, &colComm) );

   /*{
      int col_myRank, col_numProcs;
      int row_myRank, row_numProcs;

      callMPI( MPI_Comm_size( colComm, &col_numProcs) );
      callMPI( MPI_Comm_rank( colComm, &col_myRank)   );
      callMPI( MPI_Comm_size( rowComm, &row_numProcs) );
      callMPI( MPI_Comm_rank( rowComm, &row_myRank)   );

      printf("myRank= %d, %d myRowRank= %d %d %d %d %d %d\n", myRank, mpiData.numBlocks, myRow, myCol, row_numProcs, row_myRank, col_numProcs, col_myRank);

      callMPI( MPI_Barrier( mpiData.allComm ) );
   }*/

   // Initialize my block output.
   if (beta == 0.0)
      for (int j = 0; j < myColSize; ++j)
         for (int i = 0; i < myRowSize; ++i)
            C[i + j*ldc] = 0.0;
   else if (beta != 1.0)
      for (int j = 0; j < myColSize; ++j)
         for (int i = 0; i < myRowSize; ++i)
            C[i + j*ldc] *= beta;

   for (int kb = 0; kb < mpiData.numBlocks; ++kb)
   {
      // How big are the Aik and Bjk blocks?
      // # rows in Aik (m) match Cij and # of columns in Bjk (n) match Cij.
      // The # of columns (k) of Aik and Bkj must match.
      const int m = myRowSize;
      const int n = myColSize;
      const int k = mpiData.blockEnd[kb] - mpiData.blockBegin[kb];

      ValueType *Aik = new ValueType[m*k];
      ValueType *Bkj = new ValueType[k*n];

      // Copy my data into the buffer if I own the block.
      if (kb == myCol)
         for (int j = 0; j < myColSize; ++j)
            for (int i = 0; i < myRowSize; ++i)
               Aik[i + j*m] = A[i + j*lda];

      // Broadcast to my row neighbors || receive Aik.
      callMPI( MPI_Bcast( Aik, m*k, MPI_ValueType, kb, rowComm ) );

      if (kb == myRow)
         for (int j = 0; j < myColSize; ++j)
            for (int i = 0; i < myRowSize; ++i)
               Bkj[i + j*k] = B[i + j*ldb];

      callMPI( MPI_Bcast( Bkj, k*n, MPI_ValueType, kb, colComm ) );

      //matmul_block (m, n, k, alpha, Aik, m, Bkj, k, 1.0, C, ldc);
      matmul_blas (m, n, k, alpha, Aik, m, Bkj, k, 1.0, C, ldc);

      delete [] Aik;
      delete [] Bkj;
   }

   callMPI( MPI_Comm_free( &rowComm ) );
   callMPI( MPI_Comm_free( &colComm ) );

   return 0;
}

void setMatrixElement (const int i, const int j, const int n,
                       ValueType *a, ValueType *b )
{
   *a = 0.001 * i + j;
   *b = 0.1;
}
int isPowerOf2(const int x) {
   /* get all 1's in mask. Negate x so that if its a pow of 2, it has only 1 zero.
    *      XOR against all 1's mask, since a non-pow-2 would have a diff bit pattern and not
    *           come out non-zero */

   // Looking for this 100000
   // Deduct 1 to get  011111
   // And'ing the two should be zero if power of 2.
   const int minusOne = x - 1;
   const int AndValue = x & minusOne;
   return (AndValue == 0);
}

int run_matmul (const int N, MPI_Comm comm, const bool computeError)
{
   int myRank, numProcs;

   callMPI( MPI_Comm_size( comm, &numProcs) );
   callMPI( MPI_Comm_rank( comm, &myRank)   );

   // # of ranks along one direction of the matrix.
   // Split the matrix into P x P blocks.
   // numProcs is a power-of-2!!!
   //if (isPowerOf2(numProcs) == false)
   //{
   //   fprintf(stderr,"Must be a power-of-2 %d\n", numProcs);
   //   return 1;
   //}
   const int P = (int) sqrtf( float(numProcs) + 1e-10);
   if (P*P != numProcs)
   {
      fprintf(stderr,"P is not a perfect square %d %d\n", P, numProcs);
      return 1;
   }

   mpiMatMulType mpiData;

   mpiData.numBlocks  = P;
   mpiData.allComm    = comm;
   mpiData.blockBegin.resize(P,0);
   mpiData.blockEnd  .resize(P,0);

   // What is this processor's coordinates in the 2d matrix?
   mpiData.rowRank = myRank % P;
   mpiData.colRank = myRank / P;

   for (int ib = 0; ib < P; ++ib)
      partition_range( 0, N, P, ib, mpiData.blockBegin[ib], mpiData.blockEnd[ib] );

   const int myRow = mpiData.rowRank;
   const int myCol = mpiData.colRank;

   const int myRowBegin = mpiData.blockBegin[myRow];
   const int myRowEnd   = mpiData.blockEnd  [myRow];
   const int myColBegin = mpiData.blockBegin[myCol];
   const int myColEnd   = mpiData.blockEnd  [myCol];
   const int myRowSize  = myRowEnd - myRowBegin;
   const int myColSize  = myColEnd - myColBegin;

   if (computeError)
   {
      if (myRank == 0)
         for (int p = 0; p < numProcs; ++p)
         {
            const int ib = p % P;
            const int jb = p / P;
            printf("proc[%d] block %d,%d row= [%d,%d) column= [%d,%d)\n", p, ib, jb, mpiData.blockBegin[ib], mpiData.blockEnd[ib], mpiData.blockBegin[jb], mpiData.blockEnd[jb]);
         }
      //printf("myRank= %d, P= %d, myRow/Col= %d %d, myRowBegin= %d %d %d %d %d %d\n", myRank, P, myRow, myCol, myRowBegin, myRowEnd, myColBegin, myColEnd, myRowSize, myColSize);
      callMPI( MPI_Barrier(mpiData.allComm) );
   }

   // Allocate local block matrices.
   ValueType *Aij = new ValueType[myRowSize*myColSize];
   ValueType *Bij = new ValueType[myRowSize*myColSize];
   ValueType *Cij = new ValueType[myRowSize*myColSize];

   // Don't see with random #'s this time. Build a predictible matrix in parallel.
   for (int j = 0; j < myColSize; ++j)
      for (int i = 0; i < myRowSize; ++i)
      {
         const int i_global = i + myRowBegin;
         const int j_global = j + myColBegin;
         setMatrixElement( i_global, j_global, N, &Aij[i + j*myRowSize], &Bij[i + j*myRowSize] );
         Cij[i + j*myRowSize] = 0;
      }

   // Interesting scaling factors (not 0 or 1).
   ValueType alpha = 0.1, beta = 0.21;

   int niters = 1;

   const bool doReference = (N < 500) and computeError;

   // Run a few iterations to warm up the system.
   for (int iter = 0; iter < std::min(2,niters); iter++)
      matmul_block_mpi(N, N, N, alpha, Aij, myRowSize, Bij, myRowSize, beta, Cij, myRowSize, mpiData);

   // Run the test for 'a long time.'
   double tCalc = 0;
   for(;;)
   {
      myTimer_t t0 = getTimeStamp();

      for (int iter = 0; iter < niters; iter++)
      {
         if (doReference)
            for (int j = 0; j < myColSize; ++j)
               for (int i = 0; i < myRowSize; ++i)
                  Cij[i + j*myRowSize] = 0.0;

         matmul_block_mpi(N, N, N, alpha, Aij, myRowSize, Bij, myRowSize, beta, Cij, myRowSize, mpiData);
      }

      myTimer_t t1 = getTimeStamp();

      // Trick the compiler to keep the results.
      srand((int)Cij[1]);

      //tCalc += getElapsedTime(t0,t1);
      tCalc = getElapsedTime(t0,t1);

      double tCalcMax;
      callMPI( MPI_Allreduce( &tCalc, &tCalcMax, 1, MPI_DOUBLE, MPI_MAX, comm) );
      tCalc = tCalcMax;

      if (tCalc < 0.01 and computeError == false)
         niters *= 2;
      else
         break;
   }
   tCalc /= niters;

   bool passed = true;
   double error = 0.0;
   if (doReference)
   {
      // Allocate arrays with a little padding.
      ValueType *Aref = new ValueType[N*N];
      ValueType *Bref = new ValueType[N*N];
      ValueType *Cref = new ValueType[N*N];

      for (int j = 0; j < N; ++j)
         for (int i = 0; i < N; ++i)
         {
            setMatrixElement( i, j, N, &Aref[i + j*N], &Bref[i + j*N] );
            Cref[i + j*N] = 0;
         }

      matmul_blas(N, N, N, alpha, Aref, N, Bref, N, beta, Cref, N);

      // L2-norm values.
      double err2 = 0.0, ref2 = 0.0;

      // The reference matrix is 'global'. Just test against my local chunk.
      for (int j = 0; j < myColSize; ++j)
         for (int i = 0; i < myRowSize; ++i)
         {
            const int i_global = i + myRowBegin;
            const int j_global = j + myColBegin;

            double c_global = Cref[i_global + j_global*N];
            double c_local = Cij[i + j*myRowSize];
            double diff = c_local - c_global;
            //printf("err: %d %d %e %e %e %d\n", i_global, j_global, c_global, c_local, diff, myRank);
            ref2 += c_global * c_global;
            err2 += diff*diff;
         }

      // Accumlate the L2-norms from each rank.
      ValueType sbuf[2] = { err2, ref2 };
      ValueType rbuf[2];
      callMPI( MPI_Allreduce( sbuf, rbuf, 2, MPI_ValueType, MPI_SUM, comm) );

      err2 = rbuf[0];
      ref2 = rbuf[1];

      error = sqrt( err2 / ref2 );
      passed = error < 1e-10;

      if (myRank == 0 and not(passed))
         fprintf(stderr,"Matrix solution error = %e %e %e\n", sqrt( err2 / ref2), err2, ref2 );

      delete [] Aref;
      delete [] Bref;
      delete [] Cref;
   }

   double Mflops = 1e-6 * ((((2.0*N)*N)*N) + ((3.0*N)*N)) / tCalc;
   if (myRank == 0)
      printf("%d, %g, %f, %d, %d, %e\n", N, tCalc*1000, Mflops, niters, passed, error);

   delete [] Aij;
   delete [] Bij;
   delete [] Cij;

   return not(passed);
}

int main (int argc, char * argv[])
{
   callMPI( MPI_Init(&argc,&argv) );

   int myRank, numProcs;

   callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
   callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );

   int min_size = 50; // Minimum matrix size.
   int max_size = 5000; // Maximum matrix size.
   double step_size = 2; // Size increase factor.
   bool computeError = false;

   for (int i = 1; i < argc; ++i)
   {
#define check_index(i,str) \
   if ((i) >= argc) \
      { fprintf(stderr,"Missing 2nd argument for %s\n", str); return 1; }

      if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
      {
         if (myRank == 0)
         {
            fprintf(stderr,"mpi_matmul --help|-h --min --max --step|-s --error|-e|--test|-t\n");
         }

         callMPI( MPI_Finalize() );

         return 0;
      }
      else if (strcmp(argv[i],"--min") == 0)
      {
         check_index(i+1,"--min");
         i++;
         if (isdigit(*argv[i]))
            min_size = atoi( argv[i] );

         min_size = std::max(min_size, 1);
      }
      else if (strcmp(argv[i],"--max") == 0)
      {
         check_index(i+1,"--max");
         i++;
         if (isdigit(*argv[i]))
            max_size = atoi( argv[i] );
      }
      else if (strcmp(argv[i],"--step") == 0 || strcmp(argv[i],"--step") == 0)
      {
         check_index(i+1,"--step");
         i++;
         if (isdigit(*argv[i]))
            step_size = atof( argv[i] );
      }
      else if (strcmp(argv[i],"--error") == 0 || strcmp(argv[i],"-e") == 0)
         computeError = true;
      else if (strcmp(argv[i],"--test") == 0 || strcmp(argv[i],"-t") == 0)
         computeError = true;
   }

   max_size = std::max( max_size, min_size );

   // Only pass in a perfect-square # of processes to do the matmul.
   int color = 0, key = myRank;
   for (int P = 1; P*P <= numProcs; ++P)
   {
      if (min_size < P)
      {
         fprintf(stderr,"Min matrix size is less than the sqrt(# of processes), shrinking ... %d %d.\n", min_size, P);
         --P;
         break;
      }

      color = (myRank < P*P);
   }
   //printf("myRank= %d, color= %d, key= %d\n", myRank, color, key);

   MPI_Comm comm = MPI_COMM_NULL;
   callMPI( MPI_Comm_split (MPI_COMM_WORLD, color, key, &comm) );

   if (myRank == 0)
   {
      int np;
      callMPI( MPI_Comm_size(comm, &np) );
      printf("mpi_matmul: min_size= %d max_sisze= %d step_size= %.2f with %d active processes (%d)\n", min_size, max_size, step_size, np, numProcs);
   }

   if (color == 1)
   {
      for (int size = min_size; size <= max_size; size *= step_size)
      {
         int ierr = run_matmul( size, comm, computeError );
         if (ierr != 0)
            break;
      }
   }

   callMPI( MPI_Comm_free(&comm) );
   callMPI( MPI_Finalize() );

   return 0;
}
