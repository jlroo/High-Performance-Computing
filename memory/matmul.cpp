#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <cmath>
#include <algorithm>

#include <my_timer.h>

#ifndef RESTRICT
#define RESTRICT
#endif

typedef double ValueType;

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

/* This routine performs a matmul operation using a straightforward
 * "naive", three-loop method.
 *    C := beta*C + alpha*A * B => C_i,j = Sum_l A_i,l * B_l,j
 *    where A (mxk), B (kxn), and C (mxn) are matrices stored in column-major format.
 *    On exit, A and B maintain their input values. */
void matmul_naive (const int m, const int n, const int k, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc)
{
   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j)
      {
         ValueType cij(0);
         for (int l = 0; l < k; ++l)
            cij += A[i + l*lda] * B[l + j*ldb];

         C[i + j*ldc] = beta*C[i + j*ldc] + alpha*cij;
      }
}
void matmul_opt1 (const int m, const int n, const int k, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc)
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
void matmul_opt3 (const int m, const int n, const int k, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc)
{
   ValueType *RESTRICT _A = A;
   ValueType *RESTRICT _B = B;
   ValueType *RESTRICT _C = C;

   // Restructed to improved vectorization.
   for (int j = 0; j < n; ++j)
   {
      if (beta == 0.0)
         for (int i = 0; i < m; ++i) _C[i + j*ldc] = 0.0;
      else if (beta != 1.0)
         for (int i = 0; i < m; ++i) _C[i + j*ldc] *= beta;

      /*for (int l = 0; l < k; l += 2)
      {
         ValueType t0 = alpha * _B[l   + j*ldb];
         ValueType t1 = alpha * _B[l+1 + j*ldb];
         //ValueType t2 = alpha * _B[l+2 + j*ldb];
         //ValueType t3 = alpha * _B[l+3 + j*ldb];
         for (int i = 0; i < m; ++i)
         {
            ValueType cij  = t0 * _A[i + (l+0)*lda];
                      cij += t1 * _A[i + (l+1)*lda];
            _C[i + j*ldc] += cij;
         }
      }*/

#ifndef unroll
#define unroll (8)
#warning 'Setting unroll = 8'
#endif
/*      //const int unroll(8);
      int l = 0;
      const int k_stop = k - unroll;
      for (; l < k_stop; l += unroll)
      {
         ValueType temp[unroll];
         for (int ll = 0; ll < unroll; ++ll)
            temp[ll] = alpha * _B[l + ll   + j*ldb];

         for (int i = 0; i < m; ++i)
         {
            ValueType cij = _C[i + j*ldc];
            for (int ll = 0; ll < unroll; ++ll)
               cij += temp[ll] * _A[i + (l+ll)*lda];

            _C[i + j*ldc] = cij;
         }
      } */


      #define __unroll (4)
      int l = 0;
      const int k_stop = k - __unroll;
      for (; l < k_stop; l += __unroll)
      {
         ValueType t0 = alpha * _B[l    + j*ldb];
#if (__unroll > 1)
         ValueType t1 = alpha * _B[l+1  + j*ldb];
#endif
#if (__unroll > 2)
         ValueType t2 = alpha * _B[l+2  + j*ldb];
#endif
#if (__unroll > 3)
         ValueType t3 = alpha * _B[l+3  + j*ldb];
#endif

         for (int i = 0; i < m; ++i)
         {
            ValueType cij = _C[i + j*ldc];
            cij += t0 * _A[i + (l  )*lda];
#if (__unroll > 1)
            cij += t1 * _A[i + (l+1)*lda];
#endif
#if (__unroll > 2)
            cij += t2 * _A[i + (l+2)*lda];
#endif
#if (__unroll > 3)
            cij += t3 * _A[i + (l+3)*lda];
#endif
            _C[i + j*ldc] = cij;
         }
      }

      for (; l < k; ++l)
      {
         //if (B[l + j*ldb] != 0.0)
         {
            ValueType temp = alpha * _B[l + j*ldb];
            for (int i = 0; i < m; ++i)
               _C[i + j*ldc] += temp * _A[i + l*lda];
         }
      }
      

   }
}
void matmul_opt2 (const int m, const int n, const int k, const ValueType alpha, ValueType A[], const int lda, ValueType B[], const int ldb, const ValueType beta, ValueType C[], const int ldc)
{
#ifndef blockSize
#define blockSize (32)
#warning 'Setting blockSize = 32'
#endif
   //const int blockSize = 41;

   /*if (beta == 0.0)
      for (int j = 0; j < n; ++j)
         for (int i = 0; i < m; ++i)
            C[i + j*ldc] = 0.0;
   else if (beta != 1.0)
      for (int j = 0; j < n; ++j)
         for (int i = 0; i < m; ++i)
            C[i + j*ldc] *= beta;*/

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
            //const int m_blk = std::min( blockSize, m-i);
            //const int n_blk = std::min( blockSize, n-j);
            const int k_blk = std::min( blockSize, k-l);

            ValueType *A_blk = A + i + l*lda;
            ValueType *B_blk = B + l + j*ldb;
            //ValueType *C_blk = C + i + j*ldc;

            const ValueType one(1.0);
            //matmul_naive (m_blk, n_blk, k_blk, alpha, A_blk, lda, B_blk, ldb, one, C_blk, ldc);
            matmul_opt1 (m_blk, n_blk, k_blk, alpha, A_blk, lda, B_blk, ldb, one, C_blk, ldc);
         }
      }
}

int run_matmul (int n, int niters, const double tDelta, matmul_ptr matmul)
{
   int npad = 0; // Array padding.

   // Allocate arrays with a little padding.
   ValueType *a = new ValueType[n*n+npad];
   ValueType *b = new ValueType[n*n+npad];
   ValueType *c = new ValueType[n*n+npad];
   ValueType *cref = new ValueType[n*n+npad];

   // Seed A[] with random #'s, zero C[], and set B[] = I.
   srand(n);
   ValueType invRandMax = ValueType(1) / ValueType(RAND_MAX);
   for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
      {
         a[i+j*n] = ValueType(rand()) * invRandMax;
         //b[i+j*n] = 0;
         //if (i==j) b[i+j*n] = 1;
         b[i+j*n] = 0.1;
         c[i+j*n] = 0;
      }

   // Interesting scaling factors (not 0 or 1).
   ValueType alpha = 0.1, beta = 0.21;

   matmul_blas(n, n, n, alpha, a, n, b, n, beta, cref, n);

   // Run a few iterations to warm up the system.
   for (int iter = 0; iter < std::min(5,niters); iter++)
      matmul(n, n, n, alpha, a, n, b, n, beta, c, n);

   // Run the test for 'a long time.'
   double tCalc = 0;
   for(;;)
   {
      myTimer_t t0 = getTimeStamp();

      for (int iter = 0; iter < niters; iter++)
      {
         for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
               c[i + j*n] = 0.0;

         matmul(n, n, n, alpha, a, n, b, n, beta, c, n);
      }

      myTimer_t t1 = getTimeStamp();

      // Trick the compiler to keep the results.
      srand((int)c[1]);

      //tCalc += getElapsedTime(t0,t1);
      tCalc = getElapsedTime(t0,t1);
      if (tCalc < 0.01)
         niters *= 2;
      else
         break;
   }
   tCalc /= niters;

   double err2 = 0.0, ref2 = 0.0;
   for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
      {
         double diff = c[i + j*n] - cref[i + j*n];
         ref2 += cref[i + j*n] * cref[i + j*n];
         err2 += diff*diff;
      }

   double Gflops = 1e-9 * ((((2.0*n)*n)*n) + ((3.0*n)*n)) / tCalc;
   printf("%5d, %10.4f, %10.4f, %.2f%%, %d %e\n", n, tCalc*1000, Gflops, 100*tDelta/(niters*tCalc), niters, sqrt(err2 / ref2));
   //printf("%d, %f, %.2f%%\n", n, tCalc/niters, tDelta/tCalc);

   if (n<=5)
   {
      for (int i = 0; i < n; ++i)
      {
         for (int j = 0; j < n; ++j)
            printf("%f ", c[i + j*n]);
         printf("\n");
      }
   }

   delete [] a;
   delete [] b;
   delete [] c;

   return 0;
}

int main (int argc, char * argv[])
{
   int min_size = 10;
   int max_size = 1000;
   int niters = 50; // Number of samples for each test.
   int method = 0; // blas
   double stepSize = 2;

   // Allow the user to change the defaults.
   if (argc > 1)
      if (isdigit(*argv[1]))
         min_size = atoi(argv[1]);

   if (argc > 2)
      if (isdigit(*argv[2]))
         max_size = atoi(argv[2]);

   if (argc > 3)
      if (isdigit(*argv[3]))
         niters = atoi(argv[3]);

   if (argc > 4)
      if (isdigit(*argv[4]))
         method = atoi(argv[4]);

   if (argc > 5)
      if (isdigit(*argv[4]) || *argv[4] == '.')
         stepSize = atof(argv[5]);

   matmul_ptr methods[] = {matmul_blas, matmul_naive, matmul_opt1, matmul_opt2, matmul_opt3};
   const char *method_names[] = {"matmul_blas", "matmul_naive", "matmul_vect(opt1)", "matmul_block(opt2)", "matmul_unroll(opt3)"};

   if (method < 0 || method > sizeof(methods)/sizeof(methods[0]))
   {
      fprintf(stderr,"Invalid method selected %d\n", method);
      return 1;
   }
   else
      fprintf(stderr,"Using method[%d]=%s, min/max=%d/%d, step=%f\n", method, method_names[method], min_size, max_size, stepSize);

#ifdef unroll
   printf("unroll = %d\n", unroll);
#endif
#ifdef blockSize
   printf("blockSize = %d\n", blockSize);
#endif

   // Check the timer accuracy.
   double tDelta = 1e50;
   {
      int nmax = std::min(10000,max_size);

      ValueType *a = new ValueType[nmax];

      for (int i = 0; i < nmax; ++i)
         a[i] = 0.0;

      myTimer_t t0 = getTimeStamp();
      for (int i = 0; i < nmax; ++i)
      {
         myTimer_t t1 = getTimeStamp();
         a[i] = getElapsedTime(t0,t1);
         t0 = t1;
      }

      for (int i = 0; i < nmax-1; ++i)
         if (std::abs(a[i+1] - a[i]) > 0.0)
            tDelta = std::min(tDelta, std::abs(a[i+1]-a[i]));

      printf("Smallest detectable time = %e (ms)\n", tDelta*1000);

      printf("getTicksPerSecond = %e\n", getTicksPerSecond());

      delete [] a;
   }

   printf("size, time (ms), GFLOP/s, Confident(%), # runs, error\n");

   int size = min_size;
   for (; size <= max_size; size *= stepSize)
      run_matmul (size, niters, tDelta, methods[method]);

   return 0;
}
