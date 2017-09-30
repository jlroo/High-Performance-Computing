#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <string>
#include <vector>
#include <typeinfo>
#include <type_traits>

#include <cmath>
#include <algorithm>

#include <my_timer.h>

#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

#ifndef RESTRICT
#define RESTRICT
#endif

// Get machine 'tiny'
const double sfmin = std::numeric_limits<double>::epsilon();

// Library definitions
extern "C"
{
   void dgetrf_(int* m , int* n, double* A, int* lda, int* P, int* info);
   void dgetrf2_(int* m , int* n, double* A, int* lda, int* P, int* info);
   void dgetrs_( char* TRANS, int* N, int* NRHS, double A[], int* LDA, int IPIV[], double B[], int* LDB, int* INFO);
   void dlaswp_(int* m , double* A, int* lda, int* k1, int* k2, int* P, int* inc);
   //void dtrsm_(const char*, const char*, const char*, const char*, int* m, int* n, double *alpha, double *A, int *lda, double *B, *ldb);

   double dlamch_(char);
   int idamax_( int* n, double* x, int* incx);
   void dswap_( int* n, double* x, int* incx, double* y, int* incy);
   void dcopy_( int* n, double* x, int* incx, double* y, int* incy);
   void dscal_( int* n, double* alpha, double* x, int* incx);
   void daxpy_( int* n, double* alpha, double* x, int* incx, double* y, int* incy);
}

template <typename Matrix, typename Vector>
int ludec_naive( Matrix& A, Vector& P )
{
   const int m = A.mRows();
   const int n = A.nCols();

   for (int k = 0; k < n; ++k)
   {
      // Find the max.abs value on this column below the diagonal.
      int ipiv = k;
      double maxval = fabs( A(k,k) );
      for (int i = k+1; i < m; ++i)
      {
         const double val = std::abs( A(i,k) );
         if (val > maxval)
         {
            ipiv = i;
            maxval = val;
         }
      }

      // Test of singular value.
      if (maxval == 0.0)
         return (k+1);

      P(k) = (ipiv+1); // Convert to Base-1 (to match Lapack)

      // Swap rows if necessary.
      if (ipiv != k)
      {
         for (int j = 0; j < n; ++j)
            std::swap( A(k,j), A(ipiv,j) );
      }

      // Scale the elements below the diagonal. Be careful with round-off.
      if (maxval > sfmin)
      {
         const double a_kk_recip = 1.0 / A(k,k);
         for (int i = k+1; i < m; ++i)
            A(i,k) *= a_kk_recip;
      }
      else
      {
         const double a_kk = A(k,k);
         for (int i = k+1; i < m; ++i)
            A(i,k) /= a_kk;
      }

      // Update the trailing submatrix now.
      for (int j = k+1; j < n; ++j)
      {
         const double a_kj = A(k,j);
         if (a_kj != 0.0)
         {
            //#pragma ivdep
            for (int i = k+1; i < m; ++i)
               A(i,j) -= a_kj * A(i,k);
         }
      }
   }

   return 0;
}

template <typename Matrix, typename Vector>
int ludec_blas( Matrix& A, Vector& P )
{
   int m = A.mRows();
   int n = A.nCols();
   int ldim = A.lDim();

   int one = 1;

   for (int k = 0; k < n; ++k)
   {
      // Some useful array lengths we'll need to pass in as pointers.
      int m_minus_k = m-k;
      int m_minus_k_plus_one = m-(k+1);

      // Find the max.abs value on this column below the diagonal.
      const int ipiv = idamax_( &m_minus_k, &A(k,k), &one ) - 1 + k;

      P(k) = ipiv+1;

      // Swap rows if necessary.
      if (ipiv != k)
         dswap_( &n, &A(k,0), &ldim, &A(ipiv,0), &ldim );

      // Test of singular value.

      // Scale the elements below the diagonal. Be careful with round-off.
      const double a_kk = A(k,k);
      if (a_kk == 0.0)
         return (k+1);
      else if (a_kk > sfmin)
      {
         double a_kk_recip = 1.0 / a_kk;
         dscal_( &m_minus_k_plus_one, &a_kk_recip, &A(k+1,k), &one );
      }
      else
         for (int i = k+1; i < m; ++i)
            A(i,k) /= a_kk;

      // Update the trailing submatrix now.
      for (int j = k+1; j < n; ++j)
      {
         double neg_a_kj = -1.0 * A(k,j);
         if (neg_a_kj != 0.0)
            daxpy_(&m_minus_k_plus_one, &neg_a_kj, &A(k+1,k), &one, &A(k+1,j), &one);
      }
   }

   return 0;
}

const int NoTranspose = 0;
const int Transpose = ~NoTranspose; 
const int SideLeft = 0;
const int SideRight = ~SideLeft;
const int TriUpper = 0;
const int TriLower = ~TriUpper;
const int UnitDiagonal = 1;
const int NotUnitDiagonal = ~UnitDiagonal;

// c := alpha A*b + beta*c
template <typename Matrix, typename Vector>
int dgemv( const int Trans, const double alpha, const Matrix& A, const Vector& b, const double beta, Vector& c )
{
   if (Trans == Transpose) { fprintf(stderr,"Transpose not ready in dgemv\n"); exit(1); }

   const int m = A.mRows();
   const int n = b.len();
   const int k = A.nCols();
   if (c.len() != m) { fprintf(stderr,"Mismatch in dgemv %d %d %d %d\n", m, n, k, c.len()); exit(1); }

   if (beta == 0.0)
      for (int i = 0; i < m; ++i)
         c(i) = 0.0;
   else if (beta != 1.0)
      for (int i = 0; i < m; ++i)
         c(i) *= beta;

   for (int l = 0; l < k; ++l)
   {
      const double alpha_bl = alpha * b(l);
      for (int i = 0; i < m; ++i)
         c(i) += alpha_bl * A(i,l);
   }
}
// C := alpha A*B + beta*C
template <typename Matrix>
int dgemm( const int TransA, const int TransB, const double alpha, const Matrix& A, const Matrix& B, const double beta, Matrix& C )
{
   if (TransA == Transpose or TransB == Transpose) { fprintf(stderr,"Transpose not ready in dgemm\n"); exit(1); }

   const int m = A.mRows();
   const int n = B.nCols();
   const int k = A.nCols();
   if (C.mRows() != m or C.nCols() != n) { fprintf(stderr,"Mismatch in dgemm %d %d %d %d %d\n", m, n, k, C.mRows(), C.nCols()); exit(1); }

   for (int j = 0; j < n; ++j)
   {
      if (beta == 0.0)
         for (int i = 0; i < m; ++i)
            C(i,j) = 0.0;
      else if (beta != 1.0)
         for (int i = 0; i < m; ++i)
            C(i,j) *= beta;

      for (int l = 0; l < k; ++l)
      {
         const double alpha_blj = alpha * B(l,j);
         for (int i = 0; i < m; ++i)
            C(i,j) += alpha_blj * A(i,l);
      }
   }
}

template <typename Matrix>
int dtrsm( const int Side, const int UpperLower, const int TransA, const int Diag, const Matrix& A, Matrix& B )
{
   const int ncols = B.nCols();
   const int m = A.mRows();
   const int n = A.nCols();

   if (TransA != NoTranspose) { fprintf(stderr,"Tranpose not ready in dtrsm\n"); exit(1); }

   if (Side == SideLeft)
   {
      if (B.mRows() != A.nCols()) { fprintf(stderr,"Mismatch in dtrsm for SideLeft %d %d\n", B.mRows(), A.nCols()); exit(1); }

      // A is lower triangular; solve L X = B
      if (UpperLower == TriLower)
      {
         for (int j = 0; j < ncols; ++j)
            for (int k = 0; k < m; ++k)
            {
               double b_kj = B(k,j);
               if (b_kj != 0.0)
               {
                  if (Diag == NotUnitDiagonal) b_kj = ( B(k,j) /= A(k,k) );
                  for (int i = k+1; i < m; ++i)
                     B(i,j) -= b_kj * A(i,k);
               }
            }
      }
      // A is Upper triangular; solve U X = B
      else
      {
         for (int j = 0; j < ncols; ++j)
            for (int k = m-1; k >= 0; --k)
            {
               double b_kj = B(k,j);
               if (b_kj != 0.0)
               {
                  if (Diag == NotUnitDiagonal) { B(k,j) /= A(k,k); b_kj = B(k,j); }
                  for (int i = 0; i < k; ++i)
                      B(i,j) -= b_kj*A(i,k);
               }
            }
      }
   }
   else // Right matrix problem: X A = B for X
   {
      if (B.nCols() != A.mRows()) { fprintf(stderr,"Mismatch in dtrsm for SideRight %d %d\n", B.nCols(), A.mRows()); exit(1); }

      const int mrows_b = B.mRows();
      const int ncols_b = B.nCols();

      // A is Upper triangular; solve X U = B for X
      if (UpperLower == TriUpper)
      {
         for (int j = 0; j < ncols_b; ++j)
         {
            for (int k = 0; k < j; ++k)
            {
               const double a_kj = A(k,j);
               if (a_kj != 0.0)
               {
                  //#pragma omp simd
                  for (int i = 0; i < mrows_b; ++i)
                      B(i,j) -= a_kj*B(i,k);
               }
            }
            if (Diag == NotUnitDiagonal)
            {
               const double ajj_recip = 1.0 / A(j,j);
               for (int i = 0; i < mrows_b; ++i)
                  B(i,j) *= ajj_recip;
            }
         }
      }
      else
      {
         fprintf(stderr,"TriLower with SideRight not ready in dtrsm\n");
         exit(1);
      }
   }
}

template <typename Matrix, typename Vector>
int dlaswp( Matrix& A, const Vector& P )
{
   const int m = A.mRows();
   const int n = A.nCols();

   for (int i = 0; i < m; ++i)
   {
      const int ipiv = P(i)-1; // P() is base-1 alla Lapack Fortran.

      if (ipiv != i)
         for (int j = 0; j < n; ++j)
            std::swap( A(i,j), A(ipiv,j) );
   }
}

template <typename Matrix, typename Vector>
int ludec_block( Matrix& A, Vector& P, const int blockSize = 32 )
{
   int m = A.mRows();
   int n = A.nCols();
   int ldim = A.lDim();

   if (std::min(m,n) < blockSize)
   {
      fprintf(stderr,"Block size too small %d %d\n", m, blockSize);
      return ludec_naive( A, P );
   }

   // Loop over column blocks
   for (int k = 0; k < n; k += blockSize)
   {
      // First/last colums of this tall/skiny column block.
      const int k1 = std::min( n, k + blockSize);

      // Factorize the column block A(k0:m,k0:k1)
      {
         Matrix Ak ( &A(k,k), m-k, k1-k, ldim ); // column block including the diag block.
         Vector Pk ( &P(k), m-k );
         int ierr = 0;
         if (1) // Use the local version.
            ierr = ludec_naive( Ak, Pk );
         else
         { // Or use the Lapack recursive method.
            int m_blk = Ak.mRows();
            int n_blk = Ak.nCols();
            int l_blk = Ak.lDim();
            dgetrf2_( &m_blk, &n_blk, Ak.getRawPointer(), &l_blk, Pk.getRawPointer(), &ierr );
         }

         if (ierr != 0)
            return ierr+k;
      }

      // Adjust the pivot.
      for (int i = k; i < k1; ++i)
      {
         P(i) += k;
         const int ipiv = P(i)-1;

         if (ipiv != i)
         {
            // Apply the row swaps to the left columns.
            for (int j = 0; j < k; ++j)
               std::swap( A(i,j), A(ipiv,j) );

            // Apply the row swaps to the untouched lower-right columns.
            for (int j = k1; j < n; ++j)
               std::swap( A(i,j), A(ipiv,j) );
         }
      }

      // Update the trailing submatrix now.
      if (k1 < n)
      {
         // Solve for U(k:k1,k1:m) ... L_kk U_:k = A
         Matrix Uk ( &A(k,k1), k1-k, n-k1, ldim ); // row block to the right of the diag block.
         Matrix Lkk ( &A(k,k), k1-k, k1-k, ldim ); // lower tri matrix of the diagonal block.
         dtrsm( SideLeft, TriLower, NoTranspose, UnitDiagonal, Lkk, Uk );

         if (k1 < m)
         {
            Matrix Lk ( &A(k1,k), m-k1, k1-k, ldim ); // column block below diagonal block.
            Matrix Alr ( &A(k1,k1), m-k1, n-k1, ldim ); // submatrix untouched lower-right.
            // A(k+1,k+1) -= L(k+1,k) * U(k,k+1)
            dgemm( NoTranspose, NoTranspose, -1.0, Lk, Uk, 1.0, Alr );
         }
      }
   }

   return 0;
}

template <typename Matrix, typename Vector>
int ludec_tile( Matrix& A, Vector& P, const int blockSize = 32 )
{
   int m = A.mRows();
   int n = A.nCols();
   int ldim = A.lDim();

   if (std::min(m,n) < blockSize)
   {
      fprintf(stderr,"Block size too small %d %d\n", m, blockSize);
      return ludec_naive( A, P );
   }

   // This follows a very similar strategy as the classic Block approach.
   // But, instead of always solving the tall-skinny or short-fat column/row
   // matrices, solve the tall-skinny column block and then solve several
   // small tiles for the Upper blocks and the trailing submatrix.

   // Loop over column blocks
   for (int k0 = 0; k0 < n; k0 += blockSize)
   {
      // First/last colums of this tall/skiny column block.
      const int k1 = std::min( n, k0 + blockSize);

      // Factorize the column block A(k0:m,k0:k1)
      {
         Matrix Ak ( &A(k0,k0), m-k0, k1-k0, ldim ); // column block including the diag block.
         Vector Pk ( &P(k0), m-k0 );
         int ierr = ludec_naive( Ak, Pk );
         if (ierr != 0)
            return ierr+k0;
      }

      // Adjust the pivot.
      for (int i = k0; i < k1; ++i)
      {
         P(i) += k0;
         const int ipiv = P(i)-1;

         if (ipiv != i)
         {
            // Apply the row swaps to the left columns.
            for (int j = 0; j < k0; ++j)
               std::swap( A(i,j), A(ipiv,j) );

            // Apply the row swaps to the untouched righrward olumns.
            for (int j = k1; j < n; ++j)
               std::swap( A(i,j), A(ipiv,j) );
         }
      }

      // Update the trailing submatrix now.
      #pragma omp parallel for
      for (int j0 = k1; j0 < n; j0 += blockSize)
      {
         const int j1 = std::min(n, j0 + blockSize);

         // Solve the upper block matrices L_kk U(k0:k1,j0:j1) = A(k0:k1,j0:j1)
         Matrix Ukj ( &A(k0,j0), k1-k0, j1-j0, ldim ); // row block to the right of the diag.
         Matrix Lkk ( &A(k0,k0), k1-k0, k1-k0, ldim ); // lower tri matrix of the diagonal.
         dtrsm( SideLeft, TriLower, NoTranspose, UnitDiagonal, Lkk, Ukj );

         // Update the trailing blocks ... A_ij -= Lik Ukj
         for (int i0 = k1; i0 < m; i0 += blockSize)
         {
            const int i1 = std::min(m, i0 + blockSize);

            Matrix Lik ( &A(i0,k0), i1-i0, k1-k0, ldim ); // column block below diagonal block.
            Matrix Aij ( &A(i0,j0), i1-i0, j1-j0, ldim ); // submatrix untouched lower-right.

            dgemm( NoTranspose, NoTranspose, -1.0, Lik, Ukj, 1.0, Aij );
         }
      }
   }

   return 0;
}

// Tiled approach but factor only the diagonal block, not the 
// full tall-skiny column block. This limits the pivot to with
// the rows of each tile. Doubtful it's as stable.
template <typename Matrix, typename Vector>
int ludec_tile_v2( Matrix& A, Vector& P, const int blockSize = 32 )
{
   int m = A.mRows();
   int n = A.nCols();
   int ldim = A.lDim();

   if (std::min(m,n) < blockSize)
   {
      fprintf(stderr,"Block size too small %d %d\n", m, blockSize);
      return ludec_naive( A, P );
   }

   // Assuming a square matrix ...
   if (m != n)
   {
      fprintf(stderr,"Matrix is not square %d %d\n", m, n);
      return 1;
   }

   // This follows a very similar strategy as the classic Block approach.
   // But, instead of always solving the tall-skinny or short-fat column/row
   // matrices, solve the tall-skinny column block and then solve several
   // small tiles for the Upper blocks and the trailing submatrix.

   // Loop over the diagonal tiles.
   for (int k0 = 0; k0 < n; k0 += blockSize)
   {
      // First/last colums of this tile.
      const int k1 = std::min( n, k0 + blockSize);

      // Factorize the diagonal tile A(k0:k1,k0:k1)
      Matrix Akk ( &A(k0,k0), k1-k0, k1-k0, ldim );
      Vector Pk ( &P(k0), k1-k0 );
      int ierr = ludec_naive( Akk, Pk );
      if (ierr != 0)
         return ierr+k0;

      #pragma omp parallel
      {

      // Pivot the left tiles on this row.
      #pragma omp for
      for (int j0 = 0; j0 < k0; j0 += blockSize)
      {
         const int j1 = std::min(n, j0 + blockSize);

         Matrix  Lkj ( &A(k0,j0), k1-k0, j1-j0, ldim );

         dlaswp( Lkj, Pk );
      }

      // Solve the upper tiles on this row.
      #pragma omp for
      for (int j0 = k1; j0 < n; j0 += blockSize)
      {
         const int j1 = std::min(n, j0 + blockSize);

         Matrix  Ukj ( &A(k0,j0), k1-k0, j1-j0, ldim ); // row block to the right of the diag.

         // Apply the pivot before the tri solve.
         dlaswp( Ukj, Pk );

         // Solve the upper block matrices L_kk U(k0:k1,j0:j1) = A(k0:k1,j0:j1)
         dtrsm( SideLeft, TriLower, NoTranspose, UnitDiagonal, Akk, Ukj );
      }

      // Done with the local pivot. Reset to a global index.
      #pragma omp single
      for (int i = k0; i < k1; ++i)
         P(i) += k0;

      // implicit barrier after the single

      // Solve for the lower tiles on this column
      #pragma omp for
      for (int i0 = k1; i0 < m; i0 += blockSize)
      {
         const int i1 = std::min(m, i0 + blockSize);

         Matrix Lik ( &A(i0,k0), i1-i0, k1-k0, ldim ); // column block below diagonal block.

         // Solve the lower tile matrices L(i0:i1,k0:k1) U_kk = A(i0:i1,k0:k1)
         dtrsm( SideRight, TriUpper, NoTranspose, NotUnitDiagonal, Akk, Lik );

         // Update the trailing submatrix tiles ... A_ij -= Lik Ukj
         for (int j0 = k1; j0 < n; j0 += blockSize)
         {
            const int j1 = std::min(n, j0 + blockSize);

            Matrix Ukj ( &A(k0,j0), k1-k0, j1-j0, ldim ); // row block to the right of the diag.
            Matrix Aij ( &A(i0,j0), i1-i0, j1-j0, ldim ); // submatrix untouched lower-right.

            dgemm( NoTranspose, NoTranspose, -1.0, Lik, Ukj, 1.0, Aij );
         }
      }

      }// end parallel
   }

   return 0;
}

template <typename Matrix, typename Vector>
int ludec_lapack( Matrix& A, Vector& P )
{
   int ierr;

   int m = A.mRows();
   int n = A.nCols();
   int ldim = A.lDim();

   dgetrf_( &m, &n, A.getRawPointer(), &ldim, P.getRawPointer(), &ierr );

   return ierr;
}

template <typename ValueType>
struct MatrixType
{
   typedef ValueType value_type;
   typedef int index_type;

   value_type *_data;
   index_type _mrows, _ncols, _ldim;

   MatrixType (void) : _data(NULL), _mrows(0), _ncols(0), _ldim(0) {};
   MatrixType (value_type *data, const int mrows, const int ncols) : _data(data), _mrows(mrows), _ncols(ncols), _ldim(mrows) {};
   MatrixType (value_type *data, const int mrows, const int ncols, const int ldim) : _data(data), _mrows(mrows), _ncols(ncols), _ldim(ldim) {};

   void copy( const MatrixType<value_type> &From)
   {
      if (this->mRows() != From.mRows() or
          this->nCols() != From.nCols())
      {
         fprintf(stderr,"Error in MatrixType::copy: mismatch\n");
         exit(1);
      }

      for (int j = 0; j < this->nCols(); ++j)
         for (int i = 0; i < this->mRows(); ++i)
            this->operator()(i,j) = From(i,j);
   }

   void printInfo(FILE *fstream = stdout) const
   {
      fprintf(fstream, "MatrixType<%s>; _data= %lx _mrows= %d _ncols= %d _ldim= %d\n", typeid(value_type).name(), _data, _mrows, _ncols, _ldim);
   }

   void print(FILE *fstream = stdout) const
   {
      for (int i = 0; i < this->mRows(); i++)
      {
         for (int j = 0; j < this->nCols(); j++) fprintf(fstream, "%9.5f ", this->operator()(i,j));
         fprintf(fstream,"\n");
      }
   }

   inline value_type & operator() (const int i, const int j)
   {
      const int idx = i + j*this->_ldim;
      assert (idx >= 0 && idx < (this->_ncols * this->_ldim));
      return this->_data[i + j*this->_ldim];
   }
   inline const value_type & operator() (const int i, const int j) const
   {
      return this->_data[i + j*this->_ldim];
   }

   value_type * getRawPointer(void) { return this->_data; }

   int mRows(void) const { return this->_mrows; }
   int nCols(void) const { return this->_ncols; }
   int lDim(void) const { return this->_ldim; }
};

template <typename ValueType>
struct VectorType
{
   typedef ValueType value_type;
   typedef int index_type;

   value_type *_data;
   index_type _len;

   VectorType (void) : _data(NULL), _len(0) {};
   VectorType (value_type *data, const int len) : _data(data), _len(len) {};

   inline value_type & operator() (const int i)
   {
      return this->_data[i];
   }
   inline const value_type & operator() (const int i) const
   {
      return this->_data[i];
   }

   value_type * getRawPointer (void) { return this->_data; }

   int length(void) const { return this->_len; }
};

template <typename ValueType>
ValueType sqr(const ValueType &val) { return val*val; }

void show_usage( const char* prog )
{
   printf("Usage for %s\n", prog);
   printf("\t--size   | -s <int value> : Matrix size. (100)\n");
   printf("\t--method | -m <int value> : Method choices are ... (0)\n");
   printf("\t--test   | -t             : Only run verfication ... (false)\n");
   printf("\t--miters      <int value> : Max # of iterations to run ... (1)\n");
   printf("\t\t0: Lapack (vendor)\n");
   printf("\t\t1: naive\n");
   printf("\t\t2: BLAS\n");
   printf("\t\t3: block\n");
   printf("\t\t4: tile\n");
   printf("\t\t5: tile_v2\n");
}

int main (int argc, char * argv[])
{
   int n = 100;
   int miters = 1; // Max # of factorizations to try.
   int method = 0; // blas
   int blockSize = 32;

   // Get user inputs.
   {
      #define check_index(_i) { if ((_i) >= argc){ fprintf(stderr,"Missing value for argument %s\n", for
      for (int i = 1; i < argc; i++)
      {
         std::string arg = argv[i];
         if (arg == "--size" || arg == "-s")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for --size\n"); show_usage(argv[0]); return 1; }
            n = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "--method" || arg == "-m")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for --method\n"); show_usage(argv[0]); return 1; }
            method = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "--miters")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for --miters\n"); show_usage(argv[0]); return 1; }
            miters = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "--blocksize" || arg == "-b")
         {
            if ((i+1) >= argc) { fprintf(stderr,"Missing value for --blocksize\n"); show_usage(argv[0]); return 1; }
            blockSize = atoi( argv[i+1] );
            i++;
         }
         else if (arg == "--help" || arg == "-h")
         {
            show_usage(argv[0]); return 0;
         }
         else if (arg == "--test" || arg == "-t")
         {
            miters = 0;
         }
         else // Unknown
         {
            fprintf(stderr,"Unknown CLI option %s\n", arg.c_str()); show_usage(argv[0]); return 1;
         }
      }
   }

   std::vector< std::string > methods;

   methods.push_back( "lapack" );
   methods.push_back( "naive" );
   methods.push_back( "blas" );
   methods.push_back( "block" );
   methods.push_back( "tile" );
   methods.push_back( "tile_v2" );

   if (method < 0 || method >= methods.size())
   {
      fprintf(stderr,"Invalid method selected %d\n", method);
      return 1;
   }
   else
      fprintf(stderr,"Using method[%d]=%s size=%d\n", method, methods[method].c_str(), n);

   if (method == 3 or method == 4)
      fprintf(stderr,"blockSize= %d\n", blockSize);

   // Allocate arrays with a little padding.
   double *_A  = new double[n*n];
   double *_LU = new double[n*n];
   int    *_P  = new int[n];

   MatrixType<double> A(_A, n, n);
   MatrixType<double> LU(_LU, n, n);
   VectorType<int> P(_P, n);

   //A.printInfo();

   // Define A[] with random values
   srand(n);
   double invRandMax = 1.0 / double(RAND_MAX);
   for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
         A(i,j) = double( rand() ) * invRandMax;

   if (n <= 6)
   {
      printf("---A:\n");
      A.print(stdout);
   }

   double tCalc = 0.0;
   int niters = std::max(1,miters);
   double sum = 0.0;
   bool firstPass = true;
   int ierr = 0;
   while (true)
   {
      myTimer_t t0 = getTimeStamp();

      for (int iter = 0; iter < niters; iter++)
      {
         LU.copy( A );

         if (method == 0)
            ierr = ludec_lapack( LU, P );
         else if (method == 1)
            ierr = ludec_naive( LU, P );
         else if (method == 2)
            ierr = ludec_blas( LU, P );
         else if (method == 3)
            ierr = ludec_block( LU, P, blockSize );
         else if (method == 4)
            ierr = ludec_tile( LU, P, blockSize );
         else if (method == 5)
            ierr = ludec_tile_v2( LU, P, blockSize );

         if (ierr != 0)
            break;

         sum += LU(0,0);
      }

      myTimer_t t1 = getTimeStamp();

      tCalc = getElapsedTime(t0,t1);
      if ((tCalc > 0.1 and firstPass == false) or miters == 0) break;

      if (firstPass)
         firstPass = false;
      else
         niters *= 2;
   }

   if (ierr != 0)
   {
      fprintf(stderr,"Error calling LU decomp %d\n", ierr);
      return 1;
   }

   double tCopy = 0.0;
   if (miters > 0)
   {
      myTimer_t t0 = getTimeStamp();
      for (int iter = 0; iter < niters; iter++)
      {
         LU.copy( A );
         sum += LU(0,0);
      }
      myTimer_t t1 = getTimeStamp();
      tCopy = getElapsedTime(t0,t1);
   }

   printf("LU Factorization took (sec): %d %g %g %d\n", n, tCalc / niters, tCopy / niters, niters);

   if (miters == 0)
   {
      double *_LU_ref = new double[n*n];
      int    *_P_ref  = new int[n];
      MatrixType<double> LU_ref(_LU_ref, n, n);
      VectorType<int> P_ref(_P_ref, n);

      LU_ref.copy( A );

      // Run the Lapack version to get the correct answer.
      int ierr = ludec_lapack( LU_ref, P_ref );
      if (ierr != 0)
      {
         fprintf(stderr,"Error calling ludec_lapack: %d\n", ierr);
         return 1;
      }

      if (n <= 6)
      {
         printf("\n---LU:\n");
         LU.print(stdout);

         for (int i = 0; i < A.mRows(); i++)
            printf("%d\n", P(i));

         printf("\nRef\n");

         LU_ref.print(stdout);
         for (int i = 0; i < A.mRows(); i++)
            printf("%d\n", P_ref(i));
      }

      double err2 = 0, ref2 = 0;
      for (int i = 0; i < n; ++i)
         for (int j = 0; j < n; ++j)
         {
            const double val = LU(i,j);
            const double ref = LU_ref(i,j);
            err2 += sqr( val - ref );
            ref2 += sqr( ref );
         }

      if (ref2 < 1e-20) ref2 = 1.0;
      double error = sqrt( err2 / ref2 );
      bool passed = error < 1e-12;
      fprintf(stderr,"LU factorization error= %e %s\n", error, (passed ? "PASSED" : "FAILED"));

      // Test the solution of a problem Ax=b with the factorization.
      if (!passed)
      {
         double *_x = new double[n];
         double *_b = new double[n];

         MatrixType<double> x(_x, n, 1);
         MatrixType<double> b(_b, n, 1);

         for (int i = 0; i < n; ++i)
         {
            x(i,0) = double( rand() ) * invRandMax;
            //printf("P(%d)= %d\n", i, P(i));
         }

         // Compute Ax = b to set the RHS with a known vector.
         dgemm( NoTranspose, NoTranspose, 1.0, A, x, 0.0, b);

         //for (int i = 0; i < n; ++i)
         //   printf("P(%d)= %d %f %f\n", i, P(i), x(i,0), b(i,0));

         if (0)
         {
            const char *TRANS = "NoTranpose";
            int lenstr = strlen(TRANS);
            int one = 1;
            int info = 0;
            dgetrs_ (const_cast<char*>(TRANS), &n, &one, LU.getRawPointer(), &n, P.getRawPointer(), b.getRawPointer(), &n, &info );
         }
         else
         {
            dlaswp( b, P );

            //for (int i = 0; i < n; ++i)
            //   printf("P(%d)= %d %f %f\n", i, P(i), x(i,0), b(i,0));

            // Now solve for x' as Ax' = b ... => PLUx = Pb
            // 1. Ly = b
            // 2. Ux = y
            dtrsm( SideLeft, TriLower, NoTranspose,    UnitDiagonal, LU, b );
            dtrsm( SideLeft, TriUpper, NoTranspose, NotUnitDiagonal, LU, b );
         }

         double verr2 = 0, vref2 = 0;
         for (int i = 0; i < n; ++i)
         {
            double diff = x(i,0) - b(i,0);
            verr2 += sqr( diff );
            vref2 += sqr( x(i,0) );
            if (n <= 6)
               printf("x(%d)= %f %f\n", i, x(i,0), b(i,0));
         }

         if (vref2 < 1e-20) vref2 = 1.0;
         error = sqrt( verr2 / vref2 );
         passed = error < 1e-12;
         fprintf(stderr,"LU factorization matrix norm error= %e %s\n", error, (passed ? "PASSED" : "FAILED"));

         delete [] _x;
         delete [] _b;
      }

      delete [] _LU_ref;
      delete [] _P_ref;
   }

   delete [] _A;
   delete [] _LU;
   delete [] _P;

   return 0;
}
