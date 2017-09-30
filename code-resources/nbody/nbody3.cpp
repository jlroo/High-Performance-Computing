#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <vector>
#include <limits>
#include <cmath>

#include <my_timer.h>
#include <aligned_allocator.h>

#ifndef __RESTRICT
#  define __RESTRICT
#endif

#define NDIM (3)

#ifdef USE_FLOAT
typedef float  ValueType;
#else
typedef double ValueType;
#endif


const ValueType G = 1.0;
const ValueType TINY = std::numeric_limits<ValueType>::epsilon();
const ValueType TINY2 = TINY*TINY;

//#define Enable_ArrayOfStructures
#if defined(Enable_ArrayOfStructures) || defined(__AOS)
#  ifndef Enable_ArrayOfStructures
#    define Enable_ArrayOfStructures
#  endif
   /* Array-of-structures (like) format. */
#  define _index(i,j) (NDIM*(i) + (j))
#else
   /* Structure-of-arrays (like) format. */
#  define _index(i,j) ((i) + (j)*n)
#endif

#define acc_array(i,j) acc[ _index((i),(j)) ]
#define pos_array(i,j) pos[ _index((i),(j)) ]
#define vel_array(i,j) vel[ _index((i),(j)) ]

/* Generate a random double between 0,1. */
ValueType frand(void) { return ((ValueType) rand()) / RAND_MAX; }

void accel (ValueType pos[], ValueType vel[], ValueType mass[], ValueType acc[], const int n)
{
   for (int i = 0; i < n; ++i)
      for (int k = 0; k < NDIM; ++k)
         acc_array(i,k) = ValueType(0);

   for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
      {
         /* Position vector from i to j and the distance^2. */
         ValueType rx = pos_array(j,0) - pos_array(i,0);
         ValueType ry = pos_array(j,1) - pos_array(i,1);
         ValueType rz = pos_array(j,2) - pos_array(i,2);
         ValueType dsq = rx*rx + ry*ry + rz*rz + TINY2;
         ValueType dcu = dsq * std::sqrt(dsq);

         acc_array(i,0) += G * mass[j] * rx / dcu;
         acc_array(i,1) += G * mass[j] * ry / dcu;
         acc_array(i,2) += G * mass[j] * rz / dcu;
      }
}

// Operator strength reduction: replace equivalent math with cheaper operations.
void accel_opt1 (ValueType pos[], ValueType vel[], ValueType mass[], ValueType acc[], const int n)
{
   for (int i = 0; i < n; ++i)
      for (int k = 0; k < NDIM; ++k)
         acc_array(i,k) = ValueType(0);

   for (int i = 0; i < n; ++i)
   {
      for (int j = 0; j < n; ++j)
      {
         /* Position vector from i to j and the distance^2. */
         ValueType rx = pos_array(j,0) - pos_array(i,0);
         ValueType ry = pos_array(j,1) - pos_array(i,1);
         ValueType rz = pos_array(j,2) - pos_array(i,2);
         ValueType dsq = rx*rx + ry*ry + rz*rz + TINY2;
         // Store 1/r^3 once instead of 3x divisions.
         ValueType m_invR3 = mass[j] / (dsq * std::sqrt(dsq));

         acc_array(i,0) += rx * m_invR3;
         acc_array(i,1) += ry * m_invR3;
         acc_array(i,2) += rz * m_invR3;
      }

      // Pull G out of repeated calculation.
      for (int k = 0; k < NDIM; ++k)
         acc_array(i,k) *= G;
   }
}

// Store target data in registers: Compiler "may" do this automatically but
// it often helps with cache efficiency. This can be especially helpfule
// by avoiding repeated writes which are several times slower than reads.
void accel_opt2 (ValueType * __RESTRICT pos, ValueType * __RESTRICT vel, ValueType * __RESTRICT mass, ValueType * __RESTRICT acc, const int n)
{
#ifdef Enable_Alignment
__assume_aligned(pos, Alignment());
__assume_aligned(acc, Alignment());
__assume_aligned(mass, Alignment());
__assume(n % (Alignment()/sizeof(ValueType)) == 0 );
#endif

   for (int i = 0; i < n; ++i)
   {
      ValueType ax = 0, ay = 0, az = 0;
      const ValueType xi = pos_array(i,0);
      const ValueType yi = pos_array(i,1);
      const ValueType zi = pos_array(i,2);

      for (int j = 0; j < n; ++j)
      {
         /* Position vector from i to j and the distance^2. */
         ValueType rx = pos_array(j,0) - xi;
         ValueType ry = pos_array(j,1) - yi;
         ValueType rz = pos_array(j,2) - zi;
         ValueType dsq = rx*rx + ry*ry + rz*rz + TINY2;
         ValueType m_invR3 = mass[j] / (dsq * std::sqrt(dsq));

         ax += rx * m_invR3;
         ay += ry * m_invR3;
         az += rz * m_invR3;
      }

      acc_array(i,0) = G * ax;
      acc_array(i,1) = G * ay;
      acc_array(i,2) = G * az;
   }
}
void accel_block (ValueType * __RESTRICT pos, ValueType * __RESTRICT vel, ValueType * __RESTRICT mass, ValueType * __RESTRICT acc, const int n)
{
#ifdef Enable_Alignment
__assume_aligned(pos, Alignment());
__assume_aligned(acc, Alignment());
__assume_aligned(mass, Alignment());
__assume(n % (Alignment()/sizeof(ValueType)) == 0 );
#endif

   const int tileSize = 128;

   for (int istart = 0; istart < n; istart += tileSize)
   {
      const int istop = std::min(n, istart + tileSize);

      for (int i = istart; i < istop; ++i)
      {
         acc_array(i,0) = 0.0;
         acc_array(i,1) = 0.0;
         acc_array(i,2) = 0.0;
      }

      for (int jstart = 0; jstart < n; jstart += tileSize)
      {
         const int jstop = std::min(n, jstart + tileSize);

         for (int i = istart; i < istop; ++i)
         {
            ValueType ax = 0, ay = 0, az = 0;
            const ValueType xi = pos_array(i,0);
            const ValueType yi = pos_array(i,1);
            const ValueType zi = pos_array(i,2);

            for (int j = jstart; j < jstop; ++j)
            {
               /* Position vector from i to j and the distance^2. */
               ValueType rx = pos_array(j,0) - xi;
               ValueType ry = pos_array(j,1) - yi;
               ValueType rz = pos_array(j,2) - zi;
               ValueType dsq = rx*rx + ry*ry + rz*rz + TINY2;
               ValueType m_invR3 = mass[j] / (dsq * std::sqrt(dsq));

               ax += rx * m_invR3;
               ay += ry * m_invR3;
               az += rz * m_invR3;
            }

            acc_array(i,0) += G * ax;
            acc_array(i,1) += G * ay;
            acc_array(i,2) += G * az;
         }
      }
   }
}

// Cache blocking or tiling: Restructure code to reuse data in nearest cache.
// May be of lesser value in N-body since the Calc-to-Mem ratio is high.
void accel_opt3 (ValueType * __RESTRICT pos, ValueType * __RESTRICT vel, ValueType * __RESTRICT mass, ValueType * __RESTRICT acc, const int n)
{
#ifdef Enable_Alignment
__assume_aligned(pos, Alignment());
__assume_aligned(acc, Alignment());
__assume_aligned(mass, Alignment());
__assume(n % (Alignment()/sizeof(ValueType)) == 0 );
#endif

   const int tileSize = 1024;

   for (int j_start = 0; j_start < n; j_start += tileSize)
   {
      const int j_stop = std::min(n, j_start + tileSize);

      for (int i = 0; i < n; ++i)
      {
         ValueType ax = 0, ay = 0, az = 0;
         const ValueType xi = pos_array(i,0);
         const ValueType yi = pos_array(i,1);
         const ValueType zi = pos_array(i,2);

         for (int j = j_start; j < j_stop; ++j)
         {
            /* Position vector from i to j and the distance^2. */
            ValueType rx = pos_array(j,0) - xi;
            ValueType ry = pos_array(j,1) - yi;
            ValueType rz = pos_array(j,2) - zi;
            //ValueType rx = pos_array(j,0) - pos_array(i,0);
            //ValueType ry = pos_array(j,1) - pos_array(i,1);
            //ValueType rz = pos_array(j,2) - pos_array(i,2);
            ValueType r2 = rx*rx + ry*ry + rz*rz + TINY2;
            ValueType m_invR3 = mass[j] / (r2 * std::sqrt(r2));

            ax += rx * m_invR3;
            ay += ry * m_invR3;
            az += rz * m_invR3;
         }

         if (j_start > 0)
         {
            acc_array(i,0) += G * ax;
            acc_array(i,1) += G * ay;
            acc_array(i,2) += G * az;
         }
         else
         {
            acc_array(i,0) = G * ax;
            acc_array(i,1) = G * ay;
            acc_array(i,2) = G * az;
         }
      }
   }
}

//#pragma omp declare simd simdlen(32) uniform(xi,yi,zi) notinbranch
void accel_opt4_interaction (const ValueType xi, const ValueType yi, const ValueType zi,
                             const ValueType xj, const ValueType yj, const ValueType zj, const ValueType mj,
                             ValueType &ax, ValueType &ay, ValueType &az)
{
   /* Position vector from i to j and the distance^2. */
   const ValueType rx = xj - xi;
   const ValueType ry = yj - yi;
   const ValueType rz = zj - zi;
   const ValueType r2 = rx*rx + ry*ry + rz*rz + TINY2;
   const ValueType m_invR3 = mj / (r2 * std::sqrt(r2));

   ax += rx * m_invR3;
   ay += ry * m_invR3;
   az += rz * m_invR3;
}

// SIMD (vector) processing: Compute several interactions at once with SIMD operations.
void accel_opt4 (ValueType * __RESTRICT pos, ValueType * __RESTRICT vel, ValueType * __RESTRICT mass, ValueType * __RESTRICT acc, const int n)
{
#ifdef Enable_Alignment
#warning 'Enabled __assume_aligned statement'
__assume_aligned(pos, Alignment());
__assume_aligned(acc, Alignment());
__assume_aligned(mass, Alignment());
__assume(n % (Alignment()/sizeof(ValueType)) == 0 );
#endif

   // Outer loop vectorization.
   #pragma omp simd safelen(32) aligned(pos, mass, acc : Alignment())
   for (int i = 0; i < n; ++i)
   {
      ValueType ax = 0, ay = 0, az = 0;
      const ValueType xi = pos_array(i,0);
      const ValueType yi = pos_array(i,1);
      const ValueType zi = pos_array(i,2);

      // Inner loop vectorization.
      //#pragma omp simd safelen(32) aligned(pos, mass : Alignment())
      for (int j = 0; j < n; ++j)
         accel_opt4_interaction( xi, yi, zi,
                                 pos_array(j,0), pos_array(j,1), pos_array(j,2), mass[j],
                                 ax, ay, az );

      acc_array(i,0) = G * ax;
      acc_array(i,1) = G * ay;
      acc_array(i,2) = G * az;
   }
}

#define GCC_VERSION (__GNUC__ * 10000 \
                   + __GNUC_MINOR__ * 100 \
                   + __GNUC_PATCHLEVEL__)

#if defined(__INTEL_COMPILER) || GCC_VERSION > 40900

#include <immintrin.h>
// Explicit SIMD (vector) processing: Compute several interactions at once with SIMD operations.
void accel_opt5 (ValueType * __RESTRICT pos, ValueType * __RESTRICT vel, ValueType * __RESTRICT mass, ValueType * __RESTRICT acc, const int n)
{
#if defined(__INTEL_COMPILER) && defined(Enable_Alignment)
__assume_aligned(pos, Alignment());
__assume_aligned(acc, Alignment());
__assume_aligned(mass, Alignment());
__assume(n % (Alignment()/sizeof(ValueType)) == 0 );
#endif

#ifdef USE_FLOAT
   typedef __m256 vectype;
   const int veclen = 8;
   #define vecinstr(name) name ## _ps
#else
   typedef __m256d vectype;
   const int veclen = 4;
   #define vecinstr(name) name ## _pd
#endif

   // Must be aligned to 32 bytes.
   if (isAligned( pos, 32 ) == false ||
       isAligned( vel, 32 ) == false ||
       isAligned( acc, 32 ) == false ||
       isAligned( mass, 32 ) == false)
   {
      fprintf(stderr,"Data is not aligned to 32-byte boundary: %p %p %p %p\n", pos, vel, acc, mass);
      exit(1);
   }

   // N must be a multiple of the vector width since there's not remainder loop here.
   if (n % veclen != 0)
   {
      fprintf(stderr,"Number of particles (%d) not a multiple of the SIMD width (%d)\n", n, veclen);
      exit(1);
   }

   // Outer loop vectorization.
   for (int i = 0; i < n; i += veclen)
   {
      vectype ax = vecinstr(_mm256_setzero)();
      vectype ay = vecinstr(_mm256_setzero)();
      vectype az = vecinstr(_mm256_setzero)();

      const vectype xi = vecinstr(_mm256_load)( &pos_array(i,0) );
      const vectype yi = vecinstr(_mm256_load)( &pos_array(i,1) );
      const vectype zi = vecinstr(_mm256_load)( &pos_array(i,2) );

      for (int j = 0; j < n; ++j)
      {
         // Set all lanes to the same j values.
         const vectype xj = vecinstr(_mm256_set1)( pos_array(j,0) );
         const vectype yj = vecinstr(_mm256_set1)( pos_array(j,1) );
         const vectype zj = vecinstr(_mm256_set1)( pos_array(j,2) );
         const vectype mj = vecinstr(_mm256_set1)( mass[j] );

         // Position vector from i to j and the distance^3.
         const vectype rx = vecinstr(_mm256_sub)(xj, xi);
         const vectype ry = vecinstr(_mm256_sub)(yj, yi);
         const vectype rz = vecinstr(_mm256_sub)(zj, zi);

         vectype r2 = vecinstr(_mm256_mul)(rx,rx);
         r2 = vecinstr(_mm256_add)(r2, vecinstr(_mm256_mul)(ry,ry));
         r2 = vecinstr(_mm256_add)(r2, vecinstr(_mm256_mul)(rz,rz));
         r2 = vecinstr(_mm256_add)(r2, vecinstr(_mm256_set1)(TINY2));

#if defined(USE_FLOAT) && (0)
         // "approximate 1/sqrt"
         const vectype invR = vecinstr(_mm256_rsqrt)(r2);
         const vectype m_invR3 = vecinstr(_mm256_mul)(mj, vecinstr(_mm256_mul)(invR, vecinstr(_mm256_mul)(invR,invR)));
#else
         // "exact 1/sqrt"
         const vectype r3 = vecinstr(_mm256_mul)(r2, vecinstr(_mm256_sqrt)(r2));
         const vectype m_invR3 = vecinstr(_mm256_div)(mj, r3);
#endif

         ax = vecinstr(_mm256_add)(ax, vecinstr(_mm256_mul)(rx, m_invR3));
         ay = vecinstr(_mm256_add)(ay, vecinstr(_mm256_mul)(ry, m_invR3));
         az = vecinstr(_mm256_add)(az, vecinstr(_mm256_mul)(rz, m_invR3));
      }

      vectype vG = vecinstr(_mm256_set1)(G);
      vecinstr(_mm256_store)( &acc_array(i,0), vecinstr(_mm256_mul)(vG, ax));
      vecinstr(_mm256_store)( &acc_array(i,1), vecinstr(_mm256_mul)(vG, ay));
      vecinstr(_mm256_store)( &acc_array(i,2), vecinstr(_mm256_mul)(vG, az));
   }
}
#else
#warning 'SIMD AVX intrinsics not supported by this compiler.'
#endif

void update (ValueType pos[], ValueType vel[], ValueType mass[], ValueType acc[], const int n, ValueType h)
{
   for (int i = 0; i < n; ++i)
      for (int k = 0; k < NDIM; ++k)
      {
         pos_array(i,k) += vel_array(i,k)*h + acc_array(i,k)*h*h/2;
         vel_array(i,k) += acc_array(i,k)*h;
      }
}

void output (ValueType pos[], ValueType vel[], ValueType mass[], ValueType acc[], const int n, int flnum)
{
   char flname[20];
   sprintf (flname, "pos_%d.out", flnum);
   FILE *fp = fopen(flname,"w");
   if (!fp)
   {
      fprintf(stderr,"Error opening file %s\n", flname);
      exit(-1);
   }

   fwrite (&n, sizeof(int), 1, fp);
   for (int i = 0; i < n; ++i)
   {
      for (int k = 0; k < NDIM; ++k)
      {
         fwrite (&pos_array(i,k), sizeof(ValueType), 1, fp);
      }
      fwrite (&mass[i], sizeof(ValueType), 1, fp);
   }

   fclose(fp);
}

void search (ValueType pos[], ValueType vel[], ValueType mass[], ValueType acc[], const int n)
{
   ValueType minv = 1e10, maxv = 0, ave = 0;
   for (int i = 0; i < n; ++i)
   {
      ValueType vmag = 0;
      for (int k = 0; k < NDIM; ++k)
         vmag += (vel_array(i,k) * vel_array(i,k));

      vmag = sqrt(vmag);

      maxv = std::max(maxv, vmag);
      minv = std::min(minv, vmag);
      ave += vmag;
   }
   printf("min/max/ave velocity = %e, %e, %e\n", minv, maxv, ave / n);
}

void help()
{
   fprintf(stderr,"nbody3 --help|-h --nparticles|-n --nsteps|-s --stepsize|-t\n");
}

int main (int argc, char* argv[])
{
   /* Define the number of particles. The default is 100. */
   int n = 100;

   /* Define the number of steps to run. The default is 100. */
   int num_steps = 100;

   /* Pick the timestep size. */
   ValueType dt = 0.01;

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
      else if (strcmp(argv[i],"--nparticles") == 0 || strcmp(argv[i],"-n") == 0)
      {
         check_index(i+1,"--nparticles|-n");
         i++;
         if (isdigit(*argv[i]))
            n = atoi( argv[i] );
      }
      else if (strcmp(argv[i],"--nsteps") == 0 || strcmp(argv[i],"-s") == 0)
      {
         check_index(i+1,"--nsteps|-s");
         i++;
         if (isdigit(*argv[i]))
            num_steps = atoi( argv[i] );
      }
      else if (strcmp(argv[i],"--stepsize") == 0 || strcmp(argv[i],"-t") == 0)
      {
         check_index(i+1,"--stepsize|-t");
         i++;
         if (isdigit(*argv[i]) || *argv[i] == '.')
            dt = atof( argv[i] );
      }
      else
      {
         fprintf(stderr,"Unknown option %s\n", argv[i]);
         help();
         return 1;
      }
   }

   fprintf(stderr,"Number Objects = %d\n", n);
   fprintf(stderr,"Number Steps   = %d\n", num_steps);
   fprintf(stderr,"Timestep size  = %g\n", dt);
   fprintf(stderr,"Alignment      = %lu bytes\n", Alignment());
   fprintf(stderr,"ValueType      = %s\n", (sizeof(ValueType)==sizeof(double)) ? "double" : "float");
#ifdef Enable_ArrayOfStructures
   fprintf(stderr,"Format         = ArrayOfStructures\n");
#else
   fprintf(stderr,"Format         = StructureOfArrays\n");
#endif

#define _TOSTRING(s) #s
#define TOSTRING(s) _TOSTRING(s)
#ifndef ACC_FUNC
#  define ACC_FUNC accel_opt2
#endif
   fprintf(stderr,"Accel function = %s\n", TOSTRING(ACC_FUNC) );

   ValueType *pos = NULL;
   ValueType *vel = NULL;
   ValueType *acc = NULL;
   ValueType *mass = NULL;

   Allocate(pos, n*NDIM);
   Allocate(vel, n*NDIM);
   Allocate(acc, n*NDIM);
   Allocate(mass, n);

   if (1 && n == 2)
   {
      /* Initialize a 2-body problem with large mass ratio and tangential
       * velocity for the small body. */

      pos_array(0,0) = 0.0; pos_array(0,1) = 0.0; pos_array(0,2) = 0.0;
      vel_array(0,0) = 0.0; vel_array(0,1) = 0.0; vel_array(0,2) = 0.0;
      mass[0] = 1000.0;

      ValueType vy = std::sqrt(G*mass[0]);
      pos_array(1,0) = 1.0; pos_array(1,1) = 0.0; pos_array(1,2) = 0.0;
      vel_array(1,0) = 0.0; vel_array(1,1) =  vy; vel_array(1,2) = 0.0;
      mass[1] = 1.0;
   }
   else
   {
      /* Initialize the positions and velocities with random numbers (0,1]. */

      /* 1. Seed the pseudo-random generator. */
      srand(n);

      for (int i = 0; i < n; ++i)
      {
         /* 2. Set some random positions for each object {-1,1}. */
         for (int k = 0; k < NDIM; ++k)
            pos_array(i,k) = 2*(frand() - 0.5);

         /* 3. Set some random velocity (or zero). */
         for (int k = 0; k < NDIM; ++k)
            vel_array(i,k) = 0;
            //vel_array(i,k) = frand();

         /* 4. Set a random mass (> 0). */
         mass[i] = frand() + TINY;

         for (int k = 0; k < NDIM; ++k)
            acc_array(i,k) = 0;
      }
   }

   /* Run the step several times. */
   myTimer_t t_start = getTimeStamp();
   double t_accel = 0, t_update = 0, t_search = 0;
   int flnum = 0;
   for (int step = 0; step < num_steps; ++step)
   {
      /* 1. Compute the acceleration on each object. */
      myTimer_t t0 = getTimeStamp();

      ACC_FUNC( pos, vel, mass, acc, n );
      //accel_opt2( pos, vel, mass, acc, n );
      //accel_block( pos, vel, mass, acc, n );

      myTimer_t t1 = getTimeStamp();

      /* 2. Advance the position and velocities. */
      update( pos, vel, mass, acc, n, dt );

      myTimer_t t2 = getTimeStamp();

      /* 3. Find the faster moving object. */
      if (step % 10 == 0)
         search( pos, vel, mass, acc, n );

      myTimer_t t3 = getTimeStamp();

      t_accel += getElapsedTime(t0,t1);
      t_update += getElapsedTime(t1,t2);
      t_search += getElapsedTime(t2,t3);

      /* 4. Write positions. */
      if (false && (step % 1 == 0))
      {
         for (int i = 0; i < n; ++i)
         {
            for (int k = 0; k < NDIM; ++k)
               fprintf(stderr,"%f ", pos_array(i,k));
            fprintf(stderr,"%f ", mass[i]);
         }
         fprintf(stderr,"\n");
         //output (pos, vel, mass, acc, n, flnum); flnum++;
      }
   }
   double t_calc = getElapsedTime( t_start, getTimeStamp());

   float nkbytes = (float)((size_t)7 * sizeof(ValueType) * (size_t)n) / 1024.0f;
   //printf("Average time = %f (ms) per step with %d elements %.2f KB over %d steps %.3f%%, %.3f%%, %.3f%%\n", t_calc*1000.0/num_steps, n, nkbytes, num_steps, 100*t_accel/t_calc, 100*t_update/t_calc, 100*t_search/t_calc);
   printf("Average time = %f (ms) per step with %d elements %.2f KB over %d steps %f %f %f\n", t_calc*1000.0/num_steps, n, nkbytes, num_steps, t_accel*1000/num_steps, t_update*1000/num_steps, t_search*1000/num_steps);

   /*fclose(fp);*/

   /* Print out the positions (if not too large). */
   if (n < 50)
   {
      for (int i = 0; i < n; ++i)
      {
         for (int k = 0; k < NDIM; ++k)
            fprintf(stderr,"%f ", pos_array(i,k));
         for (int k = 0; k < NDIM; ++k)
            fprintf(stderr,"%f ", vel_array(i,k));

         fprintf(stderr,"%f\n", mass[i]);
      }
   }

   Deallocate(pos);
   Deallocate(vel);
   Deallocate(acc);
   Deallocate(mass);

   return 0;
}
