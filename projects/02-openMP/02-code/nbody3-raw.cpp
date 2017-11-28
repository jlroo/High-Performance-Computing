#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <vector>
#include <string>
#include <limits>
#include <cmath>

#include <my_timer.h>
#include <aligned_allocator.h>

#ifndef __RESTRICT
#  define __RESTRICT
#endif

#define NDIM (3)

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

template <typename ValueType>
constexpr ValueType square(const ValueType& x) { return x*x; }

#define G ( ValueType(1) )
#define TINY ( std::numeric_limits<ValueType>::epsilon() )
#define TINY2 ( TINY * TINY )

/* Generate a random double between 0,1. */
template <typename ValueType>
ValueType frand(void) { return ValueType( rand() ) / RAND_MAX; }

// Store target data in registers: Compiler "may" do this automatically but
// it often helps with cache efficiency. This can be especially helpfule
// by avoiding repeated writes which are several times slower than reads.
template <typename ValueType>
void accel_register (ValueType * __RESTRICT pos, ValueType * __RESTRICT vel, ValueType * __RESTRICT mass, ValueType * __RESTRICT acc, const int n)
{
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

template <typename ValueType>
void update (ValueType pos[], ValueType vel[], ValueType mass[], ValueType acc[], const int n, ValueType h)
{
   for (int i = 0; i < n; ++i)
      for (int k = 0; k < NDIM; ++k)
      {
         pos_array(i,k) += vel_array(i,k)*h + acc_array(i,k)*h*h/2;
         vel_array(i,k) += acc_array(i,k)*h;
      }
}

template <typename ValueType>
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

template <typename ValueType>
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

void help(const char* prg)
{
   if (prg) fprintf(stderr,"%s:\n", prg);
   fprintf(stderr,"\t--help | -h       : Print help message.\n");
   fprintf(stderr,"\t--nparticles | -n : # of particles (100).\n");
   fprintf(stderr,"\t--nsteps | -s     : # of steps to take (100).\n");
   fprintf(stderr,"\t--stepsize | -dt  : Delta-t step-size in seconds (0.01).\n");
   fprintf(stderr,"\t--float | -f      : Use 32-bit floats.\n");
   fprintf(stderr,"\t--double | -d     : Use 64-bit doubles. (default)\n");
}

template <typename ValueType>
int run_tests( const int n, const int num_steps, const ValueType dt)
{
   fprintf(stderr,"Number Objects = %d\n", n);
   fprintf(stderr,"Number Steps   = %d\n", num_steps);
   fprintf(stderr,"Timestep size  = %g\n", dt);
   fprintf(stderr,"Alignment      = %lu bytes\n", Alignment());
   fprintf(stderr,"ValueType      = %s\n", (sizeof(ValueType)==sizeof(double)) ? "double" : "float");
#ifdef __ENABLE_VCL_SIMD
   fprintf(stderr,"simdType       = %s\n", (sizeof(ValueType)==sizeof(double)) ? STRINGIFY(VCL_SIMD_DBL) : STRINGIFY(VCL_SIMD_FLT));
# ifdef __ENABLE_VCL_RSQRT
   fprintf(stderr,"simdType rsqrt = enabled\n");
# endif
#endif
#ifdef Enable_ArrayOfStructures
   fprintf(stderr,"Format         = ArrayOfStructures\n");
#else
   fprintf(stderr,"Format         = StructureOfArrays\n");
#endif

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
            pos_array(i,k) = 2*(frand<ValueType>() - 0.5);

         /* 3. Set some random velocity (or zero). */
         for (int k = 0; k < NDIM; ++k)
            vel_array(i,k) = 0;
            //vel_array(i,k) = frand();

         /* 4. Set a random mass (> 0). */
         mass[i] = frand<ValueType>() + TINY;

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

      accel_register( pos, vel, mass, acc, n );

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

int main (int argc, char* argv[])
{
   /* Define the number of particles. The default is 100. */
   int n = 100;

   /* Define the number of steps to run. The default is 100. */
   int num_steps = 100;

   /* Pick the timestep size. */
   double dt = 0.01;

   /* ValueType? (float or double) */
   bool useDouble = true;

   for (int i = 1; i < argc; ++i)
   {
#define check_index(i,str) \
   if ((i) >= argc) \
      { fprintf(stderr,"Missing 2nd argument for %s\n", str); help(argv[0]); return 1; }

      if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
      {
         help(argv[0]);
         return 0;
      }
      else if (strcmp(argv[i],"--nparticles") == 0 || strcmp(argv[i],"-n") == 0)
      {
         check_index(i+1,"--nparticles|-n");
         i++;
         if (not(isdigit(*argv[i])))
            { fprintf(stderr,"Invalid value for option \"--particles\" %s\n", argv[i]); help(argv[0]); return 1; }
         n = atoi( argv[i] );
      }
      else if (strcmp(argv[i],"--nsteps") == 0 || strcmp(argv[i],"-s") == 0)
      {
         check_index(i+1,"--nsteps|-s");
         i++;
         if (not(isdigit(*argv[i])))
            { fprintf(stderr,"Invalid value for option \"--nsteps\" %s\n", argv[i]); help(argv[0]); return 1; }
         num_steps = atoi( argv[i] );
      }
      else if (strcmp(argv[i],"--stepsize") == 0 || strcmp(argv[i],"-t") == 0)
      {
         check_index(i+1,"--stepsize|-t");
         i++;
         if (not(isdigit(*argv[i])))
            { fprintf(stderr,"Invalid value for option \"--stepsize\" %s\n", argv[i]); help(argv[0]); return 1; }
         dt = atof( argv[i] );
      }
      else if (strcmp(argv[i],"--double") == 0 || strcmp(argv[i],"-d") == 0)
      {
         useDouble = true;
      }
      else if (strcmp(argv[i],"--float") == 0 || strcmp(argv[i],"-f") == 0)
      {
         useDouble = false;
      }
      else
      {
         fprintf(stderr,"Unknown option %s\n", argv[i]);
         help(argv[0]);
         return 1;
      }
   }

   if (useDouble)
      return run_tests<double>( n, num_steps, dt );
   else
      return run_tests<float >( n, num_steps, dt );

   return 0;
}
