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

// HW3: Include the main MPI header.
//#include <mpi.h>

// HW3: Include this helper header to define the callMPI() macro.
// callMPI( MPI_Yxxx ) will call an MPI library function and test
// the return handle for errors. If errors are present, it will
// attempt to decode the error and print a useful message and then
// exit.

#include <my_mpi_header.h>
#ifndef __RESTRICT
#  define __RESTRICT
#endif

#define NDIM (3)

// HW3: MPI_FLOAT and MPI_DOUBLE are the MPI_Datatype parameters
// that you should pass to the MPI_* communication functions. The
// code is configured so that you can switch between double's and
// float's (if you ever wanted to), I created a special parameter
// MPI_VALUE_TYPE that takes the appropriate value depending on
// the desired datatype. So, use MPI_VALUE_TYPE instead of MPI_FLOAT
// or MPI_DOUBLE directly to ensure the code doesn't break if the
// datatype is changed.
#ifdef USE_FLOAT
typedef float  ValueType;
const MPI_Datatype MPI_VALUE_TYPE = MPI_FLOAT;
#else
typedef double ValueType;
const MPI_Datatype MPI_VALUE_TYPE = MPI_DOUBLE;
#endif

const ValueType G = 1.0;
const ValueType TINY = std::numeric_limits<ValueType>::epsilon();
const ValueType TINY2 = TINY*TINY;
const int ROOT = 0;
#define Enable_ArrayOfStructures
#if defined(Enable_ArrayOfStructures) || defined(__AOS)
#  ifndef Enable_ArrayOfStructures
#    define Enable_ArrayOfStructures
#  endif
/* Array-of-structures (like) format. */
#  define index(i,j) (NDIM*(i) + (j))
#else
/* Structure-of-arrays (like) format. */
#  define index(i,j) ((i) + (j)*n)
#endif

// HW3: This is a useful function for partitioning an iteration
// space uniforming w/o any overlap.
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

/* Generate a random double between 0,1. */
ValueType frand(void) { return ((ValueType) rand()) / RAND_MAX; }

void accel (ValueType * __RESTRICT pos,
            ValueType * __RESTRICT vel,
            ValueType * __RESTRICT mass,
            ValueType * __RESTRICT acc,
            const int n)
{
    int myRank, numProcs;
    callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
    callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );
    int partition_start, partition_end;
    partition_range( 0, n, numProcs, myRank, partition_start, partition_end );
    for (int i = partition_start; i < partition_end; ++i)
    {
        ValueType ax = 0, ay = 0, az = 0;
        const ValueType xi = pos[ index(i,0) ];
        const ValueType yi = pos[ index(i,1) ];
        const ValueType zi = pos[ index(i,2) ];
        
        for (int j = 0; j < n; ++j)
        {
            /* Position vector from i to j and the distance^2. */
            ValueType rx = pos[ index(j,0) ] - xi;
            ValueType ry = pos[ index(j,1) ] - yi;
            ValueType rz = pos[ index(j,2) ] - zi;
            ValueType dsq = rx*rx + ry*ry + rz*rz + TINY2;
            ValueType m_invR3 = mass[j] / (dsq * std::sqrt(dsq));
            
            ax += rx * m_invR3;
            ay += ry * m_invR3;
            az += rz * m_invR3;
        }
        
        acc[ index(i,0) ] = G * ax;
        acc[ index(i,1) ] = G * ay;
        acc[ index(i,2) ] = G * az;
    }
}

void update (ValueType pos[],
             ValueType vel[],
             ValueType mass[],
             ValueType acc[],
             const int n,
             ValueType h)
{
    int myRank, numProcs;
    callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
    callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );
    int partition_start, partition_end;
    partition_range( 0, n, numProcs, myRank, partition_start, partition_end );
    // HW3: Split up the iteration space.
    for (int i = partition_start; i < partition_end; ++i)
        for (int k = 0; k < NDIM; ++k)
        {
            pos[ index(i,k)] += vel[ index(i,k) ]*h + acc[ index(i,k) ]*h*h/2;
            vel[ index(i,k)] += acc[ index(i,k) ]*h;
        }
}

void search (ValueType pos[],
             ValueType vel[],
             ValueType mass[],
             ValueType acc[],
             const int n)
{
    
    int myRank, numProcs;
    callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
    callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );
    int partition_start, partition_end;
    partition_range( 0, n, numProcs, myRank, partition_start, partition_end );
    ValueType minv = 1e10, maxv = 0, ave = 0;
    // HW3: Split up the iteration space.
    
    for (int i = partition_start; i < partition_end; ++i)
    {
        ValueType vmag = 0;
        for (int k = 0; k < NDIM; ++k)
            vmag += (vel[ index(i,k) ] * vel[ index(i,k) ] );
        
        vmag = sqrt(vmag);
        maxv = std::max(maxv, vmag);
        minv = std::min(minv, vmag);
        ave += vmag;
    }
    
    ValueType g_maxv;
    callMPI(MPI_Reduce(&maxv, &g_maxv, 1, MPI_VALUE_TYPE, MPI_MAX, ROOT, MPI_COMM_WORLD));
    ValueType g_minv;
    callMPI(MPI_Reduce(&minv, &g_minv, 1, MPI_VALUE_TYPE, MPI_MIN, ROOT, MPI_COMM_WORLD));
    ValueType g_ave;
    callMPI(MPI_Reduce(&ave, &g_ave, 1, MPI_VALUE_TYPE, MPI_SUM, ROOT, MPI_COMM_WORLD));
    
    if (myRank==ROOT) {
        printf("min/max/ave velocity = %e, %e, %e\n", g_minv, g_maxv, g_ave/n);
    }
    callMPI(MPI_Barrier(MPI_COMM_WORLD));
}

void help()
{
    fprintf(stderr,"nbody3 --help|-h --nparticles|-n --nsteps|-s --stepsize|-t\n");
}

int main (int argc, char* argv[])
{
    // Initialize the MPI system before anything else.
    callMPI(MPI_Init(&argc,&argv)                   );
    //const int ROOT = 0;
    
    // Get the rank and size of the MPI world.
    int myRank, numProcs;
    callMPI( MPI_Comm_size(MPI_COMM_WORLD,&numProcs) );
    callMPI( MPI_Comm_rank(MPI_COMM_WORLD,&myRank)   );
    
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
    
    ValueType *pos = NULL;
    ValueType *vel = NULL;
    ValueType *acc = NULL;
    ValueType *mass = NULL;

    Allocate(pos, n*NDIM);
    Allocate(vel, n*NDIM);
    Allocate(acc, n*NDIM);
    Allocate(mass, n);
    
    if (myRank==ROOT) {
        // HW3: Only the master (0) rank should write these lines.
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
        // 1. Seed the pseudo-random generator.
        srand(n);
        
        for (int i = 0; i < n; ++i) {
            // 2. Set some random positions for each object {-1,1}.
            for (int k = 0; k < NDIM; ++k)
                pos[index(i,k)] = 2*(frand() - 0.5);
            
            // 3. Set some random velocity (or zero).
            for (int k = 0; k < NDIM; ++k)
                vel[index(i,k)] = 0;
            //vel[index(i,k)] = frand();
            
            // 4. Set a random mass (> 0).
            mass[i] = frand() + TINY;
            
            for (int k = 0; k < NDIM; ++k)
                acc[index(i,k)] = 0;
        }
    }

    callMPI(MPI_Bcast(pos, n*3, MPI_VALUE_TYPE, ROOT, MPI_COMM_WORLD));
    callMPI(MPI_Bcast(vel, n*3, MPI_VALUE_TYPE, ROOT, MPI_COMM_WORLD));
    callMPI(MPI_Bcast(mass, n, MPI_VALUE_TYPE, ROOT, MPI_COMM_WORLD));
    callMPI(MPI_Bcast(acc, n*3, MPI_VALUE_TYPE, ROOT, MPI_COMM_WORLD));


    // Run the step several times.
    myTimer_t t_start = getTimeStamp();
    
    // HW3: Add another time to measure the total time in the MPI library.
    double t_accel = 0, t_update = 0, t_search = 0, t_mpi;
    int flnum = 0;
    for (int step = 0; step < num_steps; ++step)
    {
        /* 1. Compute the acceleration on each object. */
        myTimer_t t0 = getTimeStamp();
        
        accel( pos, vel, mass, acc, n);
        
        myTimer_t t1 = getTimeStamp();
        
        /* 2. Advance the position and velocities. */
        update( pos, vel, mass, acc, n, dt);
        
        myTimer_t t2 = getTimeStamp();
        
        /* 3. Find the faster moving object. */
        if (step % 10 == 0)
            search( pos, vel, mass, acc, n);
        
        myTimer_t t3 = getTimeStamp();
        

        std::vector<int> counts(numProcs);
        std::vector<int> displ(numProcs);
        int partition_ranges;
        for (int i=0; i<numProcs; i++) {
            int pstart, pend;
            partition_range( 0, n, numProcs, myRank, pstart, pend );
            counts[i]=(pend-pstart)*3;
        }

        displ[0] = 0;
        for (int i = 1; i < numProcs; ++i)
            displ[i] = displ[i-1] + counts[i-1];
        
        
        callMPI(MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                pos, &counts[0], &displ[0], MPI_VALUE_TYPE,
                                MPI_COMM_WORLD ));
        
        myTimer_t t4 = getTimeStamp();
        
        t_accel += getElapsedTime(t0,t1);
        t_update += getElapsedTime(t1,t2);
        t_search += getElapsedTime(t2,t3);
        t_mpi += getElapsedTime(t3,t4);
    }
    double t_calc = getElapsedTime( t_start, getTimeStamp());
    
    float nkbytes = (float)((size_t)7 * sizeof(ValueType) * (size_t)n) / 1024.0f;
    
    // HW3: Only the master (0) rank should write the this line.
    // HW3: Extra Credit: Find the average / min / max function times for
    // each rank and see if there's a signficant variation.
    
    /*
    for (int i = 0; i < numProcs; ++i)
    {
        if (myRank == i)
            printf("I'm rank %d\n", i);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    */
    
    
    if (myRank==ROOT) {
        printf("Root[%d]: Average time = %f (ms) per step with %d elements, %.2f KB over %d steps\n", myRank, t_calc*1000.0/num_steps, n, nkbytes, num_steps);

    printf("Root[%d]: accel-time[%f] update-time[%f] search-time[%f] mpi-time[%f]\n", myRank,t_accel*1000/num_steps, t_update*1000/num_steps, t_search*1000/num_steps, t_mpi*1000/num_steps);
	}

//    callMPI(MPI_Barrier(MPI_COMM_WORLD));
    // reduction per rank
    
//    printf("Rank[%d]: Average time = %f (ms) per step with %d elements, %.2f KB over %d steps\n", myRank, t_calc*1000.0/num_steps, n, nkbytes, num_steps);
//    printf("Rank[%d]: accel-time[%f] update-time[%f] search-time[%f] mpi-time[%f]\n", myRank,t_accel*1000/num_steps, t_update*1000/num_steps, t_search*1000/num_steps, t_mpi*1000/num_steps);
    
    Deallocate(pos);
    Deallocate(vel);
    Deallocate(acc);
    Deallocate(mass);

    // Close the MPI library.
    callMPI( MPI_Finalize() );

    return 0;
}

