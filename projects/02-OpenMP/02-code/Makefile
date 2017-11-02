COMMON	= ../../common

DBG      ?=
CXX      ?= g++
CXXFLAGS  = -O3 -fopenmp -I. -std=c++11 -I$(COMMON) $(DBG)

ifeq ($(CXX),icpc)
  CXXFLAGS += -xHost #-no-vec
  CXXFLAGS += -qopt-report=5
  CXXFLAGS += -Wunknown-pragmas # Disable warning about OpenMP pragma no defined.
endif

ifeq ($(CXX),g++)
  CXXFLAGS += -mtune=native -march=native
endif

EXEC = nbody3

all: $(EXEC)

# Load common make options
include $(COMMON)/Makefile.common
LDFLAGS	  = $(COMMON_LIBS)

nbody3.o: nbody3.cpp
	$(CXX) $(CXXFLAGS) -c nbody3.cpp

nbody3: nbody3.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o nbody3 nbody3.o $(COMMON_OBJS) $(LDFLAGS)

clean: clean_common
	/bin/rm -fv $(EXEC) *.o *.optrpt

# for t in 2 4 6 8 16 32 64; \
# do export OMP_NUM_THREADS=$t; \
# echo "$t"; for n in 1000 2000 4000 8000; \
# do ./serial-nbody3 -n $n -s 500 2>/dev/null; \
# done  > $t-output.log ; done
# 
# # weak
# s=$( printf "%.0f" $(echo "sqrt ( 2 ) * 2000" | bc -l) )
# for i in 1 2 4 8 16 32 64; do var=$( printf "%.0f" $(echo "sqrt ( $i ) * 1000" | bc -l) ); echo $var; done
# for t in 1 2 4 8 16 32 64; \
# do export OMP_NUM_THREADS=$t; \
# n=$( printf "%.0f" $(echo "sqrt ( $t ) * 1000" | bc -l) ); \
# ./parallel-nbody3 -n $n -s 500; 2>/dev/null; \
# done > weak-outputv2.log


