COMMON    = ../../common

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

EXEC = fixanalyzer

all: $(EXEC)

include $(COMMON)/Makefile.common
LDFLAGS      = $(COMMON_LIBS)

util.o: util.cpp
	$(CXX) $(CXXFLAGS) -c util.cpp

fixanalyzer: util.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o fixanalyzer util.o $(COMMON_OBJS) $(LDFLAGS)

clean: clean_common
	/bin/rm -fv $(EXEC) *.o *.optrpt

## for t in 1 2 4 8 16 32 64; do export OMP_NUM_THREADS=$t;time ./fixanalyzer -p "/work/05191/jlroo/stampede2/data/01/XCME_MD_ES_20100104_20100108" -t 52 -n 4 -m 8; done


