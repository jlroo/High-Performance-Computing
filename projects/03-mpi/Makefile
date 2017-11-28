COMMON	= ../common

DBG      ?=
MPICXX   ?= mpicxx
CXXFLAGS  = -O3 -I. -I$(COMMON) $(DBG)

EXEC = mpi_nbody3

all: $(EXEC)

OBJS = $(EXEC:=.o)
DEPS = $(OBJS:.o=.d)

-include $(DEPS)

# Load common make options
include $(COMMON)/Makefile.common
LDFLAGS	= $(COMMON_LIBS)

%.o: %.cpp
	$(MPICXX) $(CXXFLAGS) -c $<
	$(MPICXX) -MM $(CXXFLAGS) $< > $*.d

mpi_nbody3: mpi_nbody3.o $(COMMON_OBJS)
	$(MPICXX) $(CXXFLAGS) -o mpi_nbody3 $^ $(LDFLAGS)

clean: clean_common
	/bin/rm -fv $(EXEC) *.d *.o *.optrpt


#### for i in 1 2 4 8 16 32 64;do for n in 1000 2000 4000 8000 16000;do ibrun -n $i ./mpi_nbody3 -n $n; 2>/dev/null; done > $i-output.log; done;