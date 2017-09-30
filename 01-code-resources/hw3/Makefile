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
