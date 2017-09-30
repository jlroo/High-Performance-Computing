COMMON	= ../common

DBG      ?=
MPICXX   ?= mpicxx
CXXFLAGS  = -O3 -I. -I$(COMMON) $(DBG)

EXEC = mpi_hello mpi_p2p_examples mpi_coll_examples mpi_task_example mpi_matmul

all: $(EXEC)

OBJS = $(EXEC:=.o)
DEPS = $(OBJS:.o=.d)

-include $(DEPS)

# Load common make options
include $(COMMON)/Makefile.common
LDFLAGS	= $(COMMON_LIBS)

ifneq (${TACC_MKL_DIR},)
LDFLAGS += -mkl -L$(TACC_MKL_LIB)
else
LDFLAGS += -lblas
endif

%.o: %.cpp
	$(MPICXX) $(CXXFLAGS) -c $<
	$(MPICXX) -MM $(CXXFLAGS) $< > $*.d

mpi_hello: mpi_hello.o $(COMMON_OBJS)
	$(MPICXX) $(CXXFLAGS) -o mpi_hello $^ $(LDFLAGS)

mpi_p2p_examples: mpi_p2p_examples.o $(COMMON_OBJS)
	$(MPICXX) $(CXXFLAGS) -o mpi_p2p_examples $^ $(LDFLAGS)

mpi_coll_examples: mpi_coll_examples.o $(COMMON_OBJS)
	$(MPICXX) $(CXXFLAGS) -o mpi_coll_examples $^ $(LDFLAGS)

mpi_task_example: mpi_task_example.o $(COMMON_OBJS)
	$(MPICXX) $(CXXFLAGS) -o mpi_task_example $^ $(LDFLAGS)

mpi_matmul: mpi_matmul.o $(COMMON_OBJS)
	$(MPICXX) $(CXXFLAGS) -o mpi_matmul $^ $(LDFLAGS)

clean: clean_common
	/bin/rm -fv $(EXEC) *.d *.o *.optrpt
