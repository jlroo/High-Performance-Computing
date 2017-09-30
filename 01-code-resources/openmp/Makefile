COMMON	= ../common

DBG      ?=
CXX      ?= g++
CXXFLAGS  = -O3 -fopenmp -I. -I$(COMMON) $(DBG)

ifeq ($(CXX),icpc)
CXXFLAGS += -xHost #-no-vec
CXXFLAGS += -opt-report=5
CXXFLAGS += -Wunknown-pragmas # Disable warning about OpenMP pragma no defined.
endif

ifeq ($(CXX),g++)
CXXFLAGS += -mtune=native
endif

ifneq ($(restrict),)
ifneq ($(restrict),0)
CXXFLAGS += -D__RESTRICT=restrict -restrict
endif
endif

ifneq ($(align),)
ifneq ($(align),0)
CXXFLAGS += -D__ALIGNMENT=$(align)
endif
endif


EXEC = omp_hello omp_false_sharing omp_for omp_list omp_lock omp_race

OBJS = $(EXEC:=.o)
DEPS = $(OBJS:.o=.d)
-include $(DEPS)

all: $(EXEC)

# Load common make options
include $(COMMON)/Makefile.common
LDFLAGS	= $(COMMON_LIBS)

omp_hello: omp_hello.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o omp_hello $^ $(LDFLAGS)

omp_false_sharing: omp_false_sharing.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o omp_false_sharing $^ $(LDFLAGS)

omp_for: omp_for.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o omp_for $^ $(LDFLAGS)

omp_list: omp_list.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o omp_list $^ $(LDFLAGS)

omp_lock: omp_lock.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o omp_lock $^ $(LDFLAGS)

omp_race: omp_race.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o omp_race $^ $(LDFLAGS)

clean: clean_common
	/bin/rm -fv $(EXEC) *.d *.o *.optrpt
