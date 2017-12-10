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



