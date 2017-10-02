COMMON	= ../common

CXX	?= g++
CXXFLAGS= -O2 -g -I. -std=c++11 -I$(COMMON)
LDFLAGS	= -lrt

ifneq ($(TACC_MKL_LIB),)
LDFLAGS	+= -lmkl_rt -L$(TACC_MKL_LIB)
else
LDFLAGS	+= -lblas -llapack -L/usr/lib64
endif

ifeq ($(CXX),icpc)
CXXFLAGS += -xHost #-xMIC-AVX512 #-no-vec
CXXFLAGS += -qopt-report=5
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

ifneq ($(debug),)
ifneq ($(debug),0)
CXXFLAGS += -DDEBUG=1
endif
endif

EXEC = cache matmul ludec

all: $(EXEC)

matmul: matmul.cpp my_timer.o
	$(CXX) $(CXXFLAGS) -o matmul matmul.cpp my_timer.o $(LDFLAGS)

cache: cache.cpp my_timer.o
	$(CXX) $(CXXFLAGS) -o cache cache.cpp my_timer.o $(LDFLAGS)

ludec: ludec.cpp my_timer.o
	$(CXX) $(CXXFLAGS) -o ludec ludec.cpp my_timer.o $(LDFLAGS)

my_timer.o: $(COMMON)/my_timer.h
	$(CXX) $(CXXFLAGS) -c $(COMMON)/my_timer.cpp -o my_timer.o

clean:
	/bin/rm -fv $(EXEC) *.o *.optrpt
