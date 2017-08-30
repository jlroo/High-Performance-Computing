COMMON	= ../common

#CXX	= icpc
CXX	?= g++
CXXFLAGS= -O3 -I. -I$(COMMON)
#LDFLAGS	= -lrt -lmkl_rt -L$(TACC_MKL_LIB)
LDFLAGS	= -lrt

ifneq ($(TACC_MKL_LIB),)
LDFLAGS	+= -lmkl_rt -L$(TACC_MKL_LIB)
else
LDFLAGS	+= -lblas -L/usr/lib64
endif

ifeq ($(CXX),icpc)
CXXFLAGS += -xHost #-no-vec
CXXFLAGS += -opt-report=5
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

ifneq ($(dtype),)
ifneq ($(dtype),0)
CXXFLAGS += -D__DTYPE=$(dtype)
endif
endif

EXEC = cache matmul

all: $(EXEC)

matmul: matmul.cpp my_timer.o
	$(CXX) $(CXXFLAGS) -o matmul matmul.cpp my_timer.o $(LDFLAGS)

cache: cache.cpp my_timer.o
	$(CXX) $(CXXFLAGS) -o cache cache.cpp my_timer.o $(LDFLAGS)

my_timer.o: $(COMMON)/my_timer.h
	$(CXX) $(CXXFLAGS) -c $(COMMON)/my_timer.cpp -o my_timer.o

clean:
	/bin/rm -fv $(EXEC) *.o *.optrpt
