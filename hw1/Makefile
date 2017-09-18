COMMON	= ../common

#CXX	= icpc
CXX    ?= g++
CXXFLAGS= -O3 -I. -I$(COMMON)
LDFLAGS	= -lrt

ifeq ($(CXX),icpc)
CXXFLAGS += -xHost #-no-vec
CXXFLAGS += -qopt-report=5
CXXFLAGS += -D__ALIGNMENT=32
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

ifneq ($(dtype),)
ifneq ($(dtype),0)
CXXFLAGS += -D__DTYPE=$(dtype)
endif
endif

EXEC = stream2

all: $(EXEC)

stream2: stream2.cpp my_timer.o
	$(CXX) $(CXXFLAGS) -o stream2 stream2.cpp my_timer.o $(LDFLAGS)

stream2_aligned: stream2_aligned.cpp my_timer.o
	$(CXX) $(CXXFLAGS) -o stream2_aligned stream2_aligned.cpp my_timer.o $(LDFLAGS)

my_timer.o: $(COMMON)/my_timer.h
	$(CXX) $(CXXFLAGS) -c $(COMMON)/my_timer.cpp -o my_timer.o

clean:
	/bin/rm -fv $(EXEC) *.o *.optrpt
