COMMON	= ../common

DBG      ?=
CXX      ?= g++
CXXFLAGS  = -O3 -I. -I$(COMMON) $(DBG)

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

ifneq ($(dtype),)
ifeq ($(dtype),float)
CXXFLAGS += -DUSE_FLOAT
endif
endif

ifneq ($(func),)
ifneq ($(func),0)
CXXFLAGS += -DACC_FUNC=$(func)
endif
endif

EXEC = nbody3

all: $(EXEC)

# Load common make options
include $(COMMON)/Makefile.common
LDFLAGS	  = $(COMMON_LIBS)

nbody3: nbody3.o $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o nbody3 $^ $(LDFLAGS)

clean: clean_common
	/bin/rm -fv $(EXEC) *.d *.o *.optrpt
