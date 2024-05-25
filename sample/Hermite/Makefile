## TARGET
TARGET= hermite hermite.mpfrc 
all: $(TARGET)

INSTALL_PATH=~/bin

VPATH=../../src 

SOURCE=${shell ls ../../src/*/*.h} ${shell ls ./*.h}

use_x86=yes
#use_k_computer=yes
#QD_PATH=/opt/qd-2.3.22
MPFRC_PATH=../../../mpfrc++
#GPERF_PATH = /opt/gperftools-2.6.1


## compiler-------------------------------------------#
ifeq ($(use_x86),yes)
CXX= g++
#CXXFLAGS += -O2 -Wall
CXXFLAGS += -g -fbounds-check -Wall -O0
endif

ifeq ($(use_k_computer),yes)
CXX = time mpiFCCpx
CXXFLAGS += -Kfast
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -Nfjcex
CXXFLAGS += -Krestp=all
# profiling
CXXFLAGS += -Ntl_trt
CXXFLAGS += -DFAPP_PROFILE
endif

## for debug------------------------------------------#
#CXXFLAGS += -D DEBUG 
CXXFLAGS += -D AR_DEBUG
CXXFLAGS += -D BINARY_DEBUG
CXXFLAGS += -D AR_COLLECT_DS_MODIFY_INFO
#CXXFLAGS += -D AR_DEEP_DEBUG
CXXFLAGS += -D HERMITE_DEBUG
CXXFLAGS += -D ADJUST_GROUP_DEBUG

## Flag ----------------------------------------------#
CXXFLAGS += -I../../src -I./
CXXFLAGS += -std=c++11
CXXFLAGS += -D AR_SLOWDOWN_TIMESCALE
#CXXFLAGS += -D AR_SLOWDOWN_MASSRATIO
CXXFLAGS += -D AR_TTL
CXXFLAGS += -D AR_SLOWDOWN_TREE
#CXXFLAGS += -D AR_TTL_GT_BINARY_INNER

## OpenMP flag----------------------------------------#
#CXXFLAGS += -fopenmp -D USE_OMP

## profile---------------------------------------#
CXXFLAGS += -D AR_PROFILE

## gperf----------------------------------------------#
ifneq (x$(GPERF_PATH),x)
CXXLIBS += -L$(GPERF_PATH)/lib -lprofiler -ltcmalloc
endif

## qd-lib---------------------------------------------#
ifneq (x$(QD_PATH),x)
CXXFLAGS += -D USE_QD ${shell ${QD_PATH}/bin/qd-config --cxxflags}
CXXLIBS += ${shell ${QD_PATH}/bin/qd-config --libs}
endif

## mpfrc++-lib----------------------------------------#
ifneq (x$(MPFRC_PATH),x)
MPFRC_FLAGS += -D USE_MPFRC -I $(MPFRC_PATH)
MPFRC_LIBS += -lmpfr -lgmp
endif


## Target --------------------------------------------#

hermite: hermite.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

hermite.mpfrc: hermite.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) $(MPFRC_FLAGS) $< -o $@ $(CXXLIBS) $(MPFRC_LIBS)

install: $(TARGET)
	install -m 755 $(TARGET) $(INSTALL_PATH)

clean:
	rm -f $(TARGET)
