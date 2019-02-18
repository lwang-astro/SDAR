## TARGET
TARGET=ar ptclgrouptest

VPATH=../src/AR ./include

use_x86=yes
#use_k_computer=yes
QD_PATH=/opt/qd-2.3.18
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
#CXXFLAGS += -D AR_DEBUG
CXXFLAGS += -D AR_DEEP_DEBUG

## Flag ----------------------------------------------#
CXXFLAGS += -I../src/AR -I./include
CXXFLAGS += -std=c++11
#CXXFLAGS += -D AR_TTL

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


## Target --------------------------------------------#

$(TARGET): %: %.cxx 
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

clean:
	rm -f $(TARGET)
