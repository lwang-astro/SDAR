TARGET= ar.logh ar.logh.sd ar.ttl ar.ttl.sd

all: $(TARGET)

VPATH=../../src 

SOURCE=${shell ls ../../src/AR/*.h} ${shell ls ../../src/Common/*.h} ${shell ls ./*.h}

use_x86=yes
#use_k_computer=yes
#QD_PATH=/opt/qd-2.3.22
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
#CXXFLAGS += -D BINARY_DEBUG
#CXXFLAGS += -D AR_DEBUG_PRINT
#CXXFLAGS += -D AR_WARN
#CXXFLAGS += -D AR_DEEP_DEBUG

## Flag ----------------------------------------------#
CXXFLAGS += -I../../src -I./
CXXFLAGS += -std=c++11

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
CXXFLAGS += -D USE_DD ${shell ${QD_PATH}/bin/qd-config --cxxflags}
CXXLIBS += ${shell ${QD_PATH}/bin/qd-config --libs}
endif


## Target --------------------------------------------#

ar.logh: ar.cxx  ${SOURCE}
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

ar.logh.sd: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_SLOWDOWN_INNER -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

ar.ttl: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL $< -o $@ $(CXXLIBS)

ar.ttl.m: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_TTL_GT_MULTI $< -o $@ $(CXXLIBS)

ar.ttl.sd: ar.cxx ${SOURCE}
#	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_INNER -D AR_TTL_GT_BINARY_INNER -D AR_SLOWDOWN_TIMESCALE $< -o $@ $(CXXLIBS)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_INNER -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

install: $(TARGET)
	mv $(TARGET) ~/bin

clean:
	rm -f $(TARGET)
