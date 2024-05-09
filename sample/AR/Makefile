TARGET= ar.logh ar.logh.sd.a ar.logh.sd.t ar.ttl ar.ttl.sd.a ar.ttl.sd.t ar.ttl.sd.t.mulpot ar.ttl.sd.t.maxpot ar.ttl.sd.t.addpot ar.ttl.sd.t.cm

INSTALL_PATH=~/bin

all: $(TARGET)

VPATH=../../src 

SOURCE=${shell ls ../../src/AR/*.h} ${shell ls ../../src/Common/*.h} ${shell ls ./*.h}

use_x86=yes
#use_k_computer=yes
#QD_PATH=../../../../
#MPFRC_PATH=../../../mpfrc++
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

## mpfrc++-lib----------------------------------------#
ifneq (x$(MPFRC_PATH),x)
CXXFLAGS += -D USE_MPFRC -I $(MPFRC_PATH)
CXXLIBS += -lmpfr -lgmp
endif

## Target --------------------------------------------#

ar.logh: ar.cxx  ${SOURCE}
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

ar.logh.sd.a: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_SLOWDOWN_ARRAY -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

ar.logh.sd.t: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

ar.ttl: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL $< -o $@ $(CXXLIBS)

ar.ttl.m: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_TTL_GT_MULTI $< -o $@ $(CXXLIBS)

ar.ttl.sd.a: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_ARRAY -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

ar.ttl.sd.t: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

ar.ttl.sd.t.cm: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

ar.ttl.sd.t.mulpot: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_TIME_FUNCTION_MUL_POT $< -o $@ $(CXXLIBS)

ar.ttl.sd.t.maxpot: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_TIME_FUNCTION_MAX_POT $< -o $@ $(CXXLIBS)

ar.ttl.sd.t.addpot: ar.cxx ${SOURCE}
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_TIME_FUNCTION_ADD_POT $< -o $@ $(CXXLIBS)

install: $(TARGET)
	install -m 755 $(TARGET) $(INSTALL_PATH)

clean:
	rm -f $(TARGET)
