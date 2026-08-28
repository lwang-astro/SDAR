TARGET= ar.logh ar.logh.sd ar.logh.ttl ar.logh.ttl.sd ar.blogh.ttl.sd ar.blogh.ttl.sd.cm ar.normblogh.ttl.sd ar.normblogh.ttl.sd.cm ar.mulall.ttl.sd ar.mulall.ttl.sd.cm ar.btlogh.ttl.sd ar.btlogh.ttl.sd.cm ar.maxpot.ttl.sd.cm ar.addpot.ttl.sd.cm ar.logh.sd.kdkpert ar.logh.ttl.sd.kdkpert
#TARGET= ar.logh ar.logh.ttl ar.logh.sd.kdkpert ar.logh.ttl.sd.kdkpert
TARGET_MPFRC= ar.logh.ttl.sd.cm ar.logh.mpfrc ar.logh.ttl.mpfrc ar.blogh.ttl.sd.cm.mpfrc ar.blogh.ttl.sd.mpfrc ar.btlogh.ttl.sd.cm.mpfrc ar.btlogh.ttl.sd.mpfrc 
INSTALL_PATH=~/bin
BUILD_DIR=./build

VPATH=../../src 

SOURCE=${shell ls ../../src/AR/*.h} ${shell ls ../../src/Common/*.h} ${shell ls ./*.h}

use_x86=yes
use_mpfrc=yes
#use_k_computer=yes
#QD_PATH=../../../../
#MPFRC_PATH=../../../mpreal
#GPERF_PATH = /opt/gperftools-2.6.1


## compiler-------------------------------------------#
ifeq ($(use_x86),yes)
CXX= g++
CXXFLAGS += -O2 -Wall
#CXXFLAGS += -g -fbounds-check -Wall -O0 -D AR_DEBUG
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
CXXFLAGS += -D SDAR_TIME_MEASURE

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
ifeq ($(use_mpfrc),yes)
MPFRC_FLAGS += -D USE_MPFRC# -I $(MPFRC_PATH)
MPFRC_LIBS += -lmpfr -lgmp
TARGET += ${TARGET_MPFRC}
endif


all: $(addprefix $(BUILD_DIR)/, $(TARGET))

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

## Target --------------------------------------------#

$(BUILD_DIR)/ar.logh: ar.cxx  ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.logh.mpfrc: ar.cxx  ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(MPFRC_FLAGS) $< -o $@ $(CXXLIBS) $(MPFRC_LIBS)

$(BUILD_DIR)/ar.logh.sd: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.logh.sd.kdkpert: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_KDK_PERT $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.logh.ttl: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.logh.ttl.mpfrc: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(MPFRC_FLAGS) -D AR_TTL $< -o $@ $(CXXLIBS) $(MPFRC_LIBS)

$(BUILD_DIR)/ar.logh.ttl.m: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_TTL_GT_MULTI $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.logh.ttl.sd: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.logh.ttl.sd.kdkpert: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_KDK_PERT $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.logh.ttl.sd.cm: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE  $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.blogh.ttl.sd: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BLOGH $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.blogh.ttl.sd.cm: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BLOGH $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.blogh.ttl.sd.mpfrc: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(MPFRC_FLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BLOGH $< -o $@ $(CXXLIBS) $(MPFRC_LIBS)

$(BUILD_DIR)/ar.blogh.ttl.sd.cm.mpfrc: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(MPFRC_FLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BLOGH $< -o $@ $(CXXLIBS) $(MPFRC_LIBS)

$(BUILD_DIR)/ar.normblogh.ttl.sd: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_NORM_BLOGH $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.normblogh.ttl.sd.cm: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_NORM_BLOGH $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.mulall.ttl.sd: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_MUL_ALL_POT $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.mulall.ttl.sd.cm: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_MUL_ALL_POT $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.btlogh.ttl.sd: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BTLOGH $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.btlogh.ttl.sd.cm: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BTLOGH $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.btlogh.ttl.sd.mpfrc: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(MPFRC_FLAGS) -D AR_TTL -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BTLOGH $< -o $@ $(CXXLIBS) $(MPFRC_LIBS)

$(BUILD_DIR)/ar.btlogh.ttl.sd.cm.mpfrc: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(MPFRC_FLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_BTLOGH $< -o $@ $(CXXLIBS) $(MPFRC_LIBS)

$(BUILD_DIR)/ar.maxpot.ttl.sd.cm: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_MAX_POT $< -o $@ $(CXXLIBS)

$(BUILD_DIR)/ar.addpot.ttl.sd.cm: ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -D AR_TTL -D USE_CM_FRAME -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE -D AR_G_FUNC_ADD_INNER_POT $< -o $@ $(CXXLIBS)

install: $(addprefix $(BUILD_DIR)/, $(TARGET))
	install -m 755 $(addprefix $(BUILD_DIR)/, $(TARGET)) $(INSTALL_PATH)

clean:
	rm -rf $(BUILD_DIR)
