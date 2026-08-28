# Exactly ONE executable per base method (4 total, see TARGET below); the
# feature suffix (.sd/.cm/.mpfrc/.kdkpert) is composed from the use_* flags,
# so every configuration builds only 4 binaries:
#   ar.logh$(FEATURE_SUFFIX)      ar.logh.ttl$(FEATURE_SUFFIX)
#   ar.blogh.ttl$(FEATURE_SUFFIX) ar.btlogh.ttl$(FEATURE_SUFFIX)
# Note: g-func methods require .sd (AR_SLOWDOWN_TREE); without it --g-func is
# silently inactive (plain LogH behavior).
INSTALL_PATH=~/bin
BUILD_DIR=./build

VPATH=../../src 

SOURCE=${shell ls ../../src/AR/*.h} ${shell ls ../../src/Common/*.h} ${shell ls ./*.h}

use_x86=yes
# feature combination: each enabled flag appends its suffix and compile flags
# to ALL 4 base methods (suffix order: .sd .cm .mpfrc .kdkpert;
# use_kdkpert requires use_sd and forces it on)
use_sd=yes
use_cm=yes
use_kdkpert=no
use_mpfrc=no
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

## feature suffix / flags (shared by all base methods)----#
# effective sd state: AR_KDK_PERT requires the slowdown tree, so use_kdkpert
# forces sd on (an internal variable is used because command-line use_sd=no
# would otherwise override an in-makefile reassignment)
ifeq ($(use_kdkpert),yes)
SD_ENABLED = yes
else
SD_ENABLED = $(use_sd)
endif

FEATURE_SUFFIX =
FEATURE_FLAGS  =
FEATURE_LIBS   =
ifeq ($(SD_ENABLED),yes)
FEATURE_SUFFIX := $(FEATURE_SUFFIX).sd
FEATURE_FLAGS += -D AR_SLOWDOWN_TREE -D AR_SLOWDOWN_TIMESCALE
endif
ifeq ($(use_cm),yes)
FEATURE_SUFFIX := $(FEATURE_SUFFIX).cm
FEATURE_FLAGS += -D USE_CM_FRAME
endif
ifeq ($(use_mpfrc),yes)
FEATURE_SUFFIX := $(FEATURE_SUFFIX).mpfrc
FEATURE_FLAGS += -D USE_MPFRC# -I $(MPFRC_PATH)
FEATURE_LIBS  += -lmpfr -lgmp
endif
ifeq ($(use_kdkpert),yes)
FEATURE_SUFFIX := $(FEATURE_SUFFIX).kdkpert
FEATURE_FLAGS += -D AR_KDK_PERT
endif

# 4 executables: one per base method with the common feature suffix
TARGET= ar.logh$(FEATURE_SUFFIX) ar.logh.ttl$(FEATURE_SUFFIX) ar.blogh.ttl$(FEATURE_SUFFIX) ar.btlogh.ttl$(FEATURE_SUFFIX)


all: $(addprefix $(BUILD_DIR)/, $(TARGET))

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

## Target --------------------------------------------#

$(BUILD_DIR)/ar.logh$(FEATURE_SUFFIX): ar.cxx  ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) $< -o $@ $(CXXLIBS) $(FEATURE_LIBS)

$(BUILD_DIR)/ar.logh.ttl$(FEATURE_SUFFIX): ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -D AR_TTL $< -o $@ $(CXXLIBS) $(FEATURE_LIBS)

$(BUILD_DIR)/ar.blogh.ttl$(FEATURE_SUFFIX): ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -D AR_TTL -D AR_G_FUNC_BLOGH $< -o $@ $(CXXLIBS) $(FEATURE_LIBS)

$(BUILD_DIR)/ar.btlogh.ttl$(FEATURE_SUFFIX): ar.cxx ${SOURCE} | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(FEATURE_FLAGS) -D AR_TTL -D AR_G_FUNC_BTLOGH $< -o $@ $(CXXLIBS) $(FEATURE_LIBS)

install: $(addprefix $(BUILD_DIR)/, $(TARGET))
	install -m 755 $(addprefix $(BUILD_DIR)/, $(TARGET)) $(INSTALL_PATH)

clean:
	rm -rf $(BUILD_DIR)
