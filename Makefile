CXX=g++
CXXFLAGS=-std=c++11 -pthread -O3 -msse4.1
BUILD=build
OBJ_DIR=$(BUILD)/objects
TARGET=GeFaST
INCLUDE=

# non-succinct compilation
LDFLAGS=
SRC=main.cpp src/Base.cpp src/Preprocessor.cpp src/SegmentFilter.cpp src/SIMD.cpp src/SwarmClustering.cpp \
    src/SwarmingSegmentFilter.cpp src/Utility.cpp src/Verification.cpp src/VerificationGotoh.cpp
OBJECTS=$(SRC:%.cpp=$(OBJ_DIR)/%.o)

# succinct compilation
SUCC_LDFLAGS=-lsdsl -lk2trees
SUCC_SRC=main.cpp $(wildcard src/*.cpp)
SUCC_OBJECTS=$(SUCC_SRC:%.cpp=$(OBJ_DIR)/%.o)

# library directories
SDSL_PREFIX?=/usr/local
K2TREES_PREFIX?=/usr/local

# preprocessor options
SUCCINCT?=0
SUCCINCT_FASTIDIOUS?=0
NO_QGRAM_FILTER?=0

PREP_OPTIONS=


# targets
non-succinct: build prepare target

succinct: build succinct-prepare succinct-target

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(PREP_OPTIONS) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

target: $(OBJECTS)
	@mkdir -p $(BUILD)
	$(CXX) $(PREP_OPTIONS) $(CXXFLAGS) -o $(BUILD)/$(TARGET) $(OBJECTS) $(INCLUDE) $(LDFLAGS)

succinct-target: $(SUCC_OBJECTS)
	$(if $(filter 00, $(SUCCINCT)$(SUCCINCT_FASTIDIOUS)), @echo "WARNING: No succinct option activated. Please see the documentation.")
	@mkdir -p $(BUILD)
	$(CXX) $(PREP_OPTIONS) $(CXXFLAGS) -o $(BUILD)/$(TARGET) $(SUCC_OBJECTS) $(INCLUDE) -L $(SDSL_PREFIX) -L $(K2TREES_PREFIX) $(SUCC_LDFLAGS)

prepare:
	$(if $(filter 1, $(NO_QGRAM_FILTER)), $(eval PREP_OPTIONS += -D QGRAM_FILTER=0))

succinct-prepare:
	$(if $(filter 1, $(NO_QGRAM_FILTER)), $(eval PREP_OPTIONS += -D QGRAM_FILTER=0))
	$(if $(filter 1, $(SUCCINCT)), $(eval PREP_OPTIONS += -D SUCCINCT))
	$(if $(filter 1, $(SUCCINCT_FASTIDIOUS)), $(eval PREP_OPTIONS += -D SUCCINCT_FASTIDIOUS))
	$(eval INCLUDE += -I$(K2TREES_PREFIX)/include/k2trees)

.PHONY: non-succinct succinct build clean

build:
	@mkdir -p $(OBJ_DIR)

clean:
	rm -rf build/*




