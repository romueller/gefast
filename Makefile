CXX=g++
CXXFLAGS=-std=c++11 -pthread -O3 -msse4.1
LDFLAGS=
BUILD=build
OBJ_DIR=$(BUILD)/objects
TARGET=GeFaST
INCLUDE=
SRC=main.cpp $(wildcard src/*.cpp)

OBJECTS=$(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(BUILD)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(BUILD)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET) $(OBJECTS) $(INCLUDE) $(LDFLAGS)

.PHONY: all build clean

build:
	@mkdir -p $(OBJ_DIR)

clean:
	rm -rf build/*




