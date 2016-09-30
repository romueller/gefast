cmake_minimum_required(VERSION 3.5)
project(SCT_PJ)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -O3")

set(SOURCE_FILES main.cpp include/Base.hpp src/Base.cpp include/Verification.hpp src/Verification.cpp include/Buffer.hpp include/Preprocessor.hpp src/Preprocessor.cpp include/Utility.hpp include/Relation.hpp include/SegmentFilter.hpp src/SegmentFilter.cpp include/Worker.hpp src/Worker.cpp include/SwarmClustering.hpp src/SwarmClustering.cpp src/Utility.cpp)
add_executable(SCT_PJ ${SOURCE_FILES})