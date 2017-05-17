cmake_minimum_required(VERSION 3.5)
project(GeFaST)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -O3 -msse4.1")

set(SOURCE_FILES main.cpp include/Base.hpp src/Base.cpp include/Verification.hpp src/Verification.cpp include/Buffer.hpp include/Preprocessor.hpp src/Preprocessor.cpp include/Utility.hpp include/Relation.hpp include/SegmentFilter.hpp src/SegmentFilter.cpp include/Worker.hpp src/Worker.cpp include/SwarmClustering.hpp src/SwarmClustering.cpp src/Utility.cpp include/VerificationGotoh.hpp src/VerificationGotoh.cpp include/SwarmingSegmentFilter.hpp src/SwarmingSegmentFilter.cpp include/SIMD.hpp src/SIMD.cpp src/simd-verification/matrix.cc src/simd-verification/search8.cc src/simd-verification/search16.cc src/simd-verification/ssse3.cc src/simd-verification/util.cc)
add_executable(GeFaST ${SOURCE_FILES})