cmake_minimum_required(VERSION 2.8)

project(OuterPoint)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
SET(CMAKE_BUILD_TYPE Debug)
MESSAGE(STATUS "BUILD TYPE:" ${CMAKE_BUILD_TYPE})

set(HDRS_DUTILS
  DUtils/Random.h
  DUtils/Timestamp.h)
set(SRCS_DUTILS
  DUtils/Random.cpp
  DUtils/Timestamp.cpp)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)

add_library(DUTILS SHARED ${SRCS_DUTILS})

add_executable(main ./main.cpp)
target_link_libraries(main ${Eigen_LIBS} ${OpenCV_LIBS} DUTILS)

