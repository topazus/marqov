cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(marqov CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_executable (marqov main.cpp)

find_package(MPI REQUIRED CXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(marqov       ${MPI_C_LIBRARIES})
SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_definitions(-O3)
#add_definitions(-g)
add_definitions(-fmax-errors=5)
add_definitions(-Wno-sign-compare)
add_definitions(-Wno-literal-suffix)


#add_subdirectory(testsuite)
#enable_testing()
