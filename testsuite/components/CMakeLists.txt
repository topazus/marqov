cmake_minimum_required(VERSION 3.0.0)
Project(components.tests CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(1-cachecontainer 1-cachecontainer.cpp)
add_executable(2-cachecontainer-long-sim 2-cachecontainer-long-sim.cpp)
add_executable(3-cachecontainer-long-obs 3-cachecontainer-long-obs.cpp)

set(comptests 1-cachecontainer 2-cachecontainer-long-sim 3-cachecontainer-long-obs)

foreach(test ${comptests})
target_link_libraries(${test} libMARQOV)
target_include_directories(${test} PUBLIC ${CMAKE_SOURCE_DIR}/src/libmarqov)
endforeach(test)

enable_testing()
add_test(1-cachecontainer 1-cachecontainer)
add_test(2-cachecontainer-long-sim 2-cachecontainer-long-sim)
add_test(3-cachecontainer-long-obs 3-cachecontainer-long-obs)
