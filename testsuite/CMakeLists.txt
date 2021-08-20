cmake_minimum_required(VERSION 3.0.0)
Project(MARQOV.tests CXX)

add_subdirectory(components)
# add_subdirectory(systems)

add_executable(4-eigentest 4-eigentest.cpp)

target_link_libraries(4-eigentest libMARQOV)
target_include_directories(4-eigentest PUBLIC ${CMAKE_SOURCE_DIR}/src/libmarqov)

enable_testing()
add_test(4-eigentest 4-eigentest)
