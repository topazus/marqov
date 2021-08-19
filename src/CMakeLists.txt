INCLUDE(CheckFunctionExists)

add_subdirectory(libmarqov)


# configuration for the main file
add_executable(main main.cpp libmarqov/util/registry.cpp)
target_include_directories(main PUBLIC ${MYHDF5INCLUDES})
target_link_libraries(main libMARQOV)
set_target_properties(main PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)


CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)

CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)

CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h)
