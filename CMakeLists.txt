cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(marqov CXX)
INCLUDE (CheckTypeSize)

# Add our cmake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# make the Release build the default and tell the user that none is specified.
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Some code to have git revision strings in our output files
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1) # the git revision is available in GIT_SHA1
# now we need to retrieve the branch
if (NOT GIT_FOUND)
    find_package(Git QUIET)
endif()
if(GIT_FOUND)
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        OUTPUT_VARIABLE GIT_BRANCH)
string(STRIP "${GIT_BRANCH}" GIT_BRANCH)
else()
set(GIT_BRANCH, "N/A")
endif()

# This option specifies whether additional scrambling is used on our integers
option(DEBIASINTEGERS "This controls whether to employ debiasing for small integers. Not for PGI." OFF)
check_type_size("__uint128_t" UINT128LEN LANGUAGE CXX)
if(DEBIASINTEGERS AND NOT HAVE_UINT128LEN)
    message(FATAL_ERROR "Debiasing of integers requires a 128bit data type")
endif()

# Try to use at least the native cpu instruction set if possible on the platform
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
CHECK_CXX_COMPILER_FLAG("-mtune=native" COMPILER_SUPPORTS_MTUNE_NATIVE)
if(COMPILER_SUPPORTS_MTUNE_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native")
endif()

# Add our sub folders
add_subdirectory(src)
add_subdirectory(demo)

add_subdirectory(testsuite)
enable_testing()
