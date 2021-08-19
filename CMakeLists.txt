cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(marqov CXX)
INCLUDE (CheckTypeSize)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
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

option(DEBIASINTEGERS "This controls whether to employ debiasing for small integers. Not for PGI." OFF)
check_type_size("__uint128_t" UINT128LEN LANGUAGE CXX)
if(DEBIASINTEGERS AND NOT HAVE_UINT128LEN)
    message(FATAL_ERROR "Debiasing of integers requires a 128bit data type")
endif()

add_subdirectory(src)
add_subdirectory(demo)

#add_subdirectory(testsuite)
#enable_testing()
