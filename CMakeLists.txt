cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(marqov CXX)

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

if(CMAKE_VERSION VERSION_LESS "3.10.0")
find_package(MPI)
else()
find_package(MPI COMPONENTS CXX)

endif()

add_subdirectory(src)

#add_subdirectory(testsuite)
#enable_testing()
