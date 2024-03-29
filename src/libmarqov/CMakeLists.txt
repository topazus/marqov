find_package(HDF5 COMPONENTS CXX)
find_package(Threads)
include(GNUInstallDirs)


set(MYHDF5INCLUDES ${HDF5_INCLUDE_DIRS})
set(MYHDF5LIBS ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if( (NOT HDF5_FOUND) )
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
        message(STATUS "Fetching libhdf-1.12.1")
	set(libhdf5major "1.12")
	set(libhdf5minor "1")
	set(libhdf5ver "${libhdf5major}.${libhdf5minor}")
	include(FetchContent)
        FetchContent_Declare(libhdf5
	URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${libhdf5major}/hdf5-${libhdf5ver}/src/hdf5-${libhdf5ver}.tar.gz
#	GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
#        GIT_TAG hdf5-1_12_1
	)
	FetchContent_GetProperties(libhdf5)
        if(NOT libhdf5_POPULATED)

            FetchContent_Populate(libhdf5)
            # Make subproject to use 'BUILD_SHARED_LIBS=ON' setting.
            set(HDF5_EXTERNALLY_CONFIGURED 1)
    	    set(HDF5_EXPORTED_TARGETS "hdf5_cpp-static")
            set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "Build FORTRAN support")
            set(HDF5_BUILD_CPP_LIB ON CACHE BOOL "Build HDF5 C++ Library")
    	    set(HDF5_BUILD_HL_LIB OFF CACHE BOOL "Build HIGH Level HDF5 Library")
            set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Shared Libraries")
            set(HDF5_BUILD_TOOLS OFF CACHE BOOL "Build HDF5 Tools")
    	    set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Build HDF5 Library Examples")
            set(BUILD_TESTING OFF CACHE BOOL "Build HDF5 Unit Testing")
            set(HDF5_ENABLE_PARALLEL OFF CACHE BOOL "Enable parallel build (requires MPI)")

#	        set (HDF5_EXTERNAL_LIB_PREFIX "marqovhdf5")
            message(STATUS "libhdf5 source dir: ${libhdf5_SOURCE_DIR}")
            message(STATUS "libhdf5 binary dir: ${libhdf5_BINARY_DIR}")
            add_subdirectory(${libhdf5_SOURCE_DIR} ${libhdf5_BINARY_DIR})
            set(HDF5_DIR ${libhdf5_BINARY_DIR} CACHE PATH "libhdf5 path"  FORCE)
        endif()

        set(MYHDF5INCLUDES ${HDF5_INCLUDE_DIR} ${libhdf5_BINARY_DIR})
        set(MYHDF5LIBS hdf5_cpp-static)
    else()
        message(FATAL_ERROR "HDF5 not found and environment(cmake < 3.12) unsuited for automatic retrieval!")
    endif()
endif()

add_library(libMARQOV marqov.cpp)
target_include_directories(libMARQOV PUBLIC ${MYHDF5INCLUDES})
target_include_directories(libMARQOV INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS libMARQOV LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

target_link_libraries(libMARQOV ${MYHDF5LIBS} Threads::Threads)
set_target_properties(libMARQOV PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES POSITION_INDEPENDENT_CODE ON)

if(DEBIASINTEGERS)
    target_compile_definitions(libMARQOV PUBLIC DEBIASINTEGERS)
endif()

if(NOT COMPILETIMELOGLEVEL)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        SET(COMPILETIMELOGLEVEL "DEBUGVERBOSE" CACHE STRING "The maximum available LogLevel" FORCE)
    else()
        SET(COMPILETIMELOGLEVEL "RELEASEVERBOSE" CACHE STRING "The maximum available LogLevel" FORCE)
    endif()
endif()

message(STATUS "Selected compile time loglevel: ${COMPILETIMELOGLEVEL}")
target_compile_definitions(libMARQOV PUBLIC MLOG_EXTLEVEL=${COMPILETIMELOGLEVEL})

set_target_properties(libMARQOV PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )
if (HDF5_IS_PARALLEL)
# If a parrallel HDF5 Version is employed, then the HDF5 headers that we include
# require MPI and hence we need to make libmarqov dependent on MPI
    find_package(MPI)
    if(MPI_CXX_FOUND)
	target_compile_options(libMARQOV PUBLIC -DMPIMARQOV)
	if(CMAKE_VERSION VERSION_LESS 3.10)
    	    target_include_directories(libMARQOV PUBLIC ${MPI_INCLUDE_PATH})
    	    target_link_libraries(libMARQOV ${MPI_CXX_LIBRARIES})
    	    target_compile_options(libMARQOV PUBLIC ${MPI_CXX_COMPILE_FLAGS})
	else()
    	    target_include_directories(libMARQOV PUBLIC ${MPI_CXX_INCLUDE_DIRS})
    	    target_link_libraries(libMARQOV ${MPI_CXX_LIBRARIES})
    	    set_target_properties(libMARQOV PROPERTIES LINK_FLAGS " ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    	    target_compile_options(libMARQOV PUBLIC ${MPI_CXX_COMPILE_OPTIONS})
    	    target_compile_definitions(libMARQOV PUBLIC ${MPI_CXX_COMPILE_DEFINITIONS})
	endif()
endif(MPI_CXX_FOUND)
endif()