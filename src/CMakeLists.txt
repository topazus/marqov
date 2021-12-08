INCLUDE(CheckFunctionExists)

add_subdirectory(libmarqov)

SET(mainTARGETS main main.mpi)


# configuration for the main file
foreach(target ${mainTARGETS})
    add_executable (${target} libmarqov/util/registry.cpp main.cpp)
    target_link_libraries(${target} libMARQOV)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)
    target_include_directories(${target} PUBLIC ${MYHDF5INCLUDES})
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )

endforeach(target)

find_package(MPI)
if(MPI_CXX_FOUND)
    target_compile_options(main.mpi PUBLIC -DMPIMARQOV)
    if(CMAKE_VERSION VERSION_LESS 3.10)
	target_include_directories(main.mpi PUBLIC ${MPI_INCLUDE_PATH})
	target_link_libraries(main.mpi ${MPI_CXX_LIBRARIES})
	target_compile_options(main.mpi PUBLIC ${MPI_CXX_COMPILE_FLAGS})
    else()
	target_include_directories(main.mpi PUBLIC ${MPI_CXX_INCLUDE_DIRS})
	target_link_libraries(main.mpi ${MPI_CXX_LIBRARIES})
	set_target_properties(main.mpi PROPERTIES LINK_FLAGS " ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
	target_compile_options(main.mpi PUBLIC ${MPI_CXX_COMPILE_OPTIONS})
	target_compile_definitions(main.mpi PUBLIC ${MPI_CXX_COMPILE_DEFINITIONS})

    endif()
endif(MPI_CXX_FOUND)


CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)

CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)

CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h)

# Let's check whether we have SWIG and can hence build the language bindings.

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.8.0)
    find_package(SWIG 3.0 OPTIONAL_COMPONENTS python java lua r ocaml)

    if(SWIG_FOUND)
      cmake_policy(SET CMP0086 NEW)
      cmake_policy(SET CMP0078 NEW)
      message(STATUS "SWIG found: ${SWIG_EXECUTABLE}.")
      set_property(SOURCE libshowcase.i PROPERTY CPLUSPLUS ON)
      FIND_PACKAGE(PythonLibs)
      if(PYTHONLIBS_FOUND)
          message(STATUS "Python found, building Python bindings")
	  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
          INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
          INCLUDE(${SWIG_USE_FILE})

          swig_add_library(pyShowcase TYPE SHARED LANGUAGE python SOURCES libshowcase.i libshowcase.cpp)
# message(STATUS ${MYHDF5INCLUDES})
#           include_directories(${MYHDF5INCLUDES})
	  SWIG_LINK_LIBRARIES(pyShowcase libMARQOV ${PYTHON_LIBRARIES})
      endif()
      FIND_PACKAGE(Lua)
      if(LUA_FOUND)
	message(STATUS "Lua found, building Lua bindings")
	  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
          INCLUDE_DIRECTORIES(${LUA_INCLUDE_DIR})
          INCLUDE(${SWIG_USE_FILE})

          swig_add_library(LuaShowcase TYPE SHARED LANGUAGE lua SOURCES libshowcase.i libshowcase.cpp)
          SWIG_LINK_LIBRARIES(LuaShowcase libMARQOV ${LUA_LIBRARIES})
      endif()
      FIND_PACKAGE(JNI)
      if(JNI_FOUND)
          message(STATUS "Java Development tools(JNI) found, building Java bindings")
	  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
          INCLUDE_DIRECTORIES(${JAVA_INCLUDE_PATH} ${JAVA_INCLUDE_PATH2})
          INCLUDE(${SWIG_USE_FILE})

          swig_add_library(JNIShowcase TYPE SHARED LANGUAGE java SOURCES libshowcase.i libshowcase.cpp)
          SWIG_LINK_LIBRARIES(JNIShowcase libMARQOV ${JNI_LIBRARIES})
      endif()

      if(CMAKE_VERSION VERSION_GREATER 3.18)# language specific checks from SWIG only available in newer cmake versions
          
          if(SWIG_r_FOUND)
          message(STATUS "R found, trying to build R bindings. Make sure development headers are present!")
	  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
          INCLUDE_DIRECTORIES("/usr/share/R/include/") # Debian Bullseye path, could also carry over to Ubuntu
          INCLUDE(${SWIG_USE_FILE})

          swig_add_library(RShowcase TYPE SHARED LANGUAGE R SOURCES libshowcase.i libshowcase.cpp)
          SWIG_LINK_LIBRARIES(RShowcase libMARQOV)
          endif()
          if(SWIG_ocaml_FOUND)
          message(STATUS "ocaml found, trying to build OCaml bindings. Make sure development headers are present!")
          INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
          INCLUDE(${SWIG_USE_FILE})

          swig_add_library(ocamlShowcase TYPE SHARED LANGUAGE ocaml SOURCES libshowcase.i libshowcase.cpp)
          target_compile_options(ocamlShowcase PRIVATE "-fpermissive")# could work for gcc
          SWIG_LINK_LIBRARIES(ocamlShowcase libMARQOV)
          endif()
      endif()
    endif()
endif()
