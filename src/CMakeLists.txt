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
          message(STATUS "R found, building R bindings")
	  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
          INCLUDE_DIRECTORIES("/usr/share/R/include/") # Debian Bullseye path, could also carry over to Ubuntu
          INCLUDE(${SWIG_USE_FILE})

          swig_add_library(RShowcase TYPE SHARED LANGUAGE R SOURCES libshowcase.i libshowcase.cpp)
          SWIG_LINK_LIBRARIES(RShowcase libMARQOV)
          endif()
          if(SWIG_ocaml_FOUND)
          message(STATUS "ocaml found, building OCaml bindings")
          INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
          INCLUDE(${SWIG_USE_FILE})

          swig_add_library(ocamlShowcase TYPE SHARED LANGUAGE ocaml SOURCES libshowcase.i libshowcase.cpp)
          target_compile_options(ocamlShowcase PRIVATE "-fpermissive")# could work for gcc
          SWIG_LINK_LIBRARIES(ocamlShowcase libMARQOV)
          endif()
      endif()
    endif()
endif()
