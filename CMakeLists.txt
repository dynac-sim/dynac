cmake_minimum_required(VERSION 2.8)

project(dynac Fortran)

set(DYNAC_VERSION 6.0.R15.beta2)
message(STATUS "Dynac version ${DYNAC_VERSION}")

# Figure out platform..
if(APPLE)
    add_definitions(-D__APPLE__)
elseif(UNIX) # Linux = UNIX and not APPLE..
    add_definitions(-D__linux__)
elseif(WINDOWS OR CYGWIN)
    add_definitions(-D_WIN32)
else()
    message(FATAL "We cannot compile on your system")
endif()

# Default compiler flags..
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -pedantic -fcheck=bounds -Wall -Wextra -ffpe-trap=invalid,overflow,underflow -fbacktrace")

# Build type for Coverage checks:
set(CMAKE_Fortran_FLAGS_DEBUGCOVERAGE "${CMAKE_Fortran_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
set(LINK_FLAGS_DEBUGCOVERAGE "${LINK_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")

# Tests setup..
if(CMAKE_Fortran_COMPILER_ID MATCHES "unknown")
    string(REGEX REPLACE ".*/([^ ])" "\\1" _fcomp_ "${CMAKE_Fortran_COMPILER}" )
else()
    set(_fcomp_ "${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}")
endif()
set(BUILDNAME "Dynac-${DYNAC_VERSION}_${_fcomp_}" CACHE STRING "Name of build on the dashboard")
mark_as_advanced(BUILDNAME)
enable_testing()
include(CTest)


add_subdirectory(source)
add_subdirectory(plot)
add_subdirectory(datafiles)
add_subdirectory(doc)

