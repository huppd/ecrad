# This CMakeLists.txt file is the entry point for the configuration and the
# compilation of the ECRAD.

FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

cmake_minimum_required(VERSION 3.0)

project("ECRAD - ECMWF" VERSION 1.4.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/module")

enable_language(Fortran)

set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
  include(cmake/compiler/pgi.cmake)
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  include(cmake/compiler/intel.cmake)
else()
  include(cmake/compiler/gnu.cmake)
endif()

option(SINGLE_PRECISION "Use single precision" OFF)
option(USE_OPENMP "Use OpenMP" ON)
option(USE_DR_HOOK "Use Dr Hook profiling system" OFF)

#if("${NetCDF_VERSION}" VERSION_GREATER_EQUAL "4")
  set(FC_FPPFLAGS "${FC_FPPFLAGS} -DNC_NETCDF")
#endif()

if(SINGLE_PRECISION)
  set(FC_FPPFLAGS "${FC_FPPFLAGS} -DSINGLE_PRECISION")
endif()

if(USE_DR_HOOK)
  message( FATAL_ERROR "Dr hook not supported in this build")
endif()

set(Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

set(CMAKE_Fortran_FLAGS_RELEASE "${FC_RELEASE_WARNFLAGS} ${FC_BASIC_FLAGS} ${FC_FPPFLAGS} ${FC_RELEASE_OPTFLAGS} ${FC_RELEASE_DEBUGFLAGS} -I ${Fortran_MODULE_DIRECTORY}")
set(CMAKE_Fortran_FLAGS_DEBUG "${FC_DEBUG_WARNFLAGS} ${FC_BASIC_FLAGS} ${FC_FPPFLAGS} ${FC_DEBUG_OPTFLAGS} ${FC_DEBUG_DEBUGFLAGS} -I ${Fortran_MODULE_DIRECTORY}")

set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Debug build")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
  message("Release build")
ELSE()
  message( FATAL_ERROR "Unsupported build type ${CMAKE_BUILD_TYPE}")
ENDIF()

if(USE_OPENMP)
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${FC_OMPFLAG}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${FC_OMPFLAG}")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include ${NETCDF_INCLUDES})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(ifsaux)
set_target_properties( ifsaux PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(ifsaux)

add_library(drhook)
set_target_properties( drhook PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(drhook)
add_dependencies(drhook ifsaux)

add_library(ifsrrtm)
set_target_properties( ifsrrtm PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(ifsrrtm)
add_dependencies(ifsrrtm drhook)

add_library(utilities)
set_target_properties( utilities PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(utilities)
add_dependencies(utilities drhook ifsaux)

add_library(radiation)
set_target_properties( radiation PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(radiation)
add_dependencies(radiation utilities)

add_library(radsurf)
set_target_properties( radsurf PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
add_subdirectory(radsurf)
add_dependencies(radsurf radiation)

add_subdirectory(driver)
message(
  STATUS
  "NETCDF_LIBRARIES_F90: ${NETCDF_LIBRARIES_F90}"
)
add_executable(ecrad ${ecrad_SRC})
set_target_properties( ecrad PROPERTIES Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY}" )
set(STATIC_LIBS ifsaux drhook ifsrrtm utilities radiation radsurf)
target_link_libraries(ecrad PUBLIC ${STATIC_LIBS} ${NETCDF_LIBRARIES_F90})


