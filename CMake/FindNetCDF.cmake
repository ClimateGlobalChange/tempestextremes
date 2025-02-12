# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# Try to find NetCDF includes and library
#
# This module defines
#
#   - NetCDF_FOUND                - System has NetCDF
#   - NetCDF_VERSION              - the version of NetCDF
#
# Following components are available:
#
#   - C                           - C interface to NetCDF          (netcdf)
#   - CXX                         - CXX4 interface to NetCDF       (netcdf_c++4)
#   - Fortran                     - Fortran interface to NetCDF    (netcdff)
#   - CXX_LEGACY                  - Legacy C++ interface to NetCDF (netcdf_c++)
#
# For each component the following are defined:
#
#   - NetCDF_<comp>_FOUND         - whether the component is found
#   - NetCDF::NetCDF_<comp>       - target of component to be used with target_link_libraries()
#
# Caveat: The targets might not link properly with static libraries, setting NetCDF_<comp>_EXTRA_LIBRARIES may be required.
#
# The following paths will be searched in order if set in CMake (first priority) or environment (second priority)
#
#   - NetCDF_<comp>_ROOT
#   - NetCDF_<comp>_DIR
#   - NetCDF_<comp>_PATH
#   - The same variables with a NETCDF, NetCDF4, or NETCDF4 prefix instead of NetCDF
#   - NetCDF_ROOT
#   - NetCDF_DIR
#   - NetCDF_PATH
#   - The same variables with a NETCDF, NetCDF4, or NETCDF4 prefix instead of NetCDF
#
# The following variables affect the targets and NetCDF*_LIBRARIES variables:
#
#   - NetCDF_<comp>_EXTRA_LIBRARIES    - added to NetCDF::NetCDF_<comp> INTERFACE_LINK_LIBRARIES and NetCDF_<comp>_LIBRARIES
#
# Notes:
#
#   - If no components are defined, only the C component will be searched.
#

list( APPEND _possible_components C CXX Fortran CXX_LEGACY )

## Header names for each component
set( NetCDF_C_INCLUDE_NAME          netcdf.h )
set( NetCDF_CXX_INCLUDE_NAME        netcdf )
set( NetCDF_CXX_LEGACY_INCLUDE_NAME netcdfcpp.h )
set( NetCDF_Fortran_INCLUDE_NAME    netcdf.mod NETCDF.mod )

## Library names for each component
set( NetCDF_C_LIBRARY_NAME          netcdf )
set( NetCDF_CXX_LIBRARY_NAME        netcdf_c++4 netcdf-cxx4 )
set( NetCDF_CXX_LEGACY_LIBRARY_NAME netcdf_c++ )
set( NetCDF_Fortran_LIBRARY_NAME    netcdff )

foreach( _comp ${_possible_components} )
  string( TOUPPER "${_comp}" _COMP )
  set( _arg_${_COMP} ${_comp} )
  set( _name_${_COMP} ${_comp} )
endforeach()

unset( _search_components )
foreach( _comp ${${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS} )
  string( TOUPPER "${_comp}" _COMP )
  set( _arg_${_COMP} ${_comp} )
  list( APPEND _search_components ${_name_${_COMP}} )
  if( NOT _name_${_COMP} )
    message( FATAL_ERROR "Find${CMAKE_FIND_PACKAGE_NAME}: COMPONENT ${_comp} is not a valid component. Valid components: ${_possible_components}" )
  endif()
endforeach()
if( NOT _search_components )
  set( _search_components C )
endif()

## Search hints for finding include directories and libraries
foreach( _comp IN ITEMS "" "C" "CXX" "Fortran" "CXX_LEGACY" )
  set( __comp "_${_comp}" )
  if( NOT _comp )
    set( __comp "" )
  endif()

  set( _search_hints${__comp} )
  foreach( _name IN ITEMS NetCDF NETCDF NetCDF4 NETCDF4 )
    foreach( _var IN ITEMS ROOT DIR PATH )
      list( APPEND _search_hints${__comp} ${${_name}${__comp}_${_var}} ENV ${_name}${__comp}_${_var} )
    endforeach()
  endforeach()

  ## Old-school HPC module env variable names
  foreach( _name IN ITEMS NetCDF NETCDF NetCDF4 NETCDF4 )
    list(APPEND _search_hints${__comp} ${${_name}${__comp}} ENV ${_name}${__comp})
  endforeach()
endforeach()

set( _found FALSE )
set( _req_vars )
foreach( _comp ${_search_components} )
  list( APPEND _req_vars NetCDF_${_comp}_INCLUDE_DIR NetCDF_${_comp}_LIBRARY )

  ## Find include directories
  find_path(NetCDF_${_comp}_INCLUDE_DIR
    NAMES ${NetCDF_${_comp}_INCLUDE_NAME}
    DOC "netcdf ${_comp} include directory"
    HINTS ${_search_hints_${_comp}} ${_search_hints}
    PATH_SUFFIXES include ../../include
  )
  mark_as_advanced(NetCDF_${_comp}_INCLUDE_DIR)

  ## Find libraries for each component
  string( TOUPPER "${_comp}" _COMP )

  find_library(NetCDF_${_comp}_LIBRARY
    NAMES ${NetCDF_${_comp}_LIBRARY_NAME}
    DOC "netcdf ${_comp} library"
    HINTS ${_search_hints_${_comp}} ${_search_hints}
    PATH_SUFFIXES lib ../../lib
  )
  mark_as_advanced(NetCDF_${_comp}_LIBRARY)
  if( NetCDF_${_comp}_LIBRARY AND NOT (NetCDF_${_comp}_LIBRARY MATCHES ".a$") )
    set( NetCDF_${_comp}_LIBRARY_SHARED TRUE )
  endif()
  if( NetCDF_${_comp}_LIBRARY AND NetCDF_${_comp}_INCLUDE_DIR )
    set( ${CMAKE_FIND_PACKAGE_NAME}_${_arg_${_COMP}}_FOUND TRUE )
    set( _found TRUE )

    if (NOT TARGET NetCDF::NetCDF_${_comp})
      add_library(NetCDF::NetCDF_${_comp} UNKNOWN IMPORTED)
      set_target_properties(NetCDF::NetCDF_${_comp} PROPERTIES
        IMPORTED_LOCATION "${NetCDF_${_comp}_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_${_comp}_INCLUDE_DIR}")
      if( DEFINED NetCDF_${_comp}_EXTRA_LIBRARIES )
        target_link_libraries(NetCDF::NetCDF_${_comp} INTERFACE ${NetCDF_${_comp}_EXTRA_LIBRARIES})
      endif()
    endif()
  endif()
endforeach()

## Find version
if (_found)
  set( _config_search_hints ${_search_hints} )
  set( _include_dirs )
  foreach( _comp ${_search_components} )
    if( DEFINED _search_hints_${_comp} )
      list( APPEND _config_search_hints ${_search_hints_${_comp}} )
    endif()
    if( DEFINED NetCDF_${_comp}_INCLUDE_DIR )
      list( APPEND _include_dirs ${NetCDF_${_comp}_INCLUDE_DIR} )
    endif()
  endforeach()
  if( _config_search_hints )
    list( REMOVE_DUPLICATES _config_search_hints )
  endif()
  if( _include_dirs )
    list( REMOVE_DUPLICATES _include_dirs )
  endif()

  find_program( NETCDF_CONFIG_EXECUTABLE
      NAMES nc-config
      HINTS ${_config_search_hints}
      PATH_SUFFIXES bin Bin ../../bin
      DOC "NetCDF nc-config helper" )
  mark_as_advanced( NETCDF_CONFIG_EXECUTABLE )

  find_file( NETCDF_META_H
    NAMES netcdf_meta.h
    HINTS ${_include_dirs}
    NO_DEFAULT_PATH
    DOC "NetCDF path to netcdf_meta.h" )
  mark_as_advanced( NETCDF_META_H )

  if( NETCDF_CONFIG_EXECUTABLE )
    execute_process( COMMAND ${NETCDF_CONFIG_EXECUTABLE} --version
      RESULT_VARIABLE _netcdf_config_result
      OUTPUT_VARIABLE _netcdf_config_version)

    if( _netcdf_config_result EQUAL 0 )
      string(REGEX REPLACE ".* ((([0-9]+)\\.)+([0-9]+)).*" "\\1" NetCDF_VERSION "${_netcdf_config_version}" )
    endif()

  elseif( NETCDF_META_H )

    file(STRINGS ${NETCDF_META_H} _netcdf_version_lines
      REGEX "#define[ \t]+NC_VERSION_(MAJOR|MINOR|PATCH|NOTE)")
    string(REGEX REPLACE ".*NC_VERSION_MAJOR *\([0-9]*\).*" "\\1" _netcdf_version_major "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_MINOR *\([0-9]*\).*" "\\1" _netcdf_version_minor "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_PATCH *\([0-9]*\).*" "\\1" _netcdf_version_patch "${_netcdf_version_lines}")
    string(REGEX REPLACE ".*NC_VERSION_NOTE *\"\([^\"]*\)\".*" "\\1" _netcdf_version_note "${_netcdf_version_lines}")
    set(NetCDF_VERSION "${_netcdf_version_major}.${_netcdf_version_minor}.${_netcdf_version_patch}${_netcdf_version_note}")
    unset(_netcdf_version_major)
    unset(_netcdf_version_minor)
    unset(_netcdf_version_patch)
    unset(_netcdf_version_note)
    unset(_netcdf_version_lines)
  endif()
endif ()

## Finalize find_package
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args( ${CMAKE_FIND_PACKAGE_NAME}
  REQUIRED_VARS ${_req_vars}
  VERSION_VAR NetCDF_VERSION
  HANDLE_COMPONENTS )

if( ${CMAKE_FIND_PACKAGE_NAME}_FOUND AND NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY )
  message( STATUS "Find${CMAKE_FIND_PACKAGE_NAME} defines targets:" )
  foreach( _comp ${_search_components} )
    string( TOUPPER "${_comp}" _COMP )

    if( ${CMAKE_FIND_PACKAGE_NAME}_${_arg_${_COMP}}_FOUND )
      message( STATUS "  - NetCDF::NetCDF_${_comp} [${NetCDF_${_comp}_LIBRARY}]")
    endif()
  endforeach()
endif()

## Backwards compatibility, only reachable if ECBUILD_2_COMPAT is ON
# Assumes the following internal variables are defined:
#   - _search_components
#   - _arg_<COMP>
#   - NetCDF_<comp>_INCLUDE_DIR
#   - NetCDF_<comp>_LIBRARY
include( netcdf_compat OPTIONAL )