#(C) https://github.com/capitalaslash/cmake-modules/blob/master/FindLibmesh.cmake
# (BSD - License)

# - Try to find Libmesh
# Once done this will define
#
#  LIBMESH_FOUND          - Libmesh has been successfully found
#  LIBMESH_INCLUDE_DIRS   - Libmesh include directories
#  LIBMESH_LIBRARIES      - Libmesh libraries
#  LIBMESH_DEFINITIONS    - Libmesh definitions
#  LIBMESH_FLAGS          - Libmesh flags
#  LIBMESH_VERSION_STRING - Libmesh version
#
#  Usage:
#  find_package(Libmesh)
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE NONE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_UPPER)
message("BUILD_TYPE_UPPER = " ${BUILD_TYPE_UPPER})
if(${BUILD_TYPE_UPPER} MATCHES DEBUG)
  set(METHOD dbg)
else()
  set(METHOD opt)
endif()
#set(METHOD opt)
message(STATUS "linking against ${METHOD} libmesh library")

# required for LIBMESH_DEFINITIONS
find_package(PkgConfig REQUIRED)

set(LIBMESH_DIR LIBMESH_DIR-NOTFOUND CACHE PATH "Libmesh installation directory")

if(LIBMESH_DIR)
  set(ENV{PKG_CONFIG_PATH} "${LIBMESH_DIR}/lib/pkgconfig")
#  message($ENV{PKG_CONFIG_PATH})
endif()

pkg_check_modules(PC_LIBMESH libmesh-${METHOD})
message(STATUS "PC_LIBMESH_FOUND = ${PC_LIBMESH_FOUND}")
message(STATUS "PC_LIBMESH_LIBRARIES = ${PC_LIBMESH_LIBRARIES}")
message(STATUS "PC_LIBMESH_LIBRARY_DIRS = ${PC_LIBMESH_LIBRARY_DIRS}")

# Get all required libraries from libmesh-config
# exec_program(${LIBMESH_CONFIG_EXECUTABLE}
#   ARGS --libs
#   OUTPUT_VARIABLE LIBMESH_LINK_FLAGS
# )

execute_process(
  COMMAND ${LIBMESH_CONFIG_EXECUTABLE} --libs
  OUTPUT_VARIABLE LIBMESH_LINK_FLAGS
  RESULT_VARIABLE LIBMESH_LINK_FLAGS_RETURN
)

# Convert the flags into a list
string(REPLACE " " ";" LIBMESH_LINK_LIST "${LIBMESH_LINK_FLAGS}")

# Initialize empty lists for libraries and library paths
set(LIBMESH_LIBRARIES "")
set(LIBRARY_DIRS "")

# Process each flag
foreach(FLAG ${LIBMESH_LINK_LIST})
  if(${FLAG} MATCHES "^-L(.*)")
    # Extract library path from -L flag
    string(REGEX REPLACE "^-L" "" LIB_PATH ${FLAG})
    list(APPEND LIBRARY_DIRS ${LIB_PATH})
  elseif(${FLAG} MATCHES "^-l(.*)")
    # Extract library name from -l flag
    string(REGEX REPLACE "^-l" "" LIB_NAME ${FLAG})
    
    # Find the actual library file
    find_library(FOUND_LIB_${LIB_NAME}
      NAMES ${LIB_NAME}
      PATHS ${LIBRARY_DIRS} ${PC_LIBMESH_LIBRARY_DIRS} ${LIBMESH_DIR}/lib
      NO_DEFAULT_PATH
    )
    
    if(FOUND_LIB_${LIB_NAME})
      list(APPEND LIBMESH_LIBRARIES ${FOUND_LIB_${LIB_NAME}})
      mark_as_advanced(FOUND_LIB_${LIB_NAME})
    endif()
  endif()
endforeach()

# Get library paths from pkg-config
set(LIBMESH_LIBRARIES "")

# Process each library from pkg-config
foreach(LIB ${PC_LIBMESH_LIBRARIES})
    # First try to find in libmesh directory
    find_library(FOUND_LIB_${LIB}
        NAMES lib${LIB}.dylib lib${LIB}.so ${LIB}
        PATHS ${LIBMESH_DIR}/lib
        NO_DEFAULT_PATH
    )
    
    # If not found in libmesh dir, try petsc dir
    if(NOT FOUND_LIB_${LIB})
        find_library(FOUND_LIB_${LIB}
            NAMES lib${LIB}.dylib lib${LIB}.so ${LIB}
            PATHS ${PETSC_DIR}/lib
            NO_DEFAULT_PATH
        )
    endif()
    
    # If still not found, try system paths
    if(NOT FOUND_LIB_${LIB})
        find_library(FOUND_LIB_${LIB}
            NAMES lib${LIB}.dylib lib${LIB}.so ${LIB}
        )
    endif()
    
    if(FOUND_LIB_${LIB})
        list(APPEND LIBMESH_LIBRARIES ${FOUND_LIB_${LIB}})
        message(STATUS "Found library ${LIB}: ${FOUND_LIB_${LIB}}")
    else()
        message(WARNING "Could not find library ${LIB}")
    endif()
endforeach()

# distinguish flags and definitions (-D...)
#message(STATUS "PC_LIBMESH_CFLAGS_OTHER = ${PC_LIBMESH_CFLAGS_OTHER}")
foreach(FLAG ${PC_LIBMESH_CFLAGS_OTHER})
  if(${FLAG} MATCHES "^[-][D].+")
    set(PC_LIBMESH_CFLAGS_DEFS "${PC_LIBMESH_CFLAGS_DEFS} ${FLAG}")
  else()
    set(PC_LIBMESH_CFLAGS_FLAGS "${PC_LIBMESH_CFLAGS_FLAGS} ${FLAG}")
  endif()
endforeach()
set(LIBMESH_DEFINITIONS ${PC_LIBMESH_CFLAGS_DEFS})
set(LIBMESH_FLAGS ${PC_LIBMESH_CFLAGS_FLAGS})

find_path(LIBMESH_INCLUDE_DIR libmesh/libmesh.h
  HINTS ${PC_LIBMESH_INCLUDEDIR} ${PC_LIBMESH_INCLUDE_DIRS}
  PATH_SUFFIXES libmesh
)

# Get include paths from libmesh-config
# exec_program(${LIBMESH_CONFIG_EXECUTABLE}
#   ARGS --include
#   OUTPUT_VARIABLE LIBMESH_INC_FLAGS
# )

execute_process(
  COMMAND ${LIBMESH_CONFIG_EXECUTABLE} --include
  OUTPUT_VARIABLE LIBMESH_INC_FLAGS
  RESULT_VARIABLE LIBMESH_INC_FLAGS_RETURN
)

# Convert to list and process include paths
string(REPLACE " " ";" LIBMESH_INC_LIST "${LIBMESH_INC_FLAGS}")
set(LIBMESH_INCLUDE_DIRS "")

foreach(FLAG ${LIBMESH_INC_LIST})
  if(${FLAG} MATCHES "^-I(.*)")
    # Extract path from -I flag
    string(REGEX REPLACE "^-I" "" INC_PATH ${FLAG})
    # Clean up path
    string(REPLACE "//" "/" INC_PATH ${INC_PATH})
    list(APPEND LIBMESH_INCLUDE_DIRS ${INC_PATH})
  endif()
endforeach()

# Also add the main libmesh include directory
if(LIBMESH_DIR)
  list(APPEND LIBMESH_INCLUDE_DIRS "${LIBMESH_DIR}/include")
endif()

message(STATUS "LIBMESH_INCLUDE_DIRS = ${LIBMESH_INCLUDE_DIRS}")

if(PC_LIBMESH_VERSION)
  set(LIBMESH_VERSION_STRING ${PC_LIBMESH_VERSION})
endif()

# handle the QUIETLY and REQUIRED arguments and set LIBMESH_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibMesh
    REQUIRED_VARS 
        LIBMESH_INCLUDE_DIR
        LIBMESH_LIBRARIES
    VERSION_VAR PC_LIBMESH_VERSION
)

mark_as_advanced(
    LIBMESH_INCLUDE_DIR
    LIBMESH_LIBRARIES
    LIBMESH_FLAGS
    LIBMESH_DEFINITIONS
)
