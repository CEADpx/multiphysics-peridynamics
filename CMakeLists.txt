# -------------------------------------------
# Copyright (c) 2021 - 2025 Prashant K. Jha
# -------------------------------------------
# https://github.com/CEADpx/multiphysics-peridynamics
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE)

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(NodalFEM_Multiphysics LANGUAGES CXX)

# ****************************************************************************
# Project setting
# ****************************************************************************
# Set default build type "Release"
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif()

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add cmake scripts path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# set flag for information dump during cmake build
set(Enable_CMAKE_Debug_Build TRUE CACHE BOOL "Output important information during cmake for isolating build issues")

# ****************************************************************************
# Package search
# ****************************************************************************
# **** bash ****#
find_program (BASH_PROGRAM bash REQUIRED)

# **** MPI ****#
find_package(MPI REQUIRED)
# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)

# **** pthreads ****#
find_package(Threads REQUIRED)

# **** libmesh ****#
find_package(LibMesh REQUIRED)

# **** petsc ****#
# Note: -DPETSC_DIR should point to the same directory that libmesh was compiled with
find_package(Petsc REQUIRED)

# **** nlohmann_json ****#
find_package(nlohmann_json REQUIRED)

# debug output flags
if (${Enable_CMAKE_Debug_Build})
message("LIBMESH_DEFINITIONS = ${PC_LIBMESH_CFLAGS_DEFS}")
message("LIBMESH_FLAGS = ${PC_LIBMESH_CFLAGS_FLAGS}")
message("LIBMESH_INCLUDE_DIRS = ${LIBMESH_INCLUDE_DIRS}")
message("LIBMESH_LIBRARIES = ${LIBMESH_LIBRARIES}")
message("PETSC_LIBRARIES = ${PETSC_LIBRARIES}")
endif()

# issue error and exit
if ((NOT MPI_FOUND) AND (NOT LIBMESH_FOUND))
    message(FATAL_ERROR "ERROR: MPI and LIBMESH are required")
endif ()

# ****************************************************************************
# Build, compiler, linker flags
# ****************************************************************************

# Set compiler flags based on compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_FLAGS} -xHost")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
    # Handle pthread separately through target_link_libraries
endif()

# RPATH settings for macOS
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${LIBMESH_DIR}/lib;${PETSC_DIR}/lib")

if (${Enable_CMAKE_Debug_Build})
    message("Build type: ${CMAKE_BUILD_TYPE}")
    message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
    message("CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
    message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
    message("CMAKE_CXX_STANDARD = ${CMAKE_CXX_STANDARD}")
endif()

# ****************************************************************************
# Include directories and library setup
# ****************************************************************************
# Collect source files
file(GLOB_RECURSE LIB_SOURCES "src/*.cpp")
file(GLOB_RECURSE LIB_HEADERS "src/*.h")

# Create our library
add_library(nodalfem STATIC ${LIB_SOURCES} ${LIB_HEADERS})

# Add include directories to the library target
target_include_directories(nodalfem PUBLIC
    ${LIBMESH_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH}
    ${PETSC_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/src
)

# Handle compiler flags and libraries
string(REPLACE " -pthread" "" CLEANED_LIBMESH_FLAGS "${LIBMESH_FLAGS}")
target_compile_options(nodalfem PUBLIC ${CLEANED_LIBMESH_FLAGS})
target_link_libraries(nodalfem PUBLIC
    ${LIBMESH_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    nlohmann_json::nlohmann_json
    Threads::Threads
)

# Function to add example executables
function(add_example EXAMPLE_NAME)
    add_executable(${EXAMPLE_NAME} examples/${EXAMPLE_NAME}/main.cpp)
    target_compile_options(${EXAMPLE_NAME} PUBLIC ${CLEANED_LIBMESH_FLAGS})
    target_link_libraries(${EXAMPLE_NAME} PUBLIC
        nodalfem
        Threads::Threads
    )
    # Set the RPATH for this executable
    set_target_properties(${EXAMPLE_NAME} PROPERTIES
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH "${LIBMESH_DIR}/lib;${PETSC_DIR}/lib"
    )
endfunction()

# Add all examples
file(GLOB EXAMPLE_DIRS "examples/*")
foreach(EXAMPLE_DIR ${EXAMPLE_DIRS})
    if(IS_DIRECTORY ${EXAMPLE_DIR})
        get_filename_component(EXAMPLE_NAME ${EXAMPLE_DIR} NAME)
        if(EXISTS ${EXAMPLE_DIR}/main.cpp)
            add_example(${EXAMPLE_NAME})
        endif()
    endif()
endforeach()

# Print debug information
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "LIBMESH_INCLUDE_DIRS = ${LIBMESH_INCLUDE_DIRS}")
message(STATUS "LIBMESH_LIBRARIES = ${LIBMESH_LIBRARIES}")
message(STATUS "LIBMESH_FLAGS = ${LIBMESH_FLAGS}")
message(STATUS "PETSC_LIBRARIES = ${PETSC_LIBRARIES}")
message(STATUS "MPI_LIBRARIES = ${MPI_LIBRARIES}")
message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}") 