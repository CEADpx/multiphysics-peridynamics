# -------------------------------------------
# Copyright (c) 2021 - 2025 Prashant K. Jha
# -------------------------------------------
# PeriDEM https://github.com/prashjha/PeriDEM
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE)
find_package(PkgConfig)

find_library(PETSC_LIBRARIES
        NAMES libpetsc.so libpetsc.dylib
        HINTS "${PETSC_LIB}")
#        HINTS /usr/lib64 /usr/local/lib64 /usr/lib/ /usr/local/lib
#"${PETSC_DIR}/lib/")


if (NOT PETSC_LIBRARIES)
    message(FATAL_ERROR "PETSC_LIBRARIES Library not found: Specify the PETSC_DIR where petsc is located")
endif ()
