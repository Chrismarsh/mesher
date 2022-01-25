###############################################################################
# CMake module to search for func library
#
# METIS_ROOT = install prefix to search

# On success, the macro sets the following variables:
# METIS_FOUND       = if the library found
# METIS_LIBRARY     = full path to the library
# METIS_INCLUDE_DIR = where to find the library headers
# also defined, but not for general use are
# METIS_LIBRARY, where to find the func library.
#
# Copyright (c) 2009 Mateusz Loskot <mateusz@loskot.net>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
###############################################################################
include (FindPackageHandleStandardArgs)

if( DEFINED ENV{Metis_DIR} )
    set( Metis_DIR "$ENV{Metis_DIR}" )
endif()

set(METIS_FOUND ON)

find_path(METIS_INCLUDE_DIR
        include/metis.h
        PATHS ${Metis_DIR}/include
        DOC "Include for metis"
        )
find_library(METIS_LIBRARY
        NAMES metis
        PATHS ${Metis_DIR}/lib
        )

find_package_handle_standard_args(Metis DEFAULT_MSG
        METIS_INCLUDE_DIR METIS_LIBRARY )


if(Metis_FOUND)
    set( Metis_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
    set( Metis_LIBRARIES ${METIS_LIBRARY} )

    mark_as_advanced(
            METIS_INCLUDE_DIR
            METIS_LIBRARY
    )

    add_library(metis::metis INTERFACE IMPORTED)
    set_target_properties(metis::metis PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Metis_INCLUDE_DIRS}")
    set_property(TARGET metis::metis PROPERTY INTERFACE_LINK_LIBRARIES "${Metis_LIBRARIES}")
endif()
