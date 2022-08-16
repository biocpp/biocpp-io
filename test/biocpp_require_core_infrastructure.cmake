# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/biocpp/biocpp-core/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Make sure that ${BIOCPP_CORE_CLONE_DIR} is available so we can access core's cmake-stuff

macro (biocpp_require_core_infrastructure)

    if (NOT IS_DIRECTORY ${BIOCPP_CORE_CLONE_DIR})
        if (IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../../submodule/biocpp-core/build_system OR
            IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/../../../biocpp-core/build_system)
            find_package (biocpp_core REQUIRED
                          HINTS ${CMAKE_CURRENT_LIST_DIR}/../../submodule/biocpp-core/build_system
                          HINTS ${CMAKE_CURRENT_LIST_DIR}/../../../biocpp-core/build_system)
        else ()
            include (FetchContent)
            FetchContent_Declare(
                biocpp_core-lib
                GIT_REPOSITORY https://github.com/biocpp/biocpp-core
                GIT_TAG main
            )

            FetchContent_Populate(biocpp_core-lib)
            find_package (biocpp_core REQUIRED
                          HINTS ${CMAKE_CURRENT_BINARY_DIR}/_deps/biocpp_core-lib-src/build_system)

        endif ()
    endif ()

    if (NOT IS_DIRECTORY ${BIOCPP_CORE_CLONE_DIR})
        message (FATAL_ERROR "BioC++ core library found, but infrastructure setup failed.")
    endif ()

endmacro ()
