# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)
project (biocpp_io_test_unit CXX)

include(../biocpp_require_core_infrastructure.cmake)
biocpp_require_core_infrastructure()

find_package (biocpp_io REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build_system)

include (${BIOCPP_CORE_CLONE_DIR}/test/biocpp-test.cmake)
include (GoogleTest OPTIONAL)

option (BIOCPP_IO_VERBOSE_TESTS "Run each test case individually" OFF)
option (BIOCPP_IO_USE_INCLUDE_DEPENDENCIES "Build tests in an hierarchical order (by an include graph, i.e. tests with less dependencies are build first)" OFF)

macro (bio_test unit_test_cpp)
    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")
    biocpp_test_component (target "${unit_test}" TARGET_NAME)
    biocpp_test_component (test_name "${unit_test}" TEST_NAME)

    add_executable (${target} ${unit_test_cpp})
    target_link_libraries (${target} biocpp::test::unit biocpp::io)
    if (BIOCPP_IO_VERBOSE_TESTS)
        gtest_discover_tests(${target} TEST_PREFIX "${test_name}::" PROPERTIES TIMEOUT "30")
    else ()
        add_test (NAME "${test_name}" COMMAND ${target})
    endif ()

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

biocpp_require_ccache ()
biocpp_require_test ()

add_subdirectories ()

# Top-level tests

bio_test(genomic_region_test.cpp)
