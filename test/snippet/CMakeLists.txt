# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)
project (biocpp_io_test_snippet CXX)

include(../biocpp_require_core_infrastructure.cmake)
biocpp_require_core_infrastructure()

find_package (biocpp_io REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build_system)

include (${BIOCPP_CORE_CLONE_DIR}/test/biocpp-test.cmake)

add_library (snippet_main snippet_main.cpp)
target_link_libraries (snippet_main PUBLIC biocpp::test gtest)

#TODO: move this to core library and import
macro (biocpp_snippet test_name_prefix snippet snippet_base_path)
    biocpp_test_component (snippet_target_name "${snippet}" TARGET_NAME)
    biocpp_test_component (snippet_test_name "${snippet}" TEST_NAME)
    biocpp_test_component (snippet_target_path "${snippet}" TARGET_PATH)
    set (target "${snippet_target_name}_snippet")

    add_executable (${target} "${snippet_base_path}/${snippet}")
    target_link_libraries (${target} PUBLIC snippet_main biocpp::io fmt::fmt-header-only)
    target_compile_options (${target} PUBLIC -Wno-unused)
    set_target_properties(${target}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${snippet_target_path}"
    )

    set (snippet_test_target "${test_name_prefix}/${snippet_test_name}_snippet")
    add_test (NAME "${snippet_test_target}" COMMAND ${target})

    set (snippet_compare_test_target "${test_name_prefix}/${snippet_test_name}_snippet_cmp_output")

    add_test (NAME "${snippet_compare_test_target}"
              COMMAND ${CMAKE_COMMAND}
                      -DTARGET_FILE=$<TARGET_FILE:${target}>
                      -DSOURCE_FILE=${snippet_base_path}/${snippet}
                      -P "${BIOCPP_CORE_CLONE_DIR}/test/snippet/compare_snippet_output.cmake")

    unset (snippet_target_name)
    unset (snippet_test_name)
    unset (snippet_target_path)
    unset (target)
endmacro ()

macro (biocpp_snippets test_name_prefix snippet_base_path)
    biocpp_test_files (snippet_files "${snippet_base_path}" "*.cpp")
    foreach (snippet ${snippet_files})
        biocpp_snippet("${test_name_prefix}" "${snippet}" "${snippet_base_path}")
    endforeach ()
endmacro ()

biocpp_require_ccache ()
biocpp_require_fmt ()
biocpp_require_test ()

biocpp_snippets ("snippet" "${CMAKE_CURRENT_SOURCE_DIR}")
#bio_snippets ("doc/snippet" "${CMAKE_CURRENT_SOURCE_DIR}/../../doc")

