# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2020-2022, deCODE Genetics
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------
#
# This CMake module will try to find The BioC++ I/O library and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package (biocpp_io [REQUIRED] ...)
#
# Since this makes a difference for CMAKE, pay attention to the case.
#
# The BioC++ I/O library has the following platform requirements:
#
#   C++20
#   pthread
#
# The BioC++ I/O library requires the following libraries:
#
#   SeqAn3    -- the succinct data structure library
#
# The BioC++ I/O library has the following optional dependencies:
#
#   ZLIB      -- zlib compression library
#   BZip2     -- libbz2 compression library
#
# If you don't wish for these to be detected (and used), you may define BIOCPP_IO_NO_ZLIB,
# BIOCPP_IO_NO_BZIP2 respectively.
#
# If you wish to require the presence of ZLIB or BZip2, just check for the module before
# finding The BioC++ I/O library, e.g. "find_package (ZLIB REQUIRED)".
#
# Once the search has been performed, the following variables will be set.
#
#   BIOCPP_IO_FOUND            -- Indicate whether The BioC++ I/O library was found and requirements met.
#
#   BIOCPP_IO_VERSION          -- The version as string, e.g. "3.0.0"
#   BIOCPP_IO_VERSION_MAJOR    -- e.g. 3
#   BIOCPP_IO_VERSION_MINOR    -- e.g. 0
#   BIOCPP_IO_VERSION_PATCH    -- e.g. 0
#
#   BIOCPP_IO_INCLUDE_DIRS     -- to be passed to include_directories ()
#   BIOCPP_IO_LIBRARIES        -- to be passed to target_link_libraries ()
#   BIOCPP_IO_DEFINITIONS      -- to be passed to add_definitions ()
#   BIOCPP_IO_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#
# Additionally, the following [IMPORTED][IMPORTED] targets are defined:
#
#   biocpp::io          -- interface target where
#                                  target_link_libraries(target biocpp::io)
#                              automatically sets
#                                  target_include_directories(target $BIOCPP_IO_INCLUDE_DIRS),
#                                  target_link_libraries(target $BIOCPP_IO_LIBRARIES),
#                                  target_compile_definitions(target $BIOCPP_IO_DEFINITIONS) and
#                                  target_compile_options(target $BIOCPP_IO_CXX_FLAGS)
#                              for a target.
#
#   [IMPORTED]: https://cmake.org/cmake/help/v3.10/prop_tgt/IMPORTED.html#prop_tgt:IMPORTED
#
# ============================================================================

cmake_minimum_required (VERSION 3.4...3.12)

# ----------------------------------------------------------------------------
# Set initial variables
# ----------------------------------------------------------------------------

# make output globally quiet if required by find_package, this effects cmake functions like `check_*`
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY})

# ----------------------------------------------------------------------------
# Greeter
# ----------------------------------------------------------------------------

string (ASCII 27 Esc)
set (ColourBold "${Esc}[1m")
set (ColourReset "${Esc}[m")

if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
    message (STATUS "${ColourBold}BioC++ I/O library:${ColourReset}")
endif ()

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (FindPackageHandleStandardArgs)

# ----------------------------------------------------------------------------
# Pretty printing and error handling
# ----------------------------------------------------------------------------

macro (bio_config_print text)
    if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
        message (STATUS "  ${text}")
    endif ()
endmacro ()

macro (bio_config_error text)
    if (${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
        message (FATAL_ERROR ${text})
    else ()
        if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message (WARNING ${text})
        endif ()
        return ()
    endif ()
endmacro ()

# ----------------------------------------------------------------------------
# Find The BioC++ I/O library include path
# ----------------------------------------------------------------------------

# Note that bio-config.cmake can be standalone and thus BIOCPP_IO_CLONE_DIR might be empty.
# * `BIOCPP_IO_CLONE_DIR` was already found in bio-config-version.cmake
# * `BIOCPP_IO_INCLUDE_DIR` was already found in bio-config-version.cmake
if (BIOCPP_IO_INCLUDE_DIR)
    bio_config_print ("Include dir found:          ${BIOCPP_IO_INCLUDE_DIR}")
else ()
    bio_config_error ("The include directory could not be found (BIOCPP_IO_INCLUDE_DIR: '${BIOCPP_IO_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Force-deactivate optional dependencies
# ----------------------------------------------------------------------------

# These two are "opt-in", because detected by CMake
# If you want to force-require these, just do find_package (zlib REQUIRED) before find_package(biocpp_io)
option (BIOCPP_IO_NO_ZLIB  "Don't use ZLIB, even if present." OFF)
option (BIOCPP_IO_NO_BZIP2 "Don't use BZip2, even if present." OFF)

# ----------------------------------------------------------------------------
# thread support (pthread, windows threads)
# ----------------------------------------------------------------------------

set (THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads QUIET)

if (Threads_FOUND)
    set (BIOCPP_IO_LIBRARIES ${BIOCPP_IO_LIBRARIES} Threads::Threads)
    if ("${CMAKE_THREAD_LIBS_INIT}" STREQUAL "")
        bio_config_print ("Thread support:             builtin.")
    else ()
        bio_config_print ("Thread support:             via ${CMAKE_THREAD_LIBS_INIT}")
    endif ()
else ()
    bio_config_print ("Thread support:             not found.")
endif ()

# ----------------------------------------------------------------------------
# Require Core library
# ----------------------------------------------------------------------------

if (BIOCPP_CORE_FOUND)
    bio_config_print ("Required dependency:        already loaded (${BIOCPP_CORE_VERSION}).")
else ()
    find_package (biocpp_core REQUIRED QUIET
                  HINTS ${BIOCPP_IO_CLONE_DIR}/submodule/biocpp-core/build_system
                        ${CMAKE_CURRENT_LIST_DIR}/../../biocpp-core/build_system)

    if (BIOCPP_CORE_FOUND)
        bio_config_print ("Required dependency:        BioC++ core library found (${BIOCPP_CORE_VERSION}).")
        set (BIOCPP_IO_LIBRARIES ${BIOCPP_IO_LIBRARIES} biocpp_core)
    else ()
        bio_config_print ("The BioC++ core library is required, but wasn't found. Get it from https://github.com/biocpp/biocpp-core")
    endif ()
endif ()

# ----------------------------------------------------------------------------
# Require SeqAn3
# ----------------------------------------------------------------------------

find_package (SeqAn3 REQUIRED QUIET
              HINTS ${CMAKE_CURRENT_LIST_DIR}/../submodules/seqan3/build_system)

if (SEQAN3_FOUND)
    bio_config_print ("Required dependency:        SeqAn3 found (${SEQAN3_VERSION}).")
    set (BIOCPP_IO_LIBRARIES ${BIOCPP_IO_LIBRARIES} seqan3_seqan3)
else ()
    bio_config_print ("The SeqAn3 library is required, but wasn't found. Get it from https://github.com/seqan/seqan3")
endif ()

# ----------------------------------------------------------------------------
# ZLIB dependency
# ----------------------------------------------------------------------------

if (NOT BIOCPP_IO_NO_ZLIB)
    find_package (ZLIB QUIET)
endif ()

if (ZLIB_FOUND)
    set (BIOCPP_IO_LIBRARIES         ${BIOCPP_IO_LIBRARIES}         ${ZLIB_LIBRARIES})
    set (BIOCPP_IO_DEFINITIONS       ${BIOCPP_IO_DEFINITIONS}       "-DBIOCPP_IO_HAS_ZLIB=1")
    bio_config_print ("Optional dependency:        ZLIB found (${ZLIB_VERSION_STRING}).")
else ()
    bio_config_print ("Optional dependency:        ZLIB not found.")
endif ()

# ----------------------------------------------------------------------------
# BZip2 dependency
# ----------------------------------------------------------------------------

if (NOT BIOCPP_IO_NO_BZIP2)
    find_package (BZip2 QUIET)
endif ()

if (NOT ZLIB_FOUND AND BZIP2_FOUND)
    # NOTE (marehr): iostream_bzip2 uses the type `uInt`, which is defined by
    # `zlib`. Therefore, `bzip2` will cause a ton of errors without `zlib`.
    message (AUTHOR_WARNING "Disabling BZip2 [which was successfully found], "
                            "because ZLIB was not found. BZip2 depends on ZLIB.")
    unset (BZIP2_FOUND)
endif ()

if (BZIP2_FOUND)
    set (BIOCPP_IO_LIBRARIES         ${BIOCPP_IO_LIBRARIES}         ${BZIP2_LIBRARIES})
    set (BIOCPP_IO_DEFINITIONS       ${BIOCPP_IO_DEFINITIONS}       "-DBIOCPP_IO_HAS_BZIP2=1")
    bio_config_print ("Optional dependency:        BZip2 found (${BZIP2_VERSION_STRING}).")
else ()
    bio_config_print ("Optional dependency:        BZip2 not found.")
endif ()

# ----------------------------------------------------------------------------
# Finish find_package call
# ----------------------------------------------------------------------------

find_package_handle_standard_args (${CMAKE_FIND_PACKAGE_NAME} REQUIRED_VARS BIOCPP_IO_INCLUDE_DIR)

# Set BIOCPP_IO_* variables with the content of ${CMAKE_FIND_PACKAGE_NAME}_(FOUND|...|VERSION)
# This needs to be done, because `find_package(biocpp_io)` might be called in any case-sensitive way and we want to
# guarantee that BIOCPP_IO_* are always set.
foreach (package_var FOUND DIR ROOT CONFIG VERSION VERSION_MAJOR VERSION_MINOR VERSION_PATCH VERSION_TWEAK VERSION_COUNT)
    set (BIOCPP_IO_${package_var} "${${CMAKE_FIND_PACKAGE_NAME}_${package_var}}")
endforeach ()

# propagate BIOCPP_IO_INCLUDE_DIR into BIOCPP_IO_INCLUDE_DIRS
set (BIOCPP_IO_INCLUDE_DIRS ${BIOCPP_IO_INCLUDE_DIR} ${BIOCPP_IO_DEPENDENCY_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# Export targets
# ----------------------------------------------------------------------------

if (BIOCPP_IO_FOUND AND NOT TARGET biocpp::io)
    add_library (biocpp_io INTERFACE)
    target_compile_definitions (biocpp_io INTERFACE ${BIOCPP_IO_DEFINITIONS})
    target_link_libraries (biocpp_io INTERFACE "${BIOCPP_IO_LIBRARIES}")
    # include bio/include/ as -I, because bio should never produce warnings.
    target_include_directories (biocpp_io INTERFACE "${BIOCPP_IO_INCLUDE_DIR}")
    # include everything except bio/include/ as -isystem, i.e.
    # a system header which suppresses warnings of external libraries. CURRENTLY EMPTY
    #target_include_directories (biocpp_io SYSTEM INTERFACE "${BIOCPP_IO_DEPENDENCY_INCLUDE_DIRS}")
    add_library (biocpp::io ALIAS biocpp_io)
endif ()

set (CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if (BIOCPP_IO_FIND_DEBUG)
  message ("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
  message ("")
  message ("  CMAKE_BUILD_TYPE            ${CMAKE_BUILD_TYPE}")
  message ("  CMAKE_SOURCE_DIR            ${CMAKE_SOURCE_DIR}")
  message ("  CMAKE_INCLUDE_PATH          ${CMAKE_INCLUDE_PATH}")
  message ("  BIOCPP_IO_INCLUDE_DIR          ${BIOCPP_IO_INCLUDE_DIR}")
  message ("")
  message ("  ${CMAKE_FIND_PACKAGE_NAME}_FOUND                ${${CMAKE_FIND_PACKAGE_NAME}_FOUND}")
  message ("  BIOCPP_IO_HAS_ZLIB             ${ZLIB_FOUND}")
  message ("  BIOCPP_IO_HAS_BZIP2            ${BZIP2_FOUND}")
  message ("")
  message ("  BIOCPP_IO_INCLUDE_DIRS         ${BIOCPP_IO_INCLUDE_DIRS}")
  message ("  BIOCPP_IO_LIBRARIES            ${BIOCPP_IO_LIBRARIES}")
  message ("  BIOCPP_IO_DEFINITIONS          ${BIOCPP_IO_DEFINITIONS}")
  message ("  BIOCPP_IO_CXX_FLAGS            ${BIOCPP_IO_CXX_FLAGS}")
  message ("")
  message ("  BIOCPP_IO_VERSION              ${BIOCPP_IO_VERSION}")
  message ("  BIOCPP_IO_VERSION_MAJOR        ${BIOCPP_IO_VERSION_MAJOR}")
  message ("  BIOCPP_IO_VERSION_MINOR        ${BIOCPP_IO_VERSION_MINOR}")
  message ("  BIOCPP_IO_VERSION_PATCH        ${BIOCPP_IO_VERSION_PATCH}")
endif ()
