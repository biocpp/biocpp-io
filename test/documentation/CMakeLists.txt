# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

### Find b.i.o
include (../../build_system/bio-config-version.cmake)
set (BIO_VERSION "${PACKAGE_VERSION}")

project (bio LANGUAGES NONE VERSION "${BIO_VERSION}")

if (NOT EXISTS "${BIO_INCLUDE_DIR}/bio/version.hpp")
    message (FATAL_ERROR "Could not find b.i.o. Not building documentation.")
endif ()
message (STATUS "Found: ${BIO_INCLUDE_DIR} (${BIO_VERSION})")

set (BIO_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (bio-doxygen.cmake)

enable_testing ()

if (BIO_USER_DOC)
    add_subdirectory(doc_usr)
endif ()

if (BIO_DEV_DOC)
    add_subdirectory(doc_dev)
endif ()
