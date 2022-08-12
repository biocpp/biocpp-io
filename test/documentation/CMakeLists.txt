# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

### Find b.i.o
include (../../build_system/biocpp_io-config-version.cmake)
set (BIOCPP_IO_VERSION "${PACKAGE_VERSION}")

project (biocpp_io_doc LANGUAGES NONE VERSION "${BIOCPP_IO_VERSION}")

if (NOT EXISTS "${BIOCPP_IO_INCLUDE_DIR}/bio/io.hpp")
    message (FATAL_ERROR "Could not find the BioC++ I/O library. Not building documentation.")
endif ()
message (STATUS "Found: ${BIOCPP_IO_INCLUDE_DIR} (${BIOCPP_IO_VERSION})")

set (BIOCPP_IO_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

#TODO: use the .cmake file from the core repo and remove this one
#TODO: also just use the doxygen_cfg.in if possible
include (bio-doxygen.cmake)

enable_testing ()

if (BIOCPP_IO_USER_DOC)
    add_subdirectory(doc_usr)
endif ()

if (BIOCPP_IO_DEV_DOC)
    add_subdirectory(doc_dev)
endif ()
