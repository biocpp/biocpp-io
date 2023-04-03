# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

## Setup I/O library
include (../../build_system/biocpp_io-config-version.cmake)

project (biocpp_io_doc LANGUAGES NONE VERSION "${BIOCPP_IO_VERSION}")

if (NOT IS_DIRECTORY "${BIOCPP_IO_CLONE_DIR}/doc")
    message (FATAL_ERROR "Could not find BioC++ I/O library. Not building documentation.")
endif ()
message (STATUS "BioC++ I/O library: ${BIOCPP_IO_CLONE_DIR} (${BIOCPP_IO_VERSION})")

### Find core infrastructure which is needed and core content which is optional
include(../biocpp_require_core_infrastructure.cmake)
biocpp_require_core_infrastructure(FALSE) # false means that find_package is not called

if (NOT IS_DIRECTORY "${BIOCPP_CORE_CLONE_DIR}/doc")
    message (FATAL_ERROR "Could not find BioC++ Core. Not building documentation.")
endif ()
message (STATUS "BioC++ Core library: ${BIOCPP_CORE_CLONE_DIR} (${BIOCPP_CORE_VERSION})")

## Configure documentation build
option(BIOCPP_DOXY_CORE "Build documentation for core module." ON)
option(BIOCPP_DOXY_IO "Build documentation for core module." ON)

if (BIOCPP_DOXY_CORE)
    include (${BIOCPP_CORE_CLONE_DIR}/test/documentation/biocpp-doxygen-module.cmake)
endif()

if (BIOCPP_DOXY_IO)
    include (bio-doxygen-module.cmake)
endif()

include (${BIOCPP_CORE_CLONE_DIR}/test/documentation/biocpp-doxygen.cmake)
