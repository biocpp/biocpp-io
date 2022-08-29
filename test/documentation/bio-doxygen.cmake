cmake_minimum_required (VERSION 3.7)

### Find doxygen and dependency to DOT tool
message (STATUS "Searching for doxygen.")
find_package (Doxygen REQUIRED)

if (NOT ${DOXYGEN_FOUND})
    message (FATAL_ERROR "Could not find doxygen. Not building documentation.")
endif ()

if (NOT ${DOXYGEN_DOT_FOUND})
    message (STATUS "Could not find dot tool. Disabling dot support.")
    set (BIOCPP_IO_DOXYGEN_HAVE_DOT "NO")
else ()
    message (STATUS "Found dot tool. Enabling dot support.")
    set (BIOCPP_IO_DOXYGEN_HAVE_DOT "YES")
endif ()

### Use mathjax instead of latex to render formulas.
set (BIOCPP_IO_DOXYGEN_USE_MATHJAX "NO")

### Number of threads to use for dot. Doxygen's default is 0 (all threads).
set (BIOCPP_IO_DOXYGEN_DOT_NUM_THREADS "0")

### Configure doc/developer targets.
set (BIOCPP_IO_DOXYGEN_SOURCE_DIR "${BIOCPP_IO_CLONE_DIR}")
set (BIOCPP_IO_DOXYFILE_IN ${BIOCPP_IO_DOXYGEN_INPUT_DIR}/bio_doxygen_cfg.in)

option(BIOCPP_IO_USER_DOC "Create build target and test for user documentation." ON)
option(BIOCPP_IO_DEV_DOC "Create build target and test for developer documentation." ON)

### Download and extract cppreference-doxygen-web.tag.xml for std:: documentation links
set(BIOCPP_IO_DOXYGEN_STD_TAGFILE "${PROJECT_BINARY_DIR}/cppreference-doxygen-web.tag.xml")
include(ExternalProject)
ExternalProject_Add (
    download-cppreference-doxygen-web-tag
    URL "https://github.com/PeterFeicht/cppreference-doc/releases/download/v20201016/html-book-20201016.tar.xz"
    URL_HASH SHA256=35d67ceb114b91d9220e7db81d00cae80f74768729b21b369bf2f17b401cbdc0
    TLS_VERIFY ON
    DOWNLOAD_DIR "${PROJECT_BINARY_DIR}"
    DOWNLOAD_NAME "html-book.tar.xz"
    DOWNLOAD_NO_EXTRACT YES
    BINARY_DIR "${PROJECT_BINARY_DIR}"
    CONFIGURE_COMMAND /bin/sh -c "xzcat html-book.tar.xz | tar -xf - cppreference-doxygen-web.tag.xml"
    BUILD_COMMAND rm "html-book.tar.xz"
    INSTALL_COMMAND ""
)

### TEST HELPER

# doxygen does not show any warnings (doxygen prints warnings / errors to cerr)
set (BIOCPP_IO_TEST_DOXYGEN_FAIL_ON_WARNINGS "
    ${DOXYGEN_EXECUTABLE} > doxygen.cout 2> doxygen.cerr;
    cat \"doxygen.cerr\";
    test ! -s \"doxygen.cerr\"")

# We search the HTML output to ensure that no `requires` clauses are at wrong places.
set (BIOCPP_IO_TEST_DOXYGEN_FAIL_ON_UNCOND_REQUIRES
     "! find . -not -name \"*_source.html\" -name \"*.html\" -print0 | xargs -0 grep \"requires\" | grep \"memname\"")


### install helper

# make sure that prefix path is /usr/local/share/doc/bio/
if (NOT DEFINED CMAKE_SIZEOF_VOID_P)
    # we need this to suppress GNUInstallDirs AUTHOR_WARNING:
    #   CMake Warning (dev) at /usr/share/cmake-3.19/Modules/GNUInstallDirs.cmake:223 (message):
    #     Unable to determine default CMAKE_INSTALL_LIBDIR directory because no
    #     target architecture is known.  Please enable at least one language before
    #     including GNUInstallDirs.
    set (CMAKE_SIZEOF_VOID_P 8)
endif ()
include (GNUInstallDirs) # this is needed to prefix the install paths
