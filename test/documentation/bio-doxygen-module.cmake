cmake_minimum_required (VERSION 3.7)

LIST(APPEND BIOCPP_DOXYGEN_MODULE_LAYOUT ${BIOCPP_IO_CLONE_DIR}/test/documentation/DoxygenLayout.xml.module)

set (BIOCPP_DOXYGEN_INCLUDE_PATH "${BIOCPP_DOXYGEN_INCLUDE_PATH} ${BIOCPP_IO_INCLUDE_DIR}")

set (BIOCPP_DOXYGEN_INPUT "${BIOCPP_DOXYGEN_INPUT}                      \
                           ${BIOCPP_IO_CLONE_DIR}/include               \
                           ${BIOCPP_IO_CLONE_DIR}/doc                   \
                           ${BIOCPP_IO_CLONE_DIR}/CHANGELOG.md")

set (BIOCPP_DOXYGEN_EXAMPLE_PATH "${BIOCPP_DOXYGEN_EXAMPLE_PATH} ${BIOCPP_IO_CLONE_DIR}")

set(BIOCPP_VERSION "${BIOCPP_VERSION} io-${BIOCPP_IO_VERSION}")
