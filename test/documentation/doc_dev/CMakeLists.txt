message (STATUS "Configuring devel doc.")

set (BIOCPP_IO_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (BIOCPP_IO_DOXYGEN_EXCLUDE_SYMBOLS "")
set (BIOCPP_IO_DOXYGEN_PREDEFINED_NDEBUG "")
set (BIOCPP_IO_DOXYGEN_ENABLED_SECTIONS "DEV")
set (BIOCPP_IO_DOXYGEN_EXTRACT_PRIVATE "YES")

configure_file (${BIOCPP_IO_DOXYFILE_IN} ${BIOCPP_IO_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${BIOCPP_IO_FOOTER_HTML_IN} ${BIOCPP_IO_DOXYGEN_OUTPUT_DIR}/bio_footer.html)

add_custom_target(doc_dev ALL
                  COMMAND ${DOXYGEN_EXECUTABLE}
                  WORKING_DIRECTORY ${BIOCPP_IO_DOXYGEN_OUTPUT_DIR}
                  BYPRODUCTS html/
                  DEPENDS download-cppreference-doxygen-web-tag
                  COMMENT "Generating developer API documentation with Doxygen"
                  VERBATIM)

# Install doc_dev documentation in ./install_doc_dev folder
# cmake --install . --prefix install_doc_dev --component doc_dev
install (
    DIRECTORY "${BIOCPP_IO_DOXYGEN_OUTPUT_DIR}/html"
    DESTINATION "${CMAKE_INSTALL_DOCDIR}"
    COMPONENT doc-dev
    EXCLUDE_FROM_ALL
    PATTERN "*.md5" EXCLUDE
    PATTERN "*.map" EXCLUDE
    PATTERN "formula.repository" EXCLUDE
)

### Enable testing

enable_testing()

message (STATUS "Add dev doc test.")
add_test(NAME doc_dev_no_warnings_test
         COMMAND bash -c "${BIOCPP_IO_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME doc_dev_no_uncond_requires_test
         COMMAND bash -c "${BIOCPP_IO_TEST_DOXYGEN_FAIL_ON_UNCOND_REQUIRES}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
