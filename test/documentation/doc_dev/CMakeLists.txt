message (STATUS "Configuring devel doc.")

set (BIO_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (BIO_DOXYGEN_EXCLUDE_SYMBOLS "")
set (BIO_DOXYGEN_PREDEFINED_NDEBUG "")
set (BIO_DOXYGEN_ENABLED_SECTIONS "DEV")
set (BIO_DOXYGEN_EXTRACT_PRIVATE "YES")

configure_file (${BIO_DOXYFILE_IN} ${BIO_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${BIO_FOOTER_HTML_IN} ${BIO_DOXYGEN_OUTPUT_DIR}/bio_footer.html)

add_custom_target(doc_dev ALL
                  COMMAND ${DOXYGEN_EXECUTABLE}
                  WORKING_DIRECTORY ${BIO_DOXYGEN_OUTPUT_DIR}
                  BYPRODUCTS html/
                  DEPENDS download-cppreference-doxygen-web-tag
                  COMMENT "Generating developer API documentation with Doxygen"
                  VERBATIM)

# Install doc_dev documentation in ./install_doc_dev folder
# cmake --install . --prefix install_doc_dev --component doc_dev
install (
    DIRECTORY "${BIO_DOXYGEN_OUTPUT_DIR}/html"
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
         COMMAND bash -c "${BIO_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME doc_dev_no_uncond_requires_test
         COMMAND bash -c "${BIO_TEST_DOXYGEN_FAIL_ON_UNCOND_REQUIRES}"
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
