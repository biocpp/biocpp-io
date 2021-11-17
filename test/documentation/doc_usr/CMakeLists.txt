message (STATUS "Configuring user doc.")

set (BIO_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
#set (BIO_DOXYGEN_EXCLUDE_SYMBOLS "detail bio::simd") #/""
set (BIO_DOXYGEN_PREDEFINED_NDEBUG "-NDEBUG") #/""
set (BIO_DOXYGEN_ENABLED_SECTIONS "") #/"DEV"
set (BIO_DOXYGEN_EXTRACT_PRIVATE "NO") #/"YES":

configure_file (${BIO_DOXYFILE_IN} ${BIO_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${BIO_FOOTER_HTML_IN} ${BIO_DOXYGEN_OUTPUT_DIR}/bio_footer.html)

add_custom_target(doc_usr ALL
                  COMMAND ${DOXYGEN_EXECUTABLE}
                  WORKING_DIRECTORY ${BIO_DOXYGEN_OUTPUT_DIR}
                  BYPRODUCTS html/
                  DEPENDS download-cppreference-doxygen-web-tag
                  COMMENT "Generating user API documentation with Doxygen"
                  VERBATIM)

# Install doc_usr documentation in ./export folder
# make DESTDIR=export install
install (
    DIRECTORY "${BIO_DOXYGEN_OUTPUT_DIR}/html"
    DESTINATION "${CMAKE_INSTALL_DOCDIR}"
    COMPONENT doc
    PATTERN "*.md5" EXCLUDE
    PATTERN "*.map" EXCLUDE
    PATTERN "formula.repository" EXCLUDE
)

### Enable testing

enable_testing()

message (STATUS "Add user doc test.")
add_test(NAME doc_usr_no_warnings_test
         COMMAND bash -c "${BIO_TEST_DOXYGEN_FAIL_ON_WARNINGS}")

add_test(NAME doc_usr_no_uncond_requires_test
         COMMAND bash -c "${BIO_TEST_DOXYGEN_FAIL_ON_UNCOND_REQUIRES}")
