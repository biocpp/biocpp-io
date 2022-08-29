message (STATUS "Configuring user doc.")

set (BIOCPP_IO_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (BIOCPP_IO_DOXYGEN_EXCLUDE_SYMBOLS "bio::io::detail bio::io::contrib") #/""
set (BIOCPP_IO_DOXYGEN_PREDEFINED_NDEBUG "-NDEBUG") #/""
set (BIOCPP_IO_DOXYGEN_ENABLED_SECTIONS "") #/"DEV"
set (BIOCPP_IO_DOXYGEN_EXTRACT_PRIVATE "NO") #/"YES":

configure_file (${BIOCPP_IO_DOXYFILE_IN} ${BIOCPP_IO_DOXYGEN_OUTPUT_DIR}/Doxyfile)

add_custom_target(doc_usr ALL
                  COMMAND ${DOXYGEN_EXECUTABLE}
                  WORKING_DIRECTORY ${BIOCPP_IO_DOXYGEN_OUTPUT_DIR}
                  BYPRODUCTS html/
                  DEPENDS download-cppreference-doxygen-web-tag
                  COMMENT "Generating user API documentation with Doxygen"
                  VERBATIM)

# Install doc_usr documentation in ./export folder
# make DESTDIR=export install
install (
    DIRECTORY "${BIOCPP_IO_DOXYGEN_OUTPUT_DIR}/html"
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
         COMMAND bash -c "${BIOCPP_IO_TEST_DOXYGEN_FAIL_ON_WARNINGS}")

#add_test(NAME doc_usr_no_uncond_requires_test
         #COMMAND bash -c "${BIOCPP_IO_TEST_DOXYGEN_FAIL_ON_UNCOND_REQUIRES}")
