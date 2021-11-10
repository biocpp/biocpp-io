if (BZIP2_FOUND)
    bio_test(bz2_istream_test.cpp)
    bio_test(bz2_ostream_test.cpp)
endif ()

if (ZLIB_FOUND)
    bio_test(gz_istream_test.cpp)
    bio_test(gz_ostream_test.cpp)

    bio_test(bgzf_istream_test.cpp)
    bio_test(bgzf_ostream_test.cpp)
endif ()

bio_test(transparent_istream_test.cpp)
bio_test(transparent_ostream_test.cpp)
