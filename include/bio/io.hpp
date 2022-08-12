// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

/*!\file
 * \brief Provides The BioC++ I/O library version macros and global variables.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

//!\brief The major version as MACRO.
#define BIOCPP_IO_VERSION_MAJOR 0
//!\brief The minor version as MACRO.
#define BIOCPP_IO_VERSION_MINOR 1
//!\brief The patch version as MACRO.
#define BIOCPP_IO_VERSION_PATCH 0

//!\brief The full version as MACRO (number).
#define BIOCPP_IO_VERSION (BIOCPP_IO_VERSION_MAJOR * 10000 + BIOCPP_IO_VERSION_MINOR * 100 + BIOCPP_IO_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define BIOCPP_IO_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define BIOCPP_IO_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH)                                                     \
    BIOCPP_IO_VERSION_CSTRING_HELPER_STR(MAJOR)                                                                        \
    "." BIOCPP_IO_VERSION_CSTRING_HELPER_STR(MINOR) "." BIOCPP_IO_VERSION_CSTRING_HELPER_STR(PATCH)

//!\brief The full version as null terminated string.
#define BIOCPP_IO_VERSION_CSTRING                                                                                      \
    BIOCPP_IO_VERSION_CSTRING_HELPER_FUNC(BIOCPP_IO_VERSION_MAJOR, BIOCPP_IO_VERSION_MINOR, BIOCPP_IO_VERSION_PATCH)

namespace bio::io
{

//!\brief The major version.
constexpr uint8_t bio_version_major = BIOCPP_IO_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t bio_version_minor = BIOCPP_IO_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t bio_version_patch = BIOCPP_IO_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t bio_version = BIOCPP_IO_VERSION;

//!\brief The full version as null terminated string.
constexpr char const * bio_version_cstring = BIOCPP_IO_VERSION_CSTRING;

} // namespace bio::io

#undef BIOCPP_IO_VERSION_CSTRING_HELPER_STR
#undef BIOCPP_IO_VERSION_CSTRING_HELPER_FUNC
