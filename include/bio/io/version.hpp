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
 * \brief Provides B.I.O. version macros and global variables.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

//!\brief The major version as MACRO.
#define BIO_VERSION_MAJOR 0
//!\brief The minor version as MACRO.
#define BIO_VERSION_MINOR 1
//!\brief The patch version as MACRO.
#define BIO_VERSION_PATCH 0

//!\brief The full version as MACRO (number).
#define BIO_VERSION (BIO_VERSION_MAJOR * 10000 + BIO_VERSION_MINOR * 100 + BIO_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define BIO_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define BIO_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH)                                                           \
    BIO_VERSION_CSTRING_HELPER_STR(MAJOR)                                                                              \
    "." BIO_VERSION_CSTRING_HELPER_STR(MINOR) "." BIO_VERSION_CSTRING_HELPER_STR(PATCH)

//!\brief The full version as null terminated string.
#define BIO_VERSION_CSTRING BIO_VERSION_CSTRING_HELPER_FUNC(BIO_VERSION_MAJOR, BIO_VERSION_MINOR, BIO_VERSION_PATCH)

namespace bio
{

//!\brief The major version.
constexpr uint8_t bio_version_major = BIO_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t bio_version_minor = BIO_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t bio_version_patch = BIO_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t bio_version = BIO_VERSION;

//!\brief The full version as null terminated string.
constexpr char const * bio_version_cstring = BIO_VERSION_CSTRING;

} // namespace bio

#undef BIO_VERSION_CSTRING_HELPER_STR
#undef BIO_VERSION_CSTRING_HELPER_FUNC
