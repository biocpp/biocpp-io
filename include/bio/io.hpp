// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

/*!\file
 * \brief Provides The BioC++ I/O library version macros and global variables.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

// ============================================================================
//  Dependencies
// ============================================================================

// Self [required]
#if __has_include(<bio/io.hpp>)
#    include <bio/io.hpp>
#else
#    error The BioC++ I/O library include directory is not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// The BioC++ core library [required]
#if __has_include(<bio/core.hpp>)
#    include <bio/core.hpp>
#else
#    error Could not find the BioC++ core library. Please add its include directory!
#endif

static_assert(bio::biocpp_core_version_major == 0 && bio::biocpp_core_version_minor == 7,
              "This version of the BioC++ I/O module requires version 0.7.x of the BioC++ core library.");

// ============================================================================
//  VERSION
// ============================================================================

//!\brief The major version as MACRO.
#define BIOCPP_IO_VERSION_MAJOR 0
//!\brief The minor version as MACRO.
#define BIOCPP_IO_VERSION_MINOR 3
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

namespace bio
{

//!\brief The major version.
constexpr std::size_t io_version_major = BIOCPP_IO_VERSION_MAJOR;
//!\brief The minor version.
constexpr std::size_t io_version_minor = BIOCPP_IO_VERSION_MINOR;
//!\brief The patch version.
constexpr std::size_t io_version_patch = BIOCPP_IO_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t io_version = BIOCPP_IO_VERSION;

//!\brief The full version as null terminated string.
constexpr std::string_view bio_version_cstring = BIOCPP_IO_VERSION_CSTRING;

} // namespace bio

#undef BIOCPP_IO_VERSION_CSTRING_HELPER_STR
#undef BIOCPP_IO_VERSION_CSTRING_HELPER_FUNC

// ============================================================================
//  WORKAROUNDS
// ============================================================================

// ============================================================================
//  NAMESPACES
// ============================================================================

/*!\namespace bio::io
 * \brief Main namespace for the I/O module.
 */
namespace bio::io
{}

/*!\if DEV
 * \namespace bio::io::detail
 * \brief The internal namespace.
 * \details
 * The contents of this namespace are not visible to consumers of the library and the documentation is
 * only generated for developers.
 * \endif
 */
namespace bio::io::detail
{}

/*!\namespace bio::io::txt
 * \brief Namespace for the Plain I/O submodule.
 */
namespace bio::io::txt
{}

/*!\if DEV
 * \namespace bio::io::txt::detail
 * \brief The txt internal namespace.
 * \details
 * The contents of this namespace are not visible to consumers of the library and the documentation is
 * only generated for developers.
 * \endif
 */
namespace bio::io::txt::detail
{}

/*!\namespace bio::io::seq
 * \brief Namespace for the Seq I/O submodule.
 */
namespace bio::io::seq
{}

/*!\if DEV
 * \namespace bio::io::seq::detail
 * \brief The seq internal namespace.
 * \details
 * The contents of this namespace are not visible to consumers of the library and the documentation is
 * only generated for developers.
 * \endif
 */
namespace bio::io::seq::detail
{}

/*!\namespace bio::io::var
 * \brief Namespace for the Var I/O module.
 */
namespace bio::io::var
{}

/*!\if DEV
 * \namespace bio::io::var::detail
 * \brief The var internal namespace.
 * \details
 * The contents of this namespace are not visible to consumers of the library and the documentation is
 * only generated for developers.
 * \endif
 */
namespace bio::io::var::detail
{}
