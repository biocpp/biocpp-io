// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <cinttypes>
#include <ciso646> // makes _LIBCPP_VERSION available
#include <cstddef> // makes __GLIBCXX__ available

/*!\file
 * \brief Provides platform and dependency checks.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

// ============================================================================
//  C++ standard and features
// ============================================================================

// C++ standard [required]
#ifdef __cplusplus
static_assert(__cplusplus >= 201709, "The BioC++ I/O library requires C++20, make sure that you have set -std=c++20.");
#else
#    error "This is not a C++ compiler."
#endif

#if __has_include(<version>)
#    include <version>
#endif

// ============================================================================
//  Dependencies
// ============================================================================

// The BioC++ I/O library [required]
#if __has_include(<bio/io.hpp>)
#    include <bio/io.hpp>
#else
#    error The BioC++ I/O library include directory is not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// SeqAn3 [required]
#if __has_include(<seqan3/version.hpp>)
#    include <seqan3/version.hpp>
static_assert(seqan3::seqan3_version_major == 3, "SeqAn >= 3.1 is required by the BioC++ I/O library.");
static_assert(seqan3::seqan3_version_minor >= 1, "SeqAn >= 3.1 is required by the BioC++ I/O library.");
#else
#    error The SeqAn3 library was not included.
#endif

// ============================================================================
//  WORKAROUNDS
// ============================================================================

//!\brief Bugs in GCC lead to our "readers" not working when declared as views.
#ifndef BIOCPP_IO_NO_VIEWBASE
#    if defined(__GNUC__) && ((__GNUC__ == 10 && __GNUC_MINOR__ == 4) || (__GNUC__ == 11 && __GNUC_MINOR__ < 3) || (__GNUC__ == 12 && __GNUC_MINOR__ < 2))
#        define BIOCPP_IO_NO_VIEWBASE 1
#    else
#        define BIOCPP_IO_NO_VIEWBASE 0
#    endif
#endif

// ============================================================================
//  Documentation
// ============================================================================

// Doxygen related
// this macro is a NO-OP unless doxygen parses it, in which case it resolves to the argument
#ifndef BIOCPP_IO_DOXYGEN_ONLY
#    define BIOCPP_IO_DOXYGEN_ONLY(x)
#endif
