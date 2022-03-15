// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides exception types.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <ios>
#include <stdexcept>

#include <bio/detail/to_string.hpp>

namespace bio
{

/*!\addtogroup bio
 * \{
 */

// ----------------------------------------------------------------------------
// generic exceptions
// ----------------------------------------------------------------------------

//!\brief All other exceptions inherit from this.
struct bio_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    explicit bio_error(auto &&... s) : std::runtime_error{(detail::to_string(s) + ...)} {}
};

//!\brief All other exceptions inherit from this.
struct unreachable_code : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit unreachable_code(auto &&... s) :
      bio_error{"Unreachable code reached.\nPlease report a bug with this message.\nDetails:\n",
                (detail::to_string(s) + ...)}
    {}
    // TODO(GCC11): When GCC10 is dropped, make the constructor take std::source_location instead.
};

// ----------------------------------------------------------------------------
// file open exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if there is no format that accepts a given file extension.
struct unhandled_extension_error : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit unhandled_extension_error(auto &&... s) : bio_error{s...} {}
};

//!\brief Thrown if there is an unspecified filesystem or stream error while opening, e.g. permission problem.
struct file_open_error : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit file_open_error(auto &&... s) : bio_error{s...} {}
};

//!\brief Thrown if there is a parse error, such as reading an unexpected character from an input stream.
struct parse_error : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit parse_error(auto &&... s) : bio_error{s...} {}
};

//!\brief Thrown if there is an io error in low level io operations such as in std::basic_streambuf operations.
struct io_error : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit io_error(auto &&... s) : bio_error{s...} {}
};

// ----------------------------------------------------------------------------
// parse exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if I/O was expecting more input (e.g. a delimiter or a new line), but the end of input was reached.
struct unexpected_end_of_input : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit unexpected_end_of_input(auto &&... s) : bio_error{s...} {}
};

// ----------------------------------------------------------------------------
// write exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if information given to output format didn't match expectations.
struct format_error : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit format_error(auto &&... s) : bio_error{s...} {}
};

//!\brief Thrown if a writer requires a header but it isn't provided.
struct missing_header_error : bio_error
{
    //!\brief Constructor that forwards the exception string.
    explicit missing_header_error(auto &&... s) : bio_error{s...} {}
};

//!\}

} // namespace bio
