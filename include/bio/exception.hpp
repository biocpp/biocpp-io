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

#include <bio/platform.hpp>

namespace bio
{

/*!\addtogroup bio
 * \{
 */

// ----------------------------------------------------------------------------
// file open exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if there is no format that accepts a given file extension.
struct unhandled_extension_error : std::invalid_argument
{
    //!\brief Constructor that forwards the exception string.
    unhandled_extension_error(std::string const & s) : std::invalid_argument{s} {}
};

//!\brief Thrown if there is an unspecified filesystem or stream error while opening, e.g. permission problem.
struct file_open_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    file_open_error(std::string const & s) : std::runtime_error{s} {}
};

//!\brief Thrown if there is a parse error, such as reading an unexpected character from an input stream.
struct parse_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    parse_error(std::string const & s) : std::runtime_error{s} {}
};

//!\brief Thrown if there is an io error in low level io operations such as in std::basic_streambuf operations.
struct io_error : std::ios_base::failure
{
    //!\brief Constructor that forwards the exception string.
    explicit io_error(std::string const & s, std::error_code const & ec = std::io_errc::stream) :
      std::ios_base::failure{s, ec}
    {}
};

// ----------------------------------------------------------------------------
// parse exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if I/O was expecting more input (e.g. a delimiter or a new line), but the end of input was reached.
struct unexpected_end_of_input : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    unexpected_end_of_input(std::string const & s) : std::runtime_error{s} {}
};

// ----------------------------------------------------------------------------
// write exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if information given to output format didn't match expectations.
struct format_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    format_error(std::string const & s) : std::runtime_error{s} {}
};

//!\brief Thrown if a writer requires a header but it isn't provided.
struct missing_header_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    missing_header_error(std::string const & s) : std::runtime_error{s} {}
};

//!\}

} // namespace bio
