// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2022, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides functions for converting strings to numbers (and vice-versa).
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <cassert>
#include <charconv>
#include <concepts>
#include <ranges>
#include <string>
#include <string_view>

#include <bio/meta/concept/core_language.hpp>

#include <bio/io/detail/to_string.hpp>
#include <bio/io/exception.hpp>

namespace bio::io::detail
{

/*!\addtogroup io
 * \{
 */

/*!\brief Turn a string into a number.
 * \param[in] input The input string.
 * \param[out] number The variable holding the result.
 * \throws bio_error If there was an error during conversion.
 * \details
 *
 * Relies on std::from_chars to efficiently convert but accepts std::string_view and throws on error so
 * there is no return value that needs to be checked.
 *
 * TODO make this public (bio::io::) since it is useful for people doing plain IO
 */
void string_to_number(std::string_view const input, meta::arithmetic auto & number)
{
    std::from_chars_result res = std::from_chars(input.data(), input.data() + input.size(), number);
    if (res.ec != std::errc{} || res.ptr != input.data() + input.size())
        throw bio_error{"Could not convert \"", input, "\" into a number."};
}

// TODO write number_to_string and append_number_to_string

//!\}

} // namespace bio::io::detail
