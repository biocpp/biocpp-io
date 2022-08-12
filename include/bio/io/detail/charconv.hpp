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

#include <seqan3/utility/concept/exposition_only/core_language.hpp>

#include <bio/io/detail/to_string.hpp>
#include <bio/io/exception.hpp>

namespace bio::detail
{

/*!\addtogroup bio
 * \{
 */

//!\brief Wrapper around standard library std::to_chars with fallback for floats on GCC10.
std::to_chars_result to_chars(char * first, char * last, auto in)
{
#if defined(__cpp_lib_to_chars) && (__cpp_lib_to_chars >= 201611)
    static constexpr bool float_support = true;
#else
    static constexpr bool float_support = false;
#endif

    if constexpr (std::integral<decltype(in)> || float_support)
    {
        return std::to_chars(first, last, in);
    }
    else // very hacky version
    {
        int count = std::sprintf(first, "%f", in);
        if (count < 0)
            return {.ptr = nullptr, .ec = std::errc::io_error};

        assert(first + count <= last);

        if (std::string_view{first, static_cast<size_t>(count)}.find('.') != std::string_view::npos)
            while (count > 0 && first[count - 1] == '0')
                --count;
        if (count > 0 && first[count - 1] == '.')
            --count;
        return {.ptr = first + count, .ec{}};
    }
}

//!\brief Wrapper around standard library std::from_chars with fallback for floats on GCC10.
std::from_chars_result from_chars(char const * first, char const * last, auto & out)
{
#if defined(__cpp_lib_to_chars) && (__cpp_lib_to_chars >= 201611)
    static constexpr bool float_support = true;
#else
    static constexpr bool float_support = false;
#endif

    // TODO always use fast_float here.
    if constexpr (std::integral<std::remove_cvref_t<decltype(out)>> || float_support)
    {
        return std::from_chars(first, last, out);
    }
    else
    {
        // THIS IS KINDA CRAZY BUT MAYBE NOT TOO BAD WITH SSO:
        std::string buffer(first, last);
        size_t      count = 0;
        if constexpr (std::same_as<float &, decltype(out)>)
            out = std::stof(buffer, &count);
        else // double
            out = std::stod(buffer, &count);

        return {.ptr = first + count, .ec{}};
    }
}

/*!\brief Turn a string into a number.
 * \param[in] input The input string.
 * \param[out] number The variable holding the result.
 * \throws bio_error If there was an error during conversion.
 * \details
 *
 * Relies on std::from_chars to efficiently convert but accepts std::string_view and throws on error so
 * there is no return value that needs to be checked.
 *
 * TODO make this public (bio::) since it is useful for people doing plain IO
 */
void string_to_number(std::string_view const input, seqan3::arithmetic auto & number)
{
    std::from_chars_result res = from_chars(input.data(), input.data() + input.size(), number);
    if (res.ec != std::errc{} || res.ptr != input.data() + input.size())
        throw bio_error{"Could not convert \"", input, "\" into a number."};
}

// TODO write number_to_string and append_number_to_string

//!\}

} // namespace bio::detail
