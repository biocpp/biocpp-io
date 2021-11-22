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

#include <cassert>
#include <charconv>
#include <concepts>
#include <ranges>
#include <string_view>

#include <bio/platform.hpp>

namespace bio::detail
{

/*!\addtogroup bio
 * \{
 */

//!\brief Wrapper around standard library std::to_chars with fallback for floats on GCC10.
std::to_chars_result to_chars(char * first, char * last, auto in)
{
#if defined(__cpp_lib_to_chars) && (__cpp_lib_to_chars >= 201611)
    constexpr static bool float_support = true;
#else
    constexpr static bool float_support = false;
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
        return {.ptr = first + count, .ec{}};
    }
}
//!\}

} // namespace bio::detail
