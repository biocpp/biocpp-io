// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2022, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides bio::detail::to_string().
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <string_view>

#include <seqan3/utility/concept/exposition_only/core_language.hpp>

#include <bio/platform.hpp>

namespace bio::detail
{

/*!\brief Convert something to a string.
 * \details
 *
 * ### Attention
 *
 * This function is not efficient. Do not use it in performance-critical code!
 */
std::string to_string(auto && in)
{
    using in_t = std::remove_cvref_t<decltype(in)>;
    if constexpr (seqan3::arithmetic<in_t>)
        return std::to_string(in);
    else if constexpr (std::same_as<in_t, std::string>)
        return in;
    else if constexpr (std::constructible_from<std::string, in_t>)
        return std::string{in};
    else
        static_assert(std::constructible_from<std::string, in_t>, "Type cannot be converted to string");
}

} // namespace bio::detail
