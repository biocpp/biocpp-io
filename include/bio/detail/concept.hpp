// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides miscellaneous utilities.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <concepts>

#include <seqan3/alphabet/concept.hpp>

#include <bio/detail/utility.hpp>

namespace bio::detail
{

/*!\addtogroup bio
 * \{
 */

/*!\interface   bio::detail::one_of <>
 * \tparam t    The query type to compare.
 * \tparam ts   The reference types.
 * \brief       Checks whether the list of reference types contains the query type.
 */
//!\cond
template <typename t, typename... ts>
concept one_of = (std::same_as<t, ts> || ...);
//!\endcond

/*!\interface   bio::detail::deliberate_alphabet <>
 * \tparam t    The query type to compare.
 * \brief       A seqan3::alphabet that is **not** a character or a number (any std::integral).
 */
//!\cond
template <typename t>
concept deliberate_alphabet = seqan3::alphabet<t> && !std::integral<std::remove_cvref_t<t>>;
//!\endcond

/*!\interface   bio::detail::decays_to <>
 * \tparam t    The type to check.
 * \brief       Shortcut for `std::same_as<std::decay_t<from_t>, to_t>`.
 */
//!\cond
template <typename from_t, typename to_t>
concept decays_to = std::same_as<std::decay_t<from_t>, to_t>;
//!\endcond

/*!\brief Pass this function a constrained functor that accepts one argument and returns std::true_type.
 * \details
 *
 * See e.g. bio::seq_io::reader_options to see how this is used.
 */
constexpr bool lazy_concept_checker(auto fun)
{
    auto fallback = []<typename T = int>(auto)
    {
        return std::false_type{};
    };
    using ret_t = decltype(detail::overloaded{fallback, fun}(1));
    return ret_t::value;
}

//!\}

} // namespace bio::detail
