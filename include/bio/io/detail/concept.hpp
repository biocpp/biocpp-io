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

#include <bio/alphabet/concept.hpp>
#include <bio/meta/overloaded.hpp>

#include <bio/io/detail/utility.hpp>

namespace bio::io::detail
{

/*!\addtogroup io
 * \{
 */

/*!\interface   bio::io::detail::deliberate_alphabet <>
 * \tparam t    The query type to compare.
 * \brief       A bio::alphabet::alphabet that is **not** a character or a number (any std::integral).
 */
//!\cond
template <typename t>
concept deliberate_alphabet = alphabet::alphabet<t> && !std::integral<std::remove_cvref_t<t>>;
//!\endcond

/*!\brief Pass this function a constrained functor that accepts one argument and returns std::true_type.
 * \details
 *
 * See e.g. bio::io::seq::reader_options to see how this is used.
 */
constexpr bool lazy_concept_checker(auto fun)
{
    auto fallback = []<typename t = int>(auto)
    {
        return std::false_type{};
    };
    using ret_t = decltype(meta::overloaded{fallback, fun}(1));
    return ret_t::value;
}

//!\}

} // namespace bio::io::detail
