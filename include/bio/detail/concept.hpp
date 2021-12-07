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

#include <bio/platform.hpp>

namespace bio::detail
{

/*!\addtogroup bio
 * \{
 */

/*!\interface   bio::one_of <>
 * \tparam t    The query type to compare.
 * \tparam ts   The reference types.
 * \brief       Checks whether the list of reference types contains the query type.
 */
//!\cond
template <typename t, typename... ts>
concept one_of = (std::same_as<t, ts> || ...);
//!\endcond

/*!\interface   bio::deliberate_alphabet <>
 * \tparam t    The query type to compare.
 * \brief       A seqan3::alphabet that is **not** a character or a number (any std::integral).
 */
//!\cond
template <typename t>
concept deliberate_alphabet = seqan3::alphabet<t> && !std::integral<std::remove_cvref_t<t>>;
//!\endcond

//!\}

} // namespace bio::detail
