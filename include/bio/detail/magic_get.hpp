// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/bio/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a very simple abstraction of tuples of two and aggregates of two; inspired by the magic_get library.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <concepts>

#include <seqan3/alphabet/concept.hpp>

#include <bio/platform.hpp>

namespace bio::detail
{

/*!\interface   bio::detail::tuple_of_two <>
 * \ingroup bio
 * \tparam t    The type to check.
 * \brief       A tuple/pair with exactly two elements.
 */
//!\cond
template <typename t>
concept tuple_of_two = (requires { std::tuple_size<t>::value; }) && (std::tuple_size<t>::value == 2);

//!\endcond

/*!\brief A type that can be converted to any reference type.
 * \ingroup bio
 * \details
 *
 * **NEVER USE THIS IN ANY EVALUATED CONTEXT!**
 *
 * It contains a unsafe reimplementation of std::declval() usable in concepts. If this is used in an evaluated context
 * it will always lead to an unallowed memory read.
 */
struct converts_to_any_lref
{
    //!\brief The conversion operator.
    template <typename t>
    constexpr operator t &() const
    {
        std::remove_reference_t<t> * ptr = nullptr;
        ptr += 13117;
        return static_cast<t &>(*ptr);
    }
};

/*!\interface   bio::detail::aggregate_of_two <>
 * \ingroup bio
 * \tparam t    The type to check.
 * \brief       A aggregate type with exactly two elements.
 */
//!\cond
// clang-format off
template <typename t>
concept aggregate_of_two = std::is_aggregate_v<t> &&
    (requires { t{converts_to_any_lref{}, converts_to_any_lref{}};                          }) &&
   !(requires { t{converts_to_any_lref{}, converts_to_any_lref{}, converts_to_any_lref{}};  });
// clang-format on
//!\endcond

/*!\interface   bio::detail::decomposable_into_two <>
 * \ingroup bio
 * \tparam t    The type to check.
 * \brief       Either bio::detail::tuple_of_two or bio::detail::aggregate_of_two. CVREF is removed.
 */
//!\cond
template <typename t>
concept decomposable_into_two = tuple_of_two<std::remove_cvref_t<t>> || aggregate_of_two<std::remove_cvref_t<t>>;
//!\endcond

/*!\brief Get the first element.
 * \ingroup bio
 * \see bio::detail::decomposable_into_two
 */
auto & get_first(decomposable_into_two auto & val)
{
    auto & [first, second] = val;
    return first;
}

/*!\brief Type of the first element with CVREF removed.
 * \ingroup bio
 * \see bio::detail::decomposable_into_two
 */
template <typename t>
using first_elem_t = std::remove_cvref_t<decltype(get_first(std::declval<t &>()))>;

/*!\brief Get the second element.
 * \ingroup bio
 * \see bio::detail::decomposable_into_two
 */
auto & get_second(decomposable_into_two auto & val)
{
    auto & [first, second] = val;
    return second;
}

/*!\brief Type of the second element with CVREF removed.
 * \ingroup bio
 * \see bio::detail::decomposable_into_two
 */
template <typename t>
using second_elem_t = std::remove_cvref_t<decltype(get_second(std::declval<t &>()))>;

//!\brief Overload of bio::detail::range_or_tuple_size for aggregates of two.
//!\ingroup bio
constexpr size_t range_or_tuple_size(aggregate_of_two auto &&)
{
    return 2;
}

} // namespace bio::detail
