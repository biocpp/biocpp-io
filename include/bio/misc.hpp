// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/b.i.o./blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <concepts>
#include <tuple>

#include <seqan3/utility/type_list/type_list.hpp>

#include <bio/platform.hpp>

/*!\file
 * \brief Provides various minor utilities.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

namespace bio
{

//-----------------------------------------------------------------------------
// ownership
//-----------------------------------------------------------------------------

/*!\brief An enum used as an argument for templates that switch between owning and non-owning behaviour.
 * \ingroup bio
 * \details
 *
 * Typically used to configure a class template to have members that are vectors/strings VS members that are views.
 * The "shallow" version of such a class is typically cheap to copy (no dynamic memory) while the "deep" version
 * is expensive to copy (holds dynamic memory).
 */
enum class ownership
{
    shallow, //!< Cheap to copy.
    deep     //!< Expensive to copy.
};

//-----------------------------------------------------------------------------
// vtag
//-----------------------------------------------------------------------------

/*!\brief The type of bio::vtag. [Default "specialisation" for 0 arguments.]
 * \tparam more_vs Any number of values [only 0 arguments pick this specialisation].
 * \ingroup bio
 * \see bio::vtag
 */
template <auto... more_vs>
struct vtag_t
{
    //!\brief The number of values stored in the tag.
    static constexpr size_t size = 0;

    //!\brief The tag converted to a tuple.
    static constexpr auto as_tuple = std::tuple{};

    //!\brief A function that checks if a value is contained in the tag.
    static constexpr bool contains(auto &&) { return false; }

    //!\brief A function that returns the index of a value or ((size_t)-1) if the value is not found.
    static constexpr size_t index_of(auto &&) { return static_cast<size_t>(-1ULL); }
};

/*!\brief The type of bio::vtag. [Specialisation for 1 or more arguments]
 * \tparam v       First value.
 * \tparam more_vs More values.
 * \ingroup bio
 * \see bio::vtag
 */
template <auto v, auto... more_vs>
struct vtag_t<v, more_vs...>
{
    //!\brief The first value in the tag.
    static constexpr auto first_value = v;

    //!\copybrief bio::vtag_t::size
    static constexpr size_t size = sizeof...(more_vs) + 1;

    //!\copybrief bio::vtag_t::as_tuple
    static constexpr auto as_tuple = std::tuple{v, more_vs...};

    //!\brief Whether all values in the tag are unique.
    static constexpr bool unique_values = ((v != more_vs) && ...);

    //!\copybrief bio::vtag_t::contains
    static constexpr bool contains(auto && s) requires std::equality_comparable_with<decltype(s), decltype(v)> &&
      (std::equality_comparable_with<decltype(s), decltype(more_vs)> &&...)
    {
        return s == v || ((s == more_vs) || ...);
    }

    //!\copybrief bio::vtag_t::index_of
    static constexpr size_t index_of(auto && s) requires std::equality_comparable_with<decltype(s), decltype(v)> &&
      (std::equality_comparable_with<decltype(s), decltype(more_vs)> &&...)
    {
        size_t c = 0;
        ((v != s && ++c) && ((more_vs != s && ++c) && ...));
        return c >= size ? static_cast<size_t>(-1ULL) : c;
    }
};

/*!\brief A value-tag template.
 * \tparam vs The values to store in the tag.
 * \ingroup bio
 * \details
 *
 * Using this template, you can easily turn a value, e.g. a literal value, into a compile-time constant with a unique
 * type.
 *
 * ### Example
 *
 * \snippet test/snippet/snippet_tag.cpp vtag
 */
template <auto... vs>
inline constexpr vtag_t<vs...> vtag{};

//-----------------------------------------------------------------------------
// ttag
//-----------------------------------------------------------------------------

/*!\brief A type-tag template.
 * \tparam type The first type to store.
 * \tparam more_types More types to store (optional).
 * \ingroup bio
 * \see seqan3::type_list
 *
 * \details
 *
 * Using this template, you can easily turn a type into a compile-time constant (value).
 *
 * ### Example
 *
 * \snippet test/snippet/snippet_tag.cpp ttag
 */
template <typename type, typename... more_types>
inline constexpr seqan3::type_list<type, more_types...> ttag{};

} // namespace bio

namespace bio::detail
{

//!\brief Check whether a type is a specialisation of seqan3::type_list.
template <typename t>
inline constexpr bool is_type_list = false;

//!\brief Check whether a type is a specialisation of seqan3::type_list.
template <typename... ts>
inline constexpr bool is_type_list<seqan3::type_list<ts...>> = true;

} // namespace bio::detail
