// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various range utilities.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <ranges>
#include <span>

#include <bio/platform.hpp>

namespace bio::detail
{

/*!\addtogroup bio
 * \{
 */

// ----------------------------------------------------------------------------
// concepts
// ----------------------------------------------------------------------------

/*!\interface bio::detail::back_insertable_with <>
 * \extends std::ranges::output_range
 * \tparam rng_t The container type.
 * \tparam val_t The type to append to the container.
 * \brief Describes range types that can grow in amortised constant time by appending an element of type val_t.
 */
//!\cond
template <typename rng_t, typename val_t>
concept back_insertable_with = std::ranges::output_range<rng_t, val_t> && requires(rng_t & v)
{
    v.push_back(std::declval<val_t>());
};
//!\endcond

/*!\interface bio::detail::back_insertable <>
 * \extends std::ranges::output_range
 * \extends std::ranges::input_range
 * \tparam rng_t The container type.
 * \brief Describes range types that can grow in amortised constant time by appending an element.
 */
//!\cond
template <typename rng_t>
concept back_insertable =
  std::ranges::input_range<rng_t> && back_insertable_with<rng_t, std::ranges::range_reference_t<rng_t>>;
//!\endcond

//!\brief A range whose value type is `char`.
template <typename t>
concept char_range =
  std::ranges::forward_range<t> && std::same_as<char, std::remove_cvref_t<std::ranges::range_value_t<t>>>;

//!\brief A range whose value type is `char` or a C-String.
template <typename t>
concept char_range_or_cstring = char_range<t> || std::same_as < std::decay_t<t>,
char const * > ;

//!\brief A range whose value type is an integral type other than `char`.
template <typename t>
concept int_range =
  std::ranges::forward_range<t> && std::integral<std::remove_cvref_t<std::ranges::range_value_t<t>>> &&
  !std::same_as<char, std::remove_cvref_t<std::ranges::range_value_t<t>>>;

/*!\interface bio::detail::out_string <>
 * \tparam rng_t The container type.
 * \brief A range that `char` can be back-inserted to, or a string_view.
 */
//!\cond
template <typename rng_t>
concept out_string = back_insertable_with<rng_t, char> || std::same_as<rng_t &, std::string_view &>;
//!\endcond

template <typename rng_t>
concept vector_like = std::ranges::random_access_range<rng_t> && std::ranges::sized_range<rng_t> &&
  std::ranges::output_range<rng_t, std::ranges::range_reference_t<rng_t>> && requires(rng_t & v)
{
    v.resize(3);
    v.clear();
};

//!\brief Helper for bio::detail::transform_view_on_string_view.
template <std::regular_invocable<char> fun_t>
void transform_view_on_string_view_impl(std::ranges::transform_view<std::string_view, fun_t> &)
{}

//!\brief Helper for bio::detail::transform_view_on_string_view.
template <std::regular_invocable<char> fun1_t, std::regular_invocable<char> fun2_t>
void transform_view_on_string_view_impl(
  std::ranges::transform_view<std::ranges::transform_view<std::string_view, fun1_t>, fun2_t> &)
{}

/*!\interface   bio::detail::transform_view_on_string_view <>
 * \tparam t    The query type to check.
 * \brief       Whether a type is a transform view (possibly nested) over std::string_view.
 *
 * Types must be one of:
 *   1. std::ranges::transform_view<std::string_view, fun_t>  && std::regular_invocable<fun_t, char>
 *   2. std::ranges::transform_view<std::ranges::transform_view<std::string_view, fun1_t>, fun2_t> &&
 * std::regular_invocable<fun1_t, char> && std::regular_invocable<fun2_t, char>
 */
//!\cond
template <typename t>
concept transform_view_on_string_view = requires
{
    transform_view_on_string_view_impl(std::declval<t &>());
};
//!\endcond

// ----------------------------------------------------------------------------
// copy functions
// ----------------------------------------------------------------------------

/*!\brief Copy elements from the first range into the second range.
 * \param[in] in The range to copy from.
 * \param[out] out The range to copy to.
 * \details
 *
 * By default, this function behaves like std::ranges::copy with a std::back_insert_iterator. It assumes that
 * the target range is empty (and clears it if possible).
 *
 * If the input range is sized and the target range offers a `.resize()` member, this function uses
 * resize and assignment instead of back-insertion.
 */
void sized_range_copy(std::ranges::input_range auto &&                                           in,
                      back_insertable_with<std::ranges::range_reference_t<decltype(in)>> auto && out)
{
    using in_t  = decltype(in);
    using out_t = decltype(out);

    if constexpr (std::ranges::sized_range<in_t> && requires(out_t out) { out.resize(0); })
    {
        out.resize(std::ranges::size(in));
        std::ranges::copy(in, std::ranges::begin(out));
    }
    else
    {
        if constexpr (requires(out_t out) { out.clear(); })
            out.clear();
        std::ranges::copy(in, std::back_inserter(out));
    }
}

//!\brief Like bio::detail::sized_range_copy except that if input and output are both std::string_view, it assigns.
void string_copy(std::string_view const in, out_string auto & out)
{
    if constexpr (std::same_as<decltype(out), std::string_view &>)
        out = in;
    else
        sized_range_copy(in, out);
}

// ----------------------------------------------------------------------------
// to_string_view
// ----------------------------------------------------------------------------

//!\brief Turn something string-viewable into a std::string_view.
constexpr std::string_view to_string_view(char const * const cstring)
{
    return std::string_view{cstring};
}

//!\overload
template <char_range rng_t>
    requires std::ranges::borrowed_range<rng_t> && std::ranges::contiguous_range<rng_t> &&
      std::ranges::sized_range<rng_t>
constexpr std::string_view to_string_view(rng_t && contig_range)
{
    return std::string_view{std::ranges::data(contig_range), std::ranges::size(contig_range)};
}

//!\overload
constexpr std::string_view to_string_view(std::string_view const in)
{
    return in;
}

// ----------------------------------------------------------------------------
// range_or_tuple_size
// ----------------------------------------------------------------------------

//!\brief Returns the size of the argument, either a range or a tuple.
constexpr size_t range_or_tuple_size(std::ranges::forward_range auto && r)
{
    return std::ranges::distance(r);
}

//!\overload
template <typename tuple_t>
    requires(requires { typename std::tuple_size<std::remove_cvref_t<tuple_t>>::type; })
constexpr size_t range_or_tuple_size(tuple_t)
{
    return std::tuple_size_v<tuple_t>;
}

// there is another overload for this in magic_get.hpp

//!\}

} // namespace bio::detail
