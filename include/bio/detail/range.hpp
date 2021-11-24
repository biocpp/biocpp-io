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

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/concept/container.hpp>

#include <bio/platform.hpp>

namespace bio::detail
{

/*!\addtogroup bio
 * \{
 */

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
void sized_range_copy(std::ranges::input_range auto &&                                                   in,
                      seqan3::back_insertable_with<std::ranges::range_reference_t<decltype(in)>> auto && out)
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
void string_copy(std::string_view const in, auto & out)
{
    if constexpr (std::same_as<decltype(out), std::string_view &>)
        out = in;
    else
        sized_range_copy(in, out);
}

//!\brief A seqan3::alphabet that is **not** a character or number (any std::integral).
template <typename t>
concept deliberate_alphabet = seqan3::alphabet<t> && !std::integral<std::remove_cvref_t<t>>;

//!\brief A range whose value type is `char`.
template <typename t>
concept char_range = std::ranges::range<t> && std::same_as<char, std::remove_cvref_t<std::ranges::range_value_t<t>>>;

//!\brief A range whose value type is an integral type other than `char`.
template <typename t>
concept int_range = std::ranges::range<t> && std::integral<std::remove_cvref_t<std::ranges::range_value_t<t>>> &&
  !std::same_as<char, std::remove_cvref_t<std::ranges::range_value_t<t>>>;

//!\brief A type that is not std::span<std::byte const>.
template <typename t>
concept not_a_byte_span = !std::same_as<t, std::span<std::byte const>>;

//!\}

} // namespace bio::detail
